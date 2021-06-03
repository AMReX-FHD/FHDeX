#include "DsmcParticleContainer.H"

//#include "particle_functions_K.H"
#include "paramplane_functions_K.H"
#include <math.h>


FhdParticleContainer::FhdParticleContainer(const Geometry & geom,
                              const DistributionMapping & dmap,
                              const BoxArray            & ba,
                              int ncells)
    : NeighborParticleContainer<FHD_realData::count, FHD_intData::count> (geom, dmap, ba, ncells)
{
    BL_PROFILE_VAR("FhdParticleContainer()",FhdParticleContainer);

    realParticles = 0;
    simParticles = 0;

    for(int i=0;i<nspecies;i++) {
        
        properties[i].mass = mass[i];
        properties[i].radius = diameter[i]/2.0;
        properties[i].Neff = particle_neff; // From DSMC, this will be set to 1 for electolyte calcs
        properties[i].R = k_B/properties[i].mass; //used a lot in kinetic stats cals, bu not otherwise necessary for electrolytes
        properties[i].T = T_init[i];

        if (particle_count[i] >= 0) {

            properties[i].total = particle_count[i];
            properties[i].n0 = particle_neff*properties[i].total/domainVol;
            
            Print() << "Species " << i << " count " << properties[i].total << "\n";
            Print() << "Species " << i << " n0 " << properties[i].total << "\n";
        }
        else {
            properties[i].total = (int)amrex::Math::ceil(particle_n0[i]*domainVol/particle_neff);
            properties[i].n0 = properties[i].total/domainVol;

            Print() << "Species " << i << " count " << properties[i].total << "\n";
            Print() << "Species " << i << " n0 adjusted to " << properties[i].n0 << "\n";
        }

        Print() << "Species " << i << "\n";
        Print() << "Mass " << properties[i].mass << "\n";
        Print() << "Radius " << properties[i].radius << "\n";

        realParticles = realParticles + properties[i].total*particle_neff;
        simParticles = simParticles + properties[i].total;

    }

    totalCollisionCells = n_cells[0]*n_cells[1]*n_cells[2];
    domainVol = (prob_hi[0] - prob_lo[0])*(prob_hi[1] - prob_lo[1])*(prob_hi[2] - prob_lo[2]);

    for (int d=0; d<AMREX_SPACEDIM; ++d)
    {
        domSize[d] = prob_hi[d] - prob_lo[d];
    }

    Print() << "Total real particles: " << realParticles << "\n";
    Print() << "Total sim particles: " << simParticles << "\n";

    Print() << "Collision cells: " << totalCollisionCells << "\n";
    Print() << "Sim particles per cell: " << simParticles/totalCollisionCells << "\n";


    int lev=0;
    for(MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const int grid_id = mfi.index();
        m_cell_vectors[grid_id].resize(box.numPts());
    }


}

void FhdParticleContainer::MoveParticlesCPP(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount)
{
    BL_PROFILE_VAR("MoveParticlesCPP()", MoveParticlesCPP);

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();
    const Real* phi = Geom(lev).ProbHi();

    int        np_tile = 0 ,       np_proc = 0 ; // particle count
    Real    moves_tile = 0.,    moves_proc = 0.; // total moves
    Real  maxspeed_proc = 0.; // max speed
    Real  maxdist_proc = 0.; // max displacement (fraction of radius)

    Real adj = 0.99999;
    Real adjalt = 2.0*(1.0-0.99999);
    Real runtime, inttime;
    int intsurf, intside, push;

    Real maxspeed = 0;
    Real maxdist = 0;

    long moves = 0;
    int reDist = 0;

    for (FhdParIter pti(* this, lev); pti.isValid(); ++pti) {

        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const long np = particles.numParticles();

        Box bx  = pti.tilebox();
        IntVect myLo = bx.smallEnd();
        IntVect myHi = bx.bigEnd();

        np_proc += np;

        for (int i = 0; i < np; ++ i) {
            ParticleType & part = particles[i];

              Real speed = 0;

              for (int d=0; d<AMREX_SPACEDIM; ++d)
              {                   
                  speed += part.rdata(FHD_realData::velx + d)*part.rdata(FHD_realData::velx + d);
              }

              if(speed > maxspeed)
              {
                  maxspeed = speed;
              }

              moves++;

              runtime = dt;

              while(runtime > 0)
              {

                  find_inter_gpu(part, runtime, paramPlaneList, paramPlaneCount, &intsurf, &inttime, &intside, ZFILL(plo), ZFILL(phi));

                  for (int d=0; d<AMREX_SPACEDIM; ++d)
                  {
                      part.pos(d) += inttime * part.rdata(FHD_realData::velx + d)*adj;
                  }

                  runtime = runtime - inttime;

                  if(intsurf > 0)
                  {
                      const paramPlane& surf = paramPlaneList[intsurf-1];//find_inter indexes from 1 to maintain compatablity with fortran version
      
                      Real posAlt[3];

                      for (int d=0; d<AMREX_SPACEDIM; ++d)
                      {
                          posAlt[d] = inttime * part.rdata(FHD_realData::velx + d)*adjalt;
                      }

                      Real dummy = 1;

                      app_bc_gpu(&surf, part, intside, domSize, &push, dummy, dummy);

                      if(push == 1)
                      {
                          for (int d=0; d<AMREX_SPACEDIM; ++d)
                          {
                              part.pos(d) += part.pos(d) + posAlt[d];
                          }
                      }
                  }

              }


            int cell[3];
            cell[0] = (int)floor((part.pos(0)-plo[0])/dx[0]);
            cell[1] = (int)floor((part.pos(1)-plo[1])/dx[1]);
            cell[2] = (int)floor((part.pos(2)-plo[2])/dx[2]);

            //cout << "current cell: " << cell[0] << ", " << cell[1] << ", " << cell[2] << "\cout";
            //n << "current pos: " << part.pos(0) << ", " << part.pos(1) << ", " << part.pos(2) << "\n";

            if((cell[0] < myLo[0]) || (cell[1] < myLo[1]) || (cell[2] < myLo[2]) || (cell[0] > myHi[0]) || (cell[1] > myHi[1]) || (cell[2] > myHi[2]))
            {
                reDist++;
            }

            if((part.idata(FHD_intData::i) != cell[0]) || (part.idata(FHD_intData::j) != cell[1]) || (part.idata(FHD_intData::k) != cell[2]))
            {
                //remove particle from old cell
                IntVect iv(part.idata(FHD_intData::i), part.idata(FHD_intData::j), part.idata(FHD_intData::k));
                long imap = tile_box.index(iv);

                int lastIndex = m_cell_vectors[pti.index()][imap].size() - 1;
                int lastPart = m_cell_vectors[pti.index()][imap][lastIndex];
                int newIndex = part.idata(FHD_intData::sorted);

                m_cell_vectors[pti.index()][imap][newIndex] = m_cell_vectors[pti.index()][imap][lastIndex];

                m_cell_vectors[pti.index()][imap].pop_back();

                //cout << "Removed " << i << " from " << iv[0] << ", " << iv[1] << ", " << iv[2] << ", now contains " << m_cell_vectors[pti.index()][imap].size() << " particles\n";

                part.idata(FHD_intData::sorted) = -1;
            }
        }

        maxspeed_proc = amrex::max(maxspeed_proc, maxspeed);
        maxdist_proc  = amrex::max(maxdist_proc, maxdist);

    }

    // gather statistics
    ParallelDescriptor::ReduceIntSum(np_proc);
    ParallelDescriptor::ReduceRealSum(moves_proc);
    ParallelDescriptor::ReduceRealMax(maxspeed_proc);
    ParallelDescriptor::ReduceRealMax(maxdist_proc);
    ParallelDescriptor::ReduceIntSum(reDist);

    // write out global diagnostics
    if (ParallelDescriptor::IOProcessor()) {
        Print() << "I see " << np_proc << " particles\n";

        Print() << reDist << " particles to be redistributed.\n";
        Print() <<"Maximum observed speed: " << sqrt(maxspeed_proc) << "\n";
        Print() <<"Maximum observed displacement (fraction of radius): " << maxdist_proc << "\n";
    }

    Redistribute();
    SortParticles();

}

void FhdParticleContainer::EvaluateStats(MultiFab& particleInstant, MultiFab& particleMeans,
                                         MultiFab& particleVars, const Real delt, int steps)
{
    BL_PROFILE_VAR("EvaluateStats()",EvaluateStats);
    
    const int lev = 0;
    
    BoxArray ba = particleMeans.boxArray();
    long cellcount = ba.numPts();

    const Real* dx = Geom(lev).CellSize();
    const Real dxInv = 1.0/dx[0];
    const Real cellVolInv = 1.0/(dx[0]*dx[0]*dx[0]);

    const Real stepsInv = 1.0/steps;
    const int stepsMinusOne = steps-1;

    // zero instantaneous values
    particleInstant.setVal(0.);
    
    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) 
    {

	PairIndex index(pti.index(), pti.LocalTileIndex());
        const int np = this->GetParticles(lev)[index].numRealParticles();
	auto& plev = this->GetParticles(lev);
	auto& ptile = plev[index];
	auto& aos   = ptile.GetArrayOfStructs();
        const Box& tile_box  = pti.tilebox();
        ParticleType* particles = aos().dataPtr();

        GpuArray<int, 3> bx_lo = {tile_box.loVect()[0], tile_box.loVect()[1], tile_box.loVect()[2]};
        GpuArray<int, 3> bx_hi = {tile_box.hiVect()[0], tile_box.hiVect()[1], tile_box.hiVect()[2]};

        Array4<Real> part_inst = particleInstant[pti].array();
        Array4<Real> part_mean = particleMeans[pti].array();
        Array4<Real> part_var = particleVars[pti].array();

        AMREX_FOR_1D( np, ni,
        {
            ParticleType & part = particles[ni];

            int i = floor(part.pos(0)*dxInv);
            int j = floor(part.pos(1)*dxInv);
            int k = floor(part.pos(2)*dxInv);
            
            amrex::Gpu::Atomic::Add(&part_inst(i,j,k,0), 1.0);

            for(int l=0;l<nspecies;l++)
            {

            }

        });

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            part_inst(i,j,k,1) = part_inst(i,j,k,1)*cellVolInv;
             
        });

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            part_mean(i,j,k,1)  = (part_mean(i,j,k,1)*stepsMinusOne + part_inst(i,j,k,1))*stepsInv;

            for(int l=0;l<nspecies;l++)
            {

            }
             
        });

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            Real del1 = part_inst(i,j,k,1) - part_mean(i,j,k,1);

            part_var(i,j,k,1)  = (part_var(i,j,k,1)*stepsMinusOne + del1*del1)*stepsInv;            
        });

    }

}

void FhdParticleContainer::SortParticles()
{
    int lev = 0;
    for (FhdParIter pti(* this, lev); pti.isValid(); ++pti) {

        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const long np = particles.numParticles();

        for (int i = 0; i < np; ++ i)
        {
            ParticleType & part = particles[i];
            if(part.idata(FHD_intData::sorted) == -1)
            {
                const IntVect& iv = this->Index(part, lev);

                part.idata(FHD_intData::i) = iv[0];
                part.idata(FHD_intData::j) = iv[1];
                part.idata(FHD_intData::k) = iv[2];

                long imap = tile_box.index(iv);
                //cout << "part " << i << " is in cell " << iv[0] << ", " << iv[1] << ", " << iv[2] << ", adding to element " << m_cell_vectors[pti.index()][imap].size() << "\n";


                part.idata(FHD_intData::sorted) = m_cell_vectors[pti.index()][imap].size();

                m_cell_vectors[pti.index()][imap].push_back(i);

            }

        }
    }
}

