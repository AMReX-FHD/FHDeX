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
        for(int i=0;i<nspecies;i++)
        {
            m_cell_vectors[i][grid_id].resize(box.numPts());
        }
    }

}

void FhdParticleContainer::MoveParticlesCPP(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount)
{
    BL_PROFILE_VAR("MoveParticlesCPP()", MoveParticlesCPP);

    const int lev = 0;
    const GpuArray<Real, 3> dx = Geom(lev).CellSizeArray();
    const GpuArray<Real, 3> plo = Geom(lev).ProbLoArray();
    const GpuArray<Real, 3> phi = Geom(lev).ProbHiArray();


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

    }

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

        for (int i = 0; i < np; i++) {
       
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

              runtime = dt*part.rdata(FHD_realData::timeFrac);

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
                       //Print() << "surf: " << intsurf-1 << "\n";
                      //app_bc_gpu(&surf, part, intside, domSize, &push, &runtime, dummy);
                      amrex::RandomEngine engine;
                      app_bc_gpu(&surf, part, intside, domSize, &push, &runtime, dummy, engine);
                       //Print() << "rt: " << runtime << "\n";

                      if(push == 1)
                      {
                          for (int d=0; d<AMREX_SPACEDIM; ++d)
                          {
                              part.pos(d) += part.pos(d) + posAlt[d];
                          }
                      }
                  }

              }
 
              part.rdata(FHD_realData::timeFrac) = 1;


            int cell[3];
            cell[0] = (int)floor((part.pos(0)-plo[0])/dx[0]);
            cell[1] = (int)floor((part.pos(1)-plo[1])/dx[1]);
            cell[2] = (int)floor((part.pos(2)-plo[2])/dx[2]);


            if((part.idata(FHD_intData::i) != cell[0]) || (part.idata(FHD_intData::j) != cell[1]) || (part.idata(FHD_intData::k) != cell[2]) || part.id() < 0)
            {
                //remove particle from old cell

                IntVect iv(part.idata(FHD_intData::i), part.idata(FHD_intData::j), part.idata(FHD_intData::k));
                long imap = tile_box.index(iv);


                int lastIndex = m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].size() - 1;

                int lastPart = m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap][lastIndex];


                int newIndex = part.idata(FHD_intData::sorted);


                m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap][newIndex] = lastPart;
                m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].pop_back();

                particles[lastPart].idata(FHD_intData::sorted) = newIndex;

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

    }

    Redistribute();
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


    }

    SortParticles();


}

void FhdParticleContainer::SortParticles()
{
    int lev = 0;

    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();
    const Real* phi = Geom(lev).ProbHi();

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

            IntVect iv ={0,0,0};


            if(part.idata(FHD_intData::sorted) == -1)
            {
                iv[0] = (int)floor((part.pos(0)-plo[0])/dx[0]);
                iv[1] = (int)floor((part.pos(1)-plo[1])/dx[1]);
                iv[2] = (int)floor((part.pos(2)-plo[2])/dx[2]);

                part.idata(FHD_intData::i) = iv[0];
                part.idata(FHD_intData::j) = iv[1];
                part.idata(FHD_intData::k) = iv[2];

                long imap = tile_box.index(iv);

                part.idata(FHD_intData::sorted) = m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].size();

                m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].push_back(i);

            }else
            {

                iv[0] = part.idata(FHD_intData::i);
                iv[1] = part.idata(FHD_intData::j);
                iv[2] = part.idata(FHD_intData::k);

                long imap = tile_box.index(iv);

                m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap][part.idata(FHD_intData::sorted)] = i;
            }

        }
    }
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

void FhdParticleContainer::Source(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount)
{


    int lev = 0;
    bool proc0_enter = true;
    //Do this all on rank 0 for now
    for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi)
    {
        if(ParallelDescriptor::MyProc() == 0 && proc0_enter)
        {
            proc0_enter = false;

            const int grid_id = mfi.index();
            const int tile_id = mfi.LocalTileIndex();
            auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

            for(int i = 0; i< paramPlaneCount; i++)
            {
                if(paramPlaneList[i].sourceLeft == 1)
                {
                    for(int j = 0; j< nspecies; j++)
                    {
                        Real density = paramPlaneList[i].densityLeft[j];
                        Real temp = paramPlaneList[i].temperatureLeft;
                        Real area = paramPlaneList[i].area;

                        Real fluxMean = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;
                        Real fluxVar = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;

                        Real totalFlux = dt*fluxMean + sqrt(dt*fluxVar)*amrex::RandomNormal(0.,1.);

                        //Print() << "Flux mean " << dt*fluxMean << ", flux sd " << sqrt(dt*fluxVar) << "\n";

                        int totalFluxInt =  (int)floor(totalFlux);
                        Real totalFluxLeftOver = totalFlux - totalFluxInt;

                        if(amrex::Random() < totalFluxLeftOver)
                        {
                            totalFluxInt++;
                        }

                        Print() << "Surface " << i << " generating " << totalFluxInt << " of species " << j << "\n";

                        for(int k=0;k<totalFluxInt;k++)
                        {
                            Real uCoord = amrex::Random()*paramPlaneList[i].uTop;
                            Real vCoord = amrex::Random()*paramPlaneList[i].vTop;

                            ParticleType p;
                            p.id() = ParticleType::NextID();

                            p.cpu() = ParallelDescriptor::MyProc();
                            p.idata(FHD_intData::sorted) = -1;

                            p.idata(FHD_intData::species) = j;

                            p.pos(0) = paramPlaneList[i].x0 + paramPlaneList[i].ux*uCoord + paramPlaneList[i].vx*vCoord;
                            p.pos(1) = paramPlaneList[i].y0 + paramPlaneList[i].uy*uCoord + paramPlaneList[i].vy*vCoord;
                            p.pos(2) = paramPlaneList[i].z0 + paramPlaneList[i].uz*uCoord + paramPlaneList[i].vz*vCoord;

                            //Print() << "origin: " << paramPlaneList[i].x0 << ", " << paramPlaneList[i].y0 << ", " << paramPlaneList[i].z0 << ", uCoord: " << uCoord << "\n";

                            //move the particle slightly off the surface so it doesn't intersect it when it moves
                            p.pos(0) = p.pos(0) + uCoord*0.00000001*paramPlaneList[i].lnx;
                            p.pos(1) = p.pos(1) + uCoord*0.00000001*paramPlaneList[i].lny;
                            p.pos(2) = p.pos(2) + uCoord*0.00000001*paramPlaneList[i].lnz;

                            p.rdata(FHD_realData::boostx) = 0;
                            p.rdata(FHD_realData::boosty) = 0;
                            p.rdata(FHD_realData::boostz) = 0;

                            p.idata(FHD_intData::i) = -100;
                            p.idata(FHD_intData::j) = -100;
                            p.idata(FHD_intData::k) = -100;

                            p.rdata(FHD_realData::R) = properties[j].R;
                            p.rdata(FHD_realData::timeFrac) = amrex::Random();

                            Real srt = sqrt(p.rdata(FHD_realData::R)*temp);

                            p.rdata(FHD_realData::velx) = srt*amrex::RandomNormal(0.,1.);
                            p.rdata(FHD_realData::vely) = srt*amrex::RandomNormal(0.,1.);
                            p.rdata(FHD_realData::velz) = sqrt(2)*srt*sqrt(-log(amrex::Random()));

                            const paramPlane surf = paramPlaneList[i];

                            rotation(surf.cosThetaLeft, surf.sinThetaLeft, surf.cosPhiLeft, surf.sinPhiLeft, &p.rdata(FHD_realData::velx), &p.rdata(FHD_realData::vely), &p.rdata(FHD_realData::velz));

                            particle_tile.push_back(p);

                        }


                    }

                }

            }

        }
    }

    Redistribute();
    SortParticles();

}

void FhdParticleContainer::PrintCellList(int i, int j, int k)
{

    int lev = 0;

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

        if((i >= myLo[0]) && (j >= myLo[1]) && (k >= myLo[2]) && (i <= myHi[0]) && (j <= myHi[1]) && (k <= myHi[2]))
        {
            IntVect iv = {i,j,k};
            long imap = tile_box.index(iv);
            
            for(int ii = 0; ii<nspecies; ii++)
            {
                int size = m_cell_vectors[ii][pti.index()][imap].size();
                cout << "Species " << ii << ":\n";
			    for(int jj = 0; jj<size; jj++)
		    	{                
		            cout << m_cell_vectors[ii][pti.index()][imap][jj] << " ";
		        }
				cout << "\n";
			}
        }
    }

}

void FhdParticleContainer::PrintCellListInternal(int i, int j, int k, int index, long imap)
{
        IntVect iv = {i,j,k};
        for(int ii = 0; ii<nspecies; ii++)
        {
            int size = m_cell_vectors[ii][index][imap].size();
            cout << "Species " << ii << ":\n";
		    for(int jj = 0; jj<size; jj++)
	    	{                
	            cout << m_cell_vectors[ii][index][imap][jj] << " ";
	        }
			cout << "\n";
		}
}

