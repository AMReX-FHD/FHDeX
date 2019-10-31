#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)

#include <common_functions.H>
#include <common_namespace.H>

#include <DsmcParticleContainer.H>
#include <ib_functions_F.H>

#include <iostream>

#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>
#include <iostream>
#include <fstream>


using namespace common;
using namespace amrex;


DsmcParticleContainer::DsmcParticleContainer(const Geometry & geom,
                              const DistributionMapping & dmap,
                              const BoxArray            & ba,
                              int ncells)
    : NeighborParticleContainer<DSMC_realData::count, DSMC_intData::count> (geom, dmap, ba, ncells)
{}

void DsmcParticleContainer::MoveParticlesDSMC(const Real dt, const surface* surfaceList, const int surfaceCount, Real time, int* flux)
{

  // Print() << "HERE MoveParticlesDSMC" << std::endline
  
    UpdateCellVectors();
    int i; int lFlux = 0; int rFlux = 0;
    const int  lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();
    const Real* phi = Geom(lev).ProbHi();

BL_PROFILE_VAR_NS("particle_move", particle_move);

BL_PROFILE_VAR_START(particle_move);

#ifdef _OPENMP
#pragma omp parallel
#endif
 
    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();
        
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const int np = particles.numParticles();

	// Print() << "FHD\n";

        move_particles_dsmc(particles.data(), &np,
                       ARLIM_3D(tile_box.loVect()), 
                       ARLIM_3D(tile_box.hiVect()),
                       m_vector_ptrs[grid_id].dataPtr(),
                       m_vector_size[grid_id].dataPtr(),
                       ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                       ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),
                       ZFILL(plo),ZFILL(phi),ZFILL(dx), &dt,
                       surfaceList, &surfaceCount, &time, flux);

        lFlux += flux[0]; rFlux += flux[1];

        // resize particle vectors after call to move_particles
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
        {
            const auto new_size = m_vector_size[grid_id](iv);
            auto& pvec = m_cell_vectors[grid_id](iv);
            pvec.resize(new_size);
        }
    }

    //reduce flux
    ParallelDescriptor::ReduceIntSum(lFlux);
    ParallelDescriptor::ReduceIntSum(rFlux);

    //reset flux
    flux[0] = lFlux; flux[1] = rFlux;

    std::ofstream outfile;
    
    if(graphene_tog==1)
      {
		char num[21];
		std::string txt=".txt";
		std::string text="test";
		sprintf(num, "%f", domega);
		outfile.open(text+num + txt, std::ios_base::app);
        //      outfile.open("out.csv", std::ios_base::app);
  for (i=0;i<1;i++)
	  {
	    outfile << surfaceList[5].dbesslist[i] << ", ";
	  }
	outfile<<"\n";
	  outfile.close();
      }
	    

BL_PROFILE_VAR_STOP(particle_move);

}


void DsmcParticleContainer::InitCollisionCells(
                              MultiFab& collisionPairs,
                              MultiFab& collisionFactor, 
                              MultiFab& cellVols, const species particleInfo, const Real delt)
{
    BL_PROFILE_VAR("InitCollisionCells()",InitCollisionCells);

    UpdateCellVectors();
    const int lev = 0;

    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) 
    {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();
        
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const int Np = particles.numParticles();

        init_cells(particles.data(),
                         ARLIM_3D(tile_box.loVect()), 
                         ARLIM_3D(tile_box.hiVect()),
                         m_vector_ptrs[grid_id].dataPtr(),
                         m_vector_size[grid_id].dataPtr(),
                         ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                         ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),
                         BL_TO_FORTRAN_3D(collisionPairs[pti]),
                         BL_TO_FORTRAN_3D(collisionFactor[pti]),
                         BL_TO_FORTRAN_3D(cellVols[pti]),&Np,&particleInfo.Neff,&particleInfo.cp,&particleInfo.d,&delt
                        );
    }
}

void DsmcParticleContainer::CollideParticles(MultiFab& collisionPairs,
                                            MultiFab& collisionFactor, 
                                            MultiFab& cellVols, const species particleInfo, const Real delt)
{
    BL_PROFILE_VAR("CollideParticles()",CollideParticles);
    
    const int lev = 0;

    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) 
    {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();
        
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const int Np = particles.numParticles();

        collide_cells(particles.data(),
                         ARLIM_3D(tile_box.loVect()), 
                         ARLIM_3D(tile_box.hiVect()),
                         m_vector_ptrs[grid_id].dataPtr(),
                         m_vector_size[grid_id].dataPtr(),
                         ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                         ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),
                         BL_TO_FORTRAN_3D(collisionPairs[pti]),
                         BL_TO_FORTRAN_3D(collisionFactor[pti]),
                         BL_TO_FORTRAN_3D(cellVols[pti]),&Np,&particleInfo.Neff,&particleInfo.cp,&particleInfo.d,&delt
                        );
    }
}

void DsmcParticleContainer::InitializeFields(MultiFab& particleInstant,
                                            MultiFab& cellVols, const species particleInfo)
{
    BL_PROFILE_VAR("InitializeFields()",InitializeFields);  

    UpdateCellVectors();
    const int lev = 0;

    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) 
    {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& parts = particle_tile.GetArrayOfStructs();
        const int Np = parts.numParticles();

        initialize_fields(parts.data(),
                         ARLIM_3D(tile_box.loVect()),
                         ARLIM_3D(tile_box.hiVect()),
                         m_vector_ptrs[grid_id].dataPtr(),
                         m_vector_size[grid_id].dataPtr(),
                         ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                         ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),   
                         BL_TO_FORTRAN_3D(particleInstant[pti]),
                         BL_TO_FORTRAN_3D(cellVols[pti]),&particleInfo.Neff, &Np, &particleInfo.R, &particleInfo.T
                        );

    }
}

void DsmcParticleContainer::EvaluateStats(MultiFab& particleInstant,
                                          MultiFab& particleMeans,
                                          MultiFab& particleVars,
                                          MultiFab& cellVols, species particleInfo,
                                          const Real delt, int steps)
{
    const int lev = 0;
    const double Neff = particleInfo.Neff;
    const double n0 = particleInfo.n0;
    const double T0 = particleInfo.T;

    double tp = 0;
    double te = 0;
    double tm = 0;

    double totalMass;

    int cellcount_tile = 0, cellcount_proc = 0;
    RealVector avcurrent_tile(3), avcurrent_proc(3);
    RealVector varcurrent_tile(3), varcurrent_proc(3);

    std::fill(avcurrent_tile.begin(), avcurrent_tile.end(), 0.);
    std::fill(avcurrent_proc.begin(), avcurrent_proc.end(), 0.);
    std::fill(varcurrent_tile.begin(), varcurrent_tile.end(), 0.);
    std::fill(varcurrent_proc.begin(), varcurrent_proc.end(), 0.);

    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& parts = particle_tile.GetArrayOfStructs();
        const int Np = parts.numParticles();

        tp = tp + Np;

        evaluate_fields(parts.data(),
                         ARLIM_3D(tile_box.loVect()),
                         ARLIM_3D(tile_box.hiVect()),
                         m_vector_ptrs[grid_id].dataPtr(),
                         m_vector_size[grid_id].dataPtr(),
                         ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                         ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),   
                         BL_TO_FORTRAN_3D(particleInstant[pti]),
                         BL_TO_FORTRAN_3D(cellVols[pti]),&Neff, &Np);
    }

    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) {
        
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& parts = particle_tile.GetArrayOfStructs();
        const int Np = parts.numParticles();

        evaluate_means(parts.data(),
                       ARLIM_3D(tile_box.loVect()),
                       ARLIM_3D(tile_box.hiVect()),
                       m_vector_ptrs[grid_id].dataPtr(),
                       m_vector_size[grid_id].dataPtr(),
                       ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                       ARLIM_3D(m_vector_ptrs[grid_id].hiVect()), 
                       BL_TO_FORTRAN_3D(particleInstant[pti]),
                       BL_TO_FORTRAN_3D(particleMeans[pti]),
                       BL_TO_FORTRAN_3D(particleVars[pti]),
                       BL_TO_FORTRAN_3D(cellVols[pti]), &Np,&Neff,&n0,&T0,&delt, &steps,
                       &cellcount_tile,avcurrent_tile.dataPtr());

        // gather statistics
        cellcount_proc += cellcount_tile;
        for (int i=0; i<3; ++i) {
            avcurrent_proc[i] += avcurrent_tile[i];
        }
    }

    // gather statistics
    ParallelDescriptor::ReduceIntSum(cellcount_proc);
    ParallelDescriptor::ReduceRealSum(avcurrent_proc.dataPtr(),3);

    // print statistics
    Print() << "Current density mean: "
            << avcurrent_proc[0]/cellcount_proc << "  "
            << avcurrent_proc[1]/cellcount_proc << "  "
            << avcurrent_proc[2]/cellcount_proc << "\n";

    // reset cell count
    cellcount_tile = 0;
    cellcount_proc = 0;
    
    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) 
    {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& parts = particle_tile.GetArrayOfStructs();
        const int Np = parts.numParticles();

        evaluate_corrs(parts.data(),
                       ARLIM_3D(tile_box.loVect()),
                       ARLIM_3D(tile_box.hiVect()),
                       m_vector_ptrs[grid_id].dataPtr(),
                       m_vector_size[grid_id].dataPtr(),
                       ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                       ARLIM_3D(m_vector_ptrs[grid_id].hiVect()), 
                       BL_TO_FORTRAN_3D(particleInstant[pti]),
                       BL_TO_FORTRAN_3D(particleMeans[pti]),
                       BL_TO_FORTRAN_3D(particleVars[pti]),
                       BL_TO_FORTRAN_3D(cellVols[pti]), &Np,&Neff,&n0,&T0,&delt, &steps,
                       &cellcount_tile, varcurrent_tile.dataPtr()
            );

        // gather statistics
        cellcount_proc += cellcount_tile;
        for (int i=0; i<3; ++i) {
            varcurrent_proc[i] += varcurrent_tile[i];
        }
    }

    // gather statistics
    ParallelDescriptor::ReduceIntSum(cellcount_proc);
    ParallelDescriptor::ReduceRealSum(varcurrent_proc.dataPtr(),3);

    // print statistics
    Print() << "Current density variance: "
            << varcurrent_proc[0]/cellcount_proc << "  "
            << varcurrent_proc[1]/cellcount_proc << "  "
            << varcurrent_proc[2]/cellcount_proc << "\n";
}

void DsmcParticleContainer::WriteParticlesAscii(std::string asciiName)
{
    WriteAsciiFile(asciiName);
}

void
DsmcParticleContainer::UpdateCellVectors()
{
    BL_PROFILE("CellSortedParticleContainer::UpdateCellVectors");
    
    const int lev = 0;

    bool needs_update = true;
    if (not m_vectors_initialized)
    {
        // this is the first call, so we must update
        m_vectors_initialized = true;
        needs_update = true;
    }
    else if ((m_BARef != this->ParticleBoxArray(lev).getRefID()) or 
             (m_DMRef != this->ParticleDistributionMap(lev).getRefID()))
    {
        // the grids have changed, so we must update
        m_BARef = this->ParticleBoxArray(lev).getRefID();
        m_DMRef = this->ParticleDistributionMap(lev).getRefID();
        needs_update = true;
    }
    
    if (not needs_update) return;

    // clear old data
    m_cell_vectors.clear();
    m_vector_size.clear();
    m_vector_ptrs.clear();
    
    // allocate storage for cell vectors. NOTE - do not tile this loop
    for(MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const int grid_id = mfi.index();
        m_cell_vectors[grid_id].resize(box);
        m_vector_size[grid_id].resize(box);
        m_vector_ptrs[grid_id].resize(box);
    }

    // insert particles into vectors - this can be tiled
#ifdef _OPENMP
#pragma omp parallel
#endif    
    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& particles = pti.GetArrayOfStructs();
        const int np    = pti.numParticles();
        for(int pindex = 0; pindex < np; ++pindex) {
            ParticleType& p = particles[pindex];
            const IntVect& iv = this->Index(p, lev);
            p.idata(DSMC_intData::sorted) = 1;
            p.idata(DSMC_intData::i) = iv[0];
            p.idata(DSMC_intData::j) = iv[1];
            p.idata(DSMC_intData::k) = iv[2];
            // note - use 1-based indexing for convenience with Fortran
            m_cell_vectors[pti.index()](iv).push_back(pindex + 1);
        }
    }
    
    UpdateFortranStructures();
}


void
DsmcParticleContainer::UpdateFortranStructures()
{
    BL_PROFILE("CellSortedParticleContainer::UpdateFortranStructures");
    
    const int lev = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const int grid_id = mfi.index();
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
        {
            m_vector_size[grid_id](iv) = m_cell_vectors[grid_id](iv).size();
            m_vector_ptrs[grid_id](iv) = m_cell_vectors[grid_id](iv).data();
        }
    }
}

void
DsmcParticleContainer::ReBin()
{
    BL_PROFILE("CellSortedParticleContainer::ReBin()");
    
    const int lev = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const int np = particles.numParticles();
        for(int pindex = 0; pindex < np; ++pindex)
        {
	    ParticleType& p = particles[pindex];
	    if (p.idata(DSMC_intData::sorted)) continue;
            const IntVect& iv = this->Index(p, lev);
            p.idata(DSMC_intData::sorted) = 1;
            p.idata(DSMC_intData::i) = iv[0];
            p.idata(DSMC_intData::j) = iv[1];
            p.idata(DSMC_intData::k) = iv[2];
            // note - use 1-based indexing for convenience with Fortran
            m_cell_vectors[pti.index()](iv).push_back(pindex + 1);
        }
    }

    UpdateFortranStructures();
}

void
DsmcParticleContainer::correctCellVectors(int old_index, int new_index, 
						int grid, const ParticleType& p)
{
    if (not p.idata(DSMC_intData::sorted)) return;
    IntVect iv(p.idata(DSMC_intData::i), p.idata(DSMC_intData::j), p.idata(DSMC_intData::k));
    //IntVect iv(AMREX_D_DECL(p.idata(IntData::i), p.idata(IntData::j), p.idata(IntData::k)));
    auto& cell_vector = m_cell_vectors[grid](iv);
    for (int i = 0; i < static_cast<int>(cell_vector.size()); ++i) {
        if (cell_vector[i] == old_index + 1) {
            cell_vector[i] = new_index + 1;
            return;
        }
    }
}

int
DsmcParticleContainer::numWrongCell()
{
    const int lev = 0;
    int num_wrong = 0;
    
#ifdef _OPENMP
#pragma omp parallel reduction(+:num_wrong)
#endif    
    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& particles = pti.GetArrayOfStructs();
        const int np    = pti.numParticles();
        for(int pindex = 0; pindex < np; ++pindex) {
            const ParticleType& p = particles[pindex];
            const IntVect& iv = this->Index(p, lev);
            if ((iv[0] != p.idata(DSMC_intData::i)) or (iv[1] != p.idata(DSMC_intData::j)) or (iv[2] != p.idata(DSMC_intData::k))) {
                num_wrong += 1;
            }
        }
    }
    
    ParallelDescriptor::ReduceIntSum(num_wrong, ParallelDescriptor::IOProcessorNumber());
    return num_wrong;
}



