#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)

#include <common_functions.H>
#include <common_namespace.H>

#include <FhdParticleContainer.H>
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


bool FhdParticleContainer::use_neighbor_list  {true};
bool FhdParticleContainer::sort_neighbor_list {false};


FhdParticleContainer::FhdParticleContainer(const Geometry & geom,
                                         const DistributionMapping & dmap,
                                         const BoxArray & ba,
                                         int n_nbhd)
    : IBMarkerContainerBase<FHD_realData, FHD_intData>(
            geom, dmap, ba, n_nbhd
        )
      , n_list(0)
{
    InitInternals(n_nbhd);
    nghost = n_nbhd;
}


void FhdParticleContainer::DoRFD(const Real dt, const Real* dxFluid, const Real* dxE, const Geometry geomF, const std::array<MultiFab, AMREX_SPACEDIM>& umac, const std::array<MultiFab, AMREX_SPACEDIM>& efield,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                           const MultiFab& cellCenters,
                                           std::array<MultiFab, AMREX_SPACEDIM>& source,
                                           std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                                           const surface* surfaceList, const int surfaceCount, int sw)
{
    UpdateCellVectors();

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();
    const Real* phi = Geom(lev).ProbHi();

#ifdef _OPENMP
#pragma omp parallel
#endif

    source[0].setVal(0.0);
    source[1].setVal(0.0);
#if (AMREX_SPACEDIM == 3)
    source[2].setVal(0.0);
#endif

    sourceTemp[0].setVal(0.0);
    sourceTemp[1].setVal(0.0);
#if (AMREX_SPACEDIM == 3)
    sourceTemp[2].setVal(0.0);
#endif


    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();
        
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const int np = particles.numParticles();
        
        //Print() << "FHD\n"; 
        do_rfd(particles.data(), &np,
                         ARLIM_3D(tile_box.loVect()),
                         ARLIM_3D(tile_box.hiVect()),
                         m_vector_ptrs[grid_id].dataPtr(),
                         m_vector_size[grid_id].dataPtr(),
                         ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                         ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),
                         ZFILL(plo), ZFILL(phi), ZFILL(dx), &dt, ZFILL(geomF.ProbLo()), ZFILL(dxFluid), ZFILL(dxE),
                         BL_TO_FORTRAN_3D(umac[0][pti]),
                         BL_TO_FORTRAN_3D(umac[1][pti]),
#if (AMREX_SPACEDIM == 3)
                         BL_TO_FORTRAN_3D(umac[2][pti]),
#endif
                         BL_TO_FORTRAN_3D(efield[0][pti]),
                         BL_TO_FORTRAN_3D(efield[1][pti]),
#if (AMREX_SPACEDIM == 3)
                         BL_TO_FORTRAN_3D(efield[2][pti]),
#endif
                         BL_TO_FORTRAN_3D(RealFaceCoords[0][pti]),
                         BL_TO_FORTRAN_3D(RealFaceCoords[1][pti]),
#if (AMREX_SPACEDIM == 3)
                         BL_TO_FORTRAN_3D(RealFaceCoords[2][pti]),
#endif
                         BL_TO_FORTRAN_3D(cellCenters[pti]),
                         BL_TO_FORTRAN_3D(sourceTemp[0][pti]),
                         BL_TO_FORTRAN_3D(sourceTemp[1][pti])
#if (AMREX_SPACEDIM == 3)
                         , BL_TO_FORTRAN_3D(sourceTemp[2][pti])
#endif
                         , surfaceList, &surfaceCount, &sw
                         );


        // resize particle vectors after call to move_particles
//        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
//        {
//            const auto new_size = m_vector_size[grid_id](iv);
//            auto& pvec = m_cell_vectors[grid_id](iv);
//            pvec.resize(new_size);
//        }
    }
}

void FhdParticleContainer::computeForcesNL(const MultiFab& charge, const MultiFab& coords, const Real* dx) {

    BL_PROFILE("FhdParticleContainer::computeForcesNL");

    double rcount = 0;
    const int lev = 0;

    buildNeighborList(CheckPair);


#ifdef _OPENMP
#pragma omp parallel
#endif

   for (FhdParIter pti(*this, lev, MFItInfo().SetDynamic(false)); pti.isValid(); ++pti) {
      
        PairIndex index(pti.index(), pti.LocalTileIndex());
        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();
        int Nn = neighbors[lev][index].size();
        int size = neighbor_list[lev][index].size();

        const Box& tile_box  = pti.tilebox();


        if(sr_tog==1) 
        {
                amrex_compute_forces_nl(particles.data(), &Np, 
                                        neighbors[lev][index].dataPtr(), &Nn,
                                        neighbor_list[lev][index].dataPtr(), &size, &rcount);
        }
        if(es_tog==3)
        {
                amrex_compute_p3m_sr_correction_nl(particles.data(), &Np, 
                                        neighbors[lev][index].dataPtr(), &Nn,
                                        neighbor_list[lev][index].dataPtr(), &size, &rcount,
                                        BL_TO_FORTRAN_3D(charge[pti]),BL_TO_FORTRAN_3D(coords[pti]), ARLIM_3D(tile_box.loVect()), ARLIM_3D(tile_box.hiVect()), ZFILL(dx)); 
        }
    }

    if(sr_tog==1) 
    {
            ParallelDescriptor::ReduceRealSum(rcount);

            Print() << rcount/2 << " close range interactions.\n";
    }
}

void FhdParticleContainer::MoveParticlesDSMC(const Real dt, const surface* surfaceList, const int surfaceCount, Real time, int* flux)
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

void FhdParticleContainer::MoveIons(const Real dt, const Real* dxFluid, const Real* dxE, const Geometry geomF, const std::array<MultiFab, AMREX_SPACEDIM>& umac, const std::array<MultiFab, AMREX_SPACEDIM>& efield,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                           std::array<MultiFab, AMREX_SPACEDIM>& source,
                                           std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                                           const MultiFab& mobility,
                                           const surface* surfaceList, const int surfaceCount, int sw)
{
    
    UpdateCellVectors();

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();
    const Real* phi = Geom(lev).ProbHi();

    double kinetic = 0;

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

        move_ions_fhd(particles.data(), &np,
                         ARLIM_3D(tile_box.loVect()),
                         ARLIM_3D(tile_box.hiVect()),
                         m_vector_ptrs[grid_id].dataPtr(),
                         m_vector_size[grid_id].dataPtr(),
                         ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                         ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),
                         ZFILL(plo), ZFILL(phi), ZFILL(dx), &dt, ZFILL(geomF.ProbLo()), ZFILL(dxFluid), ZFILL(dxE),
                         BL_TO_FORTRAN_3D(umac[0][pti]),
                         BL_TO_FORTRAN_3D(umac[1][pti]),
#if (AMREX_SPACEDIM == 3)
                         BL_TO_FORTRAN_3D(umac[2][pti]),
#endif
                         BL_TO_FORTRAN_3D(efield[0][pti]),
                         BL_TO_FORTRAN_3D(efield[1][pti]),
#if (AMREX_SPACEDIM == 3)
                         BL_TO_FORTRAN_3D(efield[2][pti]),
#endif
                         BL_TO_FORTRAN_3D(RealFaceCoords[0][pti]),
                         BL_TO_FORTRAN_3D(RealFaceCoords[1][pti]),
#if (AMREX_SPACEDIM == 3)
                         BL_TO_FORTRAN_3D(RealFaceCoords[2][pti]),
#endif
                         BL_TO_FORTRAN_3D(sourceTemp[0][pti]),
                         BL_TO_FORTRAN_3D(sourceTemp[1][pti])
#if (AMREX_SPACEDIM == 3)
                         , BL_TO_FORTRAN_3D(sourceTemp[2][pti])
#endif
                         , BL_TO_FORTRAN_3D(mobility[pti])
                         , surfaceList, &surfaceCount, &kinetic, &sw
                         );


        // resize particle vectors after call to move_particles
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
        {
            const auto new_size = m_vector_size[grid_id](iv);
            auto& pvec = m_cell_vectors[grid_id](iv);
            pvec.resize(new_size);
        }
    }
}

void FhdParticleContainer::SpreadIons(const Real dt, const Real* dxFluid, const Real* dxE, const Geometry geomF, const std::array<MultiFab, AMREX_SPACEDIM>& umac, const std::array<MultiFab, AMREX_SPACEDIM>& efield,
                                           const MultiFab& charge,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                           const MultiFab& cellCenters,
                                           std::array<MultiFab, AMREX_SPACEDIM>& source,
                                           std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                                           const surface* surfaceList, const int surfaceCount, int sw)
{
    


    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();
    const Real* phi = Geom(lev).ProbHi();

    double potential = 0;

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
        
        //Print() << "FHD\n"; 
        spread_ions_fhd(particles.data(), &np,
                         ARLIM_3D(tile_box.loVect()),
                         ARLIM_3D(tile_box.hiVect()),
                         m_vector_ptrs[grid_id].dataPtr(),
                         m_vector_size[grid_id].dataPtr(),
                         ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                         ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),
                         ZFILL(plo), ZFILL(phi), ZFILL(dx), &dt, ZFILL(geomF.ProbLo()), ZFILL(dxFluid), ZFILL(dxE),
                         BL_TO_FORTRAN_3D(umac[0][pti]),
                         BL_TO_FORTRAN_3D(umac[1][pti]),
#if (AMREX_SPACEDIM == 3)
                         BL_TO_FORTRAN_3D(umac[2][pti]),
#endif
                         BL_TO_FORTRAN_3D(efield[0][pti]),
                         BL_TO_FORTRAN_3D(efield[1][pti]),
#if (AMREX_SPACEDIM == 3)
                         BL_TO_FORTRAN_3D(efield[2][pti]),
#endif
                         BL_TO_FORTRAN_3D(charge[pti]),
                         BL_TO_FORTRAN_3D(RealFaceCoords[0][pti]),
                         BL_TO_FORTRAN_3D(RealFaceCoords[1][pti]),
#if (AMREX_SPACEDIM == 3)
                         BL_TO_FORTRAN_3D(RealFaceCoords[2][pti]),
#endif
                         BL_TO_FORTRAN_3D(cellCenters[pti]),
                         BL_TO_FORTRAN_3D(sourceTemp[0][pti]),
                         BL_TO_FORTRAN_3D(sourceTemp[1][pti])
#if (AMREX_SPACEDIM == 3)
                         , BL_TO_FORTRAN_3D(sourceTemp[2][pti])
#endif
                         , surfaceList, &surfaceCount, &potential, &sw
                         );
    }
    

    sourceTemp[0].SumBoundary(geomF.periodicity());
    sourceTemp[1].SumBoundary(geomF.periodicity());
#if (AMREX_SPACEDIM == 3)
    sourceTemp[2].SumBoundary(geomF.periodicity());
#endif

    MultiFab::Add(source[0],sourceTemp[0],0,0,source[0].nComp(),source[0].nGrow());
    MultiFab::Add(source[1],sourceTemp[1],0,0,source[1].nComp(),source[1].nGrow());
#if (AMREX_SPACEDIM == 3)
    MultiFab::Add(source[2],sourceTemp[2],0,0,source[2].nComp(),source[2].nGrow());
#endif

    source[0].FillBoundary(geomF.periodicity());
    source[1].FillBoundary(geomF.periodicity());
#if (AMREX_SPACEDIM == 3)
    source[2].FillBoundary(geomF.periodicity());
#endif

}

void FhdParticleContainer::SyncMembrane(double* spec3xPos, double* spec3yPos, double* spec3zPos, double* spec3xForce, double* spec3yForce, double* spec3zForce, const int length, const int step, const species* particleInfo)
{
    

    const int lev = 0;
    double temp;


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

        sync_particles(spec3xPos, spec3yPos, spec3zPos, spec3xForce, spec3yForce, spec3zForce, particles.data(), &np, &length);                    

    }

    //I'm sure there is an array version of this but this will do for now.
    for(int i=0;i<length;i++)
    {
        temp = spec3xPos[i];
        ParallelDescriptor::ReduceRealSum(temp);
        spec3xPos[i] = temp;

        temp = spec3yPos[i];
        ParallelDescriptor::ReduceRealSum(temp);
        spec3yPos[i] = temp;

        temp = spec3zPos[i];
        ParallelDescriptor::ReduceRealSum(temp);
        spec3zPos[i] = temp;

        spec3xForce[i] = 0;
        spec3yForce[i] = 0;
        spec3zForce[i] = 0;
    }

    if(ParallelDescriptor::MyProc() == 0)
    {

        user_force_calc(spec3xPos, spec3yPos, spec3zPos, spec3xForce, spec3yForce, spec3zForce, &length, &step, particleInfo);

    }

    for(int i=0;i<length;i++)
    {
        temp = spec3xForce[i];
        ParallelDescriptor::ReduceRealSum(temp);
        spec3xForce[i] = temp;

        temp = spec3yForce[i];
        ParallelDescriptor::ReduceRealSum(temp);
        spec3yForce[i] = temp;

        temp = spec3zForce[i];
        ParallelDescriptor::ReduceRealSum(temp);
        spec3zForce[i] = temp;

    }

    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();
        
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const int np = particles.numParticles();

        force_particles(spec3xPos, spec3yPos, spec3zPos, spec3xForce, spec3yForce, spec3zForce, particles.data(), &np, &length);                    

    }
}

void FhdParticleContainer::collectFields(const Real dt, const Real* dxPotential, 
                                         const MultiFab& RealCenterCoords, const Geometry geomP, MultiFab& charge, MultiFab& chargeTemp,
                                         MultiFab& mass, MultiFab& massTemp)
{
    BL_PROFILE_VAR("collectFields()",collectFields);
    
    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();
    const Real* phi = Geom(lev).ProbHi();

#ifdef _OPENMP
#pragma omp parallel
#endif

    charge.setVal(0.0);
    chargeTemp.setVal(0.0);

    mass.setVal(0.0);
    massTemp.setVal(0.0);

    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();
        
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const int np = particles.numParticles();
        
        collect_charge(particles.data(), &np,
                         ARLIM_3D(tile_box.loVect()),
                         ARLIM_3D(tile_box.hiVect()),
                         m_vector_ptrs[grid_id].dataPtr(),
                         m_vector_size[grid_id].dataPtr(),
                         ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                         ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),
                         ZFILL(plo), ZFILL(phi), ZFILL(dx), &dt, ZFILL(geomP.ProbLo()), ZFILL(dxPotential),
                         BL_TO_FORTRAN_3D(RealCenterCoords[pti]),
                         BL_TO_FORTRAN_3D(chargeTemp[pti]));

    }

    chargeTemp.SumBoundary(geomP.periodicity());
    //massTemp.SumBoundary(geomP.periodicity());

    MultiFab::Add(charge,chargeTemp,0,0,charge.nComp(),charge.nGrow());
    //MultiFab::Add(mass,massTemp,0,0,mass.nComp(),mass.nGrow());

    charge.FillBoundary(geomP.periodicity());
}



void FhdParticleContainer::InitCollisionCells(
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

void FhdParticleContainer::CollideParticles(MultiFab& collisionPairs,
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

void FhdParticleContainer::InitializeFields(MultiFab& particleInstant,
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

void FhdParticleContainer::EvaluateStats(
                              MultiFab& particleInstant,
                              MultiFab& particleMeans,
                              MultiFab& particleVars,
                              MultiFab& cellVols, species particleInfo, const Real delt, int steps)
{
    const int lev = 0;
    const double Neff = particleInfo.Neff;
    const double n0 = particleInfo.n0;
    const double T0 = particleInfo.T;


    double tp = 0;
    double te = 0;
    double tm = 0;

    double totalMass;    

    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) 
    {
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
                         BL_TO_FORTRAN_3D(cellVols[pti]),&Neff, &Np
                        );

    }


    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) 
    {
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

                         BL_TO_FORTRAN_3D(cellVols[pti]), &Np,&Neff,&n0,&T0,&delt, &steps);
    }

    //Print() << "c++: " << delHolder6[0] << "\n";


    //Print() << "Total mass: " << totalMass << "\n";

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

                         BL_TO_FORTRAN_3D(cellVols[pti]), &Np,&Neff,&n0,&T0,&delt, &steps
                        );
    }

}

void FhdParticleContainer::WriteParticlesAscii(std::string asciiName)
{
    WriteAsciiFile(asciiName);
}

void
FhdParticleContainer::UpdateCellVectors()
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
            p.idata(FHD_intData::sorted) = 1;
            p.idata(FHD_intData::i) = iv[0];
            p.idata(FHD_intData::j) = iv[1];
            p.idata(FHD_intData::k) = iv[2];
            // note - use 1-based indexing for convenience with Fortran
            m_cell_vectors[pti.index()](iv).push_back(pindex + 1);
        }
    }
    
    UpdateFortranStructures();
}


void
FhdParticleContainer::UpdateFortranStructures()
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
FhdParticleContainer::ReBin()
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
	    if (p.idata(FHD_intData::sorted)) continue;
            const IntVect& iv = this->Index(p, lev);
            p.idata(FHD_intData::sorted) = 1;
            p.idata(FHD_intData::i) = iv[0];
            p.idata(FHD_intData::j) = iv[1];
            p.idata(FHD_intData::k) = iv[2];
            // note - use 1-based indexing for convenience with Fortran
            m_cell_vectors[pti.index()](iv).push_back(pindex + 1);
        }
    }

    UpdateFortranStructures();
}

void
FhdParticleContainer::correctCellVectors(int old_index, int new_index, 
						int grid, const ParticleType& p)
{
    if (not p.idata(FHD_intData::sorted)) return;
    IntVect iv(p.idata(FHD_intData::i), p.idata(FHD_intData::j), p.idata(FHD_intData::k));
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
FhdParticleContainer::numWrongCell()
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
            if ((iv[0] != p.idata(FHD_intData::i)) or (iv[1] != p.idata(FHD_intData::j)) or (iv[2] != p.idata(FHD_intData::k))) {
                num_wrong += 1;
            }
        }
    }
    
    ParallelDescriptor::ReduceIntSum(num_wrong, ParallelDescriptor::IOProcessorNumber());
    return num_wrong;
}



