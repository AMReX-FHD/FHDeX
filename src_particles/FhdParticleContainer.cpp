#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>

#include "FhdParticleContainer.H"
#include "particle_functions_F.H"
#include "rng_functions_F.H"
//#include "common_namespace.H"

using namespace amrex;
using namespace common;

//constexpr Real FhdParticleContainer::min_r;
//constexpr Real FhdParticleContainer::cutoff;


FhdParticleContainer::FhdParticleContainer(const Geometry & geom,
                              const DistributionMapping & dmap,
                              const BoxArray            & ba,
                              int ncells)
    : NeighborParticleContainer<RealData::ncomps, IntData::ncomps> (geom, dmap, ba, ncells)
{}

void FhdParticleContainer::InitParticles(species* particleInfo)
{
    
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();
    const Real* plo = geom.ProbLo();
    const Real* phi = geom.ProbHi();

    int qcount = 0;

    double cosTheta, sinTheta, cosPhi, sinPhi;    

    int pcount = 0;
        
    int ll =0;
    //double initTemp = 0;
    //double pc = 0;

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

        //Assuming tile=box for now, i.e. no tiling.
        IntVect smallEnd = tile_box.smallEnd();
        IntVect bigEnd = tile_box.bigEnd();       


//        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
 //       {

            for(int i_spec=0; i_spec < nspecies; i_spec++)
            {
            for (int i_part=0; i_part<particleInfo[i_spec].ppb;i_part++)
            {
                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                p.idata(IntData::sorted) = 0;
                
//                p.pos(0) = smallEnd[0]*dx[0] + get_uniform_func()*dx[0]*(bigEnd[0]-smallEnd[0]+1);
//                p.pos(1) = smallEnd[1]*dx[1] + get_uniform_func()*dx[1]*(bigEnd[1]-smallEnd[1]+1);
//#if (BL_SPACEDIM == 3)
//                p.pos(2) = smallEnd[2]*dx[2] + get_uniform_func()*dx[2]*(bigEnd[2]-smallEnd[2]+1);
//#endif

                p.pos(0) = 10*dx[0] + ll*particleInfo[i_spec].sigma;
                p.pos(1) = 10*dx[1];
#if (BL_SPACEDIM == 3)
                p.pos(2) = 10*dx[2];
#endif
                ll++;
                
                p.rdata(RealData::q) = particleInfo[i_spec].q;

                Print() << "Pos: " << p.pos(0) << ", " << p.pos(1) << ", " << p.pos(2) << ", " << p.rdata(RealData::q) << "\n" ;

                //original position stored for MSD calculations
                p.rdata(RealData::ox) = p.pos(0);
                p.rdata(RealData::oy) = p.pos(1);
#if (BL_SPACEDIM == 3)
                p.rdata(RealData::oz) = p.pos(2);
#endif

                //p.rdata(RealData::vx) = sqrt(particleInfo.R*particleInfo.T)*get_particle_normal_func();
                //p.rdata(RealData::vy) = sqrt(particleInfo.R*particleInfo.T)*get_particle_normal_func();
                //p.rdata(RealData::vz) = sqrt(particleInfo.R*particleInfo.T)*get_particle_normal_func();

                p.rdata(RealData::vx) = 0;
                p.rdata(RealData::vy) = 0;
                p.rdata(RealData::vz) = 0;

                p.rdata(RealData::ax) = 0;
                p.rdata(RealData::ay) = 0;
                p.rdata(RealData::az) = 0;

                p.rdata(RealData::fx) = 0;
                p.rdata(RealData::fy) = 0;
                p.rdata(RealData::fz) = 0;

                p.rdata(RealData::travelTime) = 0;
                p.rdata(RealData::diffAv) = 0;
                p.rdata(RealData::stepCount) = 0;

                p.rdata(RealData::mass) = particleInfo[i_spec].m; //mass
                p.rdata(RealData::R) = particleInfo[i_spec].R; //R
                p.rdata(RealData::radius) = particleInfo[i_spec].d/2.0; //radius
                p.rdata(RealData::accelFactor) = -6*3.14159265359*p.rdata(RealData::radius)/p.rdata(RealData::mass); //acceleration factor (replace with amrex c++ constant for pi...)
                p.rdata(RealData::dragFactor) = 6*3.14159265359*p.rdata(RealData::radius); //drag factor
                //p.rdata(RealData::dragFactor) = 6*3.14159265359*dx[0]*1.322; //drag factor

                p.rdata(RealData::wetDiff) = particleInfo[i_spec].wetDiff;
                p.rdata(RealData::dryDiff) = particleInfo[i_spec].dryDiff;
                p.rdata(RealData::totalDiff) = particleInfo[i_spec].totalDiff;

                p.rdata(RealData::sigma) = particleInfo[i_spec].sigma;
                p.rdata(RealData::eepsilon) = particleInfo[i_spec].eepsilon;
                
                particle_tile.push_back(p);

                pcount++;
            }
 //           }
        }
    }

//    std::cout << "pcount: " << pcount << "\n";

    UpdateCellVectors();
    Redistribute();
    ReBin();

        Print() << "Particles1: " << TotalNumberOfParticles() << "\n";
}

void FhdParticleContainer::computeForcesNL() {

    BL_PROFILE("FhdParticleContainer::computeForcesNL");

    const int lev = 0;

    buildNeighborList(CheckPair);


#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MyParIter pti(*this, lev, MFItInfo().SetDynamic(false)); pti.isValid(); ++pti) {
        PairIndex index(pti.index(), pti.LocalTileIndex());
        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();
        int Nn = neighbors[lev][index].size();
        int size = neighbor_list[lev][index].size();

        amrex_compute_forces_nl(particles.data(), &Np, 
                                neighbors[lev][index].dataPtr(), &Nn,
                                neighbor_list[lev][index].dataPtr(), &size); 
    }
}

void FhdParticleContainer::MoveParticlesDry(const Real dt, const Real* dxFluid, const std::array<MultiFab, AMREX_SPACEDIM>& umac,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                           std::array<MultiFab, AMREX_SPACEDIM>& source,
                                           std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                                           const surface* surfaceList, const int surfaceCount)
{
    
    UpdateCellVectors();

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();
    const Real* phi = Geom(lev).ProbHi();

BL_PROFILE_VAR_NS("particle_move", particle_move);

BL_PROFILE_VAR_START(particle_move);

    //Arg1: Source multifab to be shifted. Arg2: destination multiFab. Arg3: A cell centred multifab for reference (change this later).
//    FindNodalValues(umac[0], umacNodal[0], betaCC);
//    FindNodalValues(umac[1], umacNodal[1], betaCC);

//#if (AMREX_SPACEDIM == 3)
  //  FindNodalValues(umac[2], umacNodal[2], betaCC);

    //While beta is constant we will pass a prefilled betaNodal
    //FindNodalValues(betaCC, betaNodal, betaCC);
//#endif

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

        Print() << "parts: " << np << std::endl;
        Print() << "move particles DRY\n"; 
        move_particles_dry(particles.data(), &np,
                         ARLIM_3D(tile_box.loVect()),
                         ARLIM_3D(tile_box.hiVect()),
                         m_vector_ptrs[grid_id].dataPtr(),
                         m_vector_size[grid_id].dataPtr(),
                         ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                         ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),
                         ZFILL(plo), ZFILL(phi), ZFILL(dx), &dt,
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
                         , surfaceList, &surfaceCount
                         );


        // resize particle vectors after call to move_particles
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
        {
            const auto new_size = m_vector_size[grid_id](iv);
            auto& pvec = m_cell_vectors[grid_id](iv);
            pvec.resize(new_size);
        }
    }

BL_PROFILE_VAR_STOP(particle_move);

}

        

#ifndef DSMC
void FhdParticleContainer::MoveParticles(const Real dt, const Real* dxFluid, const Real* ploFluid, const std::array<MultiFab, AMREX_SPACEDIM>& umac,
                                           std::array<MultiFab, AMREX_SPACEDIM>& umacNodal,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                           const MultiFab& betaCC, //Not necessary but may use later
                                           MultiFab& betaNodal, //Not necessary but may use later
                                           const MultiFab& rhoCC, //Not necessary but may use later
                                           const MultiFab& rhoNodal, //Not necessary but may use later
                                           std::array<MultiFab, AMREX_SPACEDIM>& source,
                                           std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                                           const surface* surfaceList, const int surfaceCount)

#endif
#ifdef DSMC
void FhdParticleContainer::MoveParticles(const Real dt, const surface* surfaceList, const int surfaceCount)

#endif
{
    
    UpdateCellVectors();

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();
    const Real* phi = Geom(lev).ProbHi();

BL_PROFILE_VAR_NS("particle_move", particle_move);

BL_PROFILE_VAR_START(particle_move);

#ifndef DSMC

    //Arg1: Source multifab to be shifted. Arg2: destination multiFab. Arg3: A cell centred multifab for reference (change this later).
    FindNodalValues(umac[0], umacNodal[0], betaCC);
    FindNodalValues(umac[1], umacNodal[1], betaCC);

#if (AMREX_SPACEDIM == 3)
    FindNodalValues(umac[2], umacNodal[2], betaCC);

    //While beta is constant we will pass a prefilled betaNodal
    //FindNodalValues(betaCC, betaNodal, betaCC);
#endif

#endif

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

        //Print() << "parts: " << np << std::endl;
        
#ifdef DSMC
        //Print() << "DSMC\n";        
        move_particles_dsmc(particles.data(), &np,
                       ARLIM_3D(tile_box.loVect()), 
                       ARLIM_3D(tile_box.hiVect()),
                       m_vector_ptrs[grid_id].dataPtr(),
                       m_vector_size[grid_id].dataPtr(),
                       ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                       ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),
                       ZFILL(plo),ZFILL(phi),ZFILL(dx), &dt,
                       surfaceList, &surfaceCount);
   
#else
        //Print() << "FHD\n"; 
//        move_particles_fhd(particles.data(), &np,
//                         ARLIM_3D(tile_box.loVect()),
//                         ARLIM_3D(tile_box.hiVect()),
//                         m_vector_ptrs[grid_id].dataPtr(),
//                         m_vector_size[grid_id].dataPtr(),
//                         ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
//                         ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),
//                         ZFILL(plo), ZFILL(phi), ZFILL(dx), &dt, ZFILL(ploFluid), ZFILL(dxFluid),
//                         BL_TO_FORTRAN_3D(umacNodal[0][pti]),
//                         BL_TO_FORTRAN_3D(umacNodal[1][pti]),
//#if (AMREX_SPACEDIM == 3)
//                         BL_TO_FORTRAN_3D(umacNodal[2][pti]),
//#endif
//                         BL_TO_FORTRAN_3D(RealFaceCoords[0][pti]),
//                         BL_TO_FORTRAN_3D(RealFaceCoords[1][pti]),
//#if (AMREX_SPACEDIM == 3)
//                         BL_TO_FORTRAN_3D(RealFaceCoords[2][pti]),
//#endif
//                         BL_TO_FORTRAN_3D(betaNodal[pti]),
//                         BL_TO_FORTRAN_3D(rhoNodal[pti]),

//                         BL_TO_FORTRAN_3D(sourceTemp[0][pti]),
//                         BL_TO_FORTRAN_3D(sourceTemp[1][pti])
//#if (AMREX_SPACEDIM == 3)
//                         , BL_TO_FORTRAN_3D(sourceTemp[2][pti])
//#endif
//                         , surfaceList, &surfaceCount
//                         );
#endif


        // resize particle vectors after call to move_particles
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
        {
            const auto new_size = m_vector_size[grid_id](iv);
            auto& pvec = m_cell_vectors[grid_id](iv);
            pvec.resize(new_size);
        }
    }

#ifndef DSMC
    sourceTemp[0].SumBoundary(Geom(lev).periodicity());
    sourceTemp[1].SumBoundary(Geom(lev).periodicity());
#if (AMREX_SPACEDIM == 3)
    sourceTemp[2].SumBoundary(Geom(lev).periodicity());
#endif
    MultiFab::Add(source[0],sourceTemp[0],0,0,source[0].nComp(),source[0].nGrow());
    MultiFab::Add(source[1],sourceTemp[1],0,0,source[1].nComp(),source[1].nGrow());
#if (AMREX_SPACEDIM == 3)
    MultiFab::Add(source[2],sourceTemp[2],0,0,source[2].nComp(),source[2].nGrow());
#endif
    source[0].FillBoundary(Geom(lev).periodicity());
    source[1].FillBoundary(Geom(lev).periodicity());
#if (AMREX_SPACEDIM == 3)
    source[2].FillBoundary(Geom(lev).periodicity());
#endif
#endif

BL_PROFILE_VAR_STOP(particle_move);

}

void FhdParticleContainer::MoveIons(const Real dt, const Real* dxFluid, const Real* dxE, const Geometry geomF, const std::array<MultiFab, AMREX_SPACEDIM>& umac, const std::array<MultiFab, AMREX_SPACEDIM>& efield,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
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

//    source[0].setVal(0.0);
//    source[1].setVal(0.0);
//#if (AMREX_SPACEDIM == 3)
//    source[2].setVal(0.0);
//#endif

//    sourceTemp[0].setVal(0.0);
//    sourceTemp[1].setVal(0.0);
//#if (AMREX_SPACEDIM == 3)
//    sourceTemp[2].setVal(0.0);
//#endif

    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();
        
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const int np = particles.numParticles();
        
        //Print() << "FHD\n"; 
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
                         , surfaceList, &surfaceCount, &sw
                         );


        // resize particle vectors after call to move_particles
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
        {
            const auto new_size = m_vector_size[grid_id](iv);
            auto& pvec = m_cell_vectors[grid_id](iv);
            pvec.resize(new_size);
        }
    }

//    sourceTemp[0].SumBoundary(geomF.periodicity());
//    sourceTemp[1].SumBoundary(geomF.periodicity());
//#if (AMREX_SPACEDIM == 3)
//    sourceTemp[2].SumBoundary(geomF.periodicity());
//#endif

//    MultiFab::Add(source[0],sourceTemp[0],0,0,source[0].nComp(),source[0].nGrow());
//    MultiFab::Add(source[1],sourceTemp[1],0,0,source[1].nComp(),source[1].nGrow());
//#if (AMREX_SPACEDIM == 3)
//    MultiFab::Add(source[2],sourceTemp[2],0,0,source[2].nComp(),source[2].nGrow());
//#endif

//    source[0].FillBoundary(geomF.periodicity());
//    source[1].FillBoundary(geomF.periodicity());
//#if (AMREX_SPACEDIM == 3)
//    source[2].FillBoundary(geomF.periodicity());
//#endif

}

void FhdParticleContainer::SpreadIons(const Real dt, const Real* dxFluid, const Real* dxE, const Geometry geomF, const std::array<MultiFab, AMREX_SPACEDIM>& umac, const std::array<MultiFab, AMREX_SPACEDIM>& efield,
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
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
        {
            const auto new_size = m_vector_size[grid_id](iv);
            auto& pvec = m_cell_vectors[grid_id](iv);
            pvec.resize(new_size);
        }
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
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
        {
            const auto new_size = m_vector_size[grid_id](iv);
            auto& pvec = m_cell_vectors[grid_id](iv);
            pvec.resize(new_size);
        }
    }

}


void FhdParticleContainer::collectFields(const Real dt, const Real* dxPotential, 
                                         const MultiFab& RealCenterCoords, const Geometry geomP, MultiFab& charge, MultiFab& chargeTemp,
                                         MultiFab& mass, MultiFab& massTemp)
{
    

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

void FhdParticleContainer::CollideParticles(
                              MultiFab& collisionPairs,
                              MultiFab& collisionFactor, 
                              MultiFab& cellVols, const species particleInfo, const Real delt)
{
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

//void FhdParticleContainer::EvaluateStatsMem(
//                              MultiFab& particleInstant,
//                              MultiFab& particleMeans,
//                              MultiFab& particleVars, Real* delHolder1, Real* delHolder2, Real* delHolder3, Real* delHolder4, Real* delHolder5, Real* delHolder6,

//                              MultiFab& particleMembraneFlux,

//                              MultiFab& cellVols, species particleInfo, const Real delt, int steps)
//{
//    const int lev = 0;
//    const double Neff = particleInfo.Neff;
//    const double n0 = particleInfo.n0;
//    const double T0 = particleInfo.T;

//    double del1 = 0;
//    double del2 = 0;
//    double del3 = 0;
//    double del4 = 0;
//    double del5 = 0;
//    double del6 = 0;

//    double tp = 0;
//    double te = 0;
//    double tm = 0;

//    double totalMass;    

//    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) 
//    {
//        const int grid_id = pti.index();
//        const int tile_id = pti.LocalTileIndex();
//        const Box& tile_box  = pti.tilebox();

//        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
//        auto& parts = particle_tile.GetArrayOfStructs();
//        const int Np = parts.numParticles();

//        tp = tp + Np;

//        evaluate_fields(parts.data(),
//                         ARLIM_3D(tile_box.loVect()),
//                         ARLIM_3D(tile_box.hiVect()),
//                         m_vector_ptrs[grid_id].dataPtr(),
//                         m_vector_size[grid_id].dataPtr(),
//                         ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
//                         ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),   
//                         BL_TO_FORTRAN_3D(particleInstant[pti]),
//                         BL_TO_FORTRAN_3D(cellVols[pti]),&Neff, &Np, &del1, &del2, &tm, &te
//                        );

//    }


//    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) 
//    {
//        const int grid_id = pti.index();
//        const int tile_id = pti.LocalTileIndex();
//        const Box& tile_box  = pti.tilebox();

//        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
//        auto& parts = particle_tile.GetArrayOfStructs();
//        const int Np = parts.numParticles();

//        evaluate_means(parts.data(),
//                         ARLIM_3D(tile_box.loVect()),
//                         ARLIM_3D(tile_box.hiVect()),
//                         m_vector_ptrs[grid_id].dataPtr(),
//                         m_vector_size[grid_id].dataPtr(),
//                         ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
//                         ARLIM_3D(m_vector_ptrs[grid_id].hiVect()), 

//                         BL_TO_FORTRAN_3D(particleInstant[pti]),
//                         BL_TO_FORTRAN_3D(particleMeans[pti]),
//                         BL_TO_FORTRAN_3D(particleVars[pti]),

//                         BL_TO_FORTRAN_3D(particleMembraneFlux[pti]),

//                         BL_TO_FORTRAN_3D(cellVols[pti]), &Np,&Neff,&n0,&T0,&delt, &steps, delHolder1, delHolder2, delHolder3, delHolder4, delHolder5, delHolder6, &totalMass
//                        );
//    }

//    //Print() << "c++: " << delHolder6[0] << "\n";

//    //this is a bit of a hack? Reduce real sum should work with vectors
//    for(int i=0;i<(n_cells[1]*n_cells[2]);i++)
//    {
//        //Fix to directly address array elements

//        del1 = delHolder1[i];
//        del2 = delHolder2[i];
//        del3 = delHolder3[i];
//        del4 = delHolder4[i];
//        del5 = delHolder5[i];
//        del6 = delHolder6[i];

//        ParallelDescriptor::ReduceRealSum(del1);
//        ParallelDescriptor::ReduceRealSum(del2);
//        ParallelDescriptor::ReduceRealSum(del3);
//        ParallelDescriptor::ReduceRealSum(del4);
//        ParallelDescriptor::ReduceRealSum(del5);
//        ParallelDescriptor::ReduceRealSum(del6);

//        delHolder1[i] = del1;
//        delHolder2[i] = del2;
//        delHolder3[i] = del3;
//        delHolder4[i] = del4;
//        delHolder5[i] = del5;
//        delHolder6[i] = del6;
//    }

//    ParallelDescriptor::ReduceRealSum(totalMass);

//    //Print() << "Total mass: " << totalMass << "\n";

//    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) 
//    {
//        const int grid_id = pti.index();
//        const int tile_id = pti.LocalTileIndex();
//        const Box& tile_box  = pti.tilebox();

//        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
//        auto& parts = particle_tile.GetArrayOfStructs();
//        const int Np = parts.numParticles();

//        evaluate_corrs(parts.data(),
//                         ARLIM_3D(tile_box.loVect()),
//                         ARLIM_3D(tile_box.hiVect()),
//                         m_vector_ptrs[grid_id].dataPtr(),
//                         m_vector_size[grid_id].dataPtr(),
//                         ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
//                         ARLIM_3D(m_vector_ptrs[grid_id].hiVect()), 

//                         BL_TO_FORTRAN_3D(particleInstant[pti]),
//                         BL_TO_FORTRAN_3D(particleMeans[pti]),
//                         BL_TO_FORTRAN_3D(particleVars[pti]),

//                         BL_TO_FORTRAN_3D(cellVols[pti]), &Np,&Neff,&n0,&T0,&delt, &steps, delHolder1, delHolder2, delHolder3, delHolder4, delHolder5, delHolder6
//                        );
//    }

//}

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
            p.idata(IntData::sorted) = 1;
            p.idata(IntData::i) = iv[0];
            p.idata(IntData::j) = iv[1];
            p.idata(IntData::k) = iv[2];
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
	    if (p.idata(IntData::sorted)) continue;
            const IntVect& iv = this->Index(p, lev);
            p.idata(IntData::sorted) = 1;
            p.idata(IntData::i) = iv[0];
            p.idata(IntData::j) = iv[1];
            p.idata(IntData::k) = iv[2];
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
    if (not p.idata(IntData::sorted)) return;
    IntVect iv(p.idata(IntData::i), p.idata(IntData::j), p.idata(IntData::k));
    //IntVect iv(AMREX_D_DECL(p.idata(IntData::i), p.idata(IntData::j), p.idata(IntData::k)));
    auto& cell_vector = m_cell_vectors[grid](iv);
    for (int i = 0; i < static_cast<int>(cell_vector.size()); ++i) {
        if (cell_vector[i] == old_index + 1) {
            cell_vector[i] = new_index + 1;
            return;
        }
    }
}


void FhdParticleContainer::WriteParticlesAscii(std::string asciiName)
{
    WriteAsciiFile(asciiName);
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
            if ((iv[0] != p.idata(IntData::i)) or (iv[1] != p.idata(IntData::j)) or (iv[2] != p.idata(IntData::k))) {
                num_wrong += 1;
            }
        }
    }
    
    ParallelDescriptor::ReduceIntSum(num_wrong, ParallelDescriptor::IOProcessorNumber());
    return num_wrong;
}


