#include "FhdParticleContainer.H"
#include "particle_functions_F.H"
#include "rng_functions_F.H"

using namespace amrex;


FhdParticleContainer::FhdParticleContainer(const Geometry & geom,
                              const DistributionMapping & dmap,
                              const BoxArray            & ba)
    : ParticleContainer<RealData::ncomps, IntData::ncomps> (geom, dmap, ba)
{}

void FhdParticleContainer::InitParticles(const int ppc, species particleInfo)
{
    
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();
    const Real* plo = geom.ProbLo();

    //std::mt19937 mt(ParallelDescriptor::MyProc()*2);
    //std::uniform_real_distribution<double> unit(0, 1.0);
    //std::normal_distribution<double> ndist(0, sqrt(particleInfo.R*particleInfo.T));

    double totalEnergy = 0;

    double cosTheta, sinTheta, cosPhi, sinPhi;    

    //double initTemp = 0;
    //double pc = 0;

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
        {
            for (int i_part=0; i_part<ppc;i_part++) {
                
                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                p.idata(IntData::sorted) = 0;
                p.idata(IntData::species) = 0;
                
                p.pos(0) = plo[0] + (iv[0]+get_uniform_func())*dx[0];
                p.pos(1) = plo[1] + (iv[1]+get_uniform_func())*dx[1];
#if (BL_SPACEDIM == 3)
                p.pos(2) = plo[2] + (iv[2]+get_uniform_func())*dx[2];
#endif
                //p.rdata(RealData::vx) = sqrt(particleInfo.R*particleInfo.T)*get_particle_normal_func();
                //p.rdata(RealData::vy) = sqrt(particleInfo.R*particleInfo.T)*get_particle_normal_func();
                //p.rdata(RealData::vz) = sqrt(particleInfo.R*particleInfo.T)*get_particle_normal_func();

                p.rdata(RealData::vx) = 0;
                p.rdata(RealData::vy) = 0;
                p.rdata(RealData::vz) = 0;

                totalEnergy = totalEnergy + p.rdata(RealData::vx)*p.rdata(RealData::vx) + p.rdata(RealData::vy)*p.rdata(RealData::vy) + p.rdata(RealData::vz)*p.rdata(RealData::vz);

                //initTemp 

                p.rdata(RealData::mass) = particleInfo.m; //mass
                p.rdata(RealData::R) = particleInfo.R; //R
                p.rdata(RealData::radius) = particleInfo.d/2.0; //radius
                p.rdata(RealData::accelFactor) = -6*3.14159265359*p.rdata(RealData::radius)/p.rdata(RealData::mass); //acceleration factor (replace with amrex c++ constant for pi...)
                p.rdata(RealData::dragFactor) = -6*3.14159265359*p.rdata(RealData::radius); //drag factor
                p.rdata(RealData::angularVel1) = 0; //angular velocity 1
                p.rdata(RealData::angularVel2) = 0; //angular velocity 2
                p.rdata(RealData::angularVel3) = 0; //angular velocity 2

                get_angles(&cosTheta, &sinTheta, &cosPhi, &sinPhi);

                //p.rdata(RealData::dirx) = sinTheta*cosPhi; //Unit vector giving orientation
                //p.rdata(RealData::diry) = sinTheta*sinPhi; 
                //p.rdata(RealData::dirz) = cosTheta;

                p.rdata(RealData::dirx) = 1; //Unit vector giving orientation
                p.rdata(RealData::diry) = 0; 
                p.rdata(RealData::dirz) = 0;

                p.rdata(RealData::propulsion) = -p.rdata(RealData::accelFactor)*9e-4*1e-1;  //propulsive acceleration

                AMREX_ASSERT(this->Index(p, lev) == iv);
                
                particle_tile.push_back(p);
            }
        }
    }

    Print() << "Initial energy: " << totalEnergy << "\n";

}


void FhdParticleContainer::MoveParticles(const Real dt, const Real* dxFluid, const Real* ploFluid, const std::array<MultiFab, AMREX_SPACEDIM>& umac,
                                           std::array<MultiFab, AMREX_SPACEDIM>& umacNodal,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                           const MultiFab& betaCC, //Not necessary but may use later
                                           MultiFab& betaNodal, //Not necessary but may use later
                                           const MultiFab& rho, //Not necessary but may use later
                                           std::array<MultiFab, AMREX_SPACEDIM>& source,
                                           std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                                           const surface* surfaceList, const int surfaceCount)
{
    

    //UpdateCellVectors();

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();
    const Real* phi = Geom(lev).ProbHi();

#ifndef DSMC

    //Arg1: Source multifab to be shifted. Arg2: destination multiFab. Arg3: A cell centred multifab for reference (change this later).
    FindNodalValues(umac[0], umacNodal[0], betaCC);
    FindNodalValues(umac[1], umacNodal[1], betaCC);

#if (AMREX_SPACEDIM == 3)
    FindNodalValues(umac[2], umacNodal[2], betaCC);
    FindNodalValues(betaCC, betaNodal, betaCC);
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
        move_particles_fhd(particles.data(), &np,
                         ARLIM_3D(tile_box.loVect()),
                         ARLIM_3D(tile_box.hiVect()),
                         m_vector_ptrs[grid_id].dataPtr(),
                         m_vector_size[grid_id].dataPtr(),
                         ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                         ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),
                         ZFILL(plo), ZFILL(phi), ZFILL(dx), &dt, ZFILL(ploFluid), ZFILL(dxFluid),
                         BL_TO_FORTRAN_3D(umacNodal[0][pti]),
                         BL_TO_FORTRAN_3D(umacNodal[1][pti]),
#if (AMREX_SPACEDIM == 3)
                         BL_TO_FORTRAN_3D(umacNodal[2][pti]),
#endif
                         BL_TO_FORTRAN_3D(RealFaceCoords[0][pti]),
                         BL_TO_FORTRAN_3D(RealFaceCoords[1][pti]),
#if (AMREX_SPACEDIM == 3)
                         BL_TO_FORTRAN_3D(RealFaceCoords[2][pti]),
#endif
                         BL_TO_FORTRAN_3D(betaCC[pti]),
                         BL_TO_FORTRAN_3D(rho[pti]),

                         BL_TO_FORTRAN_3D(sourceTemp[0][pti]),
                         BL_TO_FORTRAN_3D(sourceTemp[1][pti])
#if (AMREX_SPACEDIM == 3)
                         , BL_TO_FORTRAN_3D(sourceTemp[2][pti])
#endif
                         , surfaceList, &surfaceCount
                         );
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
}

void FhdParticleContainer::InitCollisionCells(
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

void FhdParticleContainer::InitializeFields(MultiFab& particleMembers,
                              MultiFab& particleDensity,
                              std::array<MultiFab, 3>& particleVelocity,
                              MultiFab& particleTemperature,
                              MultiFab& cellVols, const species particleInfo)
{

    //UpdateCellVectors();
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
                         BL_TO_FORTRAN_3D(particleMembers[pti]),
                         BL_TO_FORTRAN_3D(particleDensity[pti]),
                         BL_TO_FORTRAN_3D(particleVelocity[0][pti]),
                         BL_TO_FORTRAN_3D(particleVelocity[1][pti]),
                         BL_TO_FORTRAN_3D(particleVelocity[2][pti]),
                         BL_TO_FORTRAN_3D(particleTemperature[pti]),
                         BL_TO_FORTRAN_3D(cellVols[pti]),&particleInfo.Neff, &Np, &particleInfo.R, &particleInfo.T
                        );

    }
}

void FhdParticleContainer::EvaluateStats(
                              MultiFab& particleMembers,
                              MultiFab& particleDensity,
                              std::array<MultiFab, 3>& particleVelocity,
                              MultiFab& particleTemperature,
                              std::array<MultiFab, 3>& particleMomentum,
                              MultiFab& particleEnergy,
                              MultiFab& particlePressure,

                              MultiFab& particleMembersMean,
                              MultiFab& particleDensityMean,
                              std::array<MultiFab, 3>& particleVelocityMean,
                              MultiFab& particleTemperatureMean,
                              std::array<MultiFab, 3>& particleMomentumMean,
                              MultiFab& particleEnergyMean,
                              MultiFab& particlePressureMean,

                              MultiFab& particleMembersVar,
                              MultiFab& particleDensityVar,
                              std::array<MultiFab, 3>& particleVelocityVar,
                              MultiFab& particleTemperatureVar,
                              std::array<MultiFab, 3>& particleMomentumVar,
                              MultiFab& particleEnergyVar,
                              MultiFab& particlePressureVar,

                              MultiFab& particleGVar, 
                              MultiFab& particleKGCross,
                              MultiFab& particleKRhoCross, 
                              MultiFab& particleRhoGCross,
                              MultiFab& particleSpatialCross1,
                              MultiFab& particleSpatialCross2,
                              MultiFab& particleMembraneFlux,

                              MultiFab& cellVols, species particleInfo, const Real delt, int steps)
{
    const int lev = 0;
    const double Neff = particleInfo.Neff;
    const double n0 = particleInfo.n0;
    const double T0 = particleInfo.T;

    double del1 = 0;
    double del2 = 0;
    double del3 = 0;

    double tp = 0;
    double te = 0;
    double tm = 0;

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
                         BL_TO_FORTRAN_3D(particleMembers[pti]),
                         BL_TO_FORTRAN_3D(particleDensity[pti]),
                         BL_TO_FORTRAN_3D(particleVelocity[0][pti]),
                         BL_TO_FORTRAN_3D(particleVelocity[1][pti]),
                         BL_TO_FORTRAN_3D(particleVelocity[2][pti]),
                         BL_TO_FORTRAN_3D(particleTemperature[pti]),
                         BL_TO_FORTRAN_3D(particleMomentum[0][pti]),
                         BL_TO_FORTRAN_3D(particleMomentum[1][pti]),
                         BL_TO_FORTRAN_3D(particleMomentum[2][pti]),
                         BL_TO_FORTRAN_3D(particleEnergy[pti]),
                         BL_TO_FORTRAN_3D(particlePressure[pti]),
                         BL_TO_FORTRAN_3D(cellVols[pti]),&Neff, &Np, &del1, &del2, &tm, &te
                        );

    }

    ParallelDescriptor::ReduceRealSum(del1);
    ParallelDescriptor::ReduceRealSum(tp);
    ParallelDescriptor::ReduceRealSum(te);
    ParallelDescriptor::ReduceRealSum(tm);

    //Print() << "Total particles: " << tp << "\n";
    //Print() << "Total momentum: " << tm << "\n";
    //Print() << "Total energy: " << te << "\n";

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

                         BL_TO_FORTRAN_3D(particleMembers[pti]),
                         BL_TO_FORTRAN_3D(particleDensity[pti]),
                         BL_TO_FORTRAN_3D(particleVelocity[0][pti]),
                         BL_TO_FORTRAN_3D(particleVelocity[1][pti]),
                         BL_TO_FORTRAN_3D(particleVelocity[2][pti]),
                         BL_TO_FORTRAN_3D(particleTemperature[pti]),
                         BL_TO_FORTRAN_3D(particleMomentum[0][pti]),
                         BL_TO_FORTRAN_3D(particleMomentum[1][pti]),
                         BL_TO_FORTRAN_3D(particleMomentum[2][pti]),
                         BL_TO_FORTRAN_3D(particleEnergy[pti]),
                         BL_TO_FORTRAN_3D(particlePressure[pti]),

                         BL_TO_FORTRAN_3D(particleMembersMean[pti]),
                         BL_TO_FORTRAN_3D(particleDensityMean[pti]),
                         BL_TO_FORTRAN_3D(particleVelocityMean[0][pti]),
                         BL_TO_FORTRAN_3D(particleVelocityMean[1][pti]),
                         BL_TO_FORTRAN_3D(particleVelocityMean[2][pti]),
                         BL_TO_FORTRAN_3D(particleTemperatureMean[pti]),
                         BL_TO_FORTRAN_3D(particleMomentumMean[0][pti]),
                         BL_TO_FORTRAN_3D(particleMomentumMean[1][pti]),
                         BL_TO_FORTRAN_3D(particleMomentumMean[2][pti]),
                         BL_TO_FORTRAN_3D(particleEnergyMean[pti]),
                         BL_TO_FORTRAN_3D(particlePressureMean[pti]),

                         BL_TO_FORTRAN_3D(particleMembersVar[pti]),
                         BL_TO_FORTRAN_3D(particleDensityVar[pti]),
                         BL_TO_FORTRAN_3D(particleVelocityVar[0][pti]),
                         BL_TO_FORTRAN_3D(particleVelocityVar[1][pti]),
                         BL_TO_FORTRAN_3D(particleVelocityVar[2][pti]),
                         BL_TO_FORTRAN_3D(particleTemperatureVar[pti]),
                         BL_TO_FORTRAN_3D(particleMomentumVar[0][pti]),
                         BL_TO_FORTRAN_3D(particleMomentumVar[1][pti]),
                         BL_TO_FORTRAN_3D(particleMomentumVar[2][pti]),
                         BL_TO_FORTRAN_3D(particleEnergyVar[pti]),

                         BL_TO_FORTRAN_3D(particleGVar[pti]),
                         BL_TO_FORTRAN_3D(particleKGCross[pti]),
                         BL_TO_FORTRAN_3D(particleKRhoCross[pti]),
                         BL_TO_FORTRAN_3D(particleRhoGCross[pti]),
                         BL_TO_FORTRAN_3D(particleSpatialCross1[pti]),
                         BL_TO_FORTRAN_3D(particleSpatialCross2[pti]),
                         BL_TO_FORTRAN_3D(particleMembraneFlux[pti]),

                         BL_TO_FORTRAN_3D(cellVols[pti]), &Np,&Neff,&n0,&T0,&delt, &steps, &del1, &del2
                        );
    }

    ParallelDescriptor::ReduceRealSum(del2);

    //Print() << "del1: " << del1 << "\n";

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

                         BL_TO_FORTRAN_3D(particleMembers[pti]),
                         BL_TO_FORTRAN_3D(particleDensity[pti]),
                         BL_TO_FORTRAN_3D(particleVelocity[0][pti]),
                         BL_TO_FORTRAN_3D(particleVelocity[1][pti]),
                         BL_TO_FORTRAN_3D(particleVelocity[2][pti]),
                         BL_TO_FORTRAN_3D(particleTemperature[pti]),
                         BL_TO_FORTRAN_3D(particleMomentum[0][pti]),
                         BL_TO_FORTRAN_3D(particleMomentum[1][pti]),
                         BL_TO_FORTRAN_3D(particleMomentum[2][pti]),
                         BL_TO_FORTRAN_3D(particleEnergy[pti]),

                         BL_TO_FORTRAN_3D(particleMembersMean[pti]),
                         BL_TO_FORTRAN_3D(particleDensityMean[pti]),
                         BL_TO_FORTRAN_3D(particleVelocityMean[0][pti]),
                         BL_TO_FORTRAN_3D(particleVelocityMean[1][pti]),
                         BL_TO_FORTRAN_3D(particleVelocityMean[2][pti]),
                         BL_TO_FORTRAN_3D(particleTemperatureMean[pti]),
                         BL_TO_FORTRAN_3D(particleMomentumMean[0][pti]),
                         BL_TO_FORTRAN_3D(particleMomentumMean[1][pti]),
                         BL_TO_FORTRAN_3D(particleMomentumMean[2][pti]),
                         BL_TO_FORTRAN_3D(particleEnergyMean[pti]),

                         BL_TO_FORTRAN_3D(particleMembersVar[pti]),
                         BL_TO_FORTRAN_3D(particleDensityVar[pti]),
                         BL_TO_FORTRAN_3D(particleVelocityVar[0][pti]),
                         BL_TO_FORTRAN_3D(particleVelocityVar[1][pti]),
                         BL_TO_FORTRAN_3D(particleVelocityVar[2][pti]),
                         BL_TO_FORTRAN_3D(particleTemperatureVar[pti]),
                         BL_TO_FORTRAN_3D(particleMomentumVar[0][pti]),
                         BL_TO_FORTRAN_3D(particleMomentumVar[1][pti]),
                         BL_TO_FORTRAN_3D(particleMomentumVar[2][pti]),
                         BL_TO_FORTRAN_3D(particleEnergyVar[pti]),

                         BL_TO_FORTRAN_3D(particleGVar[pti]),
                         BL_TO_FORTRAN_3D(particleKGCross[pti]),
                         BL_TO_FORTRAN_3D(particleKRhoCross[pti]),
                         BL_TO_FORTRAN_3D(particleRhoGCross[pti]),
                         BL_TO_FORTRAN_3D(particleSpatialCross1[pti]),
                         BL_TO_FORTRAN_3D(particleSpatialCross2[pti]),

                         BL_TO_FORTRAN_3D(cellVols[pti]), &Np,&Neff,&n0,&T0,&delt, &steps, &del1, &del2, &del3
                        );
    }

    ParallelDescriptor::ReduceRealSum(del3);

    //MultiFab::Copy(particleMembraneFlux, particleTemperature, 0, 0, 1, 0);
    //MultiFab::Multiply(particleMembraneFlux,particleDensity,0,0,1,0);

    //Print() << "del3: " << del3 << "\n";

}

void FhdParticleContainer::UpdateCellVectors()
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
            p.idata(0) = 1;
            // note - use 1-based indexing for convenience with Fortran
            m_cell_vectors[pti.index()](iv).push_back(pindex + 1);
        }
    }
    
    UpdateFortranStructures();
}

void FhdParticleContainer::UpdateFortranStructures()
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
            if (p.idata(0)) continue;
            const IntVect& iv = this->Index(p, lev);
            p.idata(0) = 1;
            // note - use 1-based indexing for convenience with Fortran
            m_cell_vectors[pti.index()](iv).push_back(pindex + 1);
        }
    }

    UpdateFortranStructures();
}


void FhdParticleContainer::WriteParticlesAscii(int n)
{
    const std::string& pltfile = amrex::Concatenate("particles", n, 5);
    WriteAsciiFile(pltfile);
}




