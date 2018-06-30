#include "FhdParticleContainer.H"
#include "particle_functions_F.H"

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

    std::mt19937 mt(ParallelDescriptor::MyProc());
    std::uniform_real_distribution<double> unit(0, 1.0);
    std::normal_distribution<double> ndist(0, sqrt(particleInfo.R*particleInfo.T));

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
                
                p.pos(0) = plo[0] + (iv[0]+unit(mt))*dx[0];
                p.pos(1) = plo[1] + (iv[1]+unit(mt))*dx[1];
#if (BL_SPACEDIM == 3)
                p.pos(2) = plo[2] + (iv[2]+unit(mt))*dx[2];
#endif
                p.rdata(RealData::vx) = ndist(mt);
                p.rdata(RealData::vy) = ndist(mt);
                p.rdata(RealData::vz) = ndist(mt);
                p.rdata(RealData::mass) = particleInfo.m; //mass
                p.rdata(RealData::R) = particleInfo.R; //R
                p.rdata(RealData::radius) = particleInfo.d/2.0; //radius
                p.rdata(RealData::accelFactor) = -6*3.14159265359*p.rdata(RealData::radius)/p.rdata(RealData::mass); //acceleration factor (replace with amrex c++ constant for pi...)
                p.rdata(RealData::dragFactor) = -6*3.14159265359*p.rdata(RealData::radius); //drag factor
                p.rdata(RealData::angularVel1) = ndist(mt); //angular velocity 1
                p.rdata(RealData::angularVel2) = ndist(mt); //angular velocity 2
                p.rdata(RealData::angularVel3) = ndist(mt); //angular velocity 2

                AMREX_ASSERT(this->Index(p, lev) == iv);
                
                particle_tile.push_back(p);
            }
        }
    }

}


void FhdParticleContainer::MoveParticles(const Real dt, const Real* dxFluid, const Real* ploFluid, const std::array<MultiFab, AMREX_SPACEDIM>& umac,
                                           std::array<MultiFab, AMREX_SPACEDIM>& umacNodal,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                           const MultiFab& betaCC, //Not necessary but may use later
                                           MultiFab& betaNodal, //Not necessary but may use later
                                           const MultiFab& rho, //Not necessary but may use later
                                           std::array<MultiFab, AMREX_SPACEDIM>& source,
                                           std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp)
{
    

    //UpdateCellVectors();

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();

    //Arg1: Source multifab to be shifted. Arg2: destination multiFab. Arg3: A cell centred multifab for reference (change this later).
    FindNodalValues(umac[0], umacNodal[0], betaCC);
    FindNodalValues(umac[1], umacNodal[1], betaCC);

#if (AMREX_SPACEDIM == 3)
    FindNodalValues(umac[2], umacNodal[2], betaCC);
    FindNodalValues(betaCC, betaNodal, betaCC);
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
                
        move_particles_dsmc(particles.data(), &np,
                       ARLIM_3D(tile_box.loVect()), 
                       ARLIM_3D(tile_box.hiVect()),
                       m_vector_ptrs[grid_id].dataPtr(),
                       m_vector_size[grid_id].dataPtr(),
                       ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                       ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),
                       ZFILL(plo), ZFILL(dx), &dt);

        /*move_particles_fhd(particles.data(), &np,
                         ARLIM_3D(tile_box.loVect()),
                         ARLIM_3D(tile_box.hiVect()),
                         m_vector_ptrs[grid_id].dataPtr(),
                         m_vector_size[grid_id].dataPtr(),
                         ARLIM_3D(m_vector_ptrs[grid_id].loVect()),
                         ARLIM_3D(m_vector_ptrs[grid_id].hiVect()),
                         ZFILL(plo), ZFILL(dx), &dt, ZFILL(ploFluid), ZFILL(dxFluid),
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
                         );*/


        // resize particle vectors after call to move_particles
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
        {
            const auto new_size = m_vector_size[grid_id](iv);
            auto& pvec = m_cell_vectors[grid_id](iv);
            pvec.resize(new_size);
        }
    }

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

void FhdParticleContainer::EvaluateFields(MultiFab& particleMembers,
                              MultiFab& particleDensity,
                              std::array<MultiFab, 3>& particleVelocity,
                              MultiFab& particleTemperature,
                              MultiFab& particleSpeed,
                              MultiFab& cellVols, const Real Neff)
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
                         BL_TO_FORTRAN_3D(particleSpeed[pti]),
                         BL_TO_FORTRAN_3D(cellVols[pti]),&Neff, &Np
                        );

    }
}

void FhdParticleContainer::EvaluateStats(
                              MultiFab& particleMembers,
                              MultiFab& particleDensity,
                              std::array<MultiFab, 3>& particleVelocity,
                              MultiFab& particleTemperature,

                              MultiFab& particleMembersMean,
                              MultiFab& particleDensityMean,
                              std::array<MultiFab, 3>& particleVelocityMean,
                              MultiFab& particleTemperatureMean,
                              MultiFab& particleSpeedMean,

                              MultiFab& particleMembersVar,
                              MultiFab& particleDensityVar,
                              std::array<MultiFab, 3>& particleVelocityVar,
                              MultiFab& particleTemperatureVar,
                              MultiFab& particleSpeedVar,

                              MultiFab& cellVols, const Real Neff, const Real delt, int steps)
{
    const int lev = 0;

    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) 
    {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& parts = particle_tile.GetArrayOfStructs();
        const int Np = parts.numParticles();

        evaluate_stats(parts.data(),
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

                         BL_TO_FORTRAN_3D(particleMembersMean[pti]),
                         BL_TO_FORTRAN_3D(particleDensityMean[pti]),
                         BL_TO_FORTRAN_3D(particleVelocityMean[0][pti]),
                         BL_TO_FORTRAN_3D(particleVelocityMean[1][pti]),
                         BL_TO_FORTRAN_3D(particleVelocityMean[2][pti]),
                         BL_TO_FORTRAN_3D(particleTemperatureMean[pti]),
                         BL_TO_FORTRAN_3D(particleSpeedMean[pti]),

                         BL_TO_FORTRAN_3D(particleMembersVar[pti]),
                         BL_TO_FORTRAN_3D(particleDensityVar[pti]),
                         BL_TO_FORTRAN_3D(particleVelocityVar[0][pti]),
                         BL_TO_FORTRAN_3D(particleVelocityVar[1][pti]),
                         BL_TO_FORTRAN_3D(particleVelocityVar[2][pti]),
                         BL_TO_FORTRAN_3D(particleTemperatureVar[pti]),
                         BL_TO_FORTRAN_3D(particleSpeedVar[pti]),

                         BL_TO_FORTRAN_3D(cellVols[pti]), &Np,&Neff,&delt, &steps
                        );
    }

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




