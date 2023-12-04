#include "StochasticPC.H"

#include <AMReX_TracerParticle_mod_K.H>

using namespace amrex;

void
StochasticPC:: InitParticles (MultiFab& phi_fine)
{
    BL_PROFILE("StochasticPC::InitParticles");

    const int lev = 1;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();

#if (AMREX_SPACEDIM == 2)
    const Real cell_vol = dx[0]*dx[1];
#else
    const Real cell_vol = dx[0]*dx[1]*dx[2];
#endif

    for (MFIter mfi(phi_fine); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();

        const Array4<Real const>& phi_arr = phi_fine.const_array(mfi);
        /*
        Gpu::DeviceVector<unsigned int> counts(tile_box.numPts(), 0);
        unsigned int* pcount = counts.dataPtr();

        Gpu::DeviceVector<unsigned int> offsets(tile_box.numPts()+1, 0);
        unsigned int* poffset = offsets.dataPtr();

        // launch kernel to fill counts with num parts in each cell

        // fill offets
        Gpu::exclusive_scan(counts.begin(), counts.end(), offsets.begin());

        auto num_to_add = offets[tile_box.numPts()];

        auto& particles = GetParticles(lev);
        auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];

        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + num_to_add;
        particle_tile.resize(new_size);

        ParticleType* pstruct = particle_tile.GetArrayOfStructs()().data();

        if (num_to_add == 0) continue;

        // launch kernel to fill data
        amrex::ParallelForRNG(tile_box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
        {
            int ix = i - lo.x;
            int iy = j - lo.y;
            int iz = k - lo.z;
            int nx = hi.x-lo.x+1;
            int ny = hi.y-lo.y+1;
            int nz = hi.z-lo.z+1;
            unsigned int uix = amrex::min(nx-1,amrex::max(0,ix));
            unsigned int uiy = amrex::min(ny-1,amrex::max(0,iy));
            unsigned int uiz = amrex::min(nz-1,amrex::max(0,iz));
            unsigned int cellid = (uix * ny + uiy) * nz + uiz;

            auto start = poffset[cellid];
            auto stop = poffset[cellid+1];

            for (unsigned int ip = start; ip < stop; ++ip) {
                ParticleType& p = pstruct[ip];
                p.pos(0) = ...
            }
        }
        */

        Gpu::HostVector<ParticleType> host_particles;
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) 
        {
              // Real rannum = amrex::Random();
              Real rannum = 0.;
              int npart_in_cell = int(phi_arr(iv,0)*cell_vol+rannum);

              if (npart_in_cell > 0) {
#if 1
                 for (int npart = 0; npart < npart_in_cell; npart++)
                 {
#if (AMREX_SPACEDIM == 2)
                     Real r[2] = {amrex::Random(), amrex::Random()}; 
#elif (AMREX_SPACEDIM == 3)
                     Real r[3] = {amrex::Random(), amrex::Random(), amrex::Random()}; 
#endif
                     AMREX_D_TERM( Real x = plo[0] + (iv[0] + r[0])*dx[0];,
                                   Real y = plo[1] + (iv[1] + r[1])*dx[1];,
                                   Real z = plo[2] + (iv[2] + r[2])*dx[2];);

                     ParticleType p;
                     p.id()  = ParticleType::NextID();
                     p.cpu() = ParallelDescriptor::MyProc();
   
                     AMREX_D_TERM( p.pos(0) = x;,
                                   p.pos(1) = y;,
                                   p.pos(2) = z;);

                     AMREX_D_TERM( p.rdata(RealIdx::xold) = x;,
                                   p.rdata(RealIdx::yold) = y;,
                                   p.rdata(RealIdx::zold) = z;);

                     host_particles.push_back(p);
                 }
#else
                   amrex::ParallelForRNG( npart_in_cell,
                   [=] AMREX_GPU_DEVICE (int npart, RandomEngine const& engine) noexcept
                   {
#if (AMREX_SPACEDIM == 2)
                     Real r[2] = {amrex::Random(engine), amrex::Random(engine)}; 
#elif (AMREX_SPACEDIM == 3)
                     Real r[3] = {amrex::Random(engine), amrex::Random(engine), amrex::Random(engine)}; 
#endif
                     AMREX_D_TERM( Real x = plo[0] + (iv[0] + r[0])*dx[0];,
                                   Real y = plo[1] + (iv[1] + r[1])*dx[1];,
                                   Real z = plo[2] + (iv[2] + r[2])*dx[2];);

                     ParticleType p;
                     p.id()  = ParticleType::NextID();
                     p.cpu() = ParallelDescriptor::MyProc();
   
                     AMREX_D_TERM( p.pos(0) = x;,
                                   p.pos(1) = y;,
                                   p.pos(2) = z;);

                     AMREX_D_TERM( p.rdata(RealIdx::xold) = x;,
                                   p.rdata(RealIdx::yold) = y;,
                                   p.rdata(RealIdx::zold) = z;);

                     host_particles.push_back(p);
                 });
#endif
              }
        }

        auto& particles = GetParticles(lev);
        auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        particle_tile.resize(new_size);

        amrex::Print() << "INIT: OLD SIZE OF PARTICLES " << old_size << std::endl;
        amrex::Print() << "INIT: NEW SIZE OF PARTICLES " << new_size << std::endl;

        Gpu::copy(Gpu::hostToDevice,
                  host_particles.begin(),
                  host_particles.end(),
                  particle_tile.GetArrayOfStructs().begin() + old_size);
    }
}

void
StochasticPC::AddParticles (MultiFab& phi_fine, BoxArray& ba_to_exclude)
{
    BL_PROFILE("StochasticPC::AddParticles");
    const int lev = 1;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();

    Real cell_vol = dx[0]*dx[1];

    for (MFIter mfi(phi_fine); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();

        if (!ba_to_exclude.contains(tile_box)) {

        const Array4<Real const>& phi_arr = phi_fine.const_array(mfi);

        Gpu::HostVector<ParticleType> host_particles;
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) 
        {
              // Real rannum = amrex::Random();
              Real rannum = 0.;
              int npart_in_cell = int(phi_arr(iv,0)*cell_vol+rannum);

              if (npart_in_cell > 0) {
#if 1
                 for (int npart = 0; npart < npart_in_cell; npart++)
                 {
#if (AMREX_SPACEDIM == 2)
                     Real r[2] = {amrex::Random(), amrex::Random()}; 
#elif (AMREX_SPACEDIM == 3)
                     Real r[3] = {amrex::Random(), amrex::Random(), amrex::Random()}; 
#endif
                     AMREX_D_TERM( Real x = plo[0] + (iv[0] + r[0])*dx[0];,
                                   Real y = plo[1] + (iv[1] + r[1])*dx[1];,
                                   Real z = plo[2] + (iv[2] + r[2])*dx[2];);

                     ParticleType p;
                     p.id()  = ParticleType::NextID();
                     p.cpu() = ParallelDescriptor::MyProc();
   
                     AMREX_D_TERM( p.pos(0) = x;,
                                   p.pos(1) = y;,
                                   p.pos(2) = z;);

                     AMREX_D_TERM( p.rdata(RealIdx::xold) = x;,
                                   p.rdata(RealIdx::yold) = y;,
                                   p.rdata(RealIdx::zold) = z;);

                     host_particles.push_back(p);
                 }
#else
                   amrex::ParallelForRNG( npart_in_cell,
                   [=] AMREX_GPU_DEVICE (int npart, RandomEngine const& engine) noexcept
                   {
#if (AMREX_SPACEDIM == 2)
                     Real r[2] = {amrex::Random(engine), amrex::Random(engine)}; 
#elif (AMREX_SPACEDIM == 3)
                     Real r[3] = {amrex::Random(engine), amrex::Random(engine), amrex::Random(engine)}; 
#endif
                     AMREX_D_TERM( Real x = plo[0] + (iv[0] + r[0])*dx[0];,
                                   Real y = plo[1] + (iv[1] + r[1])*dx[1];,
                                   Real z = plo[2] + (iv[2] + r[2])*dx[2];);

                     ParticleType p;
                     p.id()  = ParticleType::NextID();
                     p.cpu() = ParallelDescriptor::MyProc();
   
                     AMREX_D_TERM( p.pos(0) = x;,
                                   p.pos(1) = y;,
                                   p.pos(2) = z;);

                     AMREX_D_TERM( p.rdata(RealIdx::xold) = x;,
                                   p.rdata(RealIdx::yold) = y;,
                                   p.rdata(RealIdx::zold) = z;);

                     host_particles.push_back(p);
                 });
#endif
              } // npart_in_cell
        } // iv

        auto& particles = GetParticles(lev);
        auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        particle_tile.resize(new_size);

        Gpu::copy(Gpu::hostToDevice,
                  host_particles.begin(),
                  host_particles.end(),
                  particle_tile.GetArrayOfStructs().begin() + old_size);
        } // not in ba_to_exclude
    } // mfi
}

void
StochasticPC::RemoveParticlesNotInBA (const BoxArray& ba_to_keep)
{
    BL_PROFILE("StochasticPC::RemoveParticles");
    const int lev = 1;

    for(ParIterType pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& ptile = ParticlesAt(lev, pti);
        auto& aos  = ptile.GetArrayOfStructs();
        const int np = aos.numParticles();
        auto *pstruct = aos().data();

        if (!ba_to_keep.contains(pti.tilebox())) {
            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = pstruct[i];
                p.id() = -1;
            });
        }
    }
    Redistribute();
}

void
StochasticPC::RefluxFineToCrse (const BoxArray& ba_to_keep, MultiFab& phi_for_reflux) 
{
    BL_PROFILE("StochasticPC::RefluxFineToCrse");
    const int lev = 1;

    const auto dx = Geom(lev).CellSizeArray();

    for(ParIterType pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& ptile = ParticlesAt(lev, pti);
        auto& aos  = ptile.GetArrayOfStructs();
        const int np = aos.numParticles();
        auto *pstruct = aos().data();

        if (!ba_to_keep.contains(pti.tilebox())) {

            Array4<Real> phi_arr = phi_for_reflux.array(pti.index());

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = pstruct[i];
    
                IntVect old_pos(static_cast<int>(p.rdata(RealIdx::xold)/dx[0]),static_cast<int>(p.rdata(RealIdx::yold)/dx[0]));
                IntVect new_pos(static_cast<int>(p.pos(0)              /dx[0]),static_cast<int>(p.pos(1)              /dx[1]));

                if ( ba_to_keep.contains(old_pos) && !ba_to_keep.contains(new_pos)) {
                   //amrex::Print() <<" Particle crossed out of fine grids " << old_pos << " " << new_pos << std::endl;
                   phi_arr(new_pos,0) += 1.;
                }
    
            });
        } // if not in ba_to_keep
    } // pti
}

void
StochasticPC::RefluxCrseToFine (const BoxArray& ba_to_keep, MultiFab& phi_for_reflux) 
{
    BL_PROFILE("StochasticPC::RefluxCrseToFine");
    const int lev = 1;

    const auto dx = Geom(lev).CellSizeArray();

    amrex::Print() << " PHI FOR REFLUX " << phi_for_reflux.boxArray() << std::endl;

    for(ParIterType pti(*this, lev); pti.isValid(); ++pti)
    {
        auto const& gid   = pti.index();
        auto& ptile = ParticlesAt(lev, pti);
        auto& aos  = ptile.GetArrayOfStructs();
        const int np = aos.numParticles();
        auto *pstruct = aos().data();

        if (ba_to_keep.contains(pti.tilebox())) {

            Array4<Real> phi_arr = phi_for_reflux.array(gid);
            amrex::Print() << "BOX(Phi_arr) " << Box(phi_arr) << std::endl;

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = pstruct[i];
    
                IntVect old_pos(static_cast<int>(p.rdata(RealIdx::xold)/dx[0]),static_cast<int>(p.rdata(RealIdx::yold)/dx[0]));
                IntVect new_pos(static_cast<int>(p.pos(0)              /dx[0]),static_cast<int>(p.pos(1)              /dx[1]));
    
                if (!ba_to_keep.contains(old_pos) && ba_to_keep.contains(new_pos)) {
                      if (!Box(phi_arr).contains(old_pos)) {
                         amrex::Print() <<" Particle crossed into fine grids from " << old_pos << " to " << new_pos << std::endl;
                         amrex::Print() << "Grid index is " << gid << std::endl;
                         amrex::Print() << "BUT WE HAVE A PROBLEM " << std::endl;
                      }
                      phi_arr(old_pos,0) -= 1.;
                }
            });

        } // if in ba_to_keep
    } // pti
}

void
StochasticPC::AdvectWithRandomWalk (int lev, Real dt)
{
    BL_PROFILE("StochasticPC::AdvectWithRandomWalk");
    const auto dx = Geom(lev).CellSizeArray();

    for(ParIterType pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& ptile = ParticlesAt(lev, pti);
        auto& aos  = ptile.GetArrayOfStructs();
        const int np = aos.numParticles();
        auto *pstruct = aos().data();

        Real stddev = std::sqrt(dt);

#if 0
        amrex::ParallelForRNG( np, [=] AMREX_GPU_DEVICE (int i, RandomEngine const& engine) noexcept
        {
            ParticleType& p = pstruct[i];
            AMREX_D_TERM( p.rdata(RealIdx::xold) = p.pos(0);,
                          p.rdata(RealIdx::yold) = p.pos(1);,
                          p.rdata(RealIdx::zold) = p.pos(2););

            amrex::Real incx = amrex::RandomNormal(0.,stddev,engine);
            amrex::Real incy = amrex::RandomNormal(0.,stddev,engine);

            incx = std::max(-dx[0], std::min( dx[0], incx));
            incy = std::max(-dx[1], std::min( dx[1], incy));

            p.pos(0) += static_cast<ParticleReal> (incx);
            p.pos(1) += static_cast<ParticleReal> (incy);

#if AMREX_SPACEDIM > 2
            amrex::Real incz = amrex::RandomNormal(0.,stddev,engine);
            incz = std::max(-dx[2], std::min( dx[2], incz));
            p.pos(2) += static_cast<ParticleReal> (incz);
#endif
        });

#else
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = pstruct[i];
            AMREX_D_TERM( p.rdata(RealIdx::xold) = p.pos(0);,
                          p.rdata(RealIdx::yold) = p.pos(1);,
                          p.rdata(RealIdx::zold) = p.pos(2););

            amrex::Real incx = amrex::RandomNormal(0.,stddev);
            amrex::Real incy = amrex::RandomNormal(0.,stddev);

            incx = std::max(-dx[0], std::min( dx[0], incx));
            incy = std::max(-dx[1], std::min( dx[1], incy));

            p.pos(0) += static_cast<ParticleReal> (incx);
            p.pos(1) += static_cast<ParticleReal> (incy);

#if (AMREX_SPACEDIM > 2)
            amrex::Real incz = amrex::RandomNormal(0.,stddev,engine);
            incz = std::max(-dx[2], std::min( dx[2], incz));
            p.pos(2) += static_cast<ParticleReal> (incz);
#endif

            // p.pos(0) += static_cast<ParticleReal> (1.0 * dx[0]);
            // p.pos(1) += static_cast<ParticleReal> ((2*amrex::Random(engine)-1)*move_dir*dx[1]);
#if AMREX_SPACEDIM > 2
            // p.pos(2) += static_cast<ParticleReal> ((2*amrex::Random(engine)-1)*move_dir*dx[2]);
#endif
        });
#endif
    } // pti

    // amrex::Print() << "AT  END  OF ADVECT, NUM PART AT LEVEL 0: " << NumberOfParticlesAtLevel(0) << std::endl;;
    // amrex::Print() << "AT  END  OF ADVECT, NUM PART AT LEVEL 1: " << NumberOfParticlesAtLevel(1) << std::endl;;
}
