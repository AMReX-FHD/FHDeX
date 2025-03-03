#include "StochasticPC.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_Math.H>
#include <AMReX_TracerParticle_mod_K.H>

using namespace amrex;

void
StochasticPC::InitParticles (MultiFab& phi_fine)
{
    amrex::Print() << "calling InitParrticles" << std::endl;
    amrex::Real factor = -1.;
    AddParticles(phi_fine, BoxArray{},factor);
}

void
StochasticPC::ColorParticlesWithPhi (MultiFab const& phi)
{
    BL_PROFILE("StochasticPC::ColorParticlesWithPhi");
    const int lev = 1;
    const auto dx = Geom(lev).CellSizeArray();

    amrex::Print() << "PHIARR BOX " << phi.boxArray() << std::endl;

    for(ParIterType pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& ptile = ParticlesAt(lev, pti);
        auto& aos  = ptile.GetArrayOfStructs();
        const int np = aos.numParticles();
        auto *pstruct = aos().data();

        const Array4<Real const>& phi_arr = phi.const_array(pti.index());

        // amrex::Print() << "PART BOX " << pti.tilebox() << std::endl;
        // amrex::Print() << "PTI INDEX " << pti.index() << std::endl;

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int n)
        {
            ParticleType& p = pstruct[n];
            int i = static_cast<int>(p.pos(0) / dx[0]);
            int j = static_cast<int>(p.pos(1) / dx[1]);
            int k = 0;
            p.rdata(RealIdx::zold) = phi_arr(i,j,k);
        });
    }
}

void
StochasticPC:: AddParticles (MultiFab& phi_fine, const BoxArray& ba_to_exclude, amrex::Real factor)
{
    BL_PROFILE("StochasticPC::AddParticles");

    const int lev = 1;
    const auto dx = Geom(lev).CellSizeArray();
    const auto plo = Geom(lev).ProbLoArray();

#if (AMREX_SPACEDIM == 2)
    const Real cell_vol = dx[0]*dx[1];
#else
    const Real cell_vol = dx[0]*dx[1]*dx[2];
#endif

    // We use sum to count how much phi is gained/lost when we use phi to compute an integer
    // number of particles
    //   Gpu::DeviceVector<Real> my_sum(1, 0.);
    //   Real* sum = my_sum.dataPtr();

    // We need to allow particles to be created outside the domain in cells next
    // to the particle region
    Box gdomain(Geom(lev).Domain());
    if (Geom(lev).isPeriodic(0)) gdomain.grow(0,1);
    if (Geom(lev).isPeriodic(1)) gdomain.grow(1,1);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi_fine); mfi.isValid(); ++mfi)
    {
        Box tile_box  = mfi.tilebox() & gdomain;

        BoxArray tba(tile_box);
        if (ba_to_exclude.contains(tba)) {continue;}

        if (!m_reflux_particle_locator.isValid(ba_to_exclude)) {
            m_reflux_particle_locator.build(ba_to_exclude, Geom(lev));
        }
        m_reflux_particle_locator.setGeometry(Geom(lev));

        auto assign_grid = m_reflux_particle_locator.getGridAssignor();

        const Array4<Real const>& phi_arr = phi_fine.const_array(mfi);

        // count the number of particles to create in each cell
        auto flat_index = FlatIndex(tile_box);

        Gpu::DeviceVector<unsigned int> counts(tile_box.numPts()+1, 0);
        unsigned int* pcount = counts.dataPtr();

        amrex::ParallelForRNG(tile_box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
        {
            if (assign_grid(IntVect(AMREX_D_DECL(i, j, k))) >= 0) {return;}
            Real rannum = amrex::Random(engine);
            int npart_in_cell = int(phi_arr(i,j,k,0)*cell_vol+rannum);
            pcount[flat_index(i, j, k)] += npart_in_cell;
            // if (phi_arr(i,j,k) > 0.) {
            //     amrex::Print() << " IJK/NPART/PHI/RAN " << IntVect(i,j) << " " << npart_in_cell << " given phi " <<
            //         (phi_arr(i,j,k,0)*cell_vol) << " " << rannum << std::endl;
            // }
            // sum[0] += npart_in_cell - (phi_arr(i,j,k,0)*cell_vol);
        });

        // fill offsets
        Gpu::DeviceVector<unsigned int> offsets(tile_box.numPts()+1, 0);
        Gpu::exclusive_scan(counts.begin(), counts.end(), offsets.begin());

        // the last offset is the total number of particles to add
        unsigned int num_to_add;
        Gpu::copy(Gpu::deviceToHost, offsets.begin() + tile_box.numPts(), offsets.end(), &num_to_add);
        if (num_to_add == 0) continue;

        // Get the ids and cpu numbers to assign
        int my_cpu = ParallelDescriptor::MyProc();
        Long id_start;
#ifdef AMREX_USE_OMP
#pragma omp critical (init_particles_next_id)
#endif
        {
            id_start = ParticleType::NextID();
            ParticleType::NextID(id_start+num_to_add);
        }

        // resize particle storage
        auto& particles = GetParticles(lev);
        auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + num_to_add;
        particle_tile.resize(new_size);
        amrex::Print() << "INIT: NEW SIZE OF PARTICLES IN TILE BOX " << tile_box << " " << new_size << std::endl;


        int ext_pot = 1;
        amrex::Real alpha = .3;
        amrex::Real beta = .7;
        amrex::Real gamma = 5.e-4;


        // now fill in the data
        ParticleType* pstruct = particle_tile.GetArrayOfStructs()().data() + old_size;
        unsigned int* poffset = offsets.dataPtr();
        amrex::ParallelForRNG(tile_box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
        {
            auto cellid = flat_index(i, j, k);

            auto start = poffset[cellid];
            auto stop = poffset[cellid+1];

            for (unsigned int ip = start; ip < stop; ++ip) {
                ParticleType& p = pstruct[ip];

#if (AMREX_SPACEDIM == 2)
                Real r[2] = {amrex::Random(engine), amrex::Random(engine)};
#elif (AMREX_SPACEDIM == 3)
                Real r[3] = {amrex::Random(engine), amrex::Random(engine), amrex::Random(engine)};
#endif
                if(factor > 0.)
                {
                    Real xm = plo[0] + i*dx[0];
                    Real xp = xm + dx[0];
                    Real ym = plo[1] + j*dx[1];
                    Real yp = ym + dx[1];
                    Real vpx = (xp - beta)*(xp-beta)*(xp-alpha)*(xp-alpha);
                    Real vmx = (xm - beta)*(xm-beta)*(xm-alpha)*(xm-alpha);
                    Real vpy = (yp - .5)*(yp - .5)*(yp - .5)*(yp - .5);
                    Real vmy = (ym - .5)*(ym - .5)*(ym - .5)*(ym - .5);
                    Real vsubx = (vpx - vmx)/dx[0];
                    Real vsuby = (vpy - vmy)/dx[1];

                    Real sampx,sampy;

                    if(std::abs(vsubx) >= 1.e-12)
                    {
                       sampx = -gamma * std::log(1. - r[0]*(1. - std::exp(-2*vsubx*dx[0]/gamma)))/(2.*vsubx);
                       if( sampx < 0. || sampx > dx[0])
                       {
                          amrex::Print() << " offx " << sampx << " " << r[0] << " " << vsubx << std::endl;
                       }
                     //  if(j == 50){
                     //     amrex::Print() << "center " << i << " " << j << " " << sampx << " " << r[0] << " " << vsubx << std::endl;
                     //  }
                       r[0] = sampx / dx[0];
                    }

                    if(std::abs(vsuby) >= 1.e-12)
                    {
                       sampy = -gamma * std::log(1. - r[1]*(1. - std::exp(-2*vsuby*dx[1]/gamma)))/(2.*vsuby);
                       if( sampy < 0. || sampy > dx[1])
                       {
                          amrex::Print() << "offy " <<sampy << " " << r[1] << " " << vsuby << std::endl;
                       }
                       r[1] = sampy / dx[1];
                    }
                }

                AMREX_D_TERM( Real x = plo[0] + (i + r[0])*dx[0];,
                              Real y = plo[1] + (j + r[1])*dx[1];,
                              Real z = plo[2] + (k + r[2])*dx[2];);

                p.id()  = ip + id_start;
                p.cpu() = my_cpu;

                AMREX_D_TERM( p.pos(0) = x;,
                              p.pos(1) = y;,
                              p.pos(2) = z;);

                AMREX_D_TERM( p.rdata(RealIdx::xold) = x;,
                              p.rdata(RealIdx::yold) = y;,
                              p.rdata(RealIdx::zold) = z;);
            }
        });
    }
   //  amrex::Print() << "SUM / DENS ADDED THROUGH REGRID " << sum[0] << " " << sum[0] / cell_vol << std::endl;
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

        if (ba_to_keep.contains(pti.tilebox())) {continue;}

        if (!m_reflux_particle_locator.isValid(ba_to_keep)) {
            m_reflux_particle_locator.build(ba_to_keep, Geom(lev));
        }
        m_reflux_particle_locator.setGeometry(Geom(lev));

        auto assign_grid = m_reflux_particle_locator.getGridAssignor();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = pstruct[i];
            if (assign_grid(p) < 0) {
                p.id() = -1;
            }
        });
    }
    Redistribute();
}

void
StochasticPC::RefluxFineToCrse (const BoxArray& ba_to_keep, MultiFab& phi_fine_for_reflux)
{
    BL_PROFILE("StochasticPC::RefluxFineToCrse");
    const int lev = 1;
    const auto geom_lev   = Geom(lev);
    const auto dxi_lev    = Geom(lev).InvCellSizeArray();
    const auto plo_lev    = Geom(lev).ProbLoArray();
    const auto domain_lev = Geom(lev).Domain();

    if (!m_reflux_particle_locator.isValid(ba_to_keep)) {
        m_reflux_particle_locator.build(ba_to_keep, Geom(lev));
    }
    m_reflux_particle_locator.setGeometry(Geom(lev));
    auto assign_grid = m_reflux_particle_locator.getGridAssignor();

    for(ParIterType pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& ptile = ParticlesAt(lev, pti);
        auto& aos  = ptile.GetArrayOfStructs();
        const int np = aos.numParticles();
        auto *pstruct = aos().data();

        if (!ba_to_keep.contains(pti.tilebox())) {

            Array4<Real> phi_arr = phi_fine_for_reflux.array(pti.index());

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = pstruct[i];
                auto old_pos = getOldCell(p, plo_lev, dxi_lev, domain_lev);
                auto new_pos = getNewCell(p, plo_lev, dxi_lev, domain_lev);

                // amrex::Print() << "REFLUX:OLD " << static_cast<int>(p.rdata(RealIdx::xold)*dxi_lev[0]) << " " << old_pos[0] << std::endl;
                // amrex::Print() << "REFLUX:NEW " << static_cast<int>(p.pos(0)*dxi_lev[0])               << " " << new_pos[0] << std::endl;

                // int iold = old_pos[0];
                // int inew = new_pos[0];
                // if (inew == 31) amrex::Print() <<" PARTICLE NOW AT 31 : OLD POS " << iold << std::endl;

                if ( (assign_grid(old_pos) >= 0) && (assign_grid(new_pos) < 0)) {
                   Gpu::Atomic::AddNoRet(&phi_arr(new_pos,0), 1.0);
                }
            });
        } // if not in ba_to_keep
    } // pti
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
IntVect getNewCell (StochasticPC::ParticleType const& p,
                    GpuArray<Real,AMREX_SPACEDIM> const& plo,
                    GpuArray<Real,AMREX_SPACEDIM> const& dxi,
                    const Box& domain) noexcept
{
    return getParticleCell(p, plo, dxi, domain);;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
IntVect getOldCell (StochasticPC::ParticleType const& p,
                    GpuArray<Real,AMREX_SPACEDIM> const& plo,
                    GpuArray<Real,AMREX_SPACEDIM> const& dxi,
                    const Box& domain) noexcept
{
    IntVect iv(
               AMREX_D_DECL(int(Math::floor((p.rdata(RealIdx::xold)-plo[0])*dxi[0])),
                            int(Math::floor((p.rdata(RealIdx::yold)-plo[1])*dxi[1])),
                            int(Math::floor((p.rdata(RealIdx::zold)-plo[2])*dxi[2]))));
    iv += domain.smallEnd();
    return iv;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
IntVect periodicCorrectOldCell (const IntVect& old_pos, const IntVect& new_pos,
                                const GpuArray<int,AMREX_SPACEDIM>& is_per,
                                const Box& domain) noexcept
{
    IntVect shifted = old_pos;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (!is_per[idim]) { continue; }
        if (Real(new_pos[idim] - old_pos[idim]) > 0.5*Real(domain.length(idim))) {
            shifted[idim] += domain.length(idim);
            continue;
        }
        if (Real(new_pos[idim] - old_pos[idim]) < -0.5*Real(domain.length(idim))) {
            shifted[idim] -= domain.length(idim);
            continue;
        }
    }
    return shifted;
}

void
StochasticPC::RefluxCrseToFine (const BoxArray& ba_to_keep, MultiFab& phi_for_reflux)
{
    BL_PROFILE("StochasticPC::RefluxCrseToFine");
    const int lev = 1;
    const auto geom_lev   = Geom(lev);
    const auto dxi_lev    = Geom(lev).InvCellSizeArray();
    const auto plo_lev    = Geom(lev).ProbLoArray();
    const auto domain_lev = Geom(lev).Domain();

    if (!m_reflux_particle_locator.isValid(ba_to_keep)) {
        m_reflux_particle_locator.build(ba_to_keep, Geom(lev));
    }
    m_reflux_particle_locator.setGeometry(Geom(lev));

    auto assign_grid = m_reflux_particle_locator.getGridAssignor();

    for (ParIterType pti(*this, lev); pti.isValid(); ++pti)
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
                auto old_pos = getOldCell(p, plo_lev, dxi_lev, domain_lev);
                auto new_pos = getNewCell(p, plo_lev, dxi_lev, domain_lev);

                if ( (assign_grid(old_pos) < 0) && (assign_grid(new_pos) >= 0)) {
                   Gpu::Atomic::AddNoRet(&phi_arr(old_pos,0), -1.0);
                }
            });
        } // if not in ba_to_keep
    } // pti
}

void
StochasticPC::AdvectWithRandomWalk (int lev, Real dt)
{
    BL_PROFILE("StochasticPC::AdvectWithRandomWalk");
    const auto dx = Geom(lev).CellSizeArray();

    const auto p_lo = Geom(lev).ProbLoArray();
    const auto p_hi = Geom(lev).ProbHiArray();

    AMREX_D_TERM( bool is_periodic_in_x = Geom(lev).isPeriodic(0);,
                  bool is_periodic_in_y = Geom(lev).isPeriodic(1);,
                  bool is_periodic_in_z = Geom(lev).isPeriodic(2););

    Real stddev = std::sqrt(dt);

    int ext_pot = 1;
    amrex::Real alpha = .3;
    amrex::Real beta = .7;
    amrex::Real gamma = 5.e-4;

    for(ParIterType pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& ptile = ParticlesAt(lev, pti);
        auto& aos  = ptile.GetArrayOfStructs();
        const int np = aos.numParticles();
        auto *pstruct = aos().data();

        amrex::ParallelForRNG( np, [=] AMREX_GPU_DEVICE (int i, RandomEngine const& engine) noexcept
        {
            ParticleType& p = pstruct[i];
            AMREX_D_TERM( p.rdata(RealIdx::xold) = p.pos(0);,
                          p.rdata(RealIdx::yold) = p.pos(1);,
                          p.rdata(RealIdx::zold) = p.pos(2););

            AMREX_D_TERM( Real incx = amrex::RandomNormal(0.,stddev,engine);,
                          Real incy = amrex::RandomNormal(0.,stddev,engine);,
                          Real incz = amrex::RandomNormal(0.,stddev,engine););

             if(ext_pot ==1){
                amrex::Real xloc,yloc;
                amrex::Real Vsubx, Vsuby;

                xloc = p.pos(0);
                yloc = p.pos(1);
                Vsubx = 2.*(xloc - beta) * (xloc - alpha)* (2.*xloc - alpha - beta) / gamma;
                Vsuby = 4.*(yloc - .5)*(yloc - .5)*(yloc - .5) / gamma;

                p.pos(0) += -dt*Vsubx;
                p.pos(1) += -dt*Vsuby;

             }

            AMREX_D_TERM( incx = std::max(-dx[0], std::min( dx[0], incx));,
                          incy = std::max(-dx[1], std::min( dx[1], incy));,
                          incz = std::max(-dx[2], std::min( dx[2], incz)););

            // HACK TO DO DETERMINISTIC MOVEMENT
            // AMREX_D_TERM( incx = -dx[0];,
            //               incy = +dx[1];,
            //               incz = 0.;);


            AMREX_D_TERM( p.pos(0) += static_cast<ParticleReal> (incx);,
                          p.pos(1) += static_cast<ParticleReal> (incy);,
                          p.pos(2) += static_cast<ParticleReal> (incz););

            if (is_periodic_in_x) {
                if (p.pos(0) < p_lo[0]) p.pos(0) += (p_hi[0] - p_lo[0]);
                if (p.pos(0) > p_hi[0]) p.pos(0) -= (p_hi[0] - p_lo[0]);
            }
            if (is_periodic_in_y) {
                if (p.pos(1) < p_lo[1]) p.pos(1) += (p_hi[1] - p_lo[1]);
                if (p.pos(1) > p_hi[1]) p.pos(1) -= (p_hi[1] - p_lo[1]);
            }

#if (AMREX_SPACEDIM == 3)
            if (is_periodic_in_z) {
                if (p.pos(2) < p_lo[2]) p.pos(2) += (p_hi[2] - p_lo[2]);
                if (p.pos(2) > p_hi[2]) p.pos(2) -= (p_hi[2] - p_lo[2]);
            }
#endif
        }); // np
    } // pti
}
