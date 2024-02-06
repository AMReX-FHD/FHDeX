#include "StochasticPC.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_Math.H>
#include <AMReX_TracerParticle_mod_K.H>

using namespace amrex;

void
StochasticPC::InitParticles (MultiFab& phi_fine)
{
    AddParticles(phi_fine, BoxArray{});
}

void
StochasticPC:: AddParticles (MultiFab& phi_fine, const BoxArray& ba_to_exclude)
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

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi_fine); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();

        if (ba_to_exclude.contains(tile_box)) {continue;}

        const Array4<Real const>& phi_arr = phi_fine.const_array(mfi);

        // count the number of particles to create in each cell
        auto flat_index = FlatIndex(tile_box);
        Gpu::DeviceVector<unsigned int> counts(tile_box.numPts()+1, 0);
        unsigned int* pcount = counts.dataPtr();
        amrex::ParallelForRNG(tile_box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
        {
            Real rannum = amrex::Random(engine);
            int npart_in_cell = int(phi_arr(i,j,k,0)*cell_vol+rannum);
            pcount[flat_index(i, j, k)] += npart_in_cell;
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
        amrex::Print() << "INIT: NEW SIZE OF PARTICLES " << new_size << std::endl;

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

            Array4<Real> phi_arr = phi_for_reflux.array(pti.index());

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = pstruct[i];
                auto old_pos = getOldCell(p, plo_lev, dxi_lev, domain_lev);
                auto new_pos = getNewCell(p, plo_lev, dxi_lev, domain_lev);
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
    const auto is_per     = Geom(lev).isPeriodicArray();

    if (!m_reflux_particle_locator.isValid(ba_to_keep)) {
        m_reflux_particle_locator.build(ba_to_keep, Geom(lev));
    }
    m_reflux_particle_locator.setGeometry(Geom(lev));
    auto assign_grid = m_reflux_particle_locator.getGridAssignor();

    for(ParIterType pti(*this, lev); pti.isValid(); ++pti)
    {
        auto const& gid   = pti.index();
        auto& ptile = ParticlesAt(lev, pti);
        auto& aos  = ptile.GetArrayOfStructs();
        const int np = aos.numParticles();
        auto *pstruct = aos().data();

        if (ba_to_keep.contains(pti.tilebox()))
        {
            Array4<Real> phi_arr = phi_for_reflux.array(gid);

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = pstruct[i];

                auto old_pos = getOldCell(p, plo_lev, dxi_lev, domain_lev);
                auto new_pos = getNewCell(p, plo_lev, dxi_lev, domain_lev);

                // Make a box of the cell holding the particle in its previous position
                Box bx(old_pos,old_pos);

                if ( (assign_grid(old_pos) < 0) && (assign_grid(new_pos) >= 0))
                {
                    if (Box(phi_arr).contains(old_pos)) {
                        Gpu::Atomic::AddNoRet(&phi_arr(old_pos,0), -1.0);
                    } else {
                        auto shifted_pos = periodicCorrectOldCell(old_pos, new_pos,
                                                                  is_per, domain_lev);
                        Gpu::Atomic::AddNoRet(&phi_arr(shifted_pos,0), -1.0);
                    } // else
                } // if crossed the coarse-fine boundary
            }); // i
        } // if in ba_to_keep
    } // pti
}

void
StochasticPC::AdvectWithRandomWalk (int lev, Real dt)
{
    BL_PROFILE("StochasticPC::AdvectWithRandomWalk");
    const auto dx = Geom(lev).CellSizeArray();

    Real stddev = std::sqrt(dt);

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

            AMREX_D_TERM( incx = std::max(-dx[0], std::min( dx[0], incx));,
                          incy = std::max(-dx[1], std::min( dx[1], incy));,
                          incz = std::max(-dx[2], std::min( dx[2], incz)););

            AMREX_D_TERM( p.pos(0) += static_cast<ParticleReal> (incx);,
                          p.pos(1) += static_cast<ParticleReal> (incy);,
                          p.pos(2) += static_cast<ParticleReal> (incz););
        }); // np
    } // pti
}
