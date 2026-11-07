#include "StochasticPC.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_Math.H>
#include <AMReX_ParmParse.H>
#include <AMReX_TracerParticle_mod_K.H>

using namespace amrex;

void
StochasticPC::InitParticles (MultiFab& phi_fine, PotParams const& pot)
{
    AddParticles(phi_fine, BoxArray{}, pot);
}

void
StochasticPC::ColorParticlesWithPhi (MultiFab const& phi)
{
    BL_PROFILE("StochasticPC::ColorParticlesWithPhi");
    const int lev = 1;
    const auto dxi    = Geom(lev).InvCellSizeArray();
    const auto plo    = Geom(lev).ProbLoArray();
    const auto domain = Geom(lev).Domain();

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
            // Use the same helper getNewCell/getOldCell use: the raw divide
            // omitted prob_lo and domain.smallEnd(), and static_cast truncates
            // toward zero instead of flooring, so particles just outside the
            // low side of the domain -- which AddParticles deliberately allows
            // -- collapsed onto cell 0.
            const amrex::IntVect iv = amrex::getParticleCell(p, plo, dxi, domain);
            p.rdata(RealIdx::zold) = phi_arr(iv,0);
        });
    }
}

void
StochasticPC:: AddParticles (MultiFab& phi_fine, const BoxArray& ba_to_exclude, PotParams const& pot)
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
            if (assign_grid(IntVect(AMREX_D_DECL(i, j, k))).first >= 0) {return;}
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

        const int  tilt  = (pot.ext_pot == 1 || pot.sampler_only == 1) ? 1 : 0;
        const Real gamma = pot.gamma;

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
                if (tilt == 1)
                {
                    // Boltzmann factor exp(-2 V/gamma) linearized across the cell with the
                    // cell-averaged slope dV = (V(x+) - V(x-))/dx; sample s in [0, dx] from
                    // the density ~ exp(-2 dV s/gamma) by inverting its CDF.
                    const Real xm = plo[0] + i*dx[0];
                    const Real ym = plo[1] + j*dx[1];
                    const Real vsubx = (pot.Vx(xm + dx[0]) - pot.Vx(xm))/dx[0];
                    const Real vsuby = (pot.Vy(ym + dx[1]) - pot.Vy(ym))/dx[1];

                    if (std::abs(vsubx) >= 1.e-12)
                    {
                       Real sampx = -gamma * std::log(1. - r[0]*(1. - std::exp(-2*vsubx*dx[0]/gamma)))/(2.*vsubx);
                       sampx = amrex::min(amrex::max(sampx, Real(0.)), dx[0]*(1. - 1.e-12));
                       r[0] = sampx / dx[0];
                    }

                    if (std::abs(vsuby) >= 1.e-12)
                    {
                       Real sampy = -gamma * std::log(1. - r[1]*(1. - std::exp(-2*vsuby*dx[1]/gamma)))/(2.*vsuby);
                       sampy = amrex::min(amrex::max(sampy, Real(0.)), dx[1]*(1. - 1.e-12));
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
            if (assign_grid(p).first < 0) {
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

                if ( (assign_grid(old_pos).first >= 0) && (assign_grid(new_pos).first < 0)) {
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

    const auto is_per = Geom(lev).isPeriodicArray();

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

                if ( (assign_grid(old_pos).first < 0) && (assign_grid(new_pos).first >= 0)) {
                   // The deposit belongs in the cell the particle came from.  That cell is
                   // normally inside this FAB, but a particle that crossed a periodic
                   // boundary this step has an old cell on the far side of the domain,
                   // because AdvectWithRandomWalk wrapped its position.  Undo the wrap so
                   // the deposit lands next to the new cell instead of out of bounds.
                   if (Box(phi_arr).contains(old_pos)) {
                      Gpu::Atomic::AddNoRet(&phi_arr(old_pos,0), -1.0);
                   } else {
                      auto shifted_pos = periodicCorrectOldCell(old_pos, new_pos,
                                                                is_per, domain_lev);
                      AMREX_ASSERT(Box(phi_arr).contains(shifted_pos));
                      Gpu::Atomic::AddNoRet(&phi_arr(shifted_pos,0), -1.0);
                   }
                }
            });
        } // if not in ba_to_keep
    } // pti
}

void
StochasticPC::AdvectWithRandomWalk (int lev, Real dt, PotParams const& pot)
{
    BL_PROFILE("StochasticPC::AdvectWithRandomWalk");
    const auto dx = Geom(lev).CellSizeArray();

    const auto p_lo = Geom(lev).ProbLoArray();
    const auto p_hi = Geom(lev).ProbHiArray();

    AMREX_D_TERM( bool is_periodic_in_x = Geom(lev).isPeriodic(0);,
                  bool is_periodic_in_y = Geom(lev).isPeriodic(1);,
                  bool is_periodic_in_z = Geom(lev).isPeriodic(2););

    Real stddev = std::sqrt(dt);

    const int ext_pot = pot.ext_pot;
    // amr.on_surf = 1 applies the curved-surface drift and anisotropic diffusion below;
    // the default 0 matches the flat, isotropic level-0 SPDE in mykernel.H
    int on_surf = 0;
    {
        ParmParse pp("amr");
        pp.query("on_surf", on_surf);
    }

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

             if (ext_pot == 1) {
                // drift -grad V / gamma, folded into the increment so the total move
                // is clipped to one cell below like the random part
                incx += -dt*pot.dVx(p.pos(0))/pot.gamma;
                incy += -dt*pot.dVy(p.pos(1))/pot.gamma;
             }

             if(on_surf == 1){

                 amrex::Real amp = 0.1;
                 amrex::Real sinx,siny,cosx,cosy,det,fx,fy,sig11,sig12,sig21,sig22;

                 amrex::Real xloc = p.pos(0);
                 amrex::Real yloc = p.pos(1);

                 sinx = std::sin(xloc);
                 siny = std::sin(yloc);
                 cosx = std::cos(xloc);
                 cosy = std::cos(yloc);

                 det = 1 + amp*amp*(cosx*cosx*siny*siny+sinx*sinx*cosy*cosy);

                 fx = amp*amp*(2. + amp*amp * (cosx*cosx + cosy*cosy)*sinx*cosx*siny*siny);
                 fy = amp*amp*(2. + amp*amp * (cosx*cosx + cosy*cosy)*siny*cosy*sinx*sinx);

                 fx = fx / (2.*det*det);
                 fy = fy / (2.*det*det);

                 amrex::Real detm1 = det - 1.;
                 amrex::Real cfac = (1. - 1./std::sqrt(1+detm1))/detm1;

                 sig11 = 1. - cfac * amp*amp * cosx*cosx*siny*siny;
                 sig22 = 1. - cfac * amp*amp * cosy*cosy*sinx*sinx;
                 sig12 =  - cfac * amp*amp * cosy*cosx*sinx*siny;
                 sig21 = sig12;

                 amrex::Real updatex = fx*dt + sig11*incx + sig12*incy;
                 amrex::Real updatey = fy*dt + sig21*incx + sig22*incy;

#ifndef AMREX_USE_GPU
                 if(std::abs(updatex) > dx[0] || std::abs(updatey) > dx[1])
                 {
                    amrex::Print{} << "at " << xloc << " " << yloc << " step " << updatex << " " << updatey << " with inc " << incx << " " << incy << " mesh " << dx[0] << " " << dx[1] << std::endl;
                 }
#endif

                 updatex = std::max(-dx[0], std::min( dx[0], updatex));
                 updatey = std::max(-dx[1], std::min( dx[1], updatey));

                 // amrex::Print{} << "at " << xloc << " " << yloc << " step " << updatex << " " << updatey << " with inc " << incx << " " << incy << std::endl;

                 p.pos(0) += static_cast<ParticleReal> (updatex);
                 p.pos(1) += static_cast<ParticleReal> (updatey);


             } else {

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

            }

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
