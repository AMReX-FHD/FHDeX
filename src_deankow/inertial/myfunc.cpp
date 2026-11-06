#include "myfunc.H"
#include "mykernel.H"
#include "rng_functions.H"

using namespace amrex;

void fill_phi_ghosts (MultiFab& phi, Geometry const& geom)
{
    phi.FillBoundary(geom.periodicity());
    const Box& domain = geom.Domain();
    const int jlo = domain.smallEnd(1);
    const int jhi = domain.bigEnd(1);
    const int ng = phi.nGrow();
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        auto const& p = phi.array(mfi);
        Box gbx = mfi.growntilebox(ng);
        // even reflection: phi(jlo - 1 - m) = phi(jlo + m), phi(jhi + 1 + m) = phi(jhi - m)
        Box lo = gbx; lo.setBig(1, jlo - 1);
        Box hi = gbx; hi.setSmall(1, jhi + 1);
        if (lo.ok()) {
            amrex::ParallelFor(lo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                p(i,j,k) = p(i, 2*jlo - 1 - j, k);
            });
        }
        if (hi.ok()) {
            amrex::ParallelFor(hi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                p(i,j,k) = p(i, 2*jhi + 1 - j, k);
            });
        }
    }
}

void compute_det_fluxes (MultiFab const& phi, Array<MultiFab, AMREX_SPACEDIM>& flux,
                         Geometry const& geom, InertialParams const& ip,
                         int adv_order, int limiter)
{
    const auto dx  = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const Real dvinv = geom.InvCellSize(1);
    const int jlo = geom.Domain().smallEnd(1);
    const int jhi = geom.Domain().bigEnd(1);
    const InertialParams ip_loc = ip;

    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& p  = phi.const_array(mfi);
        auto const& fx = flux[0].array(mfi);
        auto const& fv = flux[1].array(mfi);
        const Box& xbx = mfi.nodaltilebox(0);
        const Box& vbx = mfi.nodaltilebox(1);

        amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            compute_flux_x(i, j, k, fx, p, plo[1] + (j+0.5)*dx[1], adv_order, limiter);
        });
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real x = plo[0] + (i+0.5)*dx[0];
            compute_flux_v(i, j, k, fv, p, plo[1] + j*dx[1], inertial_dU(x, ip_loc), dvinv,
                           jlo, jhi, ip_loc, adv_order, limiter);
        });
    }
}

void compute_stoch_flux (MultiFab const& phi, MultiFab& sflux_v, Geometry const& geom,
                         InertialParams const& ip, Real dt)
{
    const auto dx = geom.CellSizeArray();
    const Real variance = 2.0*ip.diff_coeff/(ip.num_part*dx[0]*dx[1]*dt);
    MultiFabFillRandom(sflux_v, 0, variance, geom);
    const int jlo = geom.Domain().smallEnd(1);
    const int jhi = geom.Domain().bigEnd(1);

    for (MFIter mfi(sflux_v, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& p  = phi.const_array(mfi);
        auto const& sv = sflux_v.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            scale_stoch_flux_v(i, j, k, sv, p, jlo, jhi);
        });
    }
}

void compute_stoch_flux_given (MultiFab const& phi, MultiFab const& w, MultiFab& sflux_v,
                               Geometry const& geom, InertialParams const& ip, Real dt)
{
    const auto dx = geom.CellSizeArray();
    const Real sd = std::sqrt(2.0*ip.diff_coeff/(ip.num_part*dx[0]*dx[1]*dt));
    MultiFab::Copy(sflux_v, w, 0, 0, 1, 0);
    sflux_v.mult(sd, 0, 1);
    const int jlo = geom.Domain().smallEnd(1);
    const int jhi = geom.Domain().bigEnd(1);

    for (MFIter mfi(sflux_v, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& p  = phi.const_array(mfi);
        auto const& sv = sflux_v.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            scale_stoch_flux_v(i, j, k, sv, p, jlo, jhi);
        });
    }
}

void apply_fluxes (MultiFab const& phi_old, MultiFab& phi_new,
                   Array<MultiFab, AMREX_SPACEDIM> const& flux, Geometry const& geom, Real dt)
{
    const auto dxinv = geom.InvCellSizeArray();
    for (MFIter mfi(phi_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& po = phi_old.const_array(mfi);
        auto const& pn = phi_new.array(mfi);
        auto const& fx = flux[0].const_array(mfi);
        auto const& fv = flux[1].const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            update_phi(i, j, k, po, pn, fx, fv, dt, dxinv);
        });
    }
}
