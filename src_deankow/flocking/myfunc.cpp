#include "myfunc.H"
#include "mykernel.H"
#include "rng_functions.H"

using namespace amrex;

void compute_det_fluxes (MultiFab const& phi, MultiFab const& rhot, Real rhobar,
                         Array<MultiFab, AMREX_SPACEDIM>& flux,
                         Geometry const& geom, FlockingParams const& fp,
                         int adv_order, int limiter)
{
    const Real dthinv = geom.InvCellSize(2);
    const Real dth    = geom.CellSize(2);
    const Real thlo   = geom.ProbLo(2);
    // turning rate from the alignment interaction; 0 until it is implemented
    const Real omega  = 0.0;
    const Real rhobar_inv = 1.0/rhobar;

    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& p  = phi.const_array(mfi);
        auto const& rt = rhot.const_array(mfi);
        auto const& fx = flux[0].array(mfi);
        auto const& fy = flux[1].array(mfi);
        auto const& fz = flux[2].array(mfi);
        const Box& xbx = mfi.nodaltilebox(0);
        const Box& ybx = mfi.nodaltilebox(1);
        const Box& zbx = mfi.nodaltilebox(2);

        amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            compute_flux_x(i, j, k, fx, p, rt, rhobar_inv, thlo + (k+0.5)*dth, fp, adv_order, limiter);
        });
        amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            compute_flux_y(i, j, k, fy, p, rt, rhobar_inv, thlo + (k+0.5)*dth, fp, adv_order, limiter);
        });
        amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            compute_flux_theta(i, j, k, fz, p, dthinv, fp, omega, adv_order, limiter);
        });
    }
}

void compute_stoch_flux (MultiFab const& phi, MultiFab& sflux_z, Geometry const& geom,
                         FlockingParams const& fp, Real dt)
{
    const auto dx = geom.CellSizeArray();
    const Real dv = dx[0]*dx[1]*dx[2];
    const Real variance = 2.0*fp.diff_coeff/(fp.num_part*dv*dt);
    MultiFabFillRandom(sflux_z, 0, variance, geom);

    for (MFIter mfi(sflux_z, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& p  = phi.const_array(mfi);
        auto const& sz = sflux_z.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            scale_stoch_flux_theta(i, j, k, sz, p);
        });
    }
}

void compute_stoch_flux_given (MultiFab const& phi, MultiFab const& w, MultiFab& sflux_z,
                               Geometry const& geom, FlockingParams const& fp, Real dt)
{
    const auto dx = geom.CellSizeArray();
    const Real dv = dx[0]*dx[1]*dx[2];
    const Real sd = std::sqrt(2.0*fp.diff_coeff/(fp.num_part*dv*dt));
    MultiFab::Copy(sflux_z, w, 0, 0, 1, 0);
    sflux_z.mult(sd, 0, 1);

    for (MFIter mfi(sflux_z, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& p  = phi.const_array(mfi);
        auto const& sz = sflux_z.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            scale_stoch_flux_theta(i, j, k, sz, p);
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
        auto const& fy = flux[1].const_array(mfi);
        auto const& fz = flux[2].const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            update_phi(i, j, k, po, pn, fx, fy, fz, dt, dxinv);
        });
    }
}
