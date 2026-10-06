#include <AmrCoreFlock.H>
#include <myfunc.H>

using namespace amrex;

// One step of the flocking SPDE.
//   time_integrator 0 (Euler-Maruyama):
//       phi^{n+1} = phi^n - dt div( F(phi^n) + S(phi^n) )
//   time_integrator 1 (Heun for the deterministic part, noise built once):
//       phi^*     = phi^n - dt div( F(phi^n) + S(phi^n) )
//       phi^{n+1} = phi^n - dt div( (F(phi^n) + F(phi^*))/2 + S(phi^n) )
// F holds the deterministic fluxes (x, y advection; theta diffusion and
// turning), S the stochastic theta flux.
void
AmrCoreFlock::AdvancePhi (Real dt)
{
    std::swap(phi_old, phi_new);
    phi_old.FillBoundary(Geom(0).periodicity());

    const BoxArray& ba = phi_old.boxArray();
    const DistributionMapping& dm = phi_old.DistributionMap();

    Array<MultiFab, AMREX_SPACEDIM> flux;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        flux[d].define(amrex::convert(ba, IntVect::TheDimensionVector(d)), dm, 1, 0);
    }
    compute_det_fluxes(phi_old, flux, Geom(0), fp, adv_order, limiter);

    MultiFab sflux;
    if (dorand) {
        sflux.define(amrex::convert(ba, IntVect::TheDimensionVector(2)), dm, 1, 0);
        compute_stoch_flux(phi_old, sflux, Geom(0), fp, dt);
    }

    if (time_integrator == 0) {
        if (dorand) { MultiFab::Add(flux[2], sflux, 0, 0, 1, 0); }
        apply_fluxes(phi_old, phi_new, flux, Geom(0), dt);
    } else {
        // predictor
        Array<MultiFab, AMREX_SPACEDIM> flux1;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            flux1[d].define(flux[d].boxArray(), dm, 1, 0);
            MultiFab::Copy(flux1[d], flux[d], 0, 0, 1, 0);
        }
        if (dorand) { MultiFab::Add(flux1[2], sflux, 0, 0, 1, 0); }
        apply_fluxes(phi_old, phi_new, flux1, Geom(0), dt);
        phi_new.FillBoundary(Geom(0).periodicity());

        // corrector: average the deterministic fluxes, reuse the noise
        compute_det_fluxes(phi_new, flux1, Geom(0), fp, adv_order, limiter);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            MultiFab::LinComb(flux[d], 0.5, flux[d], 0, 0.5, flux1[d], 0, 0, 1, 0);
        }
        if (dorand) { MultiFab::Add(flux[2], sflux, 0, 0, 1, 0); }
        apply_fluxes(phi_old, phi_new, flux, Geom(0), dt);
    }

    phi_new.FillBoundary(Geom(0).periodicity());
}
