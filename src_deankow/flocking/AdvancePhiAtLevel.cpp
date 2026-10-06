#include <cmath>

#include <AmrCoreFlock.H>
#include <myfunc.H>
#include "rng_functions.H"

using namespace amrex;

// One step of the flocking SPDE.
//   time_integrator 0 (Euler-Maruyama):
//       phi^{n+1} = phi^n - dt div( F(phi^n) + S(phi^n) )
//   time_integrator 1 (Heun for the deterministic part, noise built once):
//       phi^*     = phi^n - dt div( F(phi^n) + S(phi^n) )
//       phi^{n+1} = phi^n - dt div( (F(phi^n) + F(phi^*))/2 + S(phi^n) )
//   time_integrator 2 (stochastic SSP-RK3 of Delong, Griffith, Vanden-Eijnden
//   and Donev, Phys. Rev. E 87, 033302 (2013)), with R(u, W) = -div(F(u) + S(u, W)):
//       u1        = phi^n + dt R(phi^n, W1)
//       u2        = 3/4 phi^n + 1/4 (u1 + dt R(u1, W2))
//       phi^{n+1} = 1/3 phi^n + 2/3 (u2 + dt R(u2, W3))
//   with W_i = W_A + beta_i W_B for two independent unit-variance fields
//   W_A, W_B drawn once per step, and
//       beta_1 = (2 sqrt2 + sqrt3)/5, beta_2 = (-4 sqrt2 + 3 sqrt3)/5,
//       beta_3 = (sqrt2 - 2 sqrt3)/10,
//   so that W1/6 + W2/6 + 2 W3/3 = W_A. The noise amplitude sqrt(phi) and
//   rho_tilde are evaluated at each stage. Stable for centred advection up to
//   |lambda dt| = sqrt3 on the imaginary axis.
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
    UpdateRhoTilde(phi_old);
    compute_det_fluxes(phi_old, rhot, rhobar, flux, Geom(0), fp, adv_order, limiter);

    MultiFab sflux;
    if (dorand && time_integrator != 2) {
        sflux.define(amrex::convert(ba, IntVect::TheDimensionVector(2)), dm, 1, 0);
        compute_stoch_flux(phi_old, sflux, Geom(0), fp, dt);
    }

    if (time_integrator == 2) {
        // the step above already built F(phi^n); redo the stages uniformly
        const Real s2 = std::sqrt(2.0), s3 = std::sqrt(3.0);
        const Real beta[3] = {(2.0*s2 + s3)/5.0, (-4.0*s2 + 3.0*s3)/5.0, (s2 - 2.0*s3)/10.0};
        const BoxArray zba = amrex::convert(ba, IntVect::TheDimensionVector(2));
        MultiFab wa, wb, ws;
        if (dorand) {
            wa.define(zba, dm, 1, 0);
            wb.define(zba, dm, 1, 0);
            ws.define(zba, dm, 1, 0);
            MultiFabFillRandom(wa, 0, 1.0, Geom(0));
            MultiFabFillRandom(wb, 0, 1.0, Geom(0));
        }
        MultiFab stage(ba, dm, 1, nghost);

        // out = u + dt R(u, W_s); u must have its ghost cells filled
        auto euler_stage = [&] (MultiFab const& u, int s, MultiFab& out)
        {
            if (s > 0) {
                UpdateRhoTilde(u);
                compute_det_fluxes(u, rhot, rhobar, flux, Geom(0), fp, adv_order, limiter);
            }
            if (dorand) {
                MultiFab::LinComb(ws, 1.0, wa, 0, beta[s], wb, 0, 0, 1, 0);
                compute_stoch_flux_given(u, ws, sflux, Geom(0), fp, dt);
                MultiFab::Add(flux[2], sflux, 0, 0, 1, 0);
            }
            apply_fluxes(u, out, flux, Geom(0), dt);
        };

        if (dorand && sflux.empty()) { sflux.define(zba, dm, 1, 0); }

        euler_stage(phi_old, 0, phi_new);                       // u1
        phi_new.FillBoundary(Geom(0).periodicity());
        euler_stage(phi_new, 1, stage);                         // u1 + dt R(u1, W2)
        MultiFab::LinComb(phi_new, 0.75, phi_old, 0, 0.25, stage, 0, 0, 1, 0);   // u2
        phi_new.FillBoundary(Geom(0).periodicity());
        euler_stage(phi_new, 2, stage);                         // u2 + dt R(u2, W3)
        MultiFab::LinComb(phi_new, 1.0/3.0, phi_old, 0, 2.0/3.0, stage, 0, 0, 1, 0);
    } else if (time_integrator == 0) {
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
        UpdateRhoTilde(phi_new);
        compute_det_fluxes(phi_new, rhot, rhobar, flux1, Geom(0), fp, adv_order, limiter);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            MultiFab::LinComb(flux[d], 0.5, flux[d], 0, 0.5, flux1[d], 0, 0, 1, 0);
        }
        if (dorand) { MultiFab::Add(flux[2], sflux, 0, 0, 1, 0); }
        apply_fluxes(phi_old, phi_new, flux, Geom(0), dt);
    }

    phi_new.FillBoundary(Geom(0).periodicity());
}
