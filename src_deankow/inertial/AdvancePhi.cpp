#include <cmath>

#include <AmrCoreInertial.H>
#include <myfunc.H>
#include "rng_functions.H"

using namespace amrex;

// One step of the inertial SPDE on (x, v), as the active-matter code:
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
//   drawn once per step, so that W1/6 + W2/6 + 2 W3/3 = W_A. Stable for
//   centred advection up to |lambda dt| = sqrt3 on the imaginary axis.
// F holds the deterministic fluxes (x free streaming; v friction, force and
// diffusion), S the stochastic v flux.
//
// If fsum is given (the hybrid), it receives the time-averaged face flux that
// the step applied, phi^{n+1} = phi^n - dt div(fsum), noise included:
//   Euler-Maruyama: F(phi^n) + S;  Heun: (F(phi^n) + F(phi^*))/2 + S;
//   RK3: F_1/6 + F_2/6 + 2 F_3/3 with each stage's own noise.
// The hybrid uses it to remove the SPDE flux across the particle-region
// boundary.
void
AmrCoreInertial::AdvancePhi (Real dt, Array<MultiFab, AMREX_SPACEDIM>* fsum)
{
    std::swap(phi_old, phi_new);
    fill_phi_ghosts(phi_old, Geom(0));

    const BoxArray& ba = phi_old.boxArray();
    const DistributionMapping& dm = phi_old.DistributionMap();
    const BoxArray vba = amrex::convert(ba, IntVect::TheDimensionVector(1));

    Array<MultiFab, AMREX_SPACEDIM> flux;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        flux[d].define(amrex::convert(ba, IntVect::TheDimensionVector(d)), dm, 1, 0);
    }

    if (time_integrator == 2) {
        const Real s2 = std::sqrt(2.0), s3 = std::sqrt(3.0);
        const Real beta[3] = {(2.0*s2 + s3)/5.0, (-4.0*s2 + 3.0*s3)/5.0, (s2 - 2.0*s3)/10.0};
        MultiFab wa, wb, ws, sflux;
        if (dorand) {
            wa.define(vba, dm, 1, 0);
            wb.define(vba, dm, 1, 0);
            ws.define(vba, dm, 1, 0);
            sflux.define(vba, dm, 1, 0);
            MultiFabFillRandom(wa, 0, 1.0, Geom(0));
            MultiFabFillRandom(wb, 0, 1.0, Geom(0));
        }
        MultiFab stage(ba, dm, 1, nghost);

        const Real wstage[3] = {1.0/6.0, 1.0/6.0, 2.0/3.0};

        // out = u + dt R(u, W_s); u must have its ghost cells filled
        auto euler_stage = [&] (MultiFab const& u, int s, MultiFab& out)
        {
            compute_det_fluxes(u, flux, Geom(0), ip, adv_order, limiter);
            if (dorand) {
                MultiFab::LinComb(ws, 1.0, wa, 0, beta[s], wb, 0, 0, 1, 0);
                compute_stoch_flux_given(u, ws, sflux, Geom(0), ip, dt);
                MultiFab::Add(flux[1], sflux, 0, 0, 1, 0);
            }
            if (fsum) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    MultiFab::Saxpy((*fsum)[d], wstage[s], flux[d], 0, 0, 1, 0);
                }
            }
            apply_fluxes(u, out, flux, Geom(0), dt);
        };

        euler_stage(phi_old, 0, phi_new);                       // u1
        fill_phi_ghosts(phi_new, Geom(0));
        euler_stage(phi_new, 1, stage);                         // u1 + dt R(u1, W2)
        MultiFab::LinComb(phi_new, 0.75, phi_old, 0, 0.25, stage, 0, 0, 1, 0);   // u2
        fill_phi_ghosts(phi_new, Geom(0));
        euler_stage(phi_new, 2, stage);                         // u2 + dt R(u2, W3)
        MultiFab::LinComb(phi_new, 1.0/3.0, phi_old, 0, 2.0/3.0, stage, 0, 0, 1, 0);
    } else {
        compute_det_fluxes(phi_old, flux, Geom(0), ip, adv_order, limiter);
        MultiFab sflux;
        if (dorand) {
            sflux.define(vba, dm, 1, 0);
            compute_stoch_flux(phi_old, sflux, Geom(0), ip, dt);
        }
        if (time_integrator == 0) {
            if (dorand) { MultiFab::Add(flux[1], sflux, 0, 0, 1, 0); }
            if (fsum) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { MultiFab::Copy((*fsum)[d], flux[d], 0, 0, 1, 0); }
            }
            apply_fluxes(phi_old, phi_new, flux, Geom(0), dt);
        } else {
            // predictor
            Array<MultiFab, AMREX_SPACEDIM> flux1;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                flux1[d].define(flux[d].boxArray(), dm, 1, 0);
                MultiFab::Copy(flux1[d], flux[d], 0, 0, 1, 0);
            }
            if (dorand) { MultiFab::Add(flux1[1], sflux, 0, 0, 1, 0); }
            apply_fluxes(phi_old, phi_new, flux1, Geom(0), dt);
            fill_phi_ghosts(phi_new, Geom(0));

            // corrector: average the deterministic fluxes, reuse the noise
            compute_det_fluxes(phi_new, flux1, Geom(0), ip, adv_order, limiter);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                MultiFab::LinComb(flux[d], 0.5, flux[d], 0, 0.5, flux1[d], 0, 0, 1, 0);
            }
            if (dorand) { MultiFab::Add(flux[1], sflux, 0, 0, 1, 0); }
            if (fsum) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { MultiFab::Copy((*fsum)[d], flux[d], 0, 0, 1, 0); }
            }
            apply_fluxes(phi_old, phi_new, flux, Geom(0), dt);
        }
    }

    fill_phi_ghosts(phi_new, Geom(0));
}
