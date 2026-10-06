#ifndef FLOCK_MYFUNC_H_
#define FLOCK_MYFUNC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include "FlockingParams.H"

// Deterministic fluxes (x, y advection; theta diffusion and turning) from phi,
// whose ghost cells must be filled. rhot is the smoothed spatial density
// rho_tilde (one ghost cell) that sets the speed; rhobar is its mean.
void compute_det_fluxes (amrex::MultiFab const& phi, amrex::MultiFab const& rhot,
                         amrex::Real rhobar,
                         amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& flux,
                         amrex::Geometry const& geom, FlockingParams const& fp,
                         int adv_order, int limiter);

// Stochastic theta flux from phi: fills sflux_z with N(0, 2 D/(N dV dt)) and
// multiplies by sqrt(phi_face).
void compute_stoch_flux (amrex::MultiFab const& phi, amrex::MultiFab& sflux_z,
                         amrex::Geometry const& geom, FlockingParams const& fp,
                         amrex::Real dt);

// Stochastic theta flux from phi for a given unit-variance noise field w on
// the theta faces: sflux_z = sqrt(2 D/(N dV dt)) w sqrt(phi_face).
void compute_stoch_flux_given (amrex::MultiFab const& phi, amrex::MultiFab const& w,
                               amrex::MultiFab& sflux_z, amrex::Geometry const& geom,
                               FlockingParams const& fp, amrex::Real dt);

// phi_new = phi_old - dt * div(flux)
void apply_fluxes (amrex::MultiFab const& phi_old, amrex::MultiFab& phi_new,
                 amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> const& flux,
                 amrex::Geometry const& geom, amrex::Real dt);

#endif
