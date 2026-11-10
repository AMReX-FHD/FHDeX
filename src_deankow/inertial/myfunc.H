#ifndef INERTIAL_MYFUNC_H_
#define INERTIAL_MYFUNC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include "InertialParams.H"

// Fill the ghost cells of phi: periodic in x (FillBoundary), then even
// reflection across the v walls (only the MUSCL stencils read them; the wall
// fluxes are zero).
void fill_phi_ghosts (amrex::MultiFab& phi, amrex::Geometry const& geom);

// Deterministic fluxes (x free streaming; v friction, force and diffusion)
// from phi, whose ghost cells must be filled.
void compute_det_fluxes (amrex::MultiFab const& phi,
                         amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& flux,
                         amrex::Geometry const& geom, InertialParams const& ip,
                         int adv_order, int limiter);

// Stochastic v flux from phi: fills sflux_v with N(0, 2 D/(N dV dt)) and
// multiplies by sqrt(phi_face); zero on the walls.
void compute_stoch_flux (amrex::MultiFab const& phi, amrex::MultiFab& sflux_v,
                         amrex::Geometry const& geom, InertialParams const& ip,
                         amrex::Real dt);

// Stochastic v flux for a given unit-variance noise field w on the v faces:
// sflux_v = sqrt(2 D/(N dV dt)) w sqrt(phi_face).
void compute_stoch_flux_given (amrex::MultiFab const& phi, amrex::MultiFab const& w,
                               amrex::MultiFab& sflux_v, amrex::Geometry const& geom,
                               InertialParams const& ip, amrex::Real dt);

// Hybrid: remove the SPDE flux across the boundary of the particle region.
// For each cell outside the region (mask == 0) with a neighbour inside
// (mask == 1, one ghost cell filled periodically), undo the contribution of
// that face, -dt/dx (F_hi - F_lo), using the time-averaged face flux fsum of
// the step. The particles then carry all the mass across that boundary.
void undo_cf_fluxes (amrex::MultiFab& phi, amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> const& fsum,
                     amrex::iMultiFab const& mask, amrex::Geometry const& geom, amrex::Real dt);

// phi_new = phi_old - dt * div(flux)
void apply_fluxes (amrex::MultiFab const& phi_old, amrex::MultiFab& phi_new,
                   amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> const& flux,
                   amrex::Geometry const& geom, amrex::Real dt);

#endif
