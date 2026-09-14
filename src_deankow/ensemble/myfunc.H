#ifndef MYFUNC_H_
#define MYFUNC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

using namespace amrex;

void advance_phi (amrex::MultiFab& phi_old,
                  amrex::MultiFab& phi_new,
                  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& flux,
                  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& stochFlux,
                  amrex::Real dt, amrex::Real npts_scale,
                  amrex::Geometry const& geom,
                  Vector<amrex::BCRec> const& BoundaryCondition,
                  Vector<int> const& a_ensemble_dir,
                  int a_ext_pot, Real a_alpha, Real a_beta,
                  Real a_gamma);

#endif
