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
                  MultiFab& gmetric,
                  MultiFab& sqrgmetric,
                  MultiFab& detg,
                  MultiFab& newgmetric,
                  MultiFab& newsqrgmetric,
                  MultiFab& newdetg,
                  amrex::Real dt, amrex::Real num_part, amrex::Real dorand,
                  int num_flux,int ext_pot, amrex::Geometry const& geom,
                  Vector<amrex::BCRec> const& BoundaryCondition, amrex::Real time);
#endif
