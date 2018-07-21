#ifndef _hydro_functions_H
#define _hydro_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

void advance_timestep_inertial();

void MkAdvMFluxdiv(const std::array<amrex::MultiFab, AMREX_SPACEDIM>& umac,
		   const std::array<amrex::MultiFab, AMREX_SPACEDIM>& m,
		   std::array<amrex::MultiFab, AMREX_SPACEDIM>& m_update,
		   const amrex::Real* dx,
		   const int& increment=0);

void MkAdvSFluxdiv(const std::array<amrex::MultiFab, AMREX_SPACEDIM>& umac,
		   const amrex::MultiFab& m,
		   amrex::MultiFab& m_update,
		   const amrex::Real* dx,
		   const amrex::Geometry& geom,
		   const int& increment=0);

#endif
