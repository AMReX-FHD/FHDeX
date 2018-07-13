#ifndef _hydro_functions_H
#define _hydro_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

void advance_timestep_inertial();

void MkAdvMFluxdiv(const std::array<MultiFab, AMREX_SPACEDIM>& umac,
		   const std::array<MultiFab, AMREX_SPACEDIM>& m,
		   std::array<MultiFab, AMREX_SPACEDIM>& m_update,
		   const amrex::Real* dx,
		   const bool& increment);

#endif
