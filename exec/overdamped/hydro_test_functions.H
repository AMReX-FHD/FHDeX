#ifndef _hydro_test_functions_H_
#define _hyrdo_test_functions_H_

#include "common_functions.H"
#include "gmres_functions.H"

///////////////////////////
// in advance.cpp
void advance(std::array< amrex::MultiFab, AMREX_SPACEDIM >& umac,
	     amrex::MultiFab& pres,
	     const std::array< amrex::MultiFab, AMREX_SPACEDIM >& stochMfluxdiv,
	     std::array< amrex::MultiFab, AMREX_SPACEDIM >& alpha_fc,
	     amrex::MultiFab& beta,
             amrex::MultiFab& gamma,
	     std::array< amrex::MultiFab, NUM_EDGE >& beta_ed,
	     const amrex::Geometry geom, const amrex::Real& dt);

///////////////////////////

///////////////////////////
// in WritePlotFile.cpp
void WritePlotFile(const int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   std::array< MultiFab, AMREX_SPACEDIM >& umac,
                   const MultiFab& pres);

///////////////////////////

#endif
