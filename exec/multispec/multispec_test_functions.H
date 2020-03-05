#ifndef _multispec_test_functions_H_
#define _multispec_test_functions_H_

#include "common_functions.H"
#include "gmres_functions.H"

///////////////////////////
// in advance_timestep_inertial.cpp

void advance_timestep_inertial(std::array< amrex::MultiFab, AMREX_SPACEDIM >& umac,
                               amrex::MultiFab& pres,
                               amrex::MultiFab& rho, amrex::MultiFab& rhotot,
                               const std::array< amrex::MultiFab, AMREX_SPACEDIM >& mfluxdiv_predict,
                               const std::array< amrex::MultiFab, AMREX_SPACEDIM >& mfluxdiv_correct,
                               std::array< amrex::MultiFab, AMREX_SPACEDIM >& alpha_fc,
                               const amrex::MultiFab& beta, const amrex::MultiFab& gamma,
                               const std::array< amrex::MultiFab, NUM_EDGE >& beta_ed,
                               const amrex::Geometry geom, const amrex::Real& dt);

///////////////////////////

///////////////////////////
// in WritePlotFile.cpp

void WritePlotFile(const int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   std::array< MultiFab, AMREX_SPACEDIM >& umac,
		   const MultiFab& rho,
                   const MultiFab& pres);

///////////////////////////

#endif
