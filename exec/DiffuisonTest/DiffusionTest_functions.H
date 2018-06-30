#ifndef _DiffusionTest_functions_H_
#define _DiffusionTest_functions_H_

#include "common_functions.H"
#include "gmres_functions.H"

///////////////////////////
// in WritePlotFile.cpp

void WritePlotFile(const int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   const std::array< MultiFab, AMREX_SPACEDIM >& umac);

void StagExpSolver(const std::array< amrex::MultiFab, AMREX_SPACEDIM >& alpha_fc,
                  const amrex::MultiFab& beta_cc,
                  const std::array< amrex::MultiFab, NUM_EDGE >& beta_ed,
                  const amrex::MultiFab& gamma_cc,
                  std::array< amrex::MultiFab, AMREX_SPACEDIM >& phi_fc,
                  const std::array< amrex::MultiFab, AMREX_SPACEDIM >& phiorig_fc,
                  const amrex::Real& theta,
                  const amrex::Geometry& geom);

///////////////////////////

#endif
