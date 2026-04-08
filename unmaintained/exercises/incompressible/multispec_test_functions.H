#ifndef _multispec_test_functions_H_
#define _multispec_test_functions_H_

#include "common_functions.H"
#include "gmres_functions.H"

///////////////////////////
// in advance.cpp
void advance(std::array< amrex::MultiFab, AMREX_SPACEDIM >& umac,
    std::array< amrex::MultiFab, AMREX_SPACEDIM >& umacNew,
    amrex::MultiFab& pres, amrex::MultiFab& tracer,
    amrex::MultiFab& rho, amrex::MultiFab& rhotot,
    const std::array< amrex::MultiFab, AMREX_SPACEDIM >& mfluxdiv_predict,
    const std::array< amrex::MultiFab, AMREX_SPACEDIM >& mfluxdiv_correct,
    std::array< amrex::MultiFab, AMREX_SPACEDIM >& alpha_fc,
    const amrex::MultiFab& beta, const amrex::MultiFab& gamma,
    # if (AMREX_SPACEDIM == 2)
    const std::array< amrex::MultiFab, 1 >& beta_ed,
    # endif
    # if (AMREX_SPACEDIM == 3)
    const std::array< amrex::MultiFab, 3 >& beta_ed,
    # endif
    const amrex::Geometry geom, const amrex::Real& dt);

///////////////////////////

///////////////////////////
// in WritePlotFile.cpp

void WritePlotFile(const int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   std::array< MultiFab, AMREX_SPACEDIM >& umac,
    const MultiFab& rho,
    const MultiFab& rhoxav,
    const MultiFab& tracer,
                   const MultiFab& pres);

///////////////////////////

#endif