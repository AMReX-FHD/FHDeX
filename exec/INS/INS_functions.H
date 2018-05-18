#ifndef _INS_functions_H_
#define _INS_functions_H_

#include "common_functions.H"

///////////////////////////
// in WritePlotFile.cpp

void WritePlotFile(const int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   const MultiFab& rhotot,
                   const std::array< MultiFab, AMREX_SPACEDIM >& umac,
					const MultiFab& div);

///////////////////////////


void eulerStep(const MultiFab& betaCC, const MultiFab& gammaCC,
#if AMREX_SPACEDIM == 2
                 const amrex::MultiFab& betaNodal, const amrex::MultiFab& gammaNodal, 
#endif
#if AMREX_SPACEDIM == 3
                 std::array<amrex::MultiFab, AMREX_SPACEDIM>& betaEdge, 
                 std::array<amrex::MultiFab, AMREX_SPACEDIM>& gammaEdge, 
#endif
                 std::array<MultiFab, AMREX_SPACEDIM>& umacIn, 
                 std::array<MultiFab, AMREX_SPACEDIM>& umacOut,
                 std::array<MultiFab, AMREX_SPACEDIM>& umacNew,
                 std::array<MultiFab, AMREX_SPACEDIM>& alpha,
                 const Geometry geom,
                 Real* dt);

///////////////////////////


#endif
