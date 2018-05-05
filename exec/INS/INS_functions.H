#ifndef _INS_functions_H_
#define _INS_functions_H_

#include "common_functions.H"

///////////////////////////
// in WritePlotFile.cpp

void WritePlotFile(const int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   const MultiFab& rhotot,
                   const std::array< MultiFab, AMREX_SPACEDIM >& umac);

///////////////////////////

#endif
