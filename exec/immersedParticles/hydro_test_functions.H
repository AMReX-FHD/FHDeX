#ifndef _DiffusionTest_functions_H_
#define _DiffusionTest_functions_H_

#include "common_functions.H"
#include "gmres_functions.H"

///////////////////////////
// in WritePlotFile.cpp

void WritePlotFile(const int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   std::array< MultiFab, AMREX_SPACEDIM >& umac,
		   const MultiFab& tracer,
                   const MultiFab& pres);

///////////////////////////

#endif
