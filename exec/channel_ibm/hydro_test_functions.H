#ifndef _DiffusionTest_functions_H_
#define _DiffusionTest_functions_H_

#include "common_functions.H"
#include "gmres_functions.H"

#include "IBCore.H"


///////////////////////////
// in advance.cpp
void advance(std::array< MultiFab, AMREX_SPACEDIM > & umac,
	         std::array< MultiFab, AMREX_SPACEDIM > & umacNew,
	         MultiFab & pres, MultiFab & tracer,
	         const std::array< MultiFab, AMREX_SPACEDIM > & mfluxdiv_predict,
	         const std::array< MultiFab, AMREX_SPACEDIM > & mfluxdiv_correct,
	         const std::array< MultiFab, AMREX_SPACEDIM > & alpha_fc,
	         const MultiFab & beta, const MultiFab & gamma,
# if (AMREX_SPACEDIM == 2)
	         const std::array< MultiFab, 1 > & beta_ed,
# endif
# if (AMREX_SPACEDIM == 3)
	         const std::array< MultiFab, 3 > & beta_ed,
# endif
             IBCore & ib_core,
	         const Geometry geom, const Real & dt);

///////////////////////////




///////////////////////////
// in WritePlotFile.cpp

void WritePlotFile(const int step,
                   const Real time,
                   const Geometry geom,
                   std::array< MultiFab, AMREX_SPACEDIM > & umac,
	        	   const MultiFab & tracer,
                   const MultiFab & pres);

///////////////////////////

#endif
