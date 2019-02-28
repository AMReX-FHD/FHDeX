#ifndef _INS_functions_H_
#define _INS_functions_H_

#include "common_functions.H"
#include "gmres_functions.H"
#include "FhdParticleContainer.H"

///////////////////////////
// in WritePlotFile.cpp

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   const amrex::Geometry cgeom,
                   const MultiFab& particleInstant,
                   const MultiFab& particleMeans,
                   const MultiFab& particleVars,
                   const MultiFab& particleMembraneFlux,
                   FhdParticleContainer& particles) ;

///////////////////////////


void eulerStep(const MultiFab& betaCC, const MultiFab& gammaCC,
               std::array<amrex::MultiFab, NUM_EDGE>& betaEdge, 
               std::array<MultiFab, AMREX_SPACEDIM>& umacIn, 
               std::array<MultiFab, AMREX_SPACEDIM>& umacOut,
               std::array<MultiFab, AMREX_SPACEDIM>& umacNew,
               std::array<MultiFab, AMREX_SPACEDIM>& alpha,
               const Geometry geom,
               Real* dt);

///////////////////////////

void WritePlotFileHydro(const int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   std::array< MultiFab, AMREX_SPACEDIM >& umac,
		   const MultiFab& tracer,
                   const MultiFab& pres);

void advance(  std::array< MultiFab, AMREX_SPACEDIM >& umac, 
	       std::array< MultiFab, AMREX_SPACEDIM >& umacNew, 
	       MultiFab& pres, MultiFab& tracer, 
	       const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_predict,
	       const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_correct,
	       const std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
	       const MultiFab& beta, const MultiFab& gamma,
	       const std::array< MultiFab, NUM_EDGE >& beta_ed, 
	       const Geometry geom, const Real& dt);


#endif
