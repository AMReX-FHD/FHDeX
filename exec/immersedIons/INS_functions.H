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
                   const amrex::Geometry egeom,
                   const MultiFab& particleInstant,
                   const MultiFab& particleMeans,
                   const MultiFab& particleVars,
                   FhdParticleContainer& particles,
                   const MultiFab& charge,
                   const MultiFab& potential,const std::array< MultiFab, AMREX_SPACEDIM >& efield);


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

void WritePlotFileHydro(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   std::array< MultiFab, AMREX_SPACEDIM >& umac,
		   const MultiFab& pres);

void advance(  std::array< MultiFab, AMREX_SPACEDIM >& umac,
	       MultiFab& pres,
	       const std::array< MultiFab, AMREX_SPACEDIM >& stochMfluxdiv,
	       std::array< MultiFab, AMREX_SPACEDIM >& sourceTerms,
	       std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
	       MultiFab& beta,
               MultiFab& gamma,
	       std::array< MultiFab, NUM_EDGE >& beta_ed,
	       const Geometry geom, const Real& dt);


#endif
