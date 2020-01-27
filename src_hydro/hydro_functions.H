#ifndef _hydro_functions_H_
#define _hydro_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_MLABecLaplacian.H>
#include "common_functions.H"
#include "gmres_functions.H"


using namespace amrex;

/////////////////////////////////////////////////////////////////////////////////
// in MkAdvMFluxdiv.cpp

void MkAdvMFluxdiv(const std::array<MultiFab, AMREX_SPACEDIM> & umac,
		   const std::array<MultiFab, AMREX_SPACEDIM> & m,
		   std::array<MultiFab, AMREX_SPACEDIM> & m_update,
		   const Real * dx,
		   const int & increment=0);

void MkAdvSFluxdiv(const std::array<MultiFab, AMREX_SPACEDIM> & umac,
		   const MultiFab & m,
		   MultiFab & m_update,
		   const Real * dx,
		   const Geometry & geom,
		   const int & m_comp=0,
		   const int & increment=0);

/////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////
// in MacProj.cpp

void MacProj(std::array<MultiFab, AMREX_SPACEDIM> & umac,
             const MultiFab & rho,
             const Geometry & geom,
             const bool & full_solve=true);

void ComputeMACSolverRHS(MultiFab & solverrhs,
			 const MultiFab & macrhs,
			 const std::array<MultiFab, AMREX_SPACEDIM> & umac,
			 const Geometry & geom);

void AvgFaceBcoeffsInv(std::array<MultiFab, AMREX_SPACEDIM> & facebcoef,
		       const MultiFab & rhocc);

void SetMacSolverBCs(MLABecLaplacian & mlabec);

/////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////
// in HydroBCs.cpp

void SetPressureBC(MultiFab & pressure, const Geometry & geom);

/////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
// in advance.cpp

void advanceStokes(  std::array< MultiFab, AMREX_SPACEDIM >& umac,
	       MultiFab& pres,
	       const std::array< MultiFab, AMREX_SPACEDIM >& stochMfluxdiv,
	       std::array< MultiFab, AMREX_SPACEDIM >& sourceTerms,
	       std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
	       MultiFab& beta,
               MultiFab& gamma,
	       std::array< MultiFab, NUM_EDGE >& beta_ed,
	       const Geometry geom, const Real& dt);

void advanceLowMach(  std::array< MultiFab, AMREX_SPACEDIM >& umac,
	       std::array< MultiFab, AMREX_SPACEDIM >& umacNew,
	       MultiFab& pres, MultiFab& tracer,
	       const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_predict,
	       const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_correct,
	       std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
	       const MultiFab& beta, const MultiFab& gamma,
	       const std::array< MultiFab, NUM_EDGE >& beta_ed,
	       const Geometry geom, const Real& dt);

#endif
