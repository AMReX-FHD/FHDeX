#ifndef _hydro_functions_H
#define _hydro_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_MLABecLaplacian.H>

/////////////////////////////////////////////////////////////////////////////////
// in MkAdvMFluxdiv.cpp

void MkAdvMFluxdiv(const std::array<amrex::MultiFab, AMREX_SPACEDIM>& umac,
		   const std::array<amrex::MultiFab, AMREX_SPACEDIM>& m,
		   std::array<amrex::MultiFab, AMREX_SPACEDIM>& m_update,
		   const amrex::Real* dx,
		   const int& increment=0);

void MkAdvSFluxdiv(const std::array<amrex::MultiFab, AMREX_SPACEDIM>& umac,
		   const amrex::MultiFab& m,
		   amrex::MultiFab& m_update,
		   const amrex::Real* dx,
		   const amrex::Geometry& geom,
		   const int& increment=0);

/////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
// in MacProj.cpp

void MacProj (std::array< amrex::MultiFab, AMREX_SPACEDIM >& umac,
	      const amrex::MultiFab& rho,
	      const amrex::Geometry& geom,
	      const bool& full_solve=true);

void ComputeMACSolverRHS (amrex::MultiFab& solverrhs,
			  const amrex::MultiFab& macrhs,
			  const std::array< amrex::MultiFab, AMREX_SPACEDIM >& umac,
			  const amrex::Geometry& geom);

void AvgFaceBcoeffsInv(std::array< amrex::MultiFab, AMREX_SPACEDIM >& facebcoef,
		       const amrex::MultiFab& rhocc);

void SetMacSolverBCs(amrex::MLABecLaplacian& mlabec);

/////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
// in advance.cpp
void advance(std::array< amrex::MultiFab, AMREX_SPACEDIM >& umac,
	     std::array< amrex::MultiFab, AMREX_SPACEDIM >& umacNew,
	     amrex::MultiFab& pres, amrex::MultiFab& tracer,
	     const std::array< amrex::MultiFab, AMREX_SPACEDIM >& mfluxdiv_predict,
	     const std::array< amrex::MultiFab, AMREX_SPACEDIM >& mfluxdiv_correct,
	     const std::array< amrex::MultiFab, AMREX_SPACEDIM >& alpha_fc,
	     const amrex::MultiFab& beta, const amrex::MultiFab& gamma,
# if (AMREX_SPACEDIM == 2)
	     const std::array< amrex::MultiFab, 1 >& beta_ed,
# endif
# if (AMREX_SPACEDIM == 3)
	     const std::array< amrex::MultiFab, 3 >& beta_ed,
# endif
	     const amrex::Geometry geom, const amrex::Real& dt);
/////////////////////////////////////////////////////////////////////////////////

#endif
