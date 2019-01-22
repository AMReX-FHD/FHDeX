#ifndef _multispec_functions_Ho
#define _multispec_functions_H_

#define MAX_SPECIES 10

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_MLABecLaplacian.H>

/////////////////////////////////////////////////////////////////////////////////
// in multispec_functions.cpp

void InitializeMultispecNamespace();

/////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////
// in MassFluxUtil.cpp

void ComputeMolconcMolmtot(const amrex::MultiFab& rho,
			   const amrex::MultiFab& rhotot,
			   amrex::MultiFab& molarconc,
			   amrex::MultiFab& molmtot);
  
void ComputeRhotot(const amrex::MultiFab& rho,
		   amrex::MultiFab& rhotot,
		   const amrex::Geometry& geom);

void ComputeGamma(const amrex::MultiFab& molarconc,
		  const amrex::MultiFab& Hessian,
		  amrex::MultiFab& Gamma);

void ComputeRhoWChi(const amrex::MultiFab& rho,
		    const amrex::MultiFab& rhotot,
		    const amrex::MultiFab& molarconc,
		    amrex::MultiFab& rhoWchi,
		    const amrex::MultiFab& D_bar);

void ComputeZetaByTemp(const amrex::MultiFab& molarconc,
		       const amrex::MultiFab& D_bar,
		       const amrex::MultiFab& Temp,
		       amrex::MultiFab& zeta_by_Temp,
		       const amrex::MultiFab& D_therm);

/////////////////////////////////////////////////////////////////////////////////

#endif
