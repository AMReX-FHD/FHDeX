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
		   amrex::MultiFab& rhotot);

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


/////////////////////////////////////////////////////////////////////////////////
// in ComputeMassFluxdiv.cpp

void ComputeMassFluxdiv(amrex::MultiFab& rho,
			amrex::MultiFab& rhotot,
			amrex::MultiFab& diff_mass_fluxdiv,
			amrex::MultiFab& stoch_mass_fluxdiv,
			amrex::MultiFab& diff_mass_flux,
			amrex::MultiFab& stoch_mass_flux,
			const amrex::Real& dt,
			const amrex::Real& stage_time,
			const amrex::Geometry& geom);

/////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
// in DiffusiveMassFlux.cpp

void DiffusiveMassFluxdiv(const amrex::MultiFab& rho,
			  const amrex::MultiFab& rhotot,
			  const amrex::MultiFab& molarconc,
			  const amrex::MultiFab& rhoWchi,
			  const amrex::MultiFab& Gamma,
			  const amrex::MultiFab& diff_mass_fluxdiv,
			  std::array< amrex::MultiFab, AMREX_SPACEDIM >& diff_mass_flux,
			  const amrex::Geometry& geom);

void DiffusiveMassFlux(const amrex::MultiFab& rho,
		       const amrex::MultiFab& rhotot,
		       const amrex::MultiFab& molarconc,
		       const amrex::MultiFab& rhoWchi,
		       const amrex::MultiFab& Gamma,
		       std::array< amrex::MultiFab, AMREX_SPACEDIM >& diff_mass_flux,
		       const amrex::Geometry& geom);

/////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
// in CorrectionFlux.cpp

void CorrectionFlux(const amrex::MultiFab& rho, const amrex::MultiFab& rhotot,
		    std::array< amrex::MultiFab, AMREX_SPACEDIM >& flux);

/////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
// in MatvecMul.cpp

void MatvecMul(amrex::MultiFab& x,
	  const amrex::MultiFab& A,
	  const int& nc);

/////////////////////////////////////////////////////////////////////////////////

#endif
