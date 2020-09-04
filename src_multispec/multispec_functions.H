#ifndef _multispec_functions_H_
#define _multispec_functions_H_

// used for binary diffusion coefficients
// needs to be MAX_SPECIES*(MAX_SPECIES-1)/2, MAX_SPECIES is in common_functions.H
#define MAX_ELEMENT 6

#include <AMReX_MultiFab.H>

#include "multispec_functions_F.H"
#include "multispec_namespace.H"

#include "common_functions.H"

#include "StochMassFlux.H"
#include "StochMomFlux.H"

using namespace multispec;
using namespace amrex;

/////////////////////////////////////////////////////////////////////////////////
// in multispec_functions.cpp

void InitializeMultispecNamespace();

/////////////////////////////////////////////////////////////////////////////////
// in ComputeMassFluxdiv.cpp

void ComputeMassFluxdiv(MultiFab& rho,
                        MultiFab& rhotot,
                        const MultiFab& Temp,
			MultiFab& diff_mass_fluxdiv,
                        MultiFab& stoch_mass_fluxdiv,
			std::array<MultiFab,AMREX_SPACEDIM>& diff_mass_flux,
                        std::array<MultiFab,AMREX_SPACEDIM>& stoch_mass_flux,
                        StochMassFlux& sMassFlux,
			const Real& dt, const Real& stage_time, const Geometry& geom,
                        Vector<Real>& weights);

/////////////////////////////////////////////////////////////////////////////////
// in DiffusiveMassFlux.cpp

void DiffusiveMassFluxdiv(const MultiFab& rho,
			  const MultiFab& rhotot,
			  const MultiFab& molarconc,
			  const MultiFab& rhoWchi,
			  const MultiFab& Gamma,
			  MultiFab& diff_mass_fluxdiv,
			  std::array< MultiFab, AMREX_SPACEDIM >& diff_mass_flux,
			  const Geometry& geom);

void DiffusiveMassFlux(const MultiFab& rho,
		       const MultiFab& rhotot,
		       const MultiFab& molarconc,
		       const MultiFab& rhoWchi,
		       const MultiFab& Gamma,
		       std::array< MultiFab, AMREX_SPACEDIM >& diff_mass_flux,
		       const Geometry& geom);

/////////////////////////////////////////////////////////////////////////////////
// in ComputeMixtureProperties.cpp
void ComputeMixtureProperties(const MultiFab& rho,
			      const MultiFab& rhotot,
			      MultiFab& D_bar,
			      MultiFab& D_therm,
			      MultiFab& Hessian);

/////////////////////////////////////////////////////////////////////////////////
// in CorrectionFlux.cpp

void CorrectionFlux(const MultiFab& rho, const MultiFab& rhotot,
		    std::array< MultiFab, AMREX_SPACEDIM >& flux);

/////////////////////////////////////////////////////////////////////////////////
// in CorrectionFlux.cpp

void CorrectionFlux(const MultiFab& rho, const MultiFab& rhotot,
		    std::array< MultiFab, AMREX_SPACEDIM >& flux);

/////////////////////////////////////////////////////////////////////////////////
// in InitialProjection.cpp

void InitialProjection(std::array< MultiFab, AMREX_SPACEDIM >& umac,
                       MultiFab& rho, MultiFab& rhotot,
                       MultiFab& diff_mass_fluxdiv,
                       MultiFab& stoch_mass_fluxdiv,
                       std::array< MultiFab, AMREX_SPACEDIM >& stoch_mass_flux,
                       StochMassFlux& sMassFlux,
                       const MultiFab& Temp, const MultiFab& eta,
                       const std::array< MultiFab, NUM_EDGE >& eta_ed,
                       const Real& dt, const Real& time, const Geometry& geom);

/////////////////////////////////////////////////////////////////////////////////
// in MassFluxUtil.cpp

void ComputeMolconcMolmtot(const MultiFab& rho,
			   const MultiFab& rhotot,
			   MultiFab& molarconc,
			   MultiFab& molmtot);

void ComputeGamma(const MultiFab& molarconc,
		  const MultiFab& Hessian,
		  MultiFab& Gamma);

void ComputeRhoWChi(const MultiFab& rho,
		    const MultiFab& rhotot,
		    const MultiFab& molarconc,
		    MultiFab& rhoWchi,
		    const MultiFab& D_bar);

void ComputeZetaByTemp(const MultiFab& molarconc,
 		       const MultiFab& D_bar,
 		       const MultiFab& Temp,
 		       MultiFab& zeta_by_Temp,
 		       const MultiFab& D_therm);

void ComputeSqrtLonsagerFC(const MultiFab& rho,
                           const MultiFab& rhotot,
                           std::array< MultiFab, AMREX_SPACEDIM >& sqrtLonsager_fc,
                           const Geometry& geom);

/////////////////////////////////////////////////////////////////////////////////
// in MatvecMul.cpp

void MatvecMul(MultiFab& x,
	       const MultiFab& A);


/////////////////////////////////////////////////////////////////////////////////
// in MkDiffusiveMFluxdiv.cpp

void MkDiffusiveMFluxdiv(std::array<MultiFab, AMREX_SPACEDIM> & m_update,
                         const std::array<MultiFab, AMREX_SPACEDIM> & umac,
                         const MultiFab& eta,
                         const std::array<MultiFab, NUM_EDGE> & eta_ed,
                         const MultiFab& kappa,
                         const Geometry& geom,
                         const Real* dx,
                         const int& increment);

/////////////////////////////////////////////////////////////////////////////////
// in RhoUtil.cpp

void ComputeRhotot(const MultiFab& rho, MultiFab& rhotot, int include_ghost=0);

void ConvertRhoCToC(MultiFab& rho, const MultiFab& rhotot, MultiFab& conc, int rho_to_c);

void FillRhoRhototGhost(MultiFab& rho, MultiFab& rhotot, const Geometry& geom);

void FillRhototGhost(MultiFab& rhotot_in, const MultiFab& conc_in, const Geometry& geom);
/////////////////////////////////////////////////////////////////////////////////

#endif
