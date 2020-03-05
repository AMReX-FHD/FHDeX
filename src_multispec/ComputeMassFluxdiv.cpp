#include "multispec_functions.H"
#include "multispec_functions_F.H"

#include "common_functions.H"

#include "multispec_namespace.H"
#include "common_namespace.H"

using namespace multispec;
using namespace common;
using namespace amrex;

void ComputeMassFluxdiv(MultiFab& rho, MultiFab& rhotot,
			MultiFab& diff_mass_fluxdiv,
			std::array< MultiFab, AMREX_SPACEDIM >& diff_mass_flux,
			const Real& dt, const Real& stage_time, const Geometry& geom)

// void ComputeMassFluxdiv(MultiFab& rho, MultiFab& rhotot,
// 			MultiFab& diff_mass_fluxdiv, MultiFab& stoch_mass_fluxdiv,
// 			std::array< MultiFab, AMREX_SPACEDIM >& diff_mass_flux,
// 		        std::array< MultiFab, AMREX_SPACEDIM >& stoch_mass_flux,
// 			const Real& dt, const Real& stage_time, const Geometry& geom)
  
{

  BL_PROFILE_VAR("ComputeMassFluxdiv()",ComputeMassFluxdiv);

  BoxArray ba = rho.boxArray();
  DistributionMapping dmapp = rho.DistributionMap();
  int nspecies = rho.nComp();
  int nspecies2 = nspecies*nspecies;

  const Real* dx = geom.CellSize();
  
  MultiFab rhoWchi(         ba, dmapp, nspecies2, 1);  // rho*W*chi*Gamma
  MultiFab molarconc(       ba, dmapp, nspecies, 1);   // molar concentration
  MultiFab molmtot(         ba, dmapp, 1, 1);          // total molar mass
  MultiFab Hessian(         ba, dmapp, nspecies2, 1);  // Hessian-matrix
  MultiFab Gamma(           ba, dmapp, nspecies2, 1);  // Gamma-matrix
  MultiFab D_bar(           ba, dmapp, nspecies2, 1);  // D_bar-matrix
  MultiFab D_therm(         ba, dmapp, nspecies2, 1);  // DT-matrix
  MultiFab sqrtLonsager_fc( ba, dmapp, nspecies2, 1);  // cholesky factored Lonsager on faces
   
  rhoWchi.setVal(0.);
  molarconc.setVal(0.);
  molmtot.setVal(0.);
  Hessian.setVal(0.);
  Gamma.setVal(0.);
  D_bar.setVal(0.);
  D_therm.setVal(0.);
  sqrtLonsager_fc.setVal(0.);

  // phi.FillBoundary(geom.periodicity());
  
  ComputeRhotot(rho,rhotot);
  
  // compute molmtot, molarconc (primitive variables) for 
  // each-cell from rho(conserved) 
  ComputeMolconcMolmtot(rho,rhotot,molarconc,molmtot);

  molarconc.FillBoundary(geom.periodicity()); // hack
  molmtot.FillBoundary(geom.periodicity()); // hack
  
  // populate D_bar and Hessian matrix 
  ComputeMixtureProperties(rho,rhotot,D_bar,D_therm,Hessian);

  D_bar.FillBoundary(geom.periodicity()); // hack
  Hessian.FillBoundary(geom.periodicity()); // hack
  
  // compute Gamma from Hessian
  ComputeGamma(molarconc,Hessian,Gamma);

  Gamma.FillBoundary(geom.periodicity()); // hack
  
  // compute rho*W*chi and zeta/Temp
  ComputeRhoWChi(rho,rhotot,molarconc,rhoWchi,D_bar);

  rhoWchi.FillBoundary(geom.periodicity()); // hack

  // compute diffusive mass fluxes, "-F = rho*W*chi*Gamma*grad(x) - ..."
  DiffusiveMassFluxdiv(rho,rhotot,molarconc,rhoWchi,Gamma,diff_mass_fluxdiv,diff_mass_flux,geom);
  
  // compute external forcing for manufactured solution and add to diff_mass_fluxdiv
  // we should move this to occur before the call to compute_mass_fluxdiv and into
  // the advance_timestep routines
  // external_source(rho,diff_mass_fluxdiv,stage_time,geom);

  // compute stochastic fluxdiv 
  // if (variance_coef_mass != 0.0) {

  //   // compute face-centered cholesky-factored Lonsager^(1/2)
  //   compute_sqrtLonsager_fc(rho,rhotot,sqrtLonsager_fc,geom);

  //   stochastic_mass_fluxdiv(rho,rhotot,sqrtLonsager_fc,
  //   			       stoch_mass_fluxdiv,stoch_mass_flux,
  //   			       dt,geom);

  // }

}
