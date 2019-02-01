#include "multispec_functions.H"
#include "multispec_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "multispec_namespace.H"
#include "common_namespace.H"

using namespace multispec;
using namespace common;
using namespace amrex;

void ComputeMassFluxdiv(MultiFab& rho, MultiFab& rhotot,
			MultiFab& diff_mass_fluxdiv, MultiFab& stoch_mass_fluxdiv,
			MultiFab& diff_mass_flux, MultiFab& stoch_mass_flux,
			const Real& dt, const Real& stage_time, const Geometry& geom)
  
{

  BL_PROFILE_VAR("ComputeMassFluxdiv()",ComputeMassFluxdiv);

  BoxArray ba = rho.boxArray();
  DistributionMapping dmap = rho.DistributionMap();
  int nspecies = rho.nComp();
  int nspecies2 = nspecies*nspecies;

  const Real* dx = geom.CellSize();
  
  MultiFab rhoWchi(         ba, dmap, 1, 1);  // rho*W*chi*Gamma
  MultiFab molarconc(       ba, dmap, nspecies, 1);  // molar concentration
  MultiFab molmtot(         ba, dmap, 1, 1);  // total molar mass
  MultiFab Hessian(         ba, dmap, nspecies2, 1);  // Hessian-matrix
  MultiFab Gamma(            ba, dmap, nspecies2, 1);  // Gamma-matrix
  MultiFab D_bar(           ba, dmap, nspecies2, 1);  // D_bar-matrix
  MultiFab D_therm(         ba, dmap, nspecies2, 1);  // DT-matrix
  MultiFab sqrtLonsager_fc( ba, dmap, nspecies2, 1);  // cholesky factored Lonsager on faces
   
  rhoWchi.setVal(0.);
  molarconc.setVal(0.);
  molmtot.setVal(0.);
  Hessian.setVal(0.);
  Gamma.setVal(0.);
  D_bar.setVal(0.);
  D_therm.setVal(0.);
  sqrtLonsager_fc.setVal(0.);
  
  ComputeRhotot(rho,rhotot);

  rhotot.FillBoundary(geom.periodicity());
    
  // compute molmtot, molarconc (primitive variables) for 
  // each-cell from rho(conserved) 
  ComputeMolconcMolmtot(rho,rhotot,molarconc,molmtot);

  molarconc.FillBoundary(geom.periodicity());
  molmtot.FillBoundary(geom.periodicity());
    
  // populate D_bar and Hessian matrix 
  // compute_mixture_properties(mla,rho,rhotot,D_bar,D_therm,Hessian);

  D_bar.FillBoundary(geom.periodicity());
  D_therm.FillBoundary(geom.periodicity());
  Hessian.FillBoundary(geom.periodicity());
  
  // compute Gamma from Hessian
  ComputeGamma(molarconc,Hessian,Gamma);

  Gamma.FillBoundary(geom.periodicity());
   
  // compute rho*W*chi and zeta/Temp
  ComputeRhoWChi(rho,rhotot,molarconc,rhoWchi,D_bar);

  rhoWchi.FillBoundary(geom.periodicity());

  // compute diffusive mass fluxes, "-F = rho*W*chi*Gamma*grad(x) - ..."
  // diffusive_mass_fluxdiv(mla,rho,rhotot,molarconc,rhoWchi,Gamma,
  // 			   diff_mass_fluxdiv,Temp,zeta_by_Temp,gradp_baro,
  // 			   diff_mass_flux,dx,the_bc_tower);

  // compute external forcing for manufactured solution and add to diff_mass_fluxdiv
  // we should move this to occur before the call to compute_mass_fluxdiv and into
  // the advance_timestep routines
  // external_source(mla,rho,diff_mass_fluxdiv,dx,stage_time);

  // compute stochastic fluxdiv 
  // if (variance_coef_mass != 0.0) {

  //   // compute face-centered cholesky-factored Lonsager^(1/2)
  //   compute_sqrtLonsager_fc(rho,rhotot,sqrtLonsager_fc,dx);

  //   stochastic_mass_fluxdiv(rho,rhotot,sqrtLonsager_fc,
  //   			       stoch_mass_fluxdiv,stoch_mass_flux,
  //   			       dx,dt,weights,the_bc_tower%bc_tower_array);

  // }

}
