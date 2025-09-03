#include "multispec_functions.H"
#include "StochMassFlux.H"

void ComputeMassFluxdiv(MultiFab& rho,
                        MultiFab& rhotot,
                        const MultiFab& Temp,
  MultiFab& diff_mass_fluxdiv,
                        MultiFab& stoch_mass_fluxdiv,
  std::array<MultiFab,AMREX_SPACEDIM>& diff_mass_flux,
                        std::array<MultiFab,AMREX_SPACEDIM>& stoch_mass_flux,
                        StochMassFlux& sMassFlux,
  const Real& dt, const Real& stage_time, const Geometry& geom,
                        Vector<Real>& weights,
                        MultiFab& charge,
                        std::array<MultiFab,AMREX_SPACEDIM>& grad_Epot,
                        MultiFab& Epot,
                        MultiFab& permittivity,
                        const int& zero_initial_Epot)
{

  BL_PROFILE_VAR("ComputeMassFluxdiv()",ComputeMassFluxdiv);

  BoxArray ba = rho.boxArray();
  DistributionMapping dmap = rho.DistributionMap();

  int ng = rho.nGrow();
  int nspecies2 = nspecies*nspecies;

  MultiFab rhoWchi(         ba, dmap, nspecies2, ng);  // rho*W*chi*Gamma
  MultiFab molarconc(       ba, dmap, nspecies , ng);  // molar concentration
  MultiFab molmtot(         ba, dmap, 1        , ng);  // total molar mass
  MultiFab Hessian(         ba, dmap, nspecies2, ng);  // Hessian-matrix
  MultiFab Gamma(           ba, dmap, nspecies2, ng);  // Gamma-matrix
  MultiFab D_bar(           ba, dmap, nspecies2, ng);  // D_bar-matrix
  MultiFab D_therm(         ba, dmap, nspecies2, ng);  // DT-matrix
  MultiFab zeta_by_Temp(    ba, dmap, nspecies2, ng);  // for Thermo-diffusion
 // if( use_flory_huggins == 1 )
  MultiFab massfrac(       ba, dmap, nspecies , ng);  // molar concentration

 // rhoWchi.setVal(0.);

  std::array< MultiFab, AMREX_SPACEDIM > sqrtLonsager_fc;
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    sqrtLonsager_fc[d].define(convert(ba,nodal_flag_dir[d]), dmap, nspecies2, 0);
  }
  ComputeRhotot(rho,rhotot,1);

  // compute molmtot, molarconc (primitive variables) for
  // each-cell from rho(conserved)
  ComputeMolconcMolmtot(rho,rhotot,molarconc,molmtot);

  // populate D_bar and Hessian matrix
  ComputeMixtureProperties(rho,rhotot,D_bar,D_therm,Hessian);

  // compute Gamma from Hessian

  if (use_flory_huggins == 1) {
    ComputeMassfrac(rho,rhotot,massfrac);
    ComputeFHGamma(massfrac,Gamma);
  } else {
    ComputeGamma(molarconc,Hessian,Gamma);
  }

  // compute rho*W*chi and zeta/Temp
  ComputeRhoWChi(rho,rhotot,molarconc,rhoWchi,D_bar);
  //ComputeZetaByTemp(molarconc,D_Bar,Temp,zeta_by_Temp,D_therm);
  if (is_nonisothermal == 1) {
    Abort("ComputeMassFluxDiv: implement is_nonisothermal");
  }

  // compute diffusive mass fluxes, "-F = rho*W*chi*Gamma*grad(x) - ..."
  if (use_flory_huggins == 1) {
    DiffusiveMassFluxdiv(rho,rhotot,massfrac,rhoWchi,Gamma,diff_mass_fluxdiv,diff_mass_flux,geom);
  } else {
    DiffusiveMassFluxdiv(rho,rhotot,molarconc,rhoWchi,Gamma,diff_mass_fluxdiv,diff_mass_flux,geom);
  }

  // compute external forcing for manufactured solution and add to diff_mass_fluxdiv
  // we should move this to occur before the call to compute_mass_fluxdiv and into
  // the advance_timestep routines
  // external_source(rho,diff_mass_fluxdiv,stage_time,geom);
  //if (prob_type == 4 || prob_type == 5) {
  if (prob_type == 4 ) {
    Abort("ComputMassFluxdiv: external source not implemented yet for this prob_type");
  }

  // compute stochastic fluxdiv
  if (variance_coef_mass != 0.) {

    // compute face-centered cholesky-factored Lonsager^(1/2)
    ComputeSqrtLonsagerFC(rho,rhotot,sqrtLonsager_fc,geom);

    sMassFlux.StochMassFluxDiv(rho,rhotot,sqrtLonsager_fc,stoch_mass_fluxdiv,stoch_mass_flux,
                               dt,weights);

  }

  if (use_charged_fluid) {
    ElectroDiffusiveMassFluxdiv(rho,Temp,rhoWchi,diff_mass_flux,diff_mass_fluxdiv,
                                stoch_mass_flux,charge,grad_Epot,Epot,permittivity,
                                dt,1,geom);
  }

}