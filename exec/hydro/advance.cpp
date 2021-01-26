
#include "hydro_functions.H"
#include "hydro_test_functions.H"

#include "common_functions.H"

#include "gmres_functions.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

// argv contains the name of the inputs file entered at the command line
void advance(std::array< MultiFab, AMREX_SPACEDIM >& umac,
             std::array< MultiFab, AMREX_SPACEDIM >& umacNew,
             MultiFab& pres, MultiFab& tracer,
             const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_stoch,
             std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
             MultiFab& beta, MultiFab& gamma,
             std::array< MultiFab, NUM_EDGE >& beta_ed,
             const Geometry geom, const Real& dt,
             TurbForcing& tf)
{

  BL_PROFILE_VAR("advance()",advance);

  const Real* dx = geom.CellSize();
  const Real dtinv = 1.0/dt;
  Real theta_alpha = 1.;
  Real norm_pre_rhs;

  const BoxArray& ba = beta.boxArray();
  const DistributionMapping& dmap = beta.DistributionMap();

   // rhs_p GMRES solve
   MultiFab gmres_rhs_p(ba, dmap, 1, 0);
   gmres_rhs_p.setVal(0.);

  // rhs_u GMRES solve
  std::array< MultiFab, AMREX_SPACEDIM > gmres_rhs_u;
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
      gmres_rhs_u[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
      gmres_rhs_u[d].setVal(0.);
  }

  // laplacian of umac field
  std::array< MultiFab, AMREX_SPACEDIM > Lumac;
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
      Lumac[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
      Lumac[d].setVal(0.);
  }

  // advective terms
  std::array< MultiFab, AMREX_SPACEDIM > advFluxdiv;
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
      advFluxdiv[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
      advFluxdiv[d].setVal(0.);
  }

  ///////////////////////////////////////////

  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    umac[d].FillBoundary(geom.periodicity());
  }

  //////////////////////////
  // Advance tracer

  MultiFab tracerPred(ba,dmap,1,1);
  MultiFab advFluxdivS(ba,dmap,1,1);
  
  tracer.FillBoundary(geom.periodicity());

  // compute -div(c*u)^n
  MkAdvSFluxdiv_cc(umac,tracer,advFluxdivS,geom,0,1,0);
  
  // compute c^{*,n+1} = c^n + dt * (-div(c*u)^n)
  MultiFab::Copy(tracerPred, tracer, 0, 0, 1, 0);
  MultiFab::Saxpy(tracerPred, dt, advFluxdivS, 0, 0, 1, 0);
  tracerPred.FillBoundary(geom.periodicity());

  //////////////////////////

  //////////////////////////////////////////////////
  // ADVANCE velocity field
  //////////////////////////////////////////////////

  // increment advFluxdiv
  MkAdvMFluxdiv(umac,umac,advFluxdiv,dx,0);

  // compute t^n viscous operator, (1/2) L(u)
  // passing in theta_alpha=0 so alpha_fc doesn't matter
  // beta's contain (1/2)*mu
  // this computes the NEGATIVE operator, "(alpha - L_beta)u" so we have to multiply by -1 below
  StagApplyOp(geom,beta,gamma,beta_ed,umac,Lumac,alpha_fc,dx,0.);

  for (int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 0);
    gmres_rhs_u[d].mult(dtinv, 0);
    MultiFab::Add(gmres_rhs_u[d], mfluxdiv_stoch[d], 0, 0, 1, 0);
    // account for the negative viscous operator
    MultiFab::Subtract(gmres_rhs_u[d], Lumac[d], 0, 0, 1, 0);
    MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d], 0, 0, 1, 0);
  }

  // turbulence forcing
  if (turbForcing == 1) {
      tf.AddTurbForcing(gmres_rhs_u,dt,1);
  }

  // initial guess for new solution
  // for pressure use previous solution as initial guess
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(umacNew[d], umac[d], 0, 0, 1, 0);
  }
    
  Real gmres_abs_tol_in = gmres_abs_tol; // save this

  // call GMRES to compute predictor
  GMRES gmres(ba,dmap,geom);
  gmres.Solve(gmres_rhs_u,gmres_rhs_p,umacNew,pres,
              alpha_fc,beta,beta_ed,gamma,
              theta_alpha,geom,norm_pre_rhs);

  // for the corrector gmres solve we want the stopping criteria based on the
  // norm of the preconditioned rhs from the predictor gmres solve.  otherwise
  // for cases where du in the corrector should be small the gmres stalls
  gmres_abs_tol = amrex::max(gmres_abs_tol_in, norm_pre_rhs*gmres_rel_tol);

  // Compute predictor advective term
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    umacNew[d].FillBoundary(geom.periodicity());
  }

  //////////////////////////
  // Advance tracer
  
  // compute -div(c*u)^{*,n+1} and add to -div(c*u)^n
  MkAdvSFluxdiv_cc(umacNew,tracerPred,advFluxdivS,geom,0,1,1);
  
  // compute c^{*,n+1} = c^n + (dt/2) * (-div(c*u)^n) + (dt/2) * (-div(c*u)^{*,n+1})
  MultiFab::Saxpy(tracer, dt/2.0, advFluxdivS, 0, 0, 1, 0);
  tracer.FillBoundary(geom.periodicity());  

  //////////////////////////
  
  // increment advFluxdiv
  MkAdvMFluxdiv(umacNew,umacNew,advFluxdiv,dx,1);

  // trapezoidal advective terms
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    advFluxdiv[d].mult(0.5, 1);
  }

  // Compute gmres_rhs
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 0);
    gmres_rhs_u[d].mult(dtinv);
    MultiFab::Add(gmres_rhs_u[d], mfluxdiv_stoch[d], 0, 0, 1, 0);
    // account for the negative viscous operator
    MultiFab::Subtract(gmres_rhs_u[d], Lumac[d], 0, 0, 1, 0);
    MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d], 0, 0, 1, 0);
  }

  // turbulence forcing
  if (turbForcing == 1) {
      tf.AddTurbForcing(gmres_rhs_u,dt,0);
  }

  // initial guess for new solution
  // for pressure use previous solution as initial guess
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(umacNew[d], umac[d], 0, 0, 1, 0);
  }

  // call GMRES here
  gmres.Solve(gmres_rhs_u,gmres_rhs_p,umacNew,pres,
              alpha_fc,beta,beta_ed,gamma,
              theta_alpha,geom,norm_pre_rhs);
    
  gmres_abs_tol = gmres_abs_tol_in; // Restore the desired tolerance

  for (int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(umac[d], umacNew[d], 0, 0, 1, 0);
  }
  //////////////////////////////////////////////////

}

