
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
             MultiFab& pres,
             const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_stoch,
             std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
             MultiFab& beta, MultiFab& gamma,
             std::array< MultiFab, NUM_EDGE >& beta_ed,
             const Geometry& geom, const Real& dt,
             TurbForcing& turbforce)
{

  BL_PROFILE_VAR("advance()",advance);

  const Real* dx = geom.CellSize();
  const Real dtinv = 1.0/dt;
  Real norm_pre_rhs;

  Real theta_alpha = (algorithm_type == 1) ? 0. : 1.;

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
      MultiFabPhysBCDomainVel(umac[d], geom, d);
      int is_inhomogeneous = 1;
      MultiFabPhysBCMacVel(umac[d], geom, d, is_inhomogeneous);
      umac[d].FillBoundary(geom.periodicity());
  }

  //////////////////////////////////////////////////
  // ADVANCE velocity field
  //////////////////////////////////////////////////

  /*
     Predictor

     Inertial (algorithm_type = 0) - note algorithm_type = -1 is backward Euler

     (u^{n+1,*} - u^n) / dt + grad(p) = A^n + (1/2)(div eta grad)(u^n + u^{n+1,*}) + F

     (1/dt - (1/2) div eta grad) u^{n+1,*} + grad(p) = u^n / dt + A^n + (1/2) div eta grad u^n + F

     Overdamped (algorithm_type = 1)

     grad(p) = (1/2)(div eta grad)(u^n + u^{n+1,*}) + F

     (0 - (1/2) div eta grad) u^{n+1,*} + grad(p) = (1/2) div eta grad u^n + F

  */

  if (algorithm_type == -1 || algorithm_type == 0) {
      // compute advFluxdiv
      MkAdvMFluxdiv(umac,umac,advFluxdiv,dx,0);
  }

  for (int d=0; d<AMREX_SPACEDIM; d++) {
      if (algorithm_type == -1 || algorithm_type == 0) {

          MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 0);
          gmres_rhs_u[d].mult(dtinv, 0);
          MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d], 0, 0, 1, 0);

          if (algorithm_type == 0) {

              // compute t^n viscous operator, (1/2) L(u)
              // passing in theta_alpha=0 so alpha_fc doesn't matter
              // beta's contain (1/2)*mu
              // this computes the NEGATIVE operator, "(alpha - L_beta)u" so we have to multiply by -1 below
              StagApplyOp(geom,beta,gamma,beta_ed,umac,Lumac,alpha_fc,dx,0.);
              // account for the negative viscous operator
              MultiFab::Subtract(gmres_rhs_u[d], Lumac[d], 0, 0, 1, 0);
          }

      } else if (algorithm_type == 1) {
          gmres_rhs_u[d].setVal(0.);
      }
      MultiFab::Add(gmres_rhs_u[d], mfluxdiv_stoch[d], 0, 0, 1, 0);
  }

  ExternalForce(gmres_rhs_u,gmres_rhs_p);

  // turbulence forcing
  if (turbForcing == 1) {
      turbforce.AddTurbForcing(gmres_rhs_u,dt,1);
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

  // for deterministic overdamped, we are done with the time step
  if (algorithm_type == 1 && variance_coef_mom == 0.) {
      for (int d=0; d<AMREX_SPACEDIM; d++) {
          MultiFab::Copy(umac[d], umacNew[d], 0, 0, 1, 0);
      }
      return;
  }

  // for the corrector gmres solve we want the stopping criteria based on the
  // norm of the preconditioned rhs from the predictor gmres solve.  otherwise
  // for cases where du in the corrector should be small the gmres stalls
  gmres_abs_tol = amrex::max(gmres_abs_tol_in, norm_pre_rhs*gmres_rel_tol);

  // Compute predictor advective term
  for (int d=0; d<AMREX_SPACEDIM; d++) {
      MultiFabPhysBCDomainVel(umacNew[d], geom, d);
      int is_inhomogeneous = 1;
      MultiFabPhysBCMacVel(umacNew[d], geom, d, is_inhomogeneous);
      umacNew[d].FillBoundary(geom.periodicity());
  }

  /*
    Corrector

    Inertial (algorithm_type = 0) - note algorithm_type = -1 is backward Euler

    (u^{n+1} - u^n) / dt + grad(p) = (1/2)(A^n + A^{n+1}) + (1/2)(div eta grad)(u^n + u^{n+1}) + F

    (1/dt - (1/2) div eta grad) u^{n+1} + grad(p) = u^n / dt + (1/2)(A^n + A^{n+1}) + (1/2) div eta grad u^n + F

    Overdamped (algorithm_type = 1)

    grad(p) = (1/2)(div eta grad)(u^n + u^{n+1}) + F

    (0 - (1/2) div eta grad) u^{n+1} + grad(p) = (1/2) div eta grad u^n + F

  */

  if (algorithm_type == -1 || algorithm_type == 0) {
      // increment advFluxdiv
      MkAdvMFluxdiv(umacNew,umacNew,advFluxdiv,dx,1);

      // trapezoidal advective terms
      for (int d=0; d<AMREX_SPACEDIM; d++) {
          advFluxdiv[d].mult(0.5, 1);
      }
  }

  // Compute gmres_rhs
  for (int d=0; d<AMREX_SPACEDIM; d++) {
      if (algorithm_type == -1 || algorithm_type == 0) {
          MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 0);
          gmres_rhs_u[d].mult(dtinv);
          MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d], 0, 0, 1, 0);

          if (algorithm_type == 0) {
              // account for the negative viscous operator
              MultiFab::Subtract(gmres_rhs_u[d], Lumac[d], 0, 0, 1, 0);
          }

      } else if (algorithm_type == 1) {
          gmres_rhs_u[d].setVal(0.);
      }
      MultiFab::Add(gmres_rhs_u[d], mfluxdiv_stoch[d], 0, 0, 1, 0);
  }

  ExternalForce(gmres_rhs_u,gmres_rhs_p);

  // turbulence forcing
  if (turbForcing == 1) {
      turbforce.AddTurbForcing(gmres_rhs_u,dt,0);
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

  // these calls are redundant with the calls at the beginning of the time step
  // except for diagnostics which can rely on the ghost cells for stencils such
  // as vorticity
  for (int d=0; d<AMREX_SPACEDIM; d++) {
      MultiFabPhysBCDomainVel(umac[d], geom, d);
      int is_inhomogeneous = 1;
      MultiFabPhysBCMacVel(umac[d], geom, d, is_inhomogeneous);
      umac[d].FillBoundary(geom.periodicity());
  }

  //////////////////////////////////////////////////

}

