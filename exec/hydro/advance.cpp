
#include "hydro_functions.H"

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
             const Geometry geom, const Real& dt)
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
  //////////////////////////

  MultiFab tracerPred(ba,dmap,1,1);
  MultiFab advFluxdivS(ba,dmap,1,1);
  
  // Compute tracer:
  tracer.FillBoundary(geom.periodicity());
  MkAdvSFluxdiv_cc(umac,tracer,advFluxdivS,geom,0,1,0);
  advFluxdivS.mult(dt, 1);

  // compute predictor
  MultiFab::Copy(tracerPred, tracer, 0, 0, 1, 1);
  MultiFab::Add(tracerPred, advFluxdivS, 0, 0, 1, 0);
  tracerPred.FillBoundary(geom.periodicity());
  // FIXME need to fill physical boundary condition ghost cells for tracer
  MkAdvSFluxdiv_cc(umac,tracerPred,advFluxdivS,geom,0,1,0);
  advFluxdivS.mult(dt, 1);

  // advance in time
  MultiFab::Add(tracer, tracerPred, 0, 0, 1, 0);
  MultiFab::Add(tracer, advFluxdivS, 0, 0, 1, 0);
  tracer.mult(0.5, 1);

  //////////////////////////

  //////////////////////////////////////////////////
  // ADVANCE velocity field
  //////////////////////////////////////////////////

  // increment advFluxdiv
  MkAdvMFluxdiv(umac,umac,advFluxdiv,dx,0);

  // compute t^n viscous operator
  // passing in theta_alpha=0 so alpha_fc doesn't matter
  // this computes the NEGATIVE operator so we have to multiply by -1 below
  StagApplyOp(geom,beta,gamma,beta_ed,
	      umac,Lumac,alpha_fc,dx,0.);

  for (int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 0);
    gmres_rhs_u[d].mult(dtinv, 0);
    MultiFab::Add(gmres_rhs_u[d], mfluxdiv_stoch[d], 0, 0, 1, 0);
    // account for the negative viscous operator
    MultiFab::Subtract(gmres_rhs_u[d], Lumac[d], 0, 0, 1, 0);
    MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d], 0, 0, 1, 0);
  }

  // initial guess for new solution
  // for pressure use previous solution as initial guess
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(umacNew[d], umac[d], 0, 0, 1, 0);
  }

  // call GMRES to compute predictor
  GMRES gmres(ba,dmap,geom);
  gmres.Solve(gmres_rhs_u,gmres_rhs_p,umacNew,pres,
              alpha_fc,beta,beta_ed,gamma,
              theta_alpha,geom,norm_pre_rhs);

  // Compute predictor advective term
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    umacNew[d].FillBoundary(geom.periodicity());
  }

  // increment advFluxdiv
  MkAdvMFluxdiv(umacNew,umacNew,advFluxdiv,dx,1);

  // Compute gmres_rhs

  // trapezoidal advective terms
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    advFluxdiv[d].mult(0.5, 1);
  }

  for (int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 0);

    gmres_rhs_u[d].mult(dtinv);

    MultiFab::Add(gmres_rhs_u[d], mfluxdiv_stoch[d], 0, 0, 1, 0);
    // account for the negative viscous operator
    MultiFab::Subtract(gmres_rhs_u[d], Lumac[d],     0, 0, 1, 0);
    MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],     0, 0, 1, 0);

    // initial guess for new solution
    // for pressure use previous solution as initial guess
    MultiFab::Copy(umacNew[d], umac[d], 0, 0, 1, 0);
  }

  // call GMRES here
  gmres.Solve(gmres_rhs_u,gmres_rhs_p,umacNew,pres,
              alpha_fc,beta,beta_ed,gamma,
              theta_alpha,geom,norm_pre_rhs);

  for (int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(umac[d], umacNew[d], 0, 0, 1, 0);
  }
  //////////////////////////////////////////////////

}
