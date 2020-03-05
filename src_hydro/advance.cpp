
#include "hydro_functions.H"

#include "common_functions.H"


#include "gmres_functions.H"


#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

// argv contains the name of the inputs file entered at the command line
void advanceStokes(  std::array< MultiFab, AMREX_SPACEDIM >& umac,
	       MultiFab& pres,
	       const std::array< MultiFab, AMREX_SPACEDIM >& stochMfluxdiv,
	       std::array< MultiFab, AMREX_SPACEDIM >& sourceTerms,
	       std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
	       MultiFab& beta,
               MultiFab& gamma,
	       std::array< MultiFab, NUM_EDGE >& beta_ed,
	       const Geometry geom, const Real& dt)
{

    BL_PROFILE_VAR("advance()",advance);
  
    if(fluid_tog ==1) {
        
        Real theta_alpha = 0.;
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

        //////////////////////////////////////////////////
        // ADVANCE velocity field
        //////////////////////////////////////////////////

        // add stochastic forcing to gmres_rhs_u
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Add(gmres_rhs_u[d], stochMfluxdiv[d], 0, 0, 1, 0);
            MultiFab::Add(gmres_rhs_u[d], sourceTerms[d], 0, 0, 1, 0);
        }

        if(zero_net_force == 1)
        {
            Vector<Real> mean_stress_umac(AMREX_SPACEDIM);

            SumStag(geom,gmres_rhs_u,0,mean_stress_umac,true);

            Print() << "correcting mean force: " << mean_stress_umac[0] << "\n";

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                if (geom.isPeriodic(d)) {
                    gmres_rhs_u[d].plus(-mean_stress_umac[d],0,1,0);
                }
            }
        }
      
        // call GMRES
        GMRES(gmres_rhs_u,gmres_rhs_p,umac,pres,alpha_fc,beta,beta_ed,gamma,theta_alpha,geom,norm_pre_rhs);

        for (int i=0; i<AMREX_SPACEDIM; i++) {
            umac[i].FillBoundary(geom.periodicity());
            MultiFABPhysBCDomainVel(umac[i], geom, i);
            MultiFABPhysBCMacVel(umac[i], geom, i);
        }
    }
}


// argv contains the name of the inputs file entered at the command line
//void advanceLowMach(  std::array< MultiFab, AMREX_SPACEDIM >& umac,
//	       std::array< MultiFab, AMREX_SPACEDIM >& umacNew,
//	       MultiFab& pres, MultiFab& tracer,
//	       const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_predict,
//	       const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_correct,
//	       const std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
//	       const MultiFab& beta, const MultiFab& gamma,
//	       const std::array< MultiFab, NUM_EDGE >& beta_ed,
//	       const Geometry geom, const Real& dt)
//{

//    BL_PROFILE_VAR("advance()",advance);

//    const Real * dx  = geom.CellSize();
//    const Real dtinv = 1.0/dt;
//    Real theta_alpha = 1.;
//    Real norm_pre_rhs;

//    const BoxArray &              ba = beta.boxArray();
//    const DistributionMapping & dmap = beta.DistributionMap();


//    /****************************************************************************
//     *                                                                          *
//     * Define temporary MultiFabs                                               *
//     *                                                                          *
//     ***************************************************************************/

//    // RHS pressure in GMRES
//    MultiFab gmres_rhs_p(ba, dmap, 1, 1);
//    gmres_rhs_p.setVal(0.);

//    // RHS velocities in GMRES
//    std::array< MultiFab, AMREX_SPACEDIM > gmres_rhs_u;
//    // Velocity components updated by diffusion operator
//    std::array< MultiFab, AMREX_SPACEDIM > Lumac;
//    // Advective terms
//    std::array< MultiFab, AMREX_SPACEDIM > advFluxdiv;
//    // Advective terms (for predictor)
//    std::array< MultiFab, AMREX_SPACEDIM > advFluxdivPred;
//    // Staggered momentum
//    std::array< MultiFab, AMREX_SPACEDIM > uMom;

//    for (int i=0; i<AMREX_SPACEDIM; i++) {
//           gmres_rhs_u[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
//                 Lumac[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
//            advFluxdiv[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
//        advFluxdivPred[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
//                  uMom[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
//    }

//    // Tracer concentration field for predictor
//    MultiFab tracerPred(ba, dmap, 1, 1);

//    // Tracer advective terms
//    MultiFab advFluxdivS(ba, dmap, 1, 1);


//    //___________________________________________________________________________
//    // Scaled alpha, beta, gamma:

//    // alpha_fc_0 arrays
//    std::array< MultiFab, AMREX_SPACEDIM > alpha_fc_0;

//    for (int i=0; i<AMREX_SPACEDIM; i++){
//        alpha_fc_0[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
//        alpha_fc_0[i].setVal(0.);
//    }

//    // Scaled by 1/2:
//    // beta_wtd cell centered
//    MultiFab beta_wtd(ba, dmap, 1, 1);
//    MultiFab::Copy(beta_wtd, beta, 0, 0, 1, 1);
//    beta_wtd.mult(0.5, 1);

//    // beta_wtd on nodes in 2d, on edges in 3d
//    std::array< MultiFab, NUM_EDGE > beta_ed_wtd;
//#if (AMREX_SPACEDIM == 2)
//    beta_ed_wtd[0].define(convert(ba,nodal_flag), dmap, 1, 1);
//    MultiFab::Copy(beta_ed_wtd[0], beta_ed[0], 0, 0, 1, 1);
//    beta_ed_wtd[0].mult(0.5, 1);
//#elif (AMREX_SPACEDIM == 3)
//    for(int d=0; d<AMREX_SPACEDIM; d++) {
//        beta_ed_wtd[d].define(convert(ba, nodal_flag_edge[d]), dmap, 1, 1);

//        MultiFab::Copy(beta_ed_wtd[d], beta_ed[d], 0, 0, 1, 1);
//        beta_ed_wtd[d].mult(0.5, 1);
//    }
//#endif

//    // Scaled by 1/2:
//    // gamma_wtd cell centered
//    MultiFab gamma_wtd(ba, dmap, 1, 1);
//    MultiFab::Copy(gamma_wtd, gamma, 0, 0, 1, 1);
//    gamma_wtd.mult(-0.5, 1);

//    // Scaled by -1/2:
//    // beta_negwtd cell centered
//    MultiFab beta_negwtd(ba, dmap, 1, 1);
//    MultiFab::Copy(beta_negwtd, beta, 0, 0, 1, 1);
//    beta_negwtd.mult(-0.5, 1);

//    // beta_negwtd on nodes in 2d, on edges in 3d
//    std::array< MultiFab, NUM_EDGE > beta_ed_negwtd;
//#if (AMREX_SPACEDIM == 2)
//    beta_ed_negwtd[0].define(convert(ba,nodal_flag), dmap, 1, 1);
//    MultiFab::Copy(beta_ed_negwtd[0], beta_ed[0], 0, 0, 1, 1);
//    beta_ed_negwtd[0].mult(-0.5, 1);
//#elif (AMREX_SPACEDIM == 3)
//    for(int d=0; d<AMREX_SPACEDIM; d++) {
//        beta_ed_negwtd[d].define(convert(ba,nodal_flag_edge[d]), dmap, 1, 1);

//        MultiFab::Copy(beta_ed_negwtd[d], beta_ed[d], 0, 0, 1, 1);
//        beta_ed_negwtd[d].mult(-0.5, 1);
//    }
//#endif

//    // Scaled by -1/2:
//    // gamma cell centered
//    MultiFab gamma_negwtd(ba, dmap, 1, 1);
//    MultiFab::Copy(gamma_negwtd, gamma, 0, 0, 1, 1);
//    gamma_negwtd.mult(-0.5, 1);


//    /****************************************************************************
//     *                                                                          *
//     * Apply non-stochastic boundary conditions                                 *
//     *                                                                          *
//     ***************************************************************************/


//      //This is applied in advanceStokes() - find a more sensible way to do this.

////    for (int i=0; i<AMREX_SPACEDIM; i++) {
////        umac[i].FillBoundary(geom.periodicity());
////        MultiFABPhysBCDomainVel(umac[i], i, geom, i);
////        MultiFABPhysBCMacVel(umac[i], i, geom, i);

////        MultiFABPhysBCDomainStress(sourceTerms[i], i, geom, i);
////    }


//    /****************************************************************************
//     *                                                                          *
//     * Advance tracer                                                           *
//     *                                                                          *
//     ***************************************************************************/

//    // Compute tracer:
//    tracer.FillBoundary(geom.periodicity());
//    MultiFABPhysBCPres(tracer, geom);

//    MkAdvSFluxdiv(umac, tracer, advFluxdivS, dx, geom, 0);
//    advFluxdivS.mult(dt, 1);

//    // compute predictor
//    MultiFab::Copy(tracerPred, tracer, 0, 0, 1, 0);
//    MultiFab::Add(tracerPred, advFluxdivS, 0, 0, 1, 0);

//    tracerPred.FillBoundary(geom.periodicity());
//    MultiFABPhysBCPres(tracerPred, geom);

//    MkAdvSFluxdiv(umac, tracerPred, advFluxdivS, dx, geom, 0);
//    advFluxdivS.mult(dt, 1);

//    // advance in time
//    MultiFab::Add(tracer, tracerPred,  0, 0, 1, 0);
//    MultiFab::Add(tracer, advFluxdivS, 0, 0, 1, 0);
//    tracer.mult(0.5, 1);

//    // amrex::Print() << "tracer L0 norm = " << tracer.norm0() << "\n";
//    //////////////////////////

//    //////////////////////////////////////////////////
//    // ADVANCE velocity field
//    //////////////////////////////////////////////////

//    // PREDICTOR STEP (heun's method: part 1)
//    // compute advective term

//    for (int i=0; i<AMREX_SPACEDIM; i++)
//        MultiFab::Copy(uMom[i], umac[i], 0, 0, 1, 1);

//    // let rho = 1
//    for (int d=0; d<AMREX_SPACEDIM; d++) {
//        uMom[d].mult(1.0, 1);
//    }

//    for (int i=0; i<AMREX_SPACEDIM; i++) {
//        uMom[i].FillBoundary(geom.periodicity());
//        MultiFABPhysBCDomainVel(uMom[i], i, geom, i);
//        MultiFABPhysBCMacVel(uMom[i], i, geom, i);
//    }

//    MkAdvMFluxdiv(umac, uMom, advFluxdiv, dx, 0);

//    // crank-nicolson terms
//    StagApplyOp(geom,beta_negwtd, gamma_negwtd, beta_ed_negwtd, umac, Lumac, alpha_fc_0, dx, theta_alpha);

//    for (int d=0; d<AMREX_SPACEDIM; d++) {
//        Lumac[d].FillBoundary(geom.periodicity());

//        // Only apply these BCs to the velocity term:
//        // MultiFABPhysBCDomainVel(Lumac[d], d, geom, d);
//        // MultiFABPhysBCMacVel(Lumac[d], d, geom, d);

//        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);

//        gmres_rhs_u[d].mult(dtinv, 1);
//        MultiFab::Add(gmres_rhs_u[d], mfluxdiv_predict[d], 0, 0, 1, 0);
//        MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 0);
//        MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 0);
//    }

//    std::array< MultiFab, AMREX_SPACEDIM > pg;
//    for (int i=0; i<AMREX_SPACEDIM; i++)
//        pg[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);

//    pres.setVal(0.);  // initial guess
//    SetPressureBC(pres, geom);

//    ComputeGrad(pres, pg, 0, 0, 1, geom);

//    for (int i=0; i<AMREX_SPACEDIM; i++) {
//        pg[i].FillBoundary(geom.periodicity());
//        gmres_rhs_u[i].FillBoundary(geom.periodicity());

//        MultiFab::Subtract(gmres_rhs_u[i], pg[i], 0, 0, 1, 1);

//        // Only apply these BCs to the velocity term:
//        // MultiFABPhysBCDomainVel(gmres_rhs_u[i], i, geom, i);
//        // MultiFABPhysBCMacVel(gmres_rhs_u[i], i, geom, i);
//    }

//    // initial guess for new solution
//    for (int i=0; i<AMREX_SPACEDIM; i++)
//        MultiFab::Copy(umacNew[i], umac[i], 0, 0, 1, 1);

//    // call GMRES to compute predictor
//    GMRES(gmres_rhs_u, gmres_rhs_p, umacNew, pres,
//          alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
//          geom, norm_pre_rhs);

//    // Compute predictor advective term
//    // let rho = 1
//    for (int d=0; d<AMREX_SPACEDIM; d++) {
//        umacNew[d].FillBoundary(geom.periodicity());
//        MultiFABPhysBCDomainVel(umacNew[d], d, geom, d);
//        MultiFABPhysBCMacVel(umacNew[d], d, geom, d);

//        MultiFab::Copy(uMom[d], umacNew[d], 0, 0, 1, 0);
//        uMom[d].mult(1.0, 1);

//        uMom[d].FillBoundary(geom.periodicity());
//        MultiFABPhysBCDomainVel(uMom[d], d, geom, d);
//        MultiFABPhysBCMacVel(uMom[d], d, geom, d);
//    }

//    MkAdvMFluxdiv(umacNew,uMom,advFluxdivPred,dx,0);

//    // ADVANCE STEP (crank-nicolson + heun's method)

//    // Compute gmres_rhs

//    // trapezoidal advective terms
//    for (int d=0; d<AMREX_SPACEDIM; d++) {
//        advFluxdiv[d].mult(0.5, 1);
//        advFluxdivPred[d].mult(0.5, 1);
//    }

//    // crank-nicolson terms
//    StagApplyOp(geom,beta_negwtd, gamma_negwtd, beta_ed_negwtd, umac, Lumac, alpha_fc_0, dx, theta_alpha);

//    for (int d=0; d<AMREX_SPACEDIM; d++) {
//        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);

//        gmres_rhs_u[d].mult(dtinv, 1);
//        MultiFab::Add(gmres_rhs_u[d], mfluxdiv_correct[d], 0, 0, 1, 0);
//        MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 0);
//        MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 0);
//        MultiFab::Add(gmres_rhs_u[d], advFluxdivPred[d],   0, 0, 1, 0);

//        gmres_rhs_u[d].FillBoundary(geom.periodicity());

//        MultiFab::Subtract(gmres_rhs_u[d], pg[d], 0, 0, 1, 1);

//        // Only apply these BCs to the velocity term:
//        //MultiFABPhysBCDomainVel(gmres_rhs_u[d], d, geom, d);
//        //MultiFABPhysBCMacVel(gmres_rhs_u[d], d, geom, d);

//        MultiFab::Copy(umacNew[d], umac[d], 0, 0, 1, 0);
//    }

//    //pres.setVal(0.);  // initial guess
//    //SetPressureBC(pres, geom);

//    pres.setVal(0.);

//    // call GMRES here
//    GMRES(gmres_rhs_u, gmres_rhs_p, umacNew, pres,
//          alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
//          geom, norm_pre_rhs);

//    for (int i=0; i<AMREX_SPACEDIM; i++)
//        MultiFab::Copy(umac[i], umacNew[i], 0, 0, 1, 0);

//    for (int i=0; i<AMREX_SPACEDIM; i++) {
//        umac[i].FillBoundary(geom.periodicity());
//        MultiFABPhysBCDomainVel(umac[i], i, geom, i);
//        MultiFABPhysBCMacVel(umac[i], i, geom, i);
//    }

//    //////////////////////////////////////////////////

//}

// argv contains the name of the inputs file entered at the command line
void advanceLowMach(  std::array< MultiFab, AMREX_SPACEDIM >& umac,
	       std::array< MultiFab, AMREX_SPACEDIM >& umacNew,
	       MultiFab& pres, MultiFab& tracer,
	       const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_predict,
	       const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_correct,
               std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
	       const MultiFab& beta, const MultiFab& gamma,
	       const std::array< MultiFab, NUM_EDGE >& beta_ed,
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
      Lumac[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
      Lumac[d].setVal(0.);
  }

  // advective terms
  std::array< MultiFab, AMREX_SPACEDIM > advFluxdiv;
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
      advFluxdiv[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
      advFluxdiv[d].setVal(0.);
  }

  std::array< MultiFab, AMREX_SPACEDIM > advFluxdivPred;
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
      advFluxdivPred[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
      advFluxdivPred[d].setVal(0.);
  }

  // staggered momentum
  std::array< MultiFab, AMREX_SPACEDIM > uMom;
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
      uMom[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
      uMom[d].setVal(0.);
  } 

  MultiFab tracerPred(ba,dmap,1,1);
  MultiFab advFluxdivS(ba,dmap,1,1);

  ///////////////////////////////////////////
  // Scaled alpha, beta, gamma:
  ///////////////////////////////////////////

  // alpha_fc_0 arrays
  std::array< MultiFab, AMREX_SPACEDIM > alpha_fc_0;
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
      alpha_fc_0[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
      alpha_fc_0[d].setVal(0.);
  }

  // Scaled by 1/2:
  // beta_wtd cell centered
  MultiFab beta_wtd(ba, dmap, 1, 1);
  MultiFab::Copy(beta_wtd, beta, 0, 0, 1, 1);
  beta_wtd.mult(0.5, 1);

  // beta_wtd on nodes in 2d
  // beta_wtd on edges in 3d
  std::array< MultiFab, NUM_EDGE > beta_ed_wtd;
#if (AMREX_SPACEDIM == 2)
  beta_ed_wtd[0].define(convert(ba,nodal_flag), dmap, 1, 1);
  MultiFab::Copy(beta_ed_wtd[0], beta_ed[0], 0, 0, 1, 1);
  beta_ed_wtd[0].mult(0.5, 1);
#elif (AMREX_SPACEDIM == 3)
  beta_ed_wtd[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
  beta_ed_wtd[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
  beta_ed_wtd[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
  for(int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(beta_ed_wtd[d], beta_ed[d], 0, 0, 1, 1);
    beta_ed_wtd[d].mult(0.5, 1);
  }
#endif

  // cell-centered gamma_wtd
  MultiFab gamma_wtd(ba, dmap, 1, 1);
  MultiFab::Copy(gamma_wtd, gamma, 0, 0, 1, 1);
  gamma_wtd.mult(-0.5, 1);

  // Scaled by -1/2:
  // beta_negwtd cell centered
  MultiFab beta_negwtd(ba, dmap, 1, 1);
  MultiFab::Copy(beta_negwtd, beta, 0, 0, 1, 1);
  beta_negwtd.mult(-0.5, 1);

  // beta_negwtd on nodes in 2d
  // beta_negwtd on edges in 3d
  std::array< MultiFab, NUM_EDGE > beta_ed_negwtd;
#if (AMREX_SPACEDIM == 2)
  beta_ed_negwtd[0].define(convert(ba,nodal_flag), dmap, 1, 1);
  MultiFab::Copy(beta_ed_negwtd[0], beta_ed[0], 0, 0, 1, 1);
  beta_ed_negwtd[0].mult(-0.5, 1);
#elif (AMREX_SPACEDIM == 3)
  beta_ed_negwtd[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
  beta_ed_negwtd[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
  beta_ed_negwtd[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
  for(int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(beta_ed_negwtd[d], beta_ed[d], 0, 0, 1, 1);
    beta_ed_negwtd[d].mult(-0.5, 1);
  }
#endif

  // cell-centered gamma
  MultiFab gamma_negwtd(ba, dmap, 1, 1);
  MultiFab::Copy(gamma_negwtd, gamma, 0, 0, 1, 1);
  gamma_negwtd.mult(-0.5, 1);
  ///////////////////////////////////////////

  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    umac[d].FillBoundary(geom.periodicity());
  }

  //////////////////////////
  // Advance tracer
  //////////////////////////

  // Compute tracer:
  tracer.FillBoundary(geom.periodicity());
  MkAdvSFluxdiv(umac,tracer,advFluxdivS,dx,geom,0,0);
  advFluxdivS.mult(dt, 1);

  // compute predictor
  MultiFab::Copy(tracerPred, tracer, 0, 0, 1, 0);
  MultiFab::Add(tracerPred, advFluxdivS, 0, 0, 1, 0);
  tracerPred.FillBoundary(geom.periodicity());
  MkAdvSFluxdiv(umac,tracerPred,advFluxdivS,dx,geom,0,0);
  advFluxdivS.mult(dt, 1);

  // advance in time
  MultiFab::Add(tracer, tracerPred, 0, 0, 1, 0);
  MultiFab::Add(tracer, advFluxdivS, 0, 0, 1, 0);
  tracer.mult(0.5, 1);

  // amrex::Print() << "tracer L0 norm = " << tracer.norm0() << "\n";
  //////////////////////////

  //////////////////////////////////////////////////
  // ADVANCE velocity field
  //////////////////////////////////////////////////

  // PREDICTOR STEP (heun's method: part 1)
  // compute advective term
  AMREX_D_TERM(MultiFab::Copy(uMom[0], umac[0], 0, 0, 1, 0);,
	       MultiFab::Copy(uMom[1], umac[1], 0, 0, 1, 0);,
	       MultiFab::Copy(uMom[2], umac[2], 0, 0, 1, 0););

  // let rho = 1
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    uMom[d].mult(1.0, 1);
  }

  for (int d=0; d<AMREX_SPACEDIM; d++) {
    uMom[d].FillBoundary(geom.periodicity());
  }

  MkAdvMFluxdiv(umac,uMom,advFluxdiv,dx,0);

  // crank-nicolson terms
  StagApplyOp(geom,beta_negwtd,gamma_negwtd,beta_ed_negwtd,
	      umac,Lumac,alpha_fc_0,dx,theta_alpha);

  for (int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 0);
    
    gmres_rhs_u[d].mult(dtinv, 0);

    MultiFab::Add(gmres_rhs_u[d], mfluxdiv_predict[d], 0, 0, 1, 0);
    MultiFab::Add(gmres_rhs_u[d], Lumac[d], 0, 0, 1, 0);
    MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d], 0, 0, 1, 0);

    gmres_rhs_u[d].FillBoundary(geom.periodicity());
  }

  // initial guess for new solution
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(umacNew[d], umac[d], 0, 0, 1, 0);
  }
  pres.setVal(0.);  // initial guess

  // call GMRES to compute predictor
  GMRES(gmres_rhs_u,gmres_rhs_p,umacNew,pres,
	alpha_fc,beta_wtd,beta_ed_wtd,gamma_wtd,
	theta_alpha,geom,norm_pre_rhs);

  // Compute predictor advective term
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    umacNew[d].FillBoundary(geom.periodicity());
    MultiFab::Copy(uMom[d], umacNew[d], 0, 0, 1, 0);

    // let rho = 1
    uMom[d].mult(1.0, 1);

    uMom[d].FillBoundary(geom.periodicity());
  }

  MkAdvMFluxdiv(umacNew,uMom,advFluxdivPred,dx,0);

  // ADVANCE STEP (crank-nicolson + heun's method)

  // Compute gmres_rhs

  // trapezoidal advective terms
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    advFluxdiv[d].mult(0.5, 1);
    advFluxdivPred[d].mult(0.5, 1);
  }

  // crank-nicolson terms
  StagApplyOp(geom,beta_negwtd,gamma_negwtd,beta_ed_negwtd,
	      umac,Lumac,alpha_fc_0,dx,theta_alpha);

  for (int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 0);
  
    gmres_rhs_u[d].mult(dtinv);
  
    MultiFab::Add(gmres_rhs_u[d], mfluxdiv_correct[d],  0, 0, 1, 0);
    MultiFab::Add(gmres_rhs_u[d], Lumac[d],             0, 0, 1, 0);
    MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],        0, 0, 1, 0);
    MultiFab::Add(gmres_rhs_u[d], advFluxdivPred[d],    0, 0, 1, 0);

    gmres_rhs_u[d].FillBoundary(geom.periodicity());

    // initial guess for new solution
    MultiFab::Copy(umacNew[d], umac[d], 0, 0, 1, 0);
  }
  
  pres.setVal(0.);  // initial guess

  // call GMRES here
  GMRES(gmres_rhs_u,gmres_rhs_p,umacNew,pres,
	alpha_fc,beta_wtd,beta_ed_wtd,gamma_wtd,
	theta_alpha,geom,norm_pre_rhs);

  for (int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(umac[d], umacNew[d], 0, 0, 1, 0);
  }
  //////////////////////////////////////////////////

}
