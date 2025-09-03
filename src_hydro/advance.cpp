#include "hydro_functions.H"

#include "common_functions.H"


#include "gmres_functions.H"


#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

void advanceStokes(std::array< MultiFab, AMREX_SPACEDIM >& umac,
                   MultiFab& pres,
                   const std::array< MultiFab, AMREX_SPACEDIM >& stochMfluxdiv,
                   std::array< MultiFab, AMREX_SPACEDIM >& sourceTerms,
                   std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
                   MultiFab& beta,
                   MultiFab& gamma,
                   std::array< MultiFab, NUM_EDGE >& beta_ed,
                   const Geometry geom, const Real& dt)
{
    BL_PROFILE_VAR("advanceStokes()",advance);

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



    if (zero_net_force == 1)
    {
        Vector<Real> mean_stress_umac(AMREX_SPACEDIM);

        SumStag(gmres_rhs_u,mean_stress_umac,true);

        Print() << "correcting mean force: " << mean_stress_umac[0]*(prob_hi[0]-prob_lo[0])*(prob_hi[1]-prob_lo[1])*(prob_hi[2]-prob_lo[2]);

        for (int d=1; d<AMREX_SPACEDIM; ++d) {
            Print() << ", " << mean_stress_umac[d]*(prob_hi[0]-prob_lo[0])*(prob_hi[1]-prob_lo[1])*(prob_hi[2]-prob_lo[2]);
        }

        Print() << "\n";

        Print() << "test force: " << sourceTerms[0].sum()*(prob_hi[0]-prob_lo[0])*(prob_hi[1]-prob_lo[1])*(prob_hi[2]-prob_lo[2]);

        for (int d=1; d<AMREX_SPACEDIM; ++d) {
            Print() << ", " << sourceTerms[d].sum()*(prob_hi[0]-prob_lo[0])*(prob_hi[1]-prob_lo[1])*(prob_hi[2]-prob_lo[2]);
        }

        Print() << "\n";

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            if (geom.isPeriodic(d)) {
                gmres_rhs_u[d].plus(-mean_stress_umac[d],0,1,0);
            }
        }
    }

    // call GMRES
    GMRES gmres(ba,dmap,geom);
    gmres.Solve(gmres_rhs_u,gmres_rhs_p,umac,pres,
                alpha_fc,beta,beta_ed,gamma,theta_alpha,geom,norm_pre_rhs);

    for (int i=0; i<AMREX_SPACEDIM; i++) {
        MultiFabPhysBCDomainVel(umac[i], geom, i);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umac[i], geom, i, is_inhomogeneous);
        umac[i].FillBoundary(geom.periodicity());
    }
}

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
        MultiFabPhysBCDomainVel(umac[d], geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umac[d], geom, d, is_inhomogeneous);
        umac[d].FillBoundary(geom.periodicity());
    }

    //////////////////////////
    // Advance tracer
    //////////////////////////

    // Compute tracer:
    tracer.FillBoundary(geom.periodicity());
    MkAdvSFluxdiv_cc(umac,tracer,advFluxdivS,geom,0,1,0);
    advFluxdivS.mult(dt, 1);

    // compute predictor
    MultiFab::Copy(tracerPred, tracer, 0, 0, 1, 0);
    MultiFab::Add(tracerPred, advFluxdivS, 0, 0, 1, 0);
    tracerPred.FillBoundary(geom.periodicity());
    MkAdvSFluxdiv_cc(umac,tracerPred,advFluxdivS,geom,0,1,0);
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
        MultiFabPhysBCDomainVel(uMom[d], geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(uMom[d], geom, d, is_inhomogeneous);
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
    GMRES gmres(ba,dmap,geom);
    gmres.Solve(gmres_rhs_u,gmres_rhs_p,umacNew,pres,
                alpha_fc,beta_wtd,beta_ed_wtd,gamma_wtd,
                theta_alpha,geom,norm_pre_rhs);

    // Compute predictor advective term
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        umacNew[d].FillBoundary(geom.periodicity());
        MultiFab::Copy(uMom[d], umacNew[d], 0, 0, 1, 0);

        // let rho = 1
        uMom[d].mult(1.0, 1);

        MultiFabPhysBCDomainVel(uMom[d], geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(uMom[d], geom, d, is_inhomogeneous);
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
    gmres.Solve(gmres_rhs_u,gmres_rhs_p,umacNew,pres,
                alpha_fc,beta_wtd,beta_ed_wtd,gamma_wtd,
                theta_alpha,geom,norm_pre_rhs);

    for (int d=0; d<AMREX_SPACEDIM; d++) {
        MultiFab::Copy(umac[d], umacNew[d], 0, 0, 1, 0);
    }
    //////////////////////////////////////////////////

}