
#include "hydro_functions.H"
#include "hydro_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "common_namespace.H"

#include "gmres_functions.H"
#include "gmres_functions_F.H"

#include "gmres_namespace.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;
using namespace common;
using namespace gmres;

// argv contains the name of the inputs file entered at the command line
void advance(  std::array< MultiFab, AMREX_SPACEDIM >& umac,
	       std::array< MultiFab, AMREX_SPACEDIM >& umacNew,
	       MultiFab& pres, MultiFab& tracer,
	       const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_predict,
	       const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_correct,
	       const std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
	       const MultiFab& beta, const MultiFab& gamma,
	       const std::array< MultiFab, NUM_EDGE >& beta_ed,
	       const Geometry geom, const Real& dt)
{

    BL_PROFILE_VAR("advance()",advance);

    const Real * dx  = geom.CellSize();
    const Real dtinv = 1.0/dt;
    Real theta_alpha = 1.;
    Real norm_pre_rhs;

    const BoxArray &              ba = beta.boxArray();
    const DistributionMapping & dmap = beta.DistributionMap();


    /****************************************************************************
     *                                                                          *
     * Define temporary MultiFabs                                               *
     *                                                                          *
     ***************************************************************************/

    // RHS pressure in GMRES
    MultiFab gmres_rhs_p(ba, dmap, 1, 1);
    gmres_rhs_p.setVal(0.);
    // Initial RHS pressure in GMRES
    MultiFab gmres_rhs_p_0(ba, dmap, 1, 1);
    gmres_rhs_p_0.setVal(0.);


    // RHS velocities in GMRES
    std::array< MultiFab, AMREX_SPACEDIM > gmres_rhs_u;
    // Velocity components updated by diffusion operator
    std::array< MultiFab, AMREX_SPACEDIM > Lumac;
    // Advective terms
    std::array< MultiFab, AMREX_SPACEDIM > advFluxdiv;
    // Advective terms (for predictor)
    std::array< MultiFab, AMREX_SPACEDIM > advFluxdivPred;
    // Staggered momentum
    std::array< MultiFab, AMREX_SPACEDIM > uMom;

    for (int i=0; i<AMREX_SPACEDIM; i++) {
           gmres_rhs_u[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                 Lumac[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
            advFluxdiv[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
        advFluxdivPred[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                  uMom[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
    }

    std::array< MultiFab, AMREX_SPACEDIM > pg;
    for (int i=0; i<AMREX_SPACEDIM; i++)
        pg[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);

    // Tracer concentration field for predictor
    MultiFab tracerPred(ba, dmap, 1, 1);

    // Tracer advective terms
    MultiFab advFluxdivS(ba, dmap, 1, 1);


    //___________________________________________________________________________
    // Scaled alpha, beta, gamma:

    // alpha_fc_0 arrays
    std::array< MultiFab, AMREX_SPACEDIM > alpha_fc_0;

    for (int i=0; i<AMREX_SPACEDIM; i++){
        alpha_fc_0[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
        alpha_fc_0[i].setVal(0.);
    }

    // Scaled by 1/2:
    // beta_wtd cell centered
    MultiFab beta_wtd(ba, dmap, 1, 1);
    MultiFab::Copy(beta_wtd, beta, 0, 0, 1, 1);
    beta_wtd.mult(0.5, 1);

    // beta_wtd on nodes in 2d, on edges in 3d
    std::array< MultiFab, NUM_EDGE > beta_ed_wtd;
#if (AMREX_SPACEDIM == 2)
    beta_ed_wtd[0].define(convert(ba,nodal_flag), dmap, 1, 1);
    MultiFab::Copy(beta_ed_wtd[0], beta_ed[0], 0, 0, 1, 1);
    beta_ed_wtd[0].mult(0.5, 1);
#elif (AMREX_SPACEDIM == 3)
    for(int d=0; d<AMREX_SPACEDIM; d++) {
        beta_ed_wtd[d].define(convert(ba, nodal_flag_edge[d]), dmap, 1, 1);

        MultiFab::Copy(beta_ed_wtd[d], beta_ed[d], 0, 0, 1, 1);
        beta_ed_wtd[d].mult(0.5, 1);
    }
#endif

    // Scaled by 1/2:
    // gamma_wtd cell centered
    MultiFab gamma_wtd(ba, dmap, 1, 1);
    MultiFab::Copy(gamma_wtd, gamma, 0, 0, 1, 1);
    gamma_wtd.mult(-0.5, 1);

    // Scaled by -1/2:
    // beta_negwtd cell centered
    MultiFab beta_negwtd(ba, dmap, 1, 1);
    MultiFab::Copy(beta_negwtd, beta, 0, 0, 1, 1);
    beta_negwtd.mult(-0.5, 1);

    // beta_negwtd on nodes in 2d, on edges in 3d
    std::array< MultiFab, NUM_EDGE > beta_ed_negwtd;
#if (AMREX_SPACEDIM == 2)
    beta_ed_negwtd[0].define(convert(ba,nodal_flag), dmap, 1, 1);
    MultiFab::Copy(beta_ed_negwtd[0], beta_ed[0], 0, 0, 1, 1);
    beta_ed_negwtd[0].mult(-0.5, 1);
#elif (AMREX_SPACEDIM == 3)
    for(int d=0; d<AMREX_SPACEDIM; d++) {
        beta_ed_negwtd[d].define(convert(ba,nodal_flag_edge[d]), dmap, 1, 1);

        MultiFab::Copy(beta_ed_negwtd[d], beta_ed[d], 0, 0, 1, 1);
        beta_ed_negwtd[d].mult(-0.5, 1);
    }
#endif

    // Scaled by -1/2:
    // gamma cell centered
    MultiFab gamma_negwtd(ba, dmap, 1, 1);
    MultiFab::Copy(gamma_negwtd, gamma, 0, 0, 1, 1);
    gamma_negwtd.mult(-0.5, 1);


    /****************************************************************************
     *                                                                          *
     * Apply non-stochastic boundary conditions                                 *
     *                                                                          *
     ***************************************************************************/

    for (int i=0; i<AMREX_SPACEDIM; i++) {
        umac[i].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(umac[i], i, geom);
        MultiFABPhysBCMacVel(umac[i], i, geom);
    }


    /****************************************************************************
     *                                                                          *
     * Advance tracer                                                           *
     *                                                                          *
     ***************************************************************************/

    // Compute tracer:
    tracer.FillBoundary(geom.periodicity());
    MultiFABPhysBC(tracer, geom);

    MkAdvSFluxdiv(umac, tracer, advFluxdivS, dx, geom, 0);
    advFluxdivS.mult(dt, 1);

    // compute predictor
    MultiFab::Copy(tracerPred, tracer, 0, 0, 1, 0);
    MultiFab::Add(tracerPred, advFluxdivS, 0, 0, 1, 0);

    tracerPred.FillBoundary(geom.periodicity());
    MultiFABPhysBC(tracerPred, geom);

    MkAdvSFluxdiv(umac, tracerPred, advFluxdivS, dx, geom, 0);
    advFluxdivS.mult(dt, 1);

    // advance in time
    MultiFab::Add(tracer, tracerPred,  0, 0, 1, 0);
    MultiFab::Add(tracer, advFluxdivS, 0, 0, 1, 0);
    tracer.mult(0.5, 1);

    // amrex::Print() << "tracer L0 norm = " << tracer.norm0() << "\n";
    //////////////////////////

    //////////////////////////////////////////////////
    // ADVANCE velocity field
    //////////////////////////////////////////////////

    std::array<MultiFab, AMREX_SPACEDIM> umac_0;
    std::array<MultiFab, AMREX_SPACEDIM> umac_1;
    std::array<MultiFab, AMREX_SPACEDIM> alpha_fc_1;
    std::array<MultiFab, AMREX_SPACEDIM> alpha_f;
    std::array<MultiFab, AMREX_SPACEDIM> fcoef;
    std::array<MultiFab, AMREX_SPACEDIM> force;
    std::array<MultiFab, AMREX_SPACEDIM> force_0;
    std::array<MultiFab, AMREX_SPACEDIM> force_1;
    std::array<MultiFab, AMREX_SPACEDIM> tmp_f_0;
    std::array<MultiFab, AMREX_SPACEDIM> tmp_f_1;
    std::array<MultiFab, AMREX_SPACEDIM> r_f;
    std::array<MultiFab, AMREX_SPACEDIM> tmp_f_mask;
    std::array<MultiFab, AMREX_SPACEDIM> force_est;
    std::array<MultiFab, AMREX_SPACEDIM> force_rhs;

    for (int d=0; d<AMREX_SPACEDIM; d++) {

        umac_0[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        umac_0[d].setVal(0.);

        umac_1[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        umac_1[d].setVal(0.);

        alpha_fc_1[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        alpha_fc_1[d].setVal(0.);

        alpha_f[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        alpha_f[d].setVal(dtinv);

        fcoef[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        fcoef[d].setVal(0.);

        force[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        force[d].setVal(0.);

        force_0[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        force_0[d].setVal(0.);

        force_1[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        force_1[d].setVal(0.);

        tmp_f_0[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        tmp_f_0[d].setVal(0.);

        tmp_f_1[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        tmp_f_1[d].setVal(0.);

        tmp_f_mask[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        tmp_f_mask[d].setVal(0.);

        r_f[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        r_f[d].setVal(0.);

        force_est[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        force_est[d].setVal(0.);

        force_rhs[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        force_rhs[d].setVal(0.);
    }

    MultiFab p_f(ba, dmap, 1, 1);
    p_f.setVal(0.);

    MultiFab p_0(ba, dmap, 1, 1);
    p_0.setVal(0.);

    MultiFab p_1(ba, dmap, 1, 1);
    p_1.setVal(0.);




    int tr_min = 14;
    int tr_max = 18;
    IntVect test_region_min{AMREX_D_DECL(tr_min, tr_min, tr_min)};

    std::array<Box, AMREX_SPACEDIM> test_regions;

    test_regions[0] = Box(test_region_min,
                          IntVect{AMREX_D_DECL(tr_max+1, tr_max, tr_max)});

    test_regions[1] = Box(test_region_min,
                          IntVect{AMREX_D_DECL(tr_max, tr_max+1, tr_max)});

#if (AMREX_SPACEDIM == 3)
    test_regions[2] = Box(test_region_min,
                          IntVect{AMREX_D_DECL(tr_max, tr_max, tr_max+1)});
#endif

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        fcoef[d].setVal(1e0, test_regions[d], 0, 1);
        fcoef[d].FillBoundary(geom.periodicity());
    }


    // PREDICTOR STEP (heun's method: part 1)
    // compute advective term

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(uMom[d], umac[d], 0, 0, 1, 1);

        // let rho = 1
        uMom[d].mult(1.0, 1);

        uMom[d].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(uMom[d], d, geom);
        MultiFABPhysBCMacVel(uMom[d], d, geom);
    }

    MkAdvMFluxdiv(umac, uMom, advFluxdiv, dx, 0);

    // crank-nicolson terms
    StagApplyOp(beta_negwtd, gamma_negwtd, beta_ed_negwtd, umac, Lumac, alpha_fc_0,
                dx, theta_alpha);

    // pressure boundary conditions
    pres.setVal(0.);
    SetPressureBC(pres, geom);
    ComputeGrad(pres, pg, 0, 0, 1, geom);

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(alpha_fc_1[d], fcoef[d], 0, 0, 1, 1);
        MultiFab::Multiply(alpha_fc_1[d], alpha_f[d], 0, 0, 1, 1);
        MultiFab::Add(alpha_fc_1[d], alpha_fc[d], 0, 0, 1, 1);
    }

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(umac_0[d], umac[d], 0, 0, 1, 1);
        MultiFab::Copy(umac_1[d], umac[d], 0, 0, 1, 1);
    }

    MultiFab::Copy(p_0, pres, 0, 0, 1, 1);
    MultiFab::Copy(p_1, pres, 0, 0, 1, 1);

    for (int i=0; i<20; i++) {

        //gmres_rhs_p.setVal(0.);

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            Lumac[d].FillBoundary(geom.periodicity());

            MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);
            gmres_rhs_u[d].mult(dtinv, 1);

            MultiFab::Add(gmres_rhs_u[d], mfluxdiv_predict[d], 0, 0, 1, 0);
            MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 0);
            MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 0);
            MultiFab::Add(gmres_rhs_u[d], force_0[d],          0, 0, 1, 0);

            pg[d].FillBoundary(geom.periodicity());
            gmres_rhs_u[d].FillBoundary(geom.periodicity());

            // add pressure
            MultiFab::Subtract(gmres_rhs_u[d], pg[d], 0, 0, 1, 1);

            // initial guess for new solution
            MultiFab::Copy(umacNew[d], umac_0[d], 0, 0, 1, 1);
        }

        // call GMRES to compute predictor
        GMRES(gmres_rhs_u, gmres_rhs_p, umacNew, p_0,
              alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
              geom, norm_pre_rhs);

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Copy(umac_0[d], umacNew[d], 0, 0, 1, 1);
        }

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Copy    (r_f[d], fcoef[d],   0, 0, 1, 0);
            MultiFab::Multiply(r_f[d], umacNew[d], 0, 0, 1, 0);
            r_f[d].mult(-1., 0);

            MultiFab::Copy(force_est[d], r_f[d], 0, 0, 1, 0);
            force_est[d].mult(dtinv, 0);
        }

        // Inverse Motility Matrix
        ApplyMatrix(tmp_f_0, p_f, r_f, gmres_rhs_p_0,
                    alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha, geom);

        MultiFab::Add(gmres_rhs_p, p_f, 0, 0, 1, 0);

        std::array<Real, AMREX_SPACEDIM> tmp_f_resid;
        for (int d=0; d<AMREX_SPACEDIM; ++d){
            MultiFab::Copy(tmp_f_mask[d], tmp_f_0[d],   0, 0, 1, 1);
            MultiFab::Multiply(tmp_f_mask[d], fcoef[d], 0, 0, 1, 1);

            // MultiFab::Add(force_0[d], tmp_f_mask[d],    0, 0, 1, 0);
            MultiFab::Add(force_0[d], tmp_f_0[d],    0, 0, 1, 0);
        }


        // for (int d=0; d<AMREX_SPACEDIM; ++d)
        //     MultiFab::Add(force_0[d], force_est[d], 0, 0, 1, 0);


        // Compute predictor advective term
        // let rho = 1
        for (int d=0; d<AMREX_SPACEDIM; d++) {
            umacNew[d].FillBoundary(geom.periodicity());
            MultiFABPhysBCDomainVel(umacNew[d], d, geom);
            MultiFABPhysBCMacVel(umacNew[d], d, geom);

            MultiFab::Copy(uMom[d], umacNew[d], 0, 0, 1, 0);
            uMom[d].mult(1.0, 1);

            uMom[d].FillBoundary(geom.periodicity());
            MultiFABPhysBCDomainVel(uMom[d], d, geom);
            MultiFABPhysBCMacVel(uMom[d], d, geom);
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
        // TODO: ask Andy if we should use umacNew here?
        StagApplyOp(beta_negwtd, gamma_negwtd, beta_ed_negwtd, umac, Lumac,
                    alpha_fc_0, dx, theta_alpha);

        for (int d=0; d<AMREX_SPACEDIM; d++) {
            MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);

            gmres_rhs_u[d].mult(dtinv, 1);
            MultiFab::Add(gmres_rhs_u[d], mfluxdiv_correct[d], 0, 0, 1, 0);
            MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 0);
            MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 0);
            MultiFab::Add(gmres_rhs_u[d], advFluxdivPred[d],   0, 0, 1, 0);
            MultiFab::Add(gmres_rhs_u[d], force_1[d],          0, 0, 1, 0);

            pg[d].FillBoundary(geom.periodicity());
            gmres_rhs_u[d].FillBoundary(geom.periodicity());

            // add pressure
            MultiFab::Subtract(gmres_rhs_u[d], pg[d], 0, 0, 1, 1);

            // initial guess for new solution
            MultiFab::Copy(umacNew[d], umac_1[d], 0, 0, 1, 0);
        }

        // call GMRES here
        GMRES(gmres_rhs_u, gmres_rhs_p, umacNew, p_1,
              alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
              geom, norm_pre_rhs);

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Copy(umac_1[d], umacNew[d], 0, 0, 1, 1);
        }

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            //MultiFab::Copy    (r_f[d], fcoef[d],   0, 0, 1, 0);
            MultiFab::Copy    (r_f[d], umacNew[d], 0, 0, 1, 0);
            MultiFab::Multiply(r_f[d], fcoef[d],   0, 0, 1, 0);
            r_f[d].mult(-1., 0);

            MultiFab::Copy(force_est[d], r_f[d], 0, 0, 1, 0);
            force_est[d].mult(dtinv, 0);

            // MultiFab::Add(r_f[d], force_est[d],  0, 0, 1, 0);
        }

        // Inverse Motility Matrix
        ApplyMatrix(tmp_f_1, p_f, r_f, gmres_rhs_p_0,
                    alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha, geom);

        MultiFab::Add(gmres_rhs_p, p_f, 0, 0, 1, 0);

        for (int d=0; d<AMREX_SPACEDIM; ++d){
            MultiFab::Copy(tmp_f_mask[d], tmp_f_1[d],   0, 0, 1, 1);
            MultiFab::Multiply(tmp_f_mask[d], fcoef[d], 0, 0, 1, 1);

            // MultiFab::Add(force_1[d], tmp_f_mask[d],    0, 0, 1, 0);
            MultiFab::Add(force_1[d], tmp_f_1[d],    0, 0, 1, 0);

            // force_0[d].mult(0.5, 0);
            // force_1[d].mult(0.5, 0);
            // MultiFab::Add(force[d], force_0[d], 0, 0, 1, 0);
            // MultiFab::Add(force[d], force_1[d], 0, 0, 1, 0);
        }

        // for (int d=0; d<AMREX_SPACEDIM; ++d)
        //     MultiFab::Add(force_1[d], force_est[d], 0, 0, 1, 0);

        Print() << "corrector step: " << i
                << " force = " << force_1[0].norm0()
                << ", "        << force_1[1].norm0()
                << ", "        << force_1[2].norm0()
                << " ";
        Print() << "r_f = "    << r_f[0].sum()
                << ", "        << r_f[1].sum()
                << ", "        << r_f[2].sum()
                << std::endl;

        Real norm_r, norm_fest, norm_fmask;
        StagL2Norm(r_f, 0, norm_r);
        StagL2Norm(force_est, 0, norm_fest);
        StagL2Norm(tmp_f_mask, 0, norm_fmask);
        Print() << "corrector norm_resid = "   << norm_r
                << " norm f_est = "  << norm_fest
                << " norm f_mask = " << norm_fmask << std::endl;
    }

    for (int i=0; i<AMREX_SPACEDIM; i++)
        MultiFab::Copy(umac[i], umacNew[i], 0, 0, 1, 0);

    //////////////////////////////////////////////////

}
