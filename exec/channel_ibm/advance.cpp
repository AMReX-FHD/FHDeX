#include "main_driver.H"
#include "main_driver_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "hydro_functions.H"

#include "common_namespace.H"

#include "gmres_functions.H"
#include "gmres_functions_F.H"

#include "gmres_namespace.H"

#include "IBCore.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;
using namespace common;
using namespace gmres;

// argv contains the name of the inputs file entered at the command line
void advance(std::array<MultiFab, AMREX_SPACEDIM> & umac,
             std::array<MultiFab, AMREX_SPACEDIM> & umacNew,
             MultiFab & pres, MultiFab & tracer,
             std::array<MultiFab, AMREX_SPACEDIM> & force_ibm,
             const std::array<MultiFab, AMREX_SPACEDIM> & mfluxdiv_predict,
             const std::array<MultiFab, AMREX_SPACEDIM> & mfluxdiv_correct,
             const std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
             const MultiFab & beta, const MultiFab & gamma,
             const std::array<MultiFab, NUM_EDGE> & beta_ed,
             IBCore & ib_core,
             const Geometry geom, const Real & dt)
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



    /****************************************************************************
     *                                                                          *
     * Generate immersed boundary data                                          *
     *                                                                          *
     ***************************************************************************/

    // MultiFabs storing Immersed-Boundary data
    std::array<MultiFab, AMREX_SPACEDIM> f_ibm;
    std::array<MultiFab, AMREX_SPACEDIM> vel_d;
    std::array<MultiFab, AMREX_SPACEDIM> vel_g;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_ibm[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        vel_d[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        vel_g[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
    }


    // MultiFab tagging the alpha_fc array
    std::array<MultiFab, AMREX_SPACEDIM> alpha_fc_1;

    for (int d=0; d<AMREX_SPACEDIM; ++d)
        alpha_fc_1[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);


    // Gradient of p, where div(grad(p)) = div(r), where r is the slip velocity
    // due to the unconstrained stokes problem.
    std::array<MultiFab, AMREX_SPACEDIM> tmp_ibm_f;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        tmp_ibm_f[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        tmp_ibm_f[d].setVal(0.);
    }



    //___________________________________________________________________________
    // Build coefficients for IB spreading and interpolation operators

    int lev_ib = 0;

    ib_core.ImplicitDeposition(f_ibm[0], f_ibm[1], f_ibm[2],
                               vel_d[0], vel_d[1], vel_d[2],
                               vel_g[0], vel_g[1], vel_g[2],
                               lev_ib, dt);

    for (int d=0; d<AMREX_SPACEDIM; ++d)
        f_ibm[d].FillBoundary(geom.periodicity());


    //___________________________________________________________________________
    // For debug purposes: Write out the immersed-boundary data

    VisMF::Write(f_ibm[0], "ib_data/f_u");
    VisMF::Write(f_ibm[1], "ib_data/f_v");
    VisMF::Write(f_ibm[2], "ib_data/f_w");

    VisMF::Write(vel_d[0], "ib_data/u_d");
    VisMF::Write(vel_d[1], "ib_data/v_d");
    VisMF::Write(vel_d[2], "ib_data/w_d");

    VisMF::Write(vel_g[0], "ib_data/u_g");
    VisMF::Write(vel_g[1], "ib_data/v_g");
    VisMF::Write(vel_g[2], "ib_data/w_g");


    //___________________________________________________________________________
    // Tag alpha_fc to solve implicit immersed boundary force

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(alpha_fc_1[d], f_ibm[d], 0, 0, 1, 1);
        // 1e6 is a magic number (for now)
        // alpha_fc_1[d].mult(1.e6, 1);
        alpha_fc_1[d].mult(0, 1);
        // Tag implicit force terms
        MultiFab::Add(alpha_fc_1[d], alpha_fc[d], 0, 0, 1, 1);
    }



    /****************************************************************************
     *                                                                          *
     * ADVANCE velocity field                                                   *
     *                                                                          *
     ***************************************************************************/

    // Velocities used by GMRES to solve for predictor (0) and corrector (1)
    std::array<MultiFab, AMREX_SPACEDIM> umac_0;
    std::array<MultiFab, AMREX_SPACEDIM> umac_1;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umac_0[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        umac_1[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
    }


    // Explicit force terms used by GMRES to solve for predictor (0) and
    // corrector (1)
    std::array<MultiFab, AMREX_SPACEDIM> force_0;
    std::array<MultiFab, AMREX_SPACEDIM> force_1;

    for (int d=0; d<AMREX_SPACEDIM; d++) {
        force_0[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        force_1[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
    }


    std::array<MultiFab, AMREX_SPACEDIM> tmp_f_0;
    std::array<MultiFab, AMREX_SPACEDIM> tmp_f_1;
    std::array<MultiFab, AMREX_SPACEDIM> r_f;
    std::array<MultiFab, AMREX_SPACEDIM> tmp_r_f;
    std::array<MultiFab, AMREX_SPACEDIM> tmp_f_mask;
    std::array<MultiFab, AMREX_SPACEDIM> force_est;


    for (int d=0; d<AMREX_SPACEDIM; d++) {

        tmp_f_0[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        tmp_f_0[d].setVal(0.);

        tmp_f_1[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        tmp_f_1[d].setVal(0.);

        tmp_f_mask[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        tmp_f_mask[d].setVal(0.);

        r_f[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        r_f[d].setVal(0.);

        tmp_r_f[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        tmp_r_f[d].setVal(0.);

        force_est[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        force_est[d].setVal(0.);
    }

    MultiFab p_f(ba, dmap, 1, 1);
    p_f.setVal(0.);

    MultiFab tmp_pf(ba, dmap, 1, 1);
    tmp_pf.setVal(0.);

    MultiFab p_0(ba, dmap, 1, 1);
    p_0.setVal(0.);

    MultiFab p_1(ba, dmap, 1, 1);
    p_1.setVal(0.);

    MultiFab p_ibm_0(ba, dmap, 1, 1);
    p_ibm_0.setVal(0.);

    MultiFab p_ibm_1(ba, dmap, 1, 1);
    p_ibm_1.setVal(0.);

    MultiFab zeros_cc(ba, dmap, 1, 1);
    zeros_cc.setVal(0.);

    std::array< MultiFab, AMREX_SPACEDIM > ones_fc;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        ones_fc[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        ones_fc[d].setVal(1.);
    }



    /****************************************************************************
     *                                                                          *
     * PREDICTOR STEP (heun's method: part 1) for fluid                         *
     *                                                                          *
     ***************************************************************************/

    //___________________________________________________________________________
    // Compute advective term

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(uMom[d], umac[d], 0, 0, 1, 1);

        // let rho = 1
        uMom[d].mult(1.0, 1);

        // momentum boundary conditions
        uMom[d].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(uMom[d], d, geom);
        MultiFABPhysBCMacVel(uMom[d], d, geom);
    }

    // advective momentum flux terms
    MkAdvMFluxdiv(umac, uMom, advFluxdiv, dx, 0);

    // crank-nicolson terms
    StagApplyOp(beta_negwtd, gamma_negwtd, beta_ed_negwtd, umac, Lumac, alpha_fc_0,
                dx, theta_alpha);

    // advective term boundary conditions
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        Lumac[d].FillBoundary(geom.periodicity());
        MultiFABPhysBC(Lumac[d], d, geom);

        advFluxdiv[d].FillBoundary(geom.periodicity());
        MultiFABPhysBC(advFluxdiv[d], d, geom);
    }


    //___________________________________________________________________________
    // Pressure term boundary conditions

    pres.setVal(0.);
    SetPressureBC(pres, geom);
    ComputeGrad(pres, pg, 0, 0, 1, geom);


    //___________________________________________________________________________
    // Set up initial condtions for predictor (0) and corrector (1)

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(umac_0[d],       umac[d], 0, 0, 1, 1);
        MultiFab::Copy(umac_1[d],       umac[d], 0, 0, 1, 1);
        MultiFab::Copy(force_0[d], force_ibm[d], 0, 0, 1, 1);
        MultiFab::Copy(force_1[d], force_ibm[d], 0, 0, 1, 1);
    }

    MultiFab::Copy(p_0, pres, 0, 0, 1, 1);
    MultiFab::Copy(p_1, pres, 0, 0, 1, 1);


    //___________________________________________________________________________
    // Set up the RHS for the predictor

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        // explicit part
        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);
        gmres_rhs_u[d].mult(dtinv, 1);

        MultiFab::Add(gmres_rhs_u[d], mfluxdiv_predict[d], 0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], force_0[d],          0, 0, 1, 1);

        // fill boundary before adding pressure part to prevent it from
        // overwriding any pressure gradients in the ghost cells
        gmres_rhs_u[d].FillBoundary(geom.periodicity());

        // add pressure
        pg[d].FillBoundary(geom.periodicity());
        MultiFab::Subtract(gmres_rhs_u[d], pg[d], 0, 0, 1, 1);

        // initial guess for new solution
        MultiFab::Copy(umacNew[d], umac_0[d], 0, 0, 1, 1);
    }

    //___________________________________________________________________________
    // Call GMRES to compute predictor

    GMRES(gmres_rhs_u, gmres_rhs_p, umacNew, p_0,
          alpha_fc_1, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
          geom, norm_pre_rhs);

    for (int d=0; d<AMREX_SPACEDIM; ++d)
        MultiFab::Copy(umac_0[d], umacNew[d], 0, 0, 1, 1);



    /****************************************************************************
     *                                                                          *
     * ADVANCE (CORRECTOR) STEP (crank-nicolosn heun's method: part 2)          *
     *                                                                          *
     ***************************************************************************/

    //___________________________________________________________________________
    // Compute corrector's advective term (using predictor's fluid solution)

     for (int d=0; d<AMREX_SPACEDIM; d++) {
        umacNew[d].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(umacNew[d], d, geom);
        MultiFABPhysBCMacVel(umacNew[d], d, geom);

        MultiFab::Copy(uMom[d], umacNew[d], 0, 0, 1, 1);

        // let rho = 1
        uMom[d].mult(1.0, 1);

        // momentum boundary conditions
        uMom[d].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(uMom[d], d, geom);
        MultiFABPhysBCMacVel(uMom[d], d, geom);
     }

    // advective momentum flux terms
    MkAdvMFluxdiv(umacNew, uMom, advFluxdivPred, dx, 0);

    // trapezoidal advective terms
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        advFluxdiv[d].mult(0.5, 1);
        advFluxdivPred[d].mult(0.5, 1);
    }

    // crank-nicolson terms
    // TODO: ask Andy if we should use umacNew here?
    StagApplyOp(beta_negwtd, gamma_negwtd, beta_ed_negwtd, umac, Lumac,
                alpha_fc_0, dx, theta_alpha);


    // advective term boundary conditions
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        Lumac[d].FillBoundary(geom.periodicity());
        MultiFABPhysBC(Lumac[d], d, geom);

        advFluxdiv[d].FillBoundary(geom.periodicity());
        MultiFABPhysBC(advFluxdiv[d], d, geom);

        advFluxdivPred[d].FillBoundary(geom.periodicity());
        MultiFABPhysBC(advFluxdivPred[d], d, geom);
    }


    //___________________________________________________________________________
    // Solve IBM constraint force by iterating
    // NOTE: 100 max attempts at getting the residual down, TODO: pass this as a
    // parameter (or tie to MAX ITER)

    for (int i=0; i<100; i++) {

        //_______________________________________________________________________
        // Set up the RHS for the corrector 

        for (int d=0; d<AMREX_SPACEDIM; d++) {
            // explicit part
            MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);
            gmres_rhs_u[d].mult(dtinv, 1);

            MultiFab::Add(gmres_rhs_u[d], mfluxdiv_correct[d], 0, 0, 1, 1);
            MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 1);
            MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 1);
            MultiFab::Add(gmres_rhs_u[d], advFluxdivPred[d],   0, 0, 1, 1);
            MultiFab::Add(gmres_rhs_u[d], force_1[d],          0, 0, 1, 1);

            // fill boundary before adding pressure part to prevent it from
            // overwriding any pressure gradients in the ghost cells
            gmres_rhs_u[d].FillBoundary(geom.periodicity());

            // add pressure
            pg[d].FillBoundary(geom.periodicity());
            MultiFab::Subtract(gmres_rhs_u[d], pg[d], 0, 0, 1, 1);

            // initial guess for new solution
            MultiFab::Copy(umacNew[d], umac_1[d], 0, 0, 1, 1);
        }

        //_______________________________________________________________________
        // call GMRES to compute corrector

        GMRES(gmres_rhs_u, gmres_rhs_p, umacNew, p_1,
              alpha_fc_1, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
              geom, norm_pre_rhs);


        //_______________________________________________________________________
        // Slip velocity term

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Copy    (r_f[d], umacNew[d], 0, 0, 1, 1);
            MultiFab::Multiply(r_f[d], f_ibm[d],   0, 0, 1, 1);
            r_f[d].mult(-1., 0);
        }


        //_______________________________________________________________________
        // Invert motility matrix

        // initial guess
        p_ibm_1.setVal(0.);
        for (int d=0; d<AMREX_SPACEDIM; ++d)
            tmp_ibm_f[d].setVal(0.);

        // Inverse Motility Matrix
        ApplyMatrix(tmp_f_1, p_f, r_f, p_ibm_1,
                    alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha, geom);

        // // Remove non-divergence free parts of the residual
        // p_f.mult(-1., 0);
        // MultiFABPhysBC(p_f, geom);

        // MacProj(ones_fc, p_f, tmp_pf, geom, 1);
        // tmp_pf.FillBoundary(geom.periodicity());
        // MultiFABPhysBC(tmp_pf, geom);

        // SubtractWeightedGradP(tmp_ibm_f, ones_fc, tmp_pf, geom);
        // MultiFab::Add(p_ibm_1, tmp_pf, 0, 0, 1, 1);

        MultiFab::Add(gmres_rhs_p, p_f, 0, 0, 1, 1);
        gmres_rhs_p.FillBoundary(geom.periodicity());


        //______________________________________________________________________
        // Apply immersed-boundary force

        // Inverse-motility part
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Copy(tmp_f_mask[d],     tmp_f_1[d], 0, 0, 1, 1);
            MultiFab::Multiply(tmp_f_mask[d], f_ibm[d],   0, 0, 1, 1);

            // MultiFab::Add(force_1[d], tmp_f_mask[d],    0, 0, 1, 1);
            // Add raw force to fluid => includes pressure correction terms from
            // immersed boundary => this might not be the right thing to do...
            MultiFab::Add(force_1[d], tmp_f_1[d],    0, 0, 1, 1);
        }

        // // Divergence part
        // for (int d=0; d<AMREX_SPACEDIM; ++d) {
        //     MultiFab::Subtract(umacNew[d], tmp_ibm_f[d], 0, 0, 1, 1);
        //     MultiFab::Copy(umac_1[d], umacNew[d], 0, 0, 1, 1);
        // }

        // Implicit (move to lhs) part
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Copy(force_est[d], r_f[d], 0, 0, 1, 1);
            force_est[d].mult(1e6, 0);
 
            //MultiFab::Add(force_1[d], force_est[d], 0, 0, 1, 1);
        }


        //_______________________________________________________________________
        // Updated slip velocity as error estimate

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            r_f[d].setVal(0.);
            MultiFab::Copy    (r_f[d], umacNew[d], 0, 0, 1, 1);
            MultiFab::Multiply(r_f[d], f_ibm[d],   0, 0, 1, 1);
        }


        Real norm_r;
        StagL2Norm(r_f, 0, norm_r);
        Print() << "step: " << i << " norm_resid = " << norm_r << std::endl;


        if (norm_r < 1e-12) break;
    }

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        // Output velocity solution
        MultiFab::Copy(umac[d],      umacNew[d],   0, 0, 1, 1);

        // Output immersed-boundary forces
        MultiFab::Copy(force_ibm[d], force_1[d],   0, 0, 1, 1);
        // Include divergence part
        MultiFab::Add (force_ibm[d], tmp_ibm_f[d], 0, 0, 1, 1);

        // Output pressure solution
        MultiFab::Copy(pres,         p_1,          0, 0, 1, 1);
    }
}



namespace IBM{
void MAT_GMRES(const std::array<MultiFab, AMREX_SPACEDIM> & A_MAT,
               const std::array<MultiFab, AMREX_SPACEDIM> & b_u,
               const MultiFab &                             b_p,
                     std::array<MultiFab, AMREX_SPACEDIM> & x_u,
                     MultiFab &                             x_p,
               const std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
               const MultiFab &                             beta,
               const std::array<MultiFab, NUM_EDGE> &       beta_ed,
               const MultiFab &                             gamma,
               const Real                                   theta_alpha,
               const Geometry &                             geom,
                     Real &                                 norm_pre_rhs) {

    BL_PROFILE_VAR("MAT_GMRES()", GMRES);

    // if (gmres_verbose >= 1) {
    //     Print() << "Begin call to GMRES" << std::endl;
    // }

    Vector<Real> cs(gmres_max_inner);
    Vector<Real> sn(gmres_max_inner);
    Vector<Real>  y(gmres_max_inner);
    Vector<Real>  s(gmres_max_inner+1);

    Vector<Vector<Real>> H(gmres_max_inner + 1, Vector<Real>(gmres_max_inner));

    int outer_iter, total_iter, i_copy; // for looping iteration
    int i = 0;

    Real norm_b;            // |b|;           computed once at beginning
    Real norm_pre_b;        // |M^-1 b|;      computed once at beginning
    Real norm_resid;        // |M^-1 (b-Ax)|; computed at beginning of each outer iteration
    Real norm_init_resid;   // |M^-1 (b-Ax)|; computed once at beginning
    Real norm_resid_Stokes; // |b-Ax|;        computed at beginning of each outer iteration
    Real norm_init_Stokes;  // |b-Ax|;        computed once at beginning
    Real norm_u_noprecon;   // u component of norm_resid_Stokes
    Real norm_p_noprecon;   // p component of norm_resid_Stokes
    Real norm_resid_est;

    Real norm_u; // temporary norms used to build full-state norm
    Real norm_p; // temporary norms used to build full-state norm

    Vector<Real> inner_prod_vel(AMREX_SPACEDIM);
    Real inner_prod_pres;

    const BoxArray & ba              = b_p.boxArray();
    const DistributionMapping & dmap = b_p.DistributionMap();

    // # of ghost cells must match x_u so higher-order stencils can work
    std::array< MultiFab, AMREX_SPACEDIM > r_u;
    std::array< MultiFab, AMREX_SPACEDIM > w_u;
    std::array< MultiFab, AMREX_SPACEDIM > tmp_u;
    std::array< MultiFab, AMREX_SPACEDIM > V_u;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        BoxArray ba_fc = convert(ba, nodal_flag_dir[d]);

          r_u[d].define(ba_fc, dmap, 1, x_u[d].nGrow());
          w_u[d].define(ba_fc, dmap, 1, 0);
        tmp_u[d].define(ba_fc, dmap, 1, 0);
          V_u[d].define(ba_fc, dmap, gmres_max_inner + 1, 0);
    }

    // # of ghost cells must match x_p so higher-order stencils can work
    MultiFab r_p  (ba, dmap,                  1, x_p.nGrow());
    MultiFab w_p  (ba, dmap,                  1, 0);
    MultiFab tmp_p(ba, dmap,                  1, 0);
    MultiFab V_p  (ba, dmap,gmres_max_inner + 1, 0); // Krylov vectors

    // // apply scaling factor
    // if (scale_factor != 1.) {
    //     theta_alpha = theta_alpha*scale_factor;
    //     // we will solve for scale*x_p so we need to scale the initial guess
    //     x_p.mult(scale_factor,0,1,x_p.nGrow());
    //     // scale the rhs:
    //     for (int d=0; d<AMREX_SPACEDIM; ++d) {
    //         b_u[d].mult(scale_factor,0,1,b_u[d].nGrow());
    //     }
    //     // scale the viscosities:
    //     beta.mult(scale_factor,0,1,beta.nGrow());
    //     gamma.mult(scale_factor,0,1,gamma.nGrow());
    //     for (int d=0; d<NUM_EDGE; ++d) {
    //         beta_ed[d].mult(scale_factor,0,1,beta_ed[d].nGrow());
    //     }
    // }

    // preconditioned norm_b: norm_pre_b
    // ApplyPrecon(b_u, b_p, tmp_u, tmp_p, alpha_fc, beta, beta_ed, gamma, theta_alpha, geom);

    // Compute: tmp = b - Ax
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        // tmp = A
        MultiFab::Copy(tmp_u[d], A_MAT[d], 0, 0, 1, 0);
        // tmp = Ax
        MultiFab::Multiply(tmp_u[d], x_u[d], 0, 0, 1, 0);
        // tmp = Ax - b
        MultiFab::Subtract(tmp_u[d], b_u[d], 0, 0, 1, 0);
        // tmp = b - Ax
        tmp_u[d].mult(-1., 0, 1, 0);
    }

    StagL2Norm(tmp_u, 0, norm_u);
    // CCL2Norm(tmp_p, 0, norm_p);
    // norm_p = p_norm_weight*norm_p;
    // norm_pre_b = sqrt(norm_u*norm_u+norm_p*norm_p);
    norm_pre_b = norm_u

    norm_pre_rhs = norm_pre_b;

    // calculate the l2 norm of rhs
    StagL2Norm(b_u, 0, norm_u);
    // CCL2Norm(b_p, 0, norm_p);
    // norm_p = p_norm_weight*norm_p;
    // norm_b = sqrt(norm_u*norm_u+norm_p*norm_p);
    norm_b = norm_u;

    //! If norm_b=0 we should return zero as the solution and "return" from this routine
    // It is important to use gmres_abs_tol and not 0 since sometimes due to roundoff we
    // get a nonzero number that should really be zero
    // if (gmres_verbose >= 1) {
    //     // Useful to print out to give expected scale for gmres_abs_tol
    //     Print() << "GMRES.cpp: GMRES called with ||rhs||=" << norm_b << std::endl;
    // }
    if (norm_b <= gmres_abs_tol) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) x_u[d].setVal(0.);
        // x_p.setVal(0.);
        // if (gmres_verbose >= 1) {
        //     Print() << "GMRES.cpp: converged in 0 iterations since rhs=0" << std::endl;
        // }
        return;
    }

    ///////////////////
    // begin outer iteration
    ///////////////////

    total_iter = 0;
    outer_iter = 0;

    do {
        // // Calculate tmp = Ax
        // ApplyMatrix(tmp_u, tmp_p, x_u, x_p, alpha_fc, beta, beta_ed, gamma, theta_alpha, geom);

        // // tmp = b - Ax
        // for (int d=0; d<AMREX_SPACEDIM; ++d) {
        //     MultiFab::Subtract(tmp_u[d],b_u[d], 0, 0, 1, 0);
        //     tmp_u[d].mult(-1., 0, 1, 0);
        // }
        // MultiFab::Subtract(tmp_p, b_p, 0, 0, 1, 0);
        // tmp_p.mult(-1., 0, 1, 0);

        // Compute: tmp = b - Ax
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            // tmp = A
            MultiFab::Copy(tmp_u[d], A_MAT[d], 0, 0, 1, 0);
            // tmp = Ax
            MultiFab::Multiply(tmp_u[d], x_u[d], 0, 0, 1, 0);
            // tmp = Ax - b
            MultiFab::Subtract(tmp_u[d], b_u[d], 0, 0, 1, 0);
            // tmp = b - Ax
            tmp_u[d].mult(-1., 0, 1, 0);
        }


        // un-preconditioned residuals
        StagL2Norm(tmp_u, 0, norm_u_noprecon);
        // CCL2Norm(tmp_p, 0, norm_p_noprecon);
        // norm_p_noprecon = p_norm_weight*norm_p_noprecon;
        // norm_resid_Stokes=sqrt(norm_u_noprecon*norm_u_noprecon+norm_p_noprecon*norm_p_noprecon);
        norm_resid_Stokes = norm_u_noprecon;
        if(outer_iter==0) {
            norm_init_Stokes=norm_resid_Stokes;
        }

        // if (gmres_verbose >= 2) {
        //     Print() << "total Iters = " << total_iter << std::endl;
        //     Print() << "r/(r_0,b) = " << norm_resid_Stokes/norm_init_Stokes << "  "
        //             << norm_resid_Stokes/norm_b << std::endl;
        // }
        // if (gmres_verbose >= 3) {
        //     Print() << "un-Precond. rel. resid. (u,v,p) = " << norm_resid_Stokes/norm_init_Stokes
        //             << "  " << norm_u_noprecon/norm_init_Stokes
        //             << "  " << norm_p_noprecon/norm_init_Stokes << std::endl;
        // }

        // // solve for r = M^{-1} tmp
        // // We should not be counting these toward the number of mg cycles performed
        // ApplyPrecon(tmp_u, tmp_p, r_u, r_p, alpha_fc, beta, beta_ed, gamma, theta_alpha, geom);

        // // resid = sqrt(dot_product(r, r))
        // StagL2Norm(r_u, 0, norm_u);
        // CCL2Norm(r_p, 0, norm_p);
        // norm_p = p_norm_weight*norm_p;
        // norm_resid = sqrt(norm_u*norm_u+norm_p*norm_p);

        for (int d=0; d<AMREX_SPACEDIM; ++d)
            MultiFab::Copy(r_u[d], tmp_u[d], 0, 0, 1, 0);
        
        norm_resid = norm_resid_Stokes;

        // If first iteration, save the initial preconditioned residual
        if (outer_iter==0) {
            norm_init_resid=norm_resid;
            norm_resid_est=norm_resid;
        }

        // if (gmres_verbose >= 3) {
        //     Print() << "Precond. rel. res. (u,v,p) = " << norm_resid/norm_init_resid << "  "
        //             << norm_u/norm_init_resid << "  " << norm_p/norm_init_resid << std::endl;
        // }

        // We need to test the residual now and exit OuterLoop if converged
        if (total_iter >= gmres_max_iter) {
            // if (gmres_verbose >= 1) {
            //     Print() << "GMRES did not converge in max number of total inner iterations: Exiting"
            //             << std::endl;
            // }
            break;
        }
        else if (total_iter >= gmres_min_iter) {
            // other options
            if(norm_resid <= gmres_rel_tol*std::min(norm_pre_b, norm_init_resid)) {
                // if (gmres_verbose >= 2) {
                //     Print() << "GMRES converged: Outer = " << outer_iter << ",  Inner = " << i
                //             << " Total=" << total_iter << std::endl;
                // }

                // if (norm_resid_Stokes >= 10*gmres_rel_tol*std::min(norm_b, norm_init_Stokes)) {
                //     Print() << "GMRES.cpp: Warning: gmres may not have converged: |r|/|b|= "
                //             << norm_resid_Stokes/norm_b << " |r|/|r0|="
                //             << norm_resid_Stokes/norm_init_Stokes << std::endl;
                // }

                // Only exit if the *true* preconditioned residual is less than tolerance:
                // Do not trust the gmres estimate
                break; // exit OuterLoop
            }
            else if (norm_resid <= gmres_abs_tol) {

                // if (gmres_verbose >= 2) {
                //     Print() << "GMRES converged: Outer = " << outer_iter << ",  Inner = " << i
                //             << " Total=" << total_iter << std::endl;
                // }

                break; // exit OuterLoop
            }
        }

        if (outer_iter >= gmres_max_outer) {
            // Print() << "GMRES did not converge in max number of outer iterations: Exiting" << std::endl;
            break; // exit OuterLoop
        }
        outer_iter = outer_iter + 1;

        // if (gmres_verbose >= 3) {
        //     Print() << "Begin outer iteration " << outer_iter << std::endl;
        // }

        // create the first basis in Krylov space
        // V(1) = r / norm(r)
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Copy(V_u[d] , r_u[d], 0, 0, 1, 0);
            V_u[d].mult(1./norm_resid, 0, 1, 0);
        }
        // MultiFab::Copy(V_p,r_p,0,0,1,0);
        // V_p.mult(1./norm_resid,0,1,0);

        // s = norm(r) * e_0
        std::fill(s.begin(), s.end(), 0.);
        s[0] = norm_resid;

        ///////////////////////
        // begin inner iteration
        ///////////////////////

        // i is the inner iteration loop index
        for (i=0; i<gmres_max_inner; ++i) {

            // if (gmres_verbose >= 3) {
            //     Print() << "Begin inner iteration " << i+1 << std::endl;
            // }

            total_iter = total_iter + 1;
            i_copy = i;

            // tmp=A*V(i)
            // we use r_p and r_u as temporaries to hold ith component of V
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                MultiFab::Copy(r_u[d], V_u[d], i, 0, 1, 0);
            }
            MultiFab::Copy(r_p, V_p, i, 0, 1, 0);

            ApplyMatrix(tmp_u,tmp_p,r_u,r_p,alpha_fc,beta,beta_ed,gamma,theta_alpha,geom);

            // w = M^{-1} A*V(i)
            ApplyPrecon(tmp_u,tmp_p,w_u,w_p,alpha_fc,beta,beta_ed,gamma,theta_alpha,geom);

            for (int k=0; k<=i; ++k) {
                // form H(k,i) Hessenberg matrix
                // H(k,i) = dot_product(w, V(k))
                //        = dot_product(w_u, V_u(k))+dot_product(w_p, V_p(k))
                StagInnerProd(w_u,0,V_u,k,inner_prod_vel);
                CCInnerProd(w_p,0,V_p,k,inner_prod_pres);
                H[k][i] = std::accumulate(inner_prod_vel.begin(), inner_prod_vel.end(), 0.) +
                          pow(p_norm_weight,2.0)*inner_prod_pres;

                // w = w - H(k,i) * V(k)
                //use tmp_u and tmp_p as temporaries to hold kth component of V(k)
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    MultiFab::Copy(tmp_u[d],V_u[d],k,0,1,0);
                    tmp_u[d].mult(H[k][i],0,1,0);
                    MultiFab::Subtract(w_u[d],tmp_u[d],0,0,1,0);
                }
                MultiFab::Copy(tmp_p,V_p,k,0,1,0);
                tmp_p.mult(H[k][i],0,1,0);
                MultiFab::Subtract(w_p,tmp_p,0,0,1,0);
            }

            // H(i+1,i) = norm(w)
            StagL2Norm(w_u,0,norm_u);
            CCL2Norm(w_p,0,norm_p);
            norm_p = p_norm_weight*norm_p;
            H[i+1][i] = sqrt(norm_u*norm_u+norm_p*norm_p);

            // V(i+1) = w / H(i+1,i)
            if (H[i+1][i] != 0.) {
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    MultiFab::Copy(V_u[d],w_u[d],0,i+1,1,0);
                    V_u[d].mult(1./H[i+1][i],i+1,1,0);
                }
                MultiFab::Copy(V_p,w_p,0,i+1,1,0);
                V_p.mult(1./H[i+1][i],i+1,1,0);
            }
            else {
                Abort("GMRES.cpp: error in orthogonalization");
            }

            LeastSquares(i,H,cs,sn,s); // solve least square problem
            norm_resid_est = std::abs(s[i+1]);

            if (gmres_verbose >= 2) {
                Print() << "Total iter " << total_iter << ",  est. rel. resid. |Pr|/(Pr0,b)= "
                        << norm_resid_est/norm_init_resid << "  "
                        << norm_resid_est/norm_pre_b << std::endl;
            }

            if (total_iter >= gmres_max_iter) {
                break; // exit InnerLoop
            }
            else if (total_iter >= gmres_min_iter) {
                if ((norm_resid_est <= gmres_rel_tol*std::min(norm_pre_b, norm_init_resid))
                    || (norm_resid_est <= gmres_abs_tol)) {
                    break; // exit InnerLoop
                }
            }

        } // end of inner loop

        // update the solution
        // first, solve for y
        SolveUTriangular(i_copy-1, H, s, y);

        // then, x = x + dot(V(1:i),y(1:i))
        UpdateSol(x_u,x_p,V_u,V_p,y,i_copy);

    } while (true); // end of outer loop (do iter=1,gmres_max_outer)

    // AJN - this is here since I notice epsilon roundoff errors building up
    //       just enough to destroy the asymmetry in time-advancement codes that
    //       ultimately causes lack of convergence in subsequent gmres calls
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        x_u[d].OverrideSync(geom.periodicity());
    }

    // // apply scaling factor
    // if (scale_factor != 1.) {
    //     theta_alpha = theta_alpha/scale_factor;
    //     // the solution we got is scale*x_p
    //     x_p.mult(1./scale_factor,0,1,x_p.nGrow());
    //     // unscale the rhs
    //     for (int d=0; d<AMREX_SPACEDIM; ++d) {
    //         b_u[d].mult(1./scale_factor,0,1,b_u[d].nGrow());
    //     }
    //     // unscale the viscosities
    //     beta.mult(1./scale_factor,0,1,beta.nGrow());
    //     gamma.mult(1./scale_factor,0,1,gamma.nGrow());
    //     for (int d=0; d<NUM_EDGE; ++d) {
    //         beta_ed[d].mult(1./scale_factor,0,1,beta_ed[d].nGrow());
    //     }
    // }

    if (gmres_verbose >= 1) {
        Print() << "Done with GMRES:" << std::endl;
        Print() << "  total ITERs = " << total_iter << std::endl;
        Print() << "  residual/(norm_b,initial) = " << norm_resid/norm_b << "  "
                << norm_resid/norm_init_resid << std::endl;
    }

}

void UpdateSol(std::array<MultiFab, AMREX_SPACEDIM>& x_u,
               MultiFab& x_p,
               std::array<MultiFab, AMREX_SPACEDIM>& V_u,
               MultiFab& V_p,
               Vector<Real>& y,
               int i)
{
    // set V(i) = V(i)*y(i)
    // set x = x + V(i)
    for (int iter=0; iter<=i; ++iter) {
        V_p.mult(y[iter],iter,1,0);
        MultiFab::Add(x_p,V_p,iter,0,1,0);
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            V_u[d].mult(y[iter],iter,1,0);
            MultiFab::Add(x_u[d],V_u[d],iter,0,1,0);
        }
    }
}

void LeastSquares(int i,
                  Vector<Vector<Real>>& H,
                  Vector<Real>& cs,
                  Vector<Real>& sn,
                  Vector<Real>& s)
{
    Real temp;

    // apply Givens rotation
    for (int k=0; k<=i-1; ++k) {
        temp      =  cs[k]*H[k][i] + sn[k]*H[k+1][i];
        H[k+1][i] = -sn[k]*H[k][i] + cs[k]*H[k+1][i];
        H[k][i] = temp;
    }

    // form i-th rotation matrix
    RotMat(H[i][i], H[i+1][i], cs[i], sn[i]);

    // approximate residual norm
    temp = cs[i]*s[i];
    s[i+1] = -sn[i]*s[i];
    s[i] = temp;
    H[i][i] = cs[i]*H[i][i] + sn[i]*H[i+1][i];
    H[i+1][i] = 0.;
}

void RotMat(Real a, Real b,
            Real& cs, Real& sn)
{
    Real temp;

    if (b == 0.) {
        cs = 1.;
        sn = 0.;
    }
    else if (std::abs(b) > std::abs(a)) {
        temp = a/b;
        sn = 1./sqrt(1.+temp*temp);
        cs = temp*sn;
    }
    else {
        temp = b/a;
        cs = 1./sqrt(1.+temp*temp);
        sn = temp*cs;
    }
}

void SolveUTriangular(int k, Vector<Vector<Real>>& H, Vector<Real>& s, Vector<Real>& y)
{
    Real dot;

    y[k+1] = s[k+1]/H[k+1][k+1];
    for (int i=k; i>=0; --i) {

        dot = 0.;
        for (int j=i+1; j<= k+1; ++j) {
            dot += H[i][j]*y[j];
        }

        y[i] = (s[i] - dot) / H[i][i];
    }
}

}
