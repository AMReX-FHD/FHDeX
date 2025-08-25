#include "main_driver.H"
#include "main_driver_F.H"

#include "common_functions.H"

#include "hydro_functions.H"


#include "gmres_functions.H"


#include "ib_functions.H"

#include "IBCore.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

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
             IBParticleContainer & ib_pc,
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
        MultiFabPhysBCDomainVel(umac[i], i, geom, i);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umac[i], i, geom, i, is_inhomogeneous);
    }



    /****************************************************************************
     *                                                                          *
     * Advance tracer                                                           *
     *                                                                          *
     ***************************************************************************/

    // Compute tracer:
    tracer.FillBoundary(geom.periodicity());
    MultiFabPhysBC(tracer, geom, 0, 1, SPEC_BC_COMP);

    MkAdvSFluxdiv_cc(umac, tracer, advFluxdivS, geom, 0, 1, 0);
    advFluxdivS.mult(dt, 1);

    // compute predictor
    MultiFab::Copy(tracerPred, tracer, 0, 0, 1, 0);
    MultiFab::Add(tracerPred, advFluxdivS, 0, 0, 1, 0);

    tracerPred.FillBoundary(geom.periodicity());
    MultiFabPhysBC(tracerPred, geom, 0, 1, SPEC_BC_COMP);

    MkAdvSFluxdiv_cc(umac, tracerPred, advFluxdivS, geom, 0, 1, 0);
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
        MultiFabPhysBCDomainVel(uMom[d], d, geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(uMom[d], d, geom, d, is_inhomogeneous);
    }

    // advective momentum flux terms
    MkAdvMFluxdiv(umac, uMom, advFluxdiv, dx, 0);

    // crank-nicolson terms
    StagApplyOp(geom, beta_negwtd, gamma_negwtd, beta_ed_negwtd, umac, Lumac, alpha_fc_0,
                dx, theta_alpha);

    // advective term boundary conditions
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        Lumac[d].FillBoundary(geom.periodicity());
        // MultiFabPhysBC fixme face-centered data

        advFluxdiv[d].FillBoundary(geom.periodicity());
        // MultiFabPhysBC fixme face-centered data
    }


    //___________________________________________________________________________
    // Pressure term boundary conditions

    pres.setVal(0.);
    MultiFabPhysBC(pres, geom, 0, 1, PRES_BC_COMP, 0);
    ComputeGrad(pres, pg, 0, 0, 1, PRES_BC_COMP, geom);


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

    IBGMRES(gmres_rhs_u, gmres_rhs_p, umacNew, p_0,
            alpha_fc_1, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
            ib_pc, geom, norm_pre_rhs);

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(umac_0[d], umacNew[d], 0, 0, 1, 1);
        VisMF::Write(umac_0[d], "umac_0_"+std::to_string(d));
    }



    /****************************************************************************
     *                                                                          *
     * ADVANCE (CORRECTOR) STEP (crank-nicolosn heun's method: part 2)          *
     *                                                                          *
     ***************************************************************************/

    //___________________________________________________________________________
    // Compute corrector's advective term (using predictor's fluid solution)

     for (int d=0; d<AMREX_SPACEDIM; d++) {
        umacNew[d].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(umacNew[d], d, geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umacNew[d], d, geom, d, is_inhomogeneous);

        MultiFab::Copy(uMom[d], umacNew[d], 0, 0, 1, 1);

        // let rho = 1
        uMom[d].mult(1.0, 1);

        // momentum boundary conditions
        uMom[d].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(uMom[d], d, geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(uMom[d], d, geom, d, is_inhomogeneous);
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
    StagApplyOp(geom, beta_negwtd, gamma_negwtd, beta_ed_negwtd, umac, Lumac,
                alpha_fc_0, dx, theta_alpha);


    // advective term boundary conditions
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        Lumac[d].FillBoundary(geom.periodicity());
        // MultiFabPhysBC fixme face-centered data

        advFluxdiv[d].FillBoundary(geom.periodicity());
        // MultiFabPhysBC fixme face-centered data

        advFluxdivPred[d].FillBoundary(geom.periodicity());
        // MultiFabPhysBC fixme face-centered data
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

        IBGMRES(gmres_rhs_u, gmres_rhs_p, umacNew, p_1,
                alpha_fc_1, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
                ib_pc, geom, norm_pre_rhs);


        //_______________________________________________________________________
        // Slip velocity term

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Copy    (r_f[d], umacNew[d], 0, 0, 1, 1);
            MultiFab::Multiply(r_f[d], f_ibm[d],   0, 0, 1, 1);
            r_f[d].mult(-1., 0);
        }


        // //_______________________________________________________________________
        // // Invert motility matrix

        // // initial guess
        // p_ibm_1.setVal(0.);
        // for (int d=0; d<AMREX_SPACEDIM; ++d)
        //     tmp_ibm_f[d].setVal(0.);

        // // Inverse Motility Matrix
        // ApplyMatrix(tmp_f_1, p_f, r_f, p_ibm_1,
        //             alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha, geom);

        // // // Remove non-divergence free parts of the residual
        // // p_f.mult(-1., 0);
        // // MultiFabPhysBC(p_f, geom, 0, 1, 0);

        // // MacProj(ones_fc, p_f, tmp_pf, geom, 1);
        // // tmp_pf.FillBoundary(geom.periodicity());
        // // MultiFabPhysBC(tmp_pf, geom, 0, 1, 0);

        // // SubtractWeightedGradP(tmp_ibm_f, ones_fc, tmp_pf, geom);
        // // MultiFab::Add(p_ibm_1, tmp_pf, 0, 0, 1, 1);

        // MultiFab::Add(gmres_rhs_p, p_f, 0, 0, 1, 1);
        // gmres_rhs_p.FillBoundary(geom.periodicity());


        // //______________________________________________________________________
        // // Apply immersed-boundary force

        // // Inverse-motility part
        // for (int d=0; d<AMREX_SPACEDIM; ++d) {
        //     MultiFab::Copy(tmp_f_mask[d],     tmp_f_1[d], 0, 0, 1, 1);
        //     MultiFab::Multiply(tmp_f_mask[d], f_ibm[d],   0, 0, 1, 1);

        //     // MultiFab::Add(force_1[d], tmp_f_mask[d],    0, 0, 1, 1);
        //     // Add raw force to fluid => includes pressure correction terms from
        //     // immersed boundary => this might not be the right thing to do...
        //     MultiFab::Add(force_1[d], tmp_f_1[d],    0, 0, 1, 1);
        // }

        // // // Divergence part
        // // for (int d=0; d<AMREX_SPACEDIM; ++d) {
        // //     MultiFab::Subtract(umacNew[d], tmp_ibm_f[d], 0, 0, 1, 1);
        // //     MultiFab::Copy(umac_1[d], umacNew[d], 0, 0, 1, 1);
        // // }

        // // Implicit (move to lhs) part
        // for (int d=0; d<AMREX_SPACEDIM; ++d) {
        //     MultiFab::Copy(force_est[d], r_f[d], 0, 0, 1, 1);
        //     force_est[d].mult(1e6, 0);

        //     //MultiFab::Add(force_1[d], force_est[d], 0, 0, 1, 1);
        // }


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


        // if (norm_r < 1e-12) break;
        break;
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
