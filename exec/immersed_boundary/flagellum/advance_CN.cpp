#include <main_driver.H>

#include <hydro_functions.H>

#include <gmres_functions.H>


#include <immbdy_namespace.H>

#include <IBMarkerMD.H>


using namespace amrex;


using namespace immbdy;
using namespace immbdy_md;
using namespace ib_flagellum;



// Crank-Nicolson Advance Subroutine
void advance_CN(std::array<MultiFab, AMREX_SPACEDIM >& umac,
                std::array<MultiFab, AMREX_SPACEDIM >& umacNew,
                MultiFab& pres,
                IBMarkerContainer & ib_mc,
                const std::array<MultiFab, AMREX_SPACEDIM>& mfluxdiv_predict,
                const std::array<MultiFab, AMREX_SPACEDIM>& mfluxdiv_correct,
                      std::array<MultiFab, AMREX_SPACEDIM>& alpha_fc,
                      std::array<MultiFab, AMREX_SPACEDIM>& force_ib,
                const MultiFab& beta, const MultiFab& gamma,
                const std::array<MultiFab, NUM_EDGE> & beta_ed,
                const Geometry geom, const Real& dt, Real time)
{

    BL_PROFILE_VAR("advance()", advance);

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
    // Pressure gradient at inflow/outflow
    std::array< MultiFab, AMREX_SPACEDIM > pg;

    for (int i=0; i<AMREX_SPACEDIM; i++) {
           gmres_rhs_u[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                 Lumac[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 0);
            advFluxdiv[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
        advFluxdivPred[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                  uMom[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                    pg[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);

        // Put in to fix FPE traps
        advFluxdivPred[i].setVal(0);
        advFluxdiv[i].setVal(0);
    }



    //___________________________________________________________________________
    // Scaled alpha, beta, gamma:
    //  * scaled by 1/2 => corrector/ crank-nicolson terms
    //  * scaled by -1  => move implicit term (LHS) to explicit term (RHS)

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
     * Immersed-Marker parameters                                               *
     *                                                                          *
     ***************************************************************************/

    int ib_lev = 0;

    // Parameters for calling bending force calculation
    RealVect driv_u = {0, 0, 1};

    // Slowly ramp up driving amplitude
    Real driv_amp = amrex::min(time*100, 1.);
    Print() << "driv_amp = " << driv_amp << std::endl;

    // I'm too impatient to wait... -JPB
    // Real driv_amp = 1.;



    /****************************************************************************
     *                                                                          *
     * Apply deterministic boundary conditions                                  *
     *                                                                          *
     ***************************************************************************/

    for (int d=0; d<AMREX_SPACEDIM; d++) {
        umac[d].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(umac[d], geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umac[d], geom, d, is_inhomogeneous);
    }



    /****************************************************************************
     *                                                                          *
     * PREDICTOR (midpoint) step, compute:                                      *
     * 1. x^(n+1/2) = x^n + dt/2 J(u^n)              [x = marker positions]     *
     * 2. F^(n+1/2) = f(x^(n+1/2))                   [f - IB force model  ]     *
     * 3. f^(n+1/2) = S(F^(n+1/2))                                              *
     * 4. (u^(n+1/2)-u^n)/(dt/2) + D(uu^n) + Gp = Lu^(n+1/2) + f(n+1/2)         *
     *                               Du^(n+1/2) = 0                             *
     *                                                                          *
     ***************************************************************************/


    //___________________________________________________________________________
    // Interpolate immersed boundary predictor: J(u^n)
    std::array<MultiFab, AMREX_SPACEDIM> umac_buffer;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        umac_buffer[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
        umac_buffer[d].setVal(0.);
        MultiFab::Copy(umac_buffer[d], umac[d], 0, 0, 1, umac[d].nGrow());
        umac_buffer[d].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(umac[d], geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umac[d], geom, d, is_inhomogeneous);
    }

    ib_mc.ResetPredictor(0);
    ib_mc.InterpolatePredictor(0, umac_buffer);


    //___________________________________________________________________________
    // Move predictor using previous velocity: x^(n+1/2) = x^n + dt/2 J(u^n)
    // (constrain it to move in the z = constant plane only)
    constrain_ibm_marker(ib_mc, ib_lev, IBMReal::pred_velz);
    if(immbdy::contains_fourier)
        anchor_first_marker(ib_mc, ib_lev, IBMReal::pred_velx);
    ib_mc.MovePredictor(0, dt);

    ib_mc.clearNeighbors(); // Important: clear neighbors before Redistribute
    ib_mc.Redistribute();   // Don't forget to send particles to the right CPU


    //___________________________________________________________________________
    // Update forces between markers: F^(n+1/2) = f(x^(n+1/2)) TODO: expensive
    // => use infrequently, use updateNeighbors for most steps
    ib_mc.clearNeighbors();
    ib_mc.fillNeighbors(); // Does ghost cells
    ib_mc.buildNeighborList(ib_mc.CheckPair);

    update_ibm_marker(driv_u, driv_amp, time, ib_mc, ib_lev,
                      IBMReal::pred_forcex, true,
                      geom);
    // Constrain it to move in the z = constant plane only
    constrain_ibm_marker(ib_mc, ib_lev, IBMReal::pred_forcez);
    if(immbdy::contains_fourier)
        anchor_first_marker(ib_mc, ib_lev, IBMReal::pred_forcex);
    // Sum predictor forces added to neighbors back to the real markers
    ib_mc.sumNeighbors(IBMReal::pred_forcex, AMREX_SPACEDIM, 0, 0);


    //___________________________________________________________________________
    // Spread forces to predictor f^(n+1/2) = S(F^(n+1/2))
    // Remember: Spread, Fold Fill, Sum
    std::array<MultiFab, AMREX_SPACEDIM> fc_force_pred;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        fc_force_pred[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
        fc_force_pred[d].setVal(0.);
    }

    ib_mc.SpreadPredictor(0, fc_force_pred);
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        fc_force_pred[d].SumBoundary(geom.periodicity());
    }


    //___________________________________________________________________________
    // Compute predictor velocity field
    // (u^(n+1/2)-u^n)/(dt/2) + D(uu^n) + Gp = Lu^(n+1/2) + f(n+1/2)
    //                            Du^(n+1/2) = 0

    // Compute momentum fluxes: uMom = rho*u^n
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        MultiFab::Copy(uMom[d], umac[d], 0, 0, 1, 1);

        // let rho = 1
        uMom[d].mult(1.0, 1);

        uMom[d].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(uMom[d], geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(uMom[d], geom, d, is_inhomogeneous);
    }

    // Compute advective fluxes: advFluxdiv = - D(\rho uu^n) = - D(u^n uMom)
    MkAdvMFluxdiv(umac, uMom, advFluxdiv, dx, 0);

    // (crank-nicolson terms) Explicit part of the diffusive operator Lu^n/2.
    // Note that we are using the weighted coefficients (to deal with the 1/2
    // part)
    StagApplyOp(geom, beta_negwtd, gamma_negwtd, beta_ed_negwtd, umac, Lumac,
                alpha_fc_0, dx, theta_alpha);


    //___________________________________________________________________________
    // Compute pressure, and pressure gradient due to the BC: gp = Gp
    pres.setVal(0.); // Initial guess for pressure
    MultiFabPhysBC(pres, geom, 0, 1, PRES_BC_COMP, 0); // Apply pressure boundary conditions
    for (int d=0; d<AMREX_SPACEDIM; ++d) pg[d].setVal(0);
    ComputeGrad(pres, pg, 0, 0, 1, PRES_BC_COMP, geom);

    // Construct RHS of Navier Stokes Equation
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);
        gmres_rhs_u[d].mult(dtinv, 1); // advance by dt

        MultiFab::Add(gmres_rhs_u[d], mfluxdiv_predict[d], 0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], fc_force_pred[d],    0, 0, 1, 0);

        // FillBoundary before adding boundary conditions to prevent
        // overwriting (which are defined on ghost cells)
        pg[d].FillBoundary(geom.periodicity());
        gmres_rhs_u[d].FillBoundary(geom.periodicity());

        MultiFab::Subtract(gmres_rhs_u[d], pg[d], 0, 0, 1, 1);
    }

    // Initial guess for new solution
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        MultiFab::Copy(umacNew[d], umac[d], 0, 0, 1, 1);

    // Call GMRES to compute u^(n+1/2). Lu^(n+1/2) is computed implicitly. Note
    // that we are using the un-weighted coefficients.
    GMRES gmres(ba, dmap, geom);
    gmres.Solve(gmres_rhs_u, gmres_rhs_p, umacNew, pres, alpha_fc, beta_wtd,
                beta_ed_wtd, gamma_wtd, theta_alpha, geom, norm_pre_rhs);

    // Apply boundary conditions to the solution
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umacNew[d].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(umacNew[d], geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umacNew[d], geom, d, is_inhomogeneous);
    }



    /****************************************************************************
     *                                                                          *
     * ADVANCE (midpoint) step, compute:                                        *
     * 1. x^(n+1) = x^n + dt J(u^(n+1/2))            [x = marker positions]     *
     * 2. (u^(n+1)-u^n)/dt + D(uu^(n+1/2)) + Gp = L(u^n+u^(n+1))/2 + f(n+1/2)   *
     *                                 Du^(n+1) = 0                             *
     *                                                                          *
     ***************************************************************************/


    //___________________________________________________________________________
    // Interpolate immersed boundary: J(u^(n+1/2))
    std::array<MultiFab, AMREX_SPACEDIM> umacNew_buffer;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        umacNew_buffer[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
        umacNew_buffer[d].setVal(0.);
        MultiFab::Copy(umacNew_buffer[d], umacNew[d], 0, 0, 1, umac[d].nGrow());
        umacNew_buffer[d].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(umacNew[d], geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umacNew[d], geom, d, is_inhomogeneous);
    }

    ib_mc.ResetMarkers(0);
    ib_mc.InterpolateMarkers(0, umacNew_buffer);


    //___________________________________________________________________________
    // Move markers according to velocity: x^(n+1) = x^n + dt/2 J(u^(n+1/2))
    // (constrain it to move in the z = constant plane only)
    constrain_ibm_marker(ib_mc, ib_lev, IBMReal::velz);
    if(immbdy::contains_fourier)
        anchor_first_marker(ib_mc, ib_lev, IBMReal::velx);
    ib_mc.MoveMarkers(0, dt);

    ib_mc.clearNeighbors(); // Important: clear neighbors before Redistribute
    ib_mc.Redistribute();   // Don't forget to send particles to the right CPU


    //___________________________________________________________________________
    // Update forces between markers: F^(n+1) = f(x^(n+1)) TODO: expensive =>
    // use infrequently, use updateNeighbors for most steps
    ib_mc.clearNeighbors();
    ib_mc.fillNeighbors(); // Does ghost cells
    ib_mc.buildNeighborList(ib_mc.CheckPair);

    update_ibm_marker(driv_u, driv_amp, time, ib_mc, ib_lev,
                      IBMReal::forcex, false,
                      geom);
    // Constrain it to move in the z = constant plane only
    constrain_ibm_marker(ib_mc, ib_lev, IBMReal::forcez);
    if(immbdy::contains_fourier)
        anchor_first_marker(ib_mc, ib_lev, IBMReal::forcex);
    // Sum predictor forces added to neighbors back to the real markers
    ib_mc.sumNeighbors(IBMReal::forcex, AMREX_SPACEDIM, 0, 0);


    //___________________________________________________________________________
    // Spread forces to corrector: f^(n+1) = S(F^(n+1))
    // Remember: Spread, Fold Fill, Sum
    std::array<MultiFab, AMREX_SPACEDIM> fc_force_corr;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        fc_force_corr[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
        fc_force_corr[d].setVal(0.);
    }

    ib_mc.SpreadMarkers(0, fc_force_corr);
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        fc_force_corr[d].SumBoundary(geom.periodicity());


    //__________________________________________________________________________
    // Compute corrector velocity field
    // (u^(n+1)-u^n)/(dt) + D(uu^(n+1/2)) + Gp = L(u^n+u^(n+1))/2 + f(n+1/2)
    //                                Du^(n+1) = 0

    // Compute momentum fluxes at the midpoint: uMom = rho*u^(n+1/2)
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        MultiFab::Copy(uMom[d], umacNew[d], 0, 0, 1, 0);

        // let rho = 1
        uMom[d].mult(1.0, 1);

        uMom[d].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(uMom[d], geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(uMom[d], geom, d, is_inhomogeneous);
    }

    // Compute advective fluxes at the midpoint:
    // advFluxdivPred = - D(\rho uu^(n+1/2)) = - D(u^(n+1/2) uMom)
    MkAdvMFluxdiv(umacNew, uMom, advFluxdivPred, dx, 0);

    // trapezoidal advective terms
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        advFluxdiv[d].mult(0.5, 1);
        advFluxdivPred[d].mult(0.5, 1);
    }


    // // Explicit part of the diffusive operator Lu^n/2. Note that we are using
    // // the weighted coefficients (to deal witht he 1/2 part)
    // StagApplyOp(geom, beta_negwtd, gamma_negwtd, beta_ed_negwtd,
    //             umac, Lumac, alpha_fc_0, dx, theta_alpha);


    //___________________________________________________________________________
    // Compute pressure, and pressure gradient due to the BC: gp = Gp
    // Note that the pressure gradient due to the BC is left unchanged
    pres.setVal(0.); // Initial guess for pressure
    MultiFabPhysBC(pres, geom, 0, 1, PRES_BC_COMP, 0); // Apply pressure boundary conditions

    // Construct RHS of Navier Stokes Equation
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        Lumac[d].FillBoundary(geom.periodicity());

        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);
        gmres_rhs_u[d].mult(dtinv, 1); // advance by dt

        MultiFab::Add(gmres_rhs_u[d], mfluxdiv_correct[d], 0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], advFluxdivPred[d],   0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], fc_force_corr[d],    0, 0, 1, 0);

        // FillBoundary before adding boundary conditions to prevent
        // overwriting (which are defined on ghost cells)
        gmres_rhs_u[d].FillBoundary(geom.periodicity());

        MultiFab::Subtract(gmres_rhs_u[d], pg[d], 0, 0, 1, 1);
    }

    // Initial guess for new solution
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        MultiFab::Copy(umacNew[d], umac[d], 0, 0, 1, 1);

    // Call GMRES to compute u^(n+1). Lu^(n+1)/2 is computed implicitly. Note
    // that we are using the weighted coefficients (to deal witht he 1/2 part)
    gmres.Solve(gmres_rhs_u, gmres_rhs_p, umacNew, pres, alpha_fc, beta_wtd,
                beta_ed_wtd, gamma_wtd, theta_alpha, geom, norm_pre_rhs);

    // Apply boundary conditions to the solution
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        umacNew[d].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(umacNew[d], geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umacNew[d], geom, d, is_inhomogeneous);
    }

    // Update solution, and we're done!
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(umac[d], umacNew[d],           0, 0, 1, 0);
        MultiFab::Copy(force_ib[d], fc_force_corr[d], 0, 0, 1, 0);
    }

    BL_PROFILE_VAR_STOP(advance);
}
