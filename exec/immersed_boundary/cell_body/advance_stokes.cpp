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
void advance_stokes(std::array<MultiFab, AMREX_SPACEDIM >& umac,
                    std::array<MultiFab, AMREX_SPACEDIM >& umacNew,
                    MultiFab& pres,
                    IBMarkerContainer & ib_mc,
                    const std::map<std::tuple<int, int>, double> & bond_map,
                    const std::map<int, std::vector<int>> & bond_neighbors,
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

    // beta_wtd on nodes in 2d, on edges in 3d
    std::array< MultiFab, NUM_EDGE > beta_ed_wtd;
#if (AMREX_SPACEDIM == 2)
    beta_ed_wtd[0].define(convert(ba,nodal_flag), dmap, 1, 1);
    MultiFab::Copy(beta_ed_wtd[0], beta_ed[0], 0, 0, 1, 1);
#elif (AMREX_SPACEDIM == 3)
    for(int d=0; d<AMREX_SPACEDIM; d++) {
        beta_ed_wtd[d].define(convert(ba, nodal_flag_edge[d]), dmap, 1, 1);
        MultiFab::Copy(beta_ed_wtd[d], beta_ed[d], 0, 0, 1, 1);
    }
#endif

    // Scaled by 1/2:
    // gamma_wtd cell centered
    MultiFab gamma_wtd(ba, dmap, 1, 1);
    MultiFab::Copy(gamma_wtd, gamma, 0, 0, 1, 1);



    /****************************************************************************
     *                                                                          *
     * Immersed-Marker parameters                                               *
     *                                                                          *
     ***************************************************************************/

    int ib_lev = 0;

    // Parameters for calling bending force calculation
    RealVect driv_u = {0, 0, 1};

    // // Slowly ramp up driving amplitude
    // Real driv_amp = amrex::min(time*100, 1.);
    // Print() << "driv_amp = " << driv_amp << std::endl;

    // I'm too impatient to wait... -JPB
    Real driv_amp = 1.;



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
     * ADVANCE (stokes) step, compute:                                          *
     * 1.           x^(n+1) = x^n + dt J(u^(n))      [x = marker positions]     *
     * 2. L(u^(n+1))/2 + Gp = S(f(x^(n+1)))                                     *
     *             Du^(n+1) = 0                                                 *
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

    // anchor the first two markers on flagellum with the two markers on the cells body
    // by ensuing they share the same velocities before moving
    if(immbdy::contains_fourier)
       anchor_coupling_markers(ib_mc, ib_lev, IBMReal::velx);
//     anchor_first_marker(ib_mc, ib_lev, IBMReal::velx);
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

    update_bdy_marker(bond_map, bond_neighbors, time, ib_mc, ib_lev,
                      IBMReal::forcex, false,
                      geom);

    // Constrain it to move in the z = constant plane only
    constrain_ibm_marker(ib_mc, ib_lev, IBMReal::forcez);
    //if(immbdy::contains_fourier)
    //    anchor_first_marker(ib_mc, ib_lev, IBMReal::forcex);
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


    //___________________________________________________________________________
    // Compute pressure, and pressure gradient due to the BC: gp = Gp
    // Note that the pressure gradient due to the BC is left unchanged
    pres.setVal(0.);           // Initial guess for pressure
    MultiFabPhysBC(pres, geom, 0, 1, PRES_BC_COMP, 0); // Apply pressure boundary conditions

    // Initial guess for new solution
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        MultiFab::Copy(umacNew[d], umac[d], 0, 0, 1, 1);


    advanceStokes(
            umacNew, pres,                    /* LHS */
            mfluxdiv_correct, fc_force_corr,  /* RHS */
            alpha_fc, beta_wtd, gamma_wtd, beta_ed_wtd, geom, dt
        );


    // Apply boundary conditions to the solution
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        umacNew[d].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(umacNew[d], geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umacNew[d], geom, d, is_inhomogeneous);
    }

    // Update solution, and we're done!
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(umac[d],     umacNew[d],       0, 0, 1, 0);
        MultiFab::Copy(force_ib[d], fc_force_corr[d], 0, 0, 1, 0);
    }

    BL_PROFILE_VAR_STOP(advance);
}