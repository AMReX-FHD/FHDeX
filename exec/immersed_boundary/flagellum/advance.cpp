#include <main_driver.H>

#include <hydro_functions.H>
#include <hydro_functions_F.H>

#include <common_functions.H>
#include <common_functions_F.H>

#include <common_namespace.H>

#include <gmres_functions.H>
#include <gmres_functions_F.H>

#include <gmres_namespace.H>

#include <immbdy_namespace.H>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

#include <IBMarkerContainer.H>
#include <IBMarkerMD.H>


using namespace amrex;

using namespace common;
using namespace gmres;

using namespace immbdy;
using namespace immbdy_md;
using namespace ib_flagellum;


using ParticleVector = typename IBMarkerContainer::ParticleVector;


void constrain_ibm_marker(IBMarkerContainer & ib_mc, int ib_lev, int component) {

    BL_PROFILE_VAR("constrain_ibm_marker", Constrain);

    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        PairIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = ib_mc.GetParticles(ib_lev).at(index).GetArrayOfStructs();

        long np = markers.size();

        for (int i = 0; i < np; ++i) {

            ParticleType & mark = markers[i];
            // Zero component only
            mark.rdata(component) = 0.;

            // HAXOR, holding marker 0 fixed
            if((mark.idata(IBM_intData::id_1) == 0) || (mark.idata(IBM_intData::id_1) == 1) ) {
                std::cout << "contraining marker " << mark.idata(IBM_intData::id_1) << ":"<< std::endl;
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    // std::cout << "    " << mark.rdata(IBM_realData::pred_velx + d) << std::endl;
                    // std::cout << "    " << mark.rdata(IBM_realData::velx + d) << std::endl;
                    // std::cout << "    " << mark.rdata(IBM_realData::pred_forcex + d) << std::endl;
                    // std::cout << "    " << mark.rdata(IBM_realData::forcex + d) << std::endl;

                    mark.rdata(IBM_realData::pred_velx + d)   = 0;
                    mark.rdata(IBM_realData::velx + d)        = 0;

                    mark.rdata(IBM_realData::pred_forcex + d) = 0;
                    mark.rdata(IBM_realData::forcex + d)      = 0;
                }
            }
        }
    }

    BL_PROFILE_VAR_STOP(Constrain);
}



Real theta(Real amp_ramp, Real time, int i_ib, int index_marker) {

    if(immbdy::contains_fourier) {

        int N                 = chlamy_flagellum::N[i_ib][index_marker];
        int coef_len          = ib_flagellum::fourier_coef_len;
        Vector<Real> & a_coef = chlamy_flagellum::a_coef[i_ib][index_marker];
        Vector<Real> & b_coef = chlamy_flagellum::b_coef[i_ib][index_marker];

        Real dt = 1./N;
        Real th = 0;

        time = 0; // HAXOR for now
        Real k_fact = 2*M_PI/N;
        Real t_unit = time/dt;
        for (int i=0; i < coef_len; ++i) {
            Real k = k_fact * i;
            th += a_coef[i]*cos(k*t_unit);
            th -= b_coef[i]*sin(k*t_unit);
        }

        return amp_ramp*th/N;

    } else {
        int  N          = ib_flagellum::n_marker[i_ib];
        Real L          = ib_flagellum::length[i_ib];
        Real wavelength = ib_flagellum::wavelength[i_ib];
        Real frequency  = ib_flagellum::frequency[i_ib];
        Real amplitude  = ib_flagellum::amplitude[i_ib];
        Real l_link     = L/(N-1);

        Real theta = l_link*amp_ramp*amplitude*sin(2*M_PI*frequency*time
                     + 2*M_PI/wavelength*index_marker*l_link);

        return theta;
    }
}



void update_ibm_marker(const RealVect & driv_u, Real driv_amp, Real time,
                       IBMarkerContainer & ib_mc, int ib_lev,
                       int component, bool pred_pos) {

    BL_PROFILE_VAR("update_ibm_marker", UpdateForces);

    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        PairIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = ib_mc.GetParticles(ib_lev).at(index).GetArrayOfStructs();

        // Get neighbor marker data (from neighboring threads)
        ParticleVector & nbhd_data = ib_mc.GetNeighbors(
                    ib_lev, pti.index(), pti.LocalTileIndex()
                );

        // Get neighbor list (for collision checking)
        const Vector<int> & nbhd = ib_mc.GetNeighborList(
                    ib_lev, pti.index(), pti.LocalTileIndex()
                );

        long np        = markers.size();
        int nbhd_index = 0;

        for (int i=0; i<np; ++i) {

            ParticleType & mark = markers[i];

            int i_ib    = mark.idata(IBM_intData::cpu_1);
            int N       = ib_flagellum::n_marker[i_ib];
            Real L      = ib_flagellum::length[i_ib];
            Real l_link = L/(N-1);

            Real k_spr  = ib_flagellum::k_spring[i_ib];
            Real k_driv = ib_flagellum::k_driving[i_ib];


            // Get previous and next markers connected to current marker (if they exist)
            ParticleType * next_marker = NULL;
            ParticleType * prev_marker = NULL;

            int status = ib_mc.FindConnectedMarkers(markers, mark,
                                                    nbhd_data, nbhd,
                                                    nbhd_index,
                                                    prev_marker, next_marker);


            if (status == -1) Abort("status -1 particle detected in predictor!!! flee for your life!");

            // position vectors
            RealVect prev_pos, pos, next_pos;
            if (status == 0) {
                for(int d=0; d<AMREX_SPACEDIM; ++d) {
                    prev_pos[d] = prev_marker->pos(d);
                    pos[d]      =         mark.pos(d);
                    next_pos[d] = next_marker->pos(d);

                    if (pred_pos) {
                        prev_pos[d] += prev_marker->rdata(IBM_realData::pred_posx + d);
                        pos[d]      += mark.rdata(IBM_realData::pred_posx + d);
                        next_pos[d] += next_marker->rdata(IBM_realData::pred_posx + d);
                    }
                }
            }

            // update spring forces
            if (status == 0) { // has next (p) and prev (m)

                RealVect r_p = next_pos - pos, r_m = pos - prev_pos;

                Real lp_m = r_m.vectorLength(),         lp_p = r_p.vectorLength();
                Real fm_0 = k_spr * (lp_m-l_link)/lp_m, fp_0 = k_spr * (lp_p-l_link)/lp_p;

                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    prev_marker->rdata(component + d) += fm_0 * r_m[d];
                    mark.rdata(component + d)         -= fm_0 * r_m[d];

                    mark.rdata(component + d)         += fp_0 * r_p[d];
                    next_marker->rdata(component + d) -= fp_0 * r_p[d];
                }
            }

            // update bending forces for curent, minus/prev, and next/plus
            if (status == 0) { // has next (p) and prev (m)

                // position vectors
                const RealVect & r = pos, & r_m = prev_pos, & r_p = next_pos;

                // Set bending forces to zero
                RealVect f_p = RealVect{AMREX_D_DECL(0., 0., 0.)};
                RealVect f   = RealVect{AMREX_D_DECL(0., 0., 0.)};
                RealVect f_m = RealVect{AMREX_D_DECL(0., 0., 0.)};

                // calling the active bending force calculation
                // This a simple since wave imposed

                Real th = theta(driv_amp, time, i_ib, mark.idata(IBM_intData::id_1));
                driving_f(f, f_p, f_m, r, r_p, r_m, driv_u, th, k_driv);

                // updating the force on the minus, current, and plus particles.
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    prev_marker->rdata(component + d) += f_m[d];
                    mark.rdata(component + d)         +=   f[d];
                    next_marker->rdata(component + d) += f_p[d];
                }
            }

            // Increment neighbor list
            int nn      = nbhd[nbhd_index];
            nbhd_index += nn + 1; // +1 <= because the first field contains nn
        }
    }
    BL_PROFILE_VAR_STOP(UpdateForces);
};




void advance(std::array< MultiFab, AMREX_SPACEDIM >& umac,
             std::array< MultiFab, AMREX_SPACEDIM >& umacNew,
             MultiFab& pres,
             IBMarkerContainer & ib_mc,
             const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_predict,
             const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_correct,
                   std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
             const MultiFab& beta, const MultiFab& gamma,
             const std::array< MultiFab, NUM_EDGE >& beta_ed,
             const Geometry geom, const Real& dt, Real time)
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
    // Pressure gradient at inflow/outflow
    std::array< MultiFab, AMREX_SPACEDIM > pg;

    for (int i=0; i<AMREX_SPACEDIM; i++) {
           gmres_rhs_u[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                 Lumac[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
            advFluxdiv[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
        advFluxdivPred[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                  uMom[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                    pg[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
    }

    // Tracer advective terms
    MultiFab advFluxdivS(ba, dmap, 1, 1);


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
    MultiFab beta_old(ba, dmap, 1, 1);
    MultiFab::Copy(beta_wtd, beta, 0, 0, 1, 1);
    MultiFab::Copy(beta_old, beta, 0, 0, 1, 1);
    beta_wtd.mult(0.5, 1);

    // beta_wtd on nodes in 2d, on edges in 3d
    std::array< MultiFab, NUM_EDGE > beta_ed_wtd;
    std::array< MultiFab, NUM_EDGE > beta_ed_old;
#if (AMREX_SPACEDIM == 2)
    beta_ed_wtd[0].define(convert(ba,nodal_flag), dmap, 1, 1);
    MultiFab::Copy(beta_ed_wtd[0], beta_ed[0], 0, 0, 1, 1);
    beta_ed_wtd[0].mult(0.5, 1);
    beta_ed_old[0].define(convert(ba,nodal_flag), dmap, 1, 1);
    MultiFab::Copy(beta_ed_old[0], beta_ed[0], 0, 0, 1, 1);
#elif (AMREX_SPACEDIM == 3)
    for(int d=0; d<AMREX_SPACEDIM; d++) {
        beta_ed_wtd[d].define(convert(ba, nodal_flag_edge[d]), dmap, 1, 1);
        beta_ed_old[d].define(convert(ba, nodal_flag_edge[d]), dmap, 1, 1);

        MultiFab::Copy(beta_ed_wtd[d], beta_ed[d], 0, 0, 1, 1);
        MultiFab::Copy(beta_ed_old[d], beta_ed[d], 0, 0, 1, 1);
        beta_ed_wtd[d].mult(0.5, 1);
    }
#endif

    // Scaled by 1/2:
    // gamma_wtd cell centered
    MultiFab gamma_wtd(ba, dmap, 1, 1);
    MultiFab gamma_old(ba, dmap, 1, 1);
    MultiFab::Copy(gamma_wtd, gamma, 0, 0, 1, 1);
    MultiFab::Copy(gamma_old, gamma, 0, 0, 1, 1);
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
    Real driv_amp = std::min(time*10, 1.);
    Print() << "driv_amp = " << driv_amp << std::endl;

    // I'm too impatient to wait... -JPB
    // Real driv_amp = 1.;



    /****************************************************************************
     *                                                                          *
     * Apply deterministic boundary conditions                                  *
     *                                                                          *
     ***************************************************************************/

    for (int i=0; i<AMREX_SPACEDIM; i++) {
        umac[i].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(umac[i], i, geom, i);
        MultiFABPhysBCMacVel(umac[i], i, geom, i);
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
        MultiFab::Copy(umac_buffer[d], umac[d], 0, 0, 1, umac[d].nGrow());
        umac_buffer[d].FillBoundary(geom.periodicity());
    }

    ib_mc.ResetPredictor(0);
    ib_mc.InterpolatePredictor(0, umac_buffer);


    //___________________________________________________________________________
    // Move predictor using previous velocity: x^(n+1/2) = x^n + dt/2 J(u^n)
    // (constrain it to move in the z = constant plane only)
    constrain_ibm_marker(ib_mc, ib_lev, IBM_realData::pred_velz);
    ib_mc.MovePredictor(0, 0.5*dt);
    ib_mc.Redistribute(); // just in case (maybe can be removed)


    //___________________________________________________________________________
    // Update forces between markers: F^(n+1/2) = f(x^(n+1/2)) TODO: expensive
    // => use infrequently, use updateNeighbors for most steps
    ib_mc.clearNeighbors();
    ib_mc.fillNeighbors(); // Does ghost cells
    ib_mc.buildNeighborList(ib_mc.CheckPair);

    update_ibm_marker(driv_u, driv_amp, time, ib_mc, ib_lev, IBM_realData::pred_forcex, true);
    // Constrain it to move in the z = constant plane only
    constrain_ibm_marker(ib_mc, ib_lev, IBM_realData::pred_forcez);
    // Sum predictor forces added to neighbors back to the real markers
    ib_mc.sumNeighbors(IBM_realData::pred_forcex, AMREX_SPACEDIM, 0, 0);


    //___________________________________________________________________________
    // Spread forces to predictor f^(n+1/2) = S(F^(n+1/2))
    std::array<MultiFab, AMREX_SPACEDIM> fc_force_pred;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        fc_force_pred[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
        fc_force_pred[d].setVal(0.);
    }

    ib_mc.SpreadPredictor(0, fc_force_pred);
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        fc_force_pred[d].SumBoundary(geom.periodicity());


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
        MultiFABPhysBCDomainVel(uMom[d], d, geom, d);
        MultiFABPhysBCMacVel(uMom[d], d, geom, d);
    }

    // Compute advective fluxes: advFluxdiv = - D(\rho uu^n) = - D(u^n uMom)
    MkAdvMFluxdiv(umac, uMom, advFluxdiv, dx, 0);

    // Compute pressure gradient due to the BC: gp = Gp
    pres.setVal(0.); // Initial guess for pressure
    SetPressureBC(pres, geom); // Apply pressure boundary conditions
    ComputeGrad(pres, pg, 0, 0, 1, geom);

    // Construct RHS of Navier Stokes Equation
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);
        gmres_rhs_u[d].mult(dtinv*2, 1); // advance by dt/2

        MultiFab::Add(gmres_rhs_u[d], mfluxdiv_predict[d], 0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], fc_force_pred[d],    0, 0, 1, 1);

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
    GMRES(gmres_rhs_u, gmres_rhs_p, umacNew, pres,
          alpha_fc, beta_old, beta_ed_old, gamma_old, theta_alpha,
          geom, norm_pre_rhs);

    // Apply boundary conditions to the solution
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umacNew[d].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(umacNew[d], d, geom, d);
        MultiFABPhysBCMacVel(umacNew[d], d, geom, d);
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
        MultiFab::Copy(umacNew_buffer[d], umacNew[d], 0, 0, 1, umac[d].nGrow());
        umacNew_buffer[d].FillBoundary(geom.periodicity());
    }

    ib_mc.ResetMarkers(0);
    ib_mc.InterpolateMarkers(0, umacNew_buffer);


    //___________________________________________________________________________
    // Move markers according to velocity: x^(n+1) = x^n + dt/2 J(u^(n+1/2))
    // (constrain it to move in the z = constant plane only)
    constrain_ibm_marker(ib_mc, ib_lev, IBM_realData::velz);
    ib_mc.MoveMarkers(0, dt);
    ib_mc.Redistribute(); // Don't forget to send particles to the right CPU


    // //___________________________________________________________________________
    // // Update forces between markers: F^(n+1) = f(x^(n+1)) TODO: expensive =>
    // // use infrequently, use updateNeighbors for most steps
    // ib_mc.clearNeighbors();
    // ib_mc.fillNeighbors(); // Does ghost cells
    // ib_mc.buildNeighborList(ib_mc.CheckPair);

    // update_ibm_marker(driv_u, driv_amp, time, ib_mc, ib_lev, IBM_realData::forcex, false);
    // // Constrain it to move in the z = constant plane only
    // constrain_ibm_marker(ib_mc, ib_lev, IBM_realData::forcez);
    // // Sum predictor forces added to neighbors back to the real markers
    // ib_mc.sumNeighbors(IBM_realData::forcex, AMREX_SPACEDIM, 0, 0);


    // //___________________________________________________________________________
    // // Spread forces to corrector: f^(n+1) = S(F^(n+1))
    // std::array<MultiFab, AMREX_SPACEDIM> fc_force_corr;
    // for (int d=0; d<AMREX_SPACEDIM; ++d){
    //     fc_force_corr[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
    //     fc_force_corr[d].setVal(0.);
    // }

    // ib_mc.SpreadMarkers(0, fc_force_corr);
    // for (int d=0; d<AMREX_SPACEDIM; ++d)
    //     fc_force_corr[d].SumBoundary(geom.periodicity());


    //___________________________________________________________________________
    // Compute corrector velocity field
    // (u^(n+1)-u^n)/(dt) + D(uu^(n+1/2)) + Gp = L(u^n+u^(n+1))/2 + f(n+1/2)
    //                                Du^(n+1) = 0

    // Compute momentum fluxes at the midpoint: uMom = rho*u^(n+1/2)
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        MultiFab::Copy(uMom[d], umacNew[d], 0, 0, 1, 0);

        // let rho = 1
        uMom[d].mult(1.0, 1);

        uMom[d].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(uMom[d], d, geom, d);
        MultiFABPhysBCMacVel(uMom[d], d, geom, d);
    }

    // Compute advective fluxes at the midpoint:
    // advFluxdivPred = - D(\rho uu^(n+1/2)) = - D(u^(n+1/2) uMom)
    MkAdvMFluxdiv(umacNew, uMom, advFluxdivPred, dx, 0);

    // Explicit part of the diffusive operator Lu^n/2. Note that we are using
    // the weighted coefficients (to deal witht he 1/2 part)
    StagApplyOp(beta_negwtd, gamma_negwtd, beta_ed_negwtd,
                umac, Lumac, alpha_fc_0, dx, theta_alpha);


    // Note that the pressure gradient due to the BC is left unchanged
    pres.setVal(0.); // Initial guess for pressure
    SetPressureBC(pres, geom); // Apply pressure boundary conditions

    // Construct RHS of Navier Stokes Equation
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        Lumac[d].FillBoundary(geom.periodicity());

        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);
        gmres_rhs_u[d].mult(dtinv, 1); // advance by dt

        MultiFab::Add(gmres_rhs_u[d], mfluxdiv_correct[d], 0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], advFluxdivPred[d],   0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], fc_force_pred[d],    0, 0, 1, 1);

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
    GMRES(gmres_rhs_u, gmres_rhs_p, umacNew, pres,
          alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
          geom, norm_pre_rhs);

    // Apply boundary conditions to the solution
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        umacNew[d].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(umacNew[d], d, geom, d);
        MultiFABPhysBCMacVel(umacNew[d], d, geom, d);
    }

    // Update solution, and we're done!
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        MultiFab::Copy(umac[d], umacNew[d], 0, 0, 1, 0);

    BL_PROFILE_VAR_STOP(advance);
}



// Crank-Nicolson Advance Subroutine
void advance_CN(std::array<MultiFab, AMREX_SPACEDIM >& umac,
                std::array<MultiFab, AMREX_SPACEDIM >& umacNew,
                MultiFab& pres,
                IBMarkerContainer & ib_mc,
                const std::array<MultiFab, AMREX_SPACEDIM>& mfluxdiv_predict,
                const std::array<MultiFab, AMREX_SPACEDIM>& mfluxdiv_correct,
                      std::array<MultiFab, AMREX_SPACEDIM>& alpha_fc,
                const MultiFab& beta, const MultiFab& gamma,
                const std::array<MultiFab, NUM_EDGE> & beta_ed,
                const Geometry geom, const Real& dt, Real time)
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
    // Pressure gradient at inflow/outflow
    std::array< MultiFab, AMREX_SPACEDIM > pg;

    for (int i=0; i<AMREX_SPACEDIM; i++) {
           gmres_rhs_u[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                 Lumac[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
            advFluxdiv[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
        advFluxdivPred[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                  uMom[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                    pg[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
    }

    // Tracer advective terms
    MultiFab advFluxdivS(ba, dmap, 1, 1);


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
    Real driv_amp = std::min(time*10, 1.);
    Print() << "driv_amp = " << driv_amp << std::endl;

    // I'm too impatient to wait... -JPB
    // Real driv_amp = 1.;



    /****************************************************************************
     *                                                                          *
     * Apply deterministic boundary conditions                                  *
     *                                                                          *
     ***************************************************************************/

    for (int i=0; i<AMREX_SPACEDIM; i++) {
        umac[i].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(umac[i], i, geom, i);
        MultiFABPhysBCMacVel(umac[i], i, geom, i);
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
        MultiFab::Copy(umac_buffer[d], umac[d], 0, 0, 1, umac[d].nGrow());
        umac_buffer[d].FillBoundary(geom.periodicity());
    }

    ib_mc.ResetPredictor(0);
    ib_mc.InterpolatePredictor(0, umac_buffer);


    //___________________________________________________________________________
    // Move predictor using previous velocity: x^(n+1/2) = x^n + dt/2 J(u^n)
    // (constrain it to move in the z = constant plane only)
    constrain_ibm_marker(ib_mc, ib_lev, IBM_realData::pred_velz);
    ib_mc.MovePredictor(0, 0.5*dt);
    ib_mc.Redistribute(); // just in case (maybe can be removed)


    //___________________________________________________________________________
    // Update forces between markers: F^(n+1/2) = f(x^(n+1/2)) TODO: expensive
    // => use infrequently, use updateNeighbors for most steps
    ib_mc.clearNeighbors();
    ib_mc.fillNeighbors(); // Does ghost cells
    ib_mc.buildNeighborList(ib_mc.CheckPair);

    update_ibm_marker(driv_u, driv_amp, time, ib_mc, ib_lev, IBM_realData::pred_forcex, true);
    // Constrain it to move in the z = constant plane only
    constrain_ibm_marker(ib_mc, ib_lev, IBM_realData::pred_forcez);
    // Sum predictor forces added to neighbors back to the real markers
    ib_mc.sumNeighbors(IBM_realData::pred_forcex, AMREX_SPACEDIM, 0, 0);


    //___________________________________________________________________________
    // Spread forces to predictor f^(n+1/2) = S(F^(n+1/2))
    std::array<MultiFab, AMREX_SPACEDIM> fc_force_pred;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        fc_force_pred[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
        fc_force_pred[d].setVal(0.);
    }

    ib_mc.SpreadPredictor(0, fc_force_pred);
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        fc_force_pred[d].SumBoundary(geom.periodicity());


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
        MultiFABPhysBCDomainVel(uMom[d], d, geom, d);
        MultiFABPhysBCMacVel(uMom[d], d, geom, d);
    }

    // Compute advective fluxes: advFluxdiv = - D(\rho uu^n) = - D(u^n uMom)
    MkAdvMFluxdiv(umac, uMom, advFluxdiv, dx, 0);

    // (crank-nicolson terms) Explicit part of the diffusive operator Lu^n/2.
    // Note that we are using the weighted coefficients (to deal with the 1/2
    // part)
    StagApplyOp(beta_negwtd, gamma_negwtd, beta_ed_negwtd, umac, Lumac,
                alpha_fc_0, dx, theta_alpha);


    // Compute pressure gradient due to the BC: gp = Gp
    pres.setVal(0.); // Initial guess for pressure
    SetPressureBC(pres, geom); // Apply pressure boundary conditions
    ComputeGrad(pres, pg, 0, 0, 1, geom);

    // Construct RHS of Navier Stokes Equation
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);
        gmres_rhs_u[d].mult(dtinv, 1); // advance by dt

        MultiFab::Add(gmres_rhs_u[d], mfluxdiv_predict[d], 0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], fc_force_pred[d],    0, 0, 1, 1);

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
    GMRES(gmres_rhs_u, gmres_rhs_p, umacNew, pres,
          alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
          geom, norm_pre_rhs);

    // Apply boundary conditions to the solution
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umacNew[d].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(umacNew[d], d, geom, d);
        MultiFABPhysBCMacVel(umacNew[d], d, geom, d);
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
        MultiFab::Copy(umacNew_buffer[d], umacNew[d], 0, 0, 1, umac[d].nGrow());
        umacNew_buffer[d].FillBoundary(geom.periodicity());
    }

    ib_mc.ResetMarkers(0);
    ib_mc.InterpolateMarkers(0, umacNew_buffer);


    //___________________________________________________________________________
    // Move markers according to velocity: x^(n+1) = x^n + dt/2 J(u^(n+1/2))
    // (constrain it to move in the z = constant plane only)
    constrain_ibm_marker(ib_mc, ib_lev, IBM_realData::velz);
    ib_mc.MoveMarkers(0, dt);
    ib_mc.Redistribute(); // Don't forget to send particles to the right CPU


    //___________________________________________________________________________
    // Update forces between markers: F^(n+1) = f(x^(n+1)) TODO: expensive =>
    // use infrequently, use updateNeighbors for most steps
    ib_mc.clearNeighbors();
    ib_mc.fillNeighbors(); // Does ghost cells
    ib_mc.buildNeighborList(ib_mc.CheckPair);

    update_ibm_marker(driv_u, driv_amp, time, ib_mc, ib_lev, IBM_realData::forcex, false);
    // Constrain it to move in the z = constant plane only
    constrain_ibm_marker(ib_mc, ib_lev, IBM_realData::forcez);
    // Sum predictor forces added to neighbors back to the real markers
    ib_mc.sumNeighbors(IBM_realData::forcex, AMREX_SPACEDIM, 0, 0);


    //___________________________________________________________________________
    // Spread forces to corrector: f^(n+1) = S(F^(n+1))
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
        MultiFABPhysBCDomainVel(uMom[d], d, geom, d);
        MultiFABPhysBCMacVel(uMom[d], d, geom, d);
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
    // StagApplyOp(beta_negwtd, gamma_negwtd, beta_ed_negwtd,
    //             umac, Lumac, alpha_fc_0, dx, theta_alpha);


    // Note that the pressure gradient due to the BC is left unchanged
    pres.setVal(0.); // Initial guess for pressure
    SetPressureBC(pres, geom); // Apply pressure boundary conditions

    // Construct RHS of Navier Stokes Equation
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        Lumac[d].FillBoundary(geom.periodicity());

        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);
        gmres_rhs_u[d].mult(dtinv, 1); // advance by dt

        MultiFab::Add(gmres_rhs_u[d], mfluxdiv_correct[d], 0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], advFluxdivPred[d],   0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], fc_force_corr[d],    0, 0, 1, 1);

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
    GMRES(gmres_rhs_u, gmres_rhs_p, umacNew, pres,
          alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
          geom, norm_pre_rhs);

    // Apply boundary conditions to the solution
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        umacNew[d].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(umacNew[d], d, geom, d);
        MultiFABPhysBCMacVel(umacNew[d], d, geom, d);
    }

    // Update solution, and we're done!
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        MultiFab::Copy(umac[d], umacNew[d], 0, 0, 1, 0);

    BL_PROFILE_VAR_STOP(advance);
}
