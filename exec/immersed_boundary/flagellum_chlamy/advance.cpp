#include "main_driver.H"

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

#include <IBMarkerContainer.H>
#include <IBMarkerMD.H>


using namespace amrex;
using namespace common;
using namespace gmres;
using namespace immbdy_md;


using ParticleVector = typename IBMarkerContainer::ParticleVector;


// argv contains the name of the inputs file entered at the command line
void advance(std::array< MultiFab, AMREX_SPACEDIM >& umac,
             std::array< MultiFab, AMREX_SPACEDIM >& umacNew,
             MultiFab& pres, MultiFab& tracer,
             IBMarkerContainer & ib_mc,
             const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_predict,
             const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_correct,
             std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
             const MultiFab& beta, const MultiFab& gamma,
             const std::array< MultiFab, NUM_EDGE >& beta_ed,
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

    for (int i=0; i<AMREX_SPACEDIM; i++) {
           gmres_rhs_u[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                 Lumac[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
            advFluxdiv[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
        advFluxdivPred[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                  uMom[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
    }

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
     * Immersed-Marker parameters                                               *
     *                                                                          *
     ***************************************************************************/

    int ib_lev = 0;

    // Parameters for spring force calculation
    Real spr_k = 100000.0 ; // spring constant

    // initial distance btw markers. TODO: Need to update depending on initial
    // coordinates.
    Real l_db = 0.025;

    // Parameters for calling bending force calculation
    Real driv_k = 100000.0; //bending stiffness
    RealVect driv_u = {0, 0, 1};

    //Real driv_period = 100;  //This is actually angular frequency =  2*pi/T
    Real driv_period = 1.33/2.0;  // This is 50 times larger than the actual period
    Real length_flagellum = 0.5;
    //Real driv_amp = 15 * std::min(time*10, 1.);
    Real driv_amp = 15.0 * std::min(time*100, 1.);
    Print() << "driv_amp = " << driv_amp << std::endl;


    /****************************************************************************
     *                                                                          *
     * Apply non-stochastic boundary conditions                                 *
     *                                                                          *
     ***************************************************************************/

    for (int i=0; i<AMREX_SPACEDIM; i++) {
        umac[i].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(umac[i], i, geom, i);
        MultiFABPhysBCMacVel(umac[i], i, geom, i);
    }


    /****************************************************************************
     *                                                                          *
     * Advance tracer                                                           *
     *                                                                          *
     ***************************************************************************/

    // Compute tracer:
    BL_PROFILE_VAR("adv_compute tracer",TRACER);
    tracer.FillBoundary(geom.periodicity());
    MultiFABPhysBC(tracer, geom);

    MkAdvSFluxdiv(umac, tracer, advFluxdivS, dx, geom, 0);
    advFluxdivS.mult(dt, 1);
    BL_PROFILE_VAR_STOP(TRACER);

    // compute predictor
    BL_PROFILE_VAR("adv_compute predictor",PRED);

    MultiFab::Copy(tracerPred, tracer, 0, 0, 1, 0);
    MultiFab::Add(tracerPred, advFluxdivS, 0, 0, 1, 0);

    tracerPred.FillBoundary(geom.periodicity());
    MultiFABPhysBC(tracerPred, geom);

    MkAdvSFluxdiv(umac, tracerPred, advFluxdivS, dx, geom, 0);
    advFluxdivS.mult(dt, 1);
    BL_PROFILE_VAR_STOP(PRED);
    // advance in time
    MultiFab::Add(tracer, tracerPred,  0, 0, 1, 0);
    MultiFab::Add(tracer, advFluxdivS, 0, 0, 1, 0);
    tracer.mult(0.5, 1);

    // amrex::Print() << "tracer L0 norm = " << tracer.norm0() << "\n";
    //////////////////////////

    //////////////////////////////////////////////////
    // ADVANCE velocity field
    //////////////////////////////////////////////////

    // PREDICTOR STEP (heun's method: part 1)
    // compute advective term

    for (int i=0; i<AMREX_SPACEDIM; i++)
        MultiFab::Copy(uMom[i], umac[i], 0, 0, 1, 1);

    // let rho = 1
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        uMom[d].mult(1.0, 1);
    }

    for (int i=0; i<AMREX_SPACEDIM; i++) {
        uMom[i].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(uMom[i], i, geom, i);
        MultiFABPhysBCMacVel(uMom[i], i, geom, i);
    }

    MkAdvMFluxdiv(umac, uMom, advFluxdiv, dx, 0);

    // crank-nicolson terms
    StagApplyOp(beta_negwtd, gamma_negwtd, beta_ed_negwtd, umac, Lumac, alpha_fc_0, dx, theta_alpha);

    //___________________________________________________________________________
    // Interpolate immersed boundary predictor
    std::array<MultiFab, AMREX_SPACEDIM> umac_buffer;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        umac_buffer[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
        MultiFab::Copy(umac_buffer[d], umac[d], 0, 0, 1, umac[d].nGrow());
        umac_buffer[d].FillBoundary(geom.periodicity());
    }


    ib_mc.ResetPredictor(0);
    ib_mc.InterpolatePredictor(0, umac_buffer);


    // Constrain it to move in the z = constant plane only
    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        PairIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = ib_mc.GetParticles(ib_lev).at(index).GetArrayOfStructs();

        long np = markers.size();

        for (int i = 0; i < np; ++i) {

            ParticleType & mark = markers[i];
            // Zero z-velocity only
            mark.rdata(IBM_realData::pred_velz) = 0.;

	    // fix the first marker in x and y as well. can also use id_1==0
            if (mark.idata(IBM_intData::id_0) == -1) {
	    	mark.rdata(IBM_realData::pred_velx) = 0.;
                mark.rdata(IBM_realData::pred_vely) = 0.;
            }

        }
    }



    //___________________________________________________________________________
    // Move markers according to predictor velocity
    ib_mc.MovePredictor(0, dt);
    ib_mc.Redistribute(); // just in case (maybe can be removed)


    //___________________________________________________________________________
    // Update forces between markers
    // TODO: expensive => use infrequently, use updateNeighbors for most steps
    ib_mc.clearNeighbors();
    ib_mc.fillNeighbors(); // Does ghost cells
    ib_mc.buildNeighborList(ib_mc.CheckPair);


    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        PairIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = ib_mc.GetParticles(ib_lev).at(index).GetArrayOfStructs();

        // Get neighbor marker data (from neighboring threads)
        ParticleVector & nbhd_data = ib_mc.GetNeighbors(ib_lev, pti.index(),
                                                        pti.LocalTileIndex());


        // Get neighbor list (for collision checking)
        const Vector<int> & nbhd = ib_mc.GetNeighborList(ib_lev, pti.index(),
                                                         pti.LocalTileIndex());
        long np = markers.size();
        int nbhd_index = 0;

        for (int i=0; i<np; ++i) {

            BL_PROFILE_VAR("adv_find neighbors", PREDFINDNEIGBORS);

            ParticleType & mark = markers[i];

            // Get previous and next markers connected to current marker (if they exist)
            ParticleType * next_marker = NULL;
            ParticleType * prev_marker = NULL;

            int status = ib_mc.FindConnectedMarkers(markers, mark,
                                                    nbhd_data, nbhd,
                                                    nbhd_index,
                                                    prev_marker, next_marker);

            BL_PROFILE_VAR_STOP(PREDFINDNEIGBORS);

            if (status == -1) Abort("status -1 particle detected in predictor!!! flee for your life!");

            // update spring forces
            if (status == 0) { // has next (p) and prev (m)
                BL_PROFILE_VAR("adv_updating predictor spring forces", PREDUPDATESPRINGFORCES);

                RealVect r_p, r_m;
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    r_m[d] = mark.pos(d) + mark.rdata(IBM_realData::pred_posx + d)
                         - (prev_marker->pos(d) + prev_marker->rdata(IBM_realData::pred_posx + d));

                    r_p[d] = next_marker->pos(d) + next_marker->rdata(IBM_realData::pred_posx + d)
                         - (mark.pos(d) + mark.rdata(IBM_realData::pred_posx + d));
                }

                Real lp_m = r_m.vectorLength(),       lp_p = r_p.vectorLength();
                Real fm_0 = spr_k * (lp_m-l_db)/lp_m, fp_0 = spr_k * (lp_p-l_db)/lp_p;

                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    prev_marker->rdata(IBM_realData::pred_forcex + d) += fm_0 * r_m[d];
                    mark.rdata(IBM_realData::pred_forcex + d)         -= fm_0 * r_m[d];

                    mark.rdata(IBM_realData::pred_forcex + d)         += fp_0 * r_p[d];
                    next_marker->rdata(IBM_realData::pred_forcex + d) -= fp_0 * r_p[d];
                }

                BL_PROFILE_VAR_STOP(PREDUPDATESPRINGFORCES);
            }

            // update bending forces for curent, minus/prev, and next/plus
            if (status == 0) { // has next (p) and prev (m)

                BL_PROFILE_VAR("adv_Predictor bending forces",predictorbendingforces);

                // position vectors
                RealVect r, r_m, r_p;
                for(int d=0; d<AMREX_SPACEDIM; ++d) {
                    r[d]   = mark.pos(d) + mark.rdata(IBM_realData::pred_posx + d);
                    r_m[d] = prev_marker->pos(d) + prev_marker->rdata(IBM_realData::pred_posx + d);
                    r_p[d] = next_marker->pos(d) + next_marker->rdata(IBM_realData::pred_posx + d);
                }


                // Set bending forces to zero
                RealVect f_p = RealVect{AMREX_D_DECL(0., 0., 0.)};
                RealVect f   = RealVect{AMREX_D_DECL(0., 0., 0.)};
                RealVect f_m = RealVect{AMREX_D_DECL(0., 0., 0.)};

                // calling the active bending force calculation
                // This a simple since wave imposed
                //Real theta = l_db*driv_amp*sin(driv_period*time
                //             + 2*M_PI/length_flagellum*mark.idata(IBM_intData::id_1)*l_db);


                //Fourier series parameters based on normalized axial location (0.06-0.96)
                Real s = mark.idata(IBM_intData::id_1)/20+0.01; //note: id_1 = 0 for first marker created

		Real a0 = -0.2194*s*s*s*s + 0.5067*s*s*s - 0.3265*s*s + 0.1561*s - 0.1626;
		Real a1 = -0.1155*s*s*s*s + 1.0627*s*s*s - 1.8144*s*s + 0.4139*s + 0.1751;
		Real b1 = 0.2912*s*s*s*s - 0.712*s*s*s - 1.075*s*s + 1.6789*s - 0.3073;
		Real a2 = -0.1537*s*s*s*s - 0.1021*s*s*s + 0.1863*s*s + 0.08*s - 0.041;
		Real b2 = 0.1462*s*s*s*s - 0.6082*s*s*s + 1.0747*s*s - 0.5667*s + 0.0623;
		Real a3 = -0.2303*s*s*s*s + 0.5034*s*s*s - 0.5192*s*s + 0.268*s - 0.0419;
		Real b3 = 0.1073*s*s*s*s - 0.1903*s*s*s + 0.1145*s*s + 0.0115*s - 0.022;

                Real w  = 2*M_PI; //for normailized beat period of 1

                Real T = driv_period;  // T=13.3 ms for 75.2Hz beating

                Real theta_raw = a0 + a1*cos(w*time/T)   + b1*sin(w*time/T)
                                    + a2*cos(2*w*time/T) + b2*sin(2*w*time/T)
                                    + a3*cos(3*w*time/T) + b3*sin(3*w*time/T);

                Real theta = l_db*driv_amp*theta_raw;

		//std::cout << "predictor theta = " << theta << std::endl;

                driving_f(f, f_p, f_m, r, r_p, r_m, driv_u, theta, driv_k);

                // updating the force on the minus, current, and plus particles.
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    prev_marker->rdata(IBM_realData::pred_forcex + d) += f_m[d];
                    mark.rdata(IBM_realData::pred_forcex + d)         +=   f[d];
                    next_marker->rdata(IBM_realData::pred_forcex + d) += f_p[d];
                }

                BL_PROFILE_VAR_STOP(predictorbendingforces);
            }

            // Increment neighbor list
            int nn      = nbhd[nbhd_index];
            nbhd_index += nn + 1; // +1 <= because the first field contains nn
        }
    }

    // Constrain it to move in the z = constant plane only
    // Set the forces in the z direction to zero
    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {

        BL_PROFILE_VAR("adv_Constrain z pred", CONSTRAINZPRED);

        PairIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = ib_mc.GetParticles(ib_lev).at(index).GetArrayOfStructs();

        long np = markers.size();

        for (int i = 0; i < np; ++i) {

            ParticleType & mark = markers[i];

            // Zero z-force only
            mark.rdata(IBM_realData::pred_forcez) = 0.;

            // fix the first marker in x and y as well. can also use id_1==0
            if (mark.idata(IBM_intData::id_0) == -1) {
                mark.rdata(IBM_realData::pred_forcex) = 0.;
                mark.rdata(IBM_realData::pred_forcey) = 0.;
            }

        }

        BL_PROFILE_VAR_STOP(CONSTRAINZPRED);
    }

    // Sum predictor forces added to neighbors back to the real markers
    ib_mc.sumNeighbors(IBM_realData::pred_forcex, AMREX_SPACEDIM, 0, 0);


    //___________________________________________________________________________
    // Spread forces to predictor
    BL_PROFILE_VAR("adv_spread forces to predictor", spreadpredfor);

    std::array<MultiFab, AMREX_SPACEDIM> fc_force_pred;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        fc_force_pred[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
        fc_force_pred[d].setVal(0.);
    }

    BL_PROFILE_VAR_STOP(spreadpredfor);

	BL_PROFILE_VAR("adv_spread predictor forces", SREADPREDICTORFORCES);

    ib_mc.SpreadPredictor(0, fc_force_pred);
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        fc_force_pred[d].SumBoundary(geom.periodicity());
	}

	BL_PROFILE_VAR_STOP(SREADPREDICTORFORCES);

	BL_PROFILE_VAR("adv_fill neighbors", filltheneighbors);

    for (int d=0; d<AMREX_SPACEDIM; d++) {
        Lumac[d].FillBoundary(geom.periodicity());

        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);

        gmres_rhs_u[d].mult(dtinv, 1);
        MultiFab::Add(gmres_rhs_u[d], mfluxdiv_predict[d], 0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], fc_force_pred[d],    0, 0, 1, 0);
    }

	BL_PROFILE_VAR_STOP(filltheneighbors);


    std::array< MultiFab, AMREX_SPACEDIM > pg;
    for (int i=0; i<AMREX_SPACEDIM; i++)
        pg[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);

    pres.setVal(0.);  // initial guess
    SetPressureBC(pres, geom);

    ComputeGrad(pres, pg, 0, 0, 1, geom);

    for (int i=0; i<AMREX_SPACEDIM; i++) {
        pg[i].FillBoundary(geom.periodicity());
        gmres_rhs_u[i].FillBoundary(geom.periodicity());

        MultiFab::Subtract(gmres_rhs_u[i], pg[i], 0, 0, 1, 1);
    }

    // initial guess for new solution
    for (int i=0; i<AMREX_SPACEDIM; i++)
        MultiFab::Copy(umacNew[i], umac[i], 0, 0, 1, 1);

    // call GMRES to compute predictor
    GMRES(gmres_rhs_u, gmres_rhs_p, umacNew, pres,
          alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
          geom, norm_pre_rhs);

	BL_PROFILE_VAR("adv_compute predictor advective", computepredictorad);

    // Compute predictor advective term
    // let rho = 1
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        umacNew[d].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(umacNew[d], d, geom, d);
        MultiFABPhysBCMacVel(umacNew[d], d, geom, d);

        MultiFab::Copy(uMom[d], umacNew[d], 0, 0, 1, 0);
        uMom[d].mult(1.0, 1);

        uMom[d].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(uMom[d], d, geom, d);
        MultiFABPhysBCMacVel(uMom[d], d, geom, d);
    }

    BL_PROFILE_VAR_STOP(computepredictorad);


    MkAdvMFluxdiv(umacNew,uMom,advFluxdivPred,dx,0);

    // ADVANCE STEP (crank-nicolson + heun's method)

    // Compute gmres_rhs

    // trapezoidal advective terms
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        advFluxdiv[d].mult(0.5, 1);
        advFluxdivPred[d].mult(0.5, 1);
    }

    // crank-nicolson terms
    StagApplyOp(beta_negwtd, gamma_negwtd, beta_ed_negwtd, umac, Lumac, alpha_fc_0, dx, theta_alpha);



    //___________________________________________________________________________
    // Interpolate immersed boundary
    std::array<MultiFab, AMREX_SPACEDIM> umacNew_buffer;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        umacNew_buffer[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
        MultiFab::Copy(umacNew_buffer[d], umacNew[d], 0, 0, 1, umac[d].nGrow());
        umacNew_buffer[d].FillBoundary(geom.periodicity());
    }

    ib_mc.ResetMarkers(0);
    ib_mc.InterpolateMarkers(0, umacNew_buffer);


    // Constrain it to move in the z = constant plane only
    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        PairIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = ib_mc.GetParticles(ib_lev).at(index).GetArrayOfStructs();

        long np = markers.size();

        for (int i = 0; i < np; ++i) {

            ParticleType & mark = markers[i];
            // Zero z-velocity only
            mark.rdata(IBM_realData::velz) = 0.;

            // fix the first marker in x and y as well. can also use id_1==0
            if (mark.idata(IBM_intData::id_0) == -1) {
                mark.rdata(IBM_realData::velx) = 0.;
                mark.rdata(IBM_realData::vely) = 0.;
            }
        }
    }


    //___________________________________________________________________________
    // Move markers according to velocity
    ib_mc.MoveMarkers(0, dt);
    ib_mc.Redistribute(); // Don't forget to send particles to the right CPU


    //___________________________________________________________________________
    // Update forces between markers
    // TODO: expensive => use infrequently, use updateNeighbors for most steps
    ib_mc.clearNeighbors();
    ib_mc.fillNeighbors(); // Does ghost cells
    ib_mc.buildNeighborList(ib_mc.CheckPair);


    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {

        BL_PROFILE_VAR("adv_grab marker data", grabmarkerdata);

        // Get marker data (local to current thread)
        PairIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = ib_mc.GetParticles(ib_lev).at(index).GetArrayOfStructs();

        // Get neighbor marker data (from neighboring threads)
        ParticleVector & nbhd_data = ib_mc.GetNeighbors(ib_lev, pti.index(),
                                                        pti.LocalTileIndex());


        // Get neighbor list (for collision checking)
        const Vector<int> & nbhd = ib_mc.GetNeighborList(ib_lev, pti.index(),
                                                         pti.LocalTileIndex());

        long np = markers.size();
        int nbhd_index = 0;

        BL_PROFILE_VAR_STOP(grabmarkerdata);

        for (int i=0; i<np; ++i) {

            BL_PROFILE_VAR("adv_match markers", matchmarkerdata);

            ParticleType & mark = markers[i];

            // Get previous and next markers connected to current marker (if they exist)
            ParticleType * next_marker = NULL;
            ParticleType * prev_marker = NULL;

            int status = ib_mc.FindConnectedMarkers(markers, mark,
                                                    nbhd_data, nbhd,
                                                    nbhd_index,
                                                    prev_marker, next_marker);

            BL_PROFILE_VAR_STOP(matchmarkerdata);

            if (status == -1) Abort("status -1 particle detected in corrector!!! flee for your life!");

            // update spring forces
            if (status == 0) { // has next (p) and prev (m)

                BL_PROFILE_VAR("adv_updatingspringforces", UPDATESRPINGFORCES);

                RealVect r_p, r_m;
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    r_m[d] = mark.pos(d) - prev_marker->pos(d);

                    r_p[d] = next_marker->pos(d) - mark.pos(d);
                }

                Real lp_m = r_m.vectorLength(),       lp_p = r_p.vectorLength();
                Real fm_0 = spr_k * (lp_m-l_db)/lp_m, fp_0 = spr_k * (lp_p-l_db)/lp_p;

                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    prev_marker->rdata(IBM_realData::forcex + d) += fm_0 * r_m[d];
                    mark.rdata(IBM_realData::forcex + d)         -= fm_0 * r_m[d];

                    mark.rdata(IBM_realData::forcex + d)         += fp_0 * r_p[d];
                    next_marker->rdata(IBM_realData::forcex + d) -= fp_0 * r_p[d];
                }

                BL_PROFILE_VAR_STOP(UPDATESRPINGFORCES);
            }

            // update bending forces for curent, minus/prev, and next/plus
            if (status == 0) { // has next (p) and prev (m)

                BL_PROFILE_VAR("adv_update bending forces", UPDATEBENDINGFORCES);
                //
                // position vectors
                RealVect r, r_m, r_p;
                for(int d=0; d<AMREX_SPACEDIM; ++d) {
                    r[d]   = mark.pos(d);
                    r_m[d] = prev_marker->pos(d);
                    r_p[d] = next_marker->pos(d);
                }

                // Set bending forces to zero
                RealVect f_p = RealVect{0., 0., 0.};
                RealVect f   = RealVect{0., 0., 0.};
                RealVect f_m = RealVect{0., 0., 0.};

                // calling the active bending force calculation
                //Real theta = l_db*driv_amp*sin(driv_period*time
                //            + 2*M_PI/length_flagellum*mark.idata(IBM_intData::id_1)*l_db);

                //Fourier series parameters based on normalized axial location (0.06-0.96)
                Real s = mark.idata(IBM_intData::id_1)/20+0.01; //note: id_1 = 0 for first marker created

                Real a0 = -0.2194*s*s*s*s + 0.5067*s*s*s - 0.3265*s*s + 0.1561*s - 0.1626;
                Real a1 = -0.1155*s*s*s*s + 1.0627*s*s*s - 1.8144*s*s + 0.4139*s + 0.1751;
                Real b1 = 0.2912*s*s*s*s - 0.712*s*s*s - 1.075*s*s + 1.6789*s - 0.3073;
                Real a2 = -0.1537*s*s*s*s - 0.1021*s*s*s + 0.1863*s*s + 0.08*s - 0.041;
                Real b2 = 0.1462*s*s*s*s - 0.6082*s*s*s + 1.0747*s*s - 0.5667*s + 0.0623;
                Real a3 = -0.2303*s*s*s*s + 0.5034*s*s*s - 0.5192*s*s + 0.268*s - 0.0419;
                Real b3 = 0.1073*s*s*s*s - 0.1903*s*s*s + 0.1145*s*s + 0.0115*s - 0.022;

                Real w  = 2*M_PI; //for normailized beat period of 1

                Real T = driv_period;  // T=13.3 ms for 75.2Hz beating

                Real theta_raw = a0 + a1*cos(w*time/T)   + b1*sin(w*time/T)
                                    + a2*cos(2*w*time/T) + b2*sin(2*w*time/T)
                                    + a3*cos(3*w*time/T) + b3*sin(3*w*time/T);

                Real theta = l_db*driv_amp*theta_raw;

                driving_f(f, f_p, f_m, r, r_p, r_m, driv_u, theta, driv_k);

                // updating the force on the minus, current, and plus particles.
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    prev_marker->rdata(IBM_realData::forcex + d) += f_m[d];
                    mark.rdata(IBM_realData::forcex + d)         +=   f[d];
                    next_marker->rdata(IBM_realData::forcex + d) += f_p[d];
                }

                BL_PROFILE_VAR_STOP(UPDATEBENDINGFORCES);
            }

            // Increment neighbor list
            int nn      = nbhd[nbhd_index];
            nbhd_index += nn + 1; // +1 <= because the first field contains nn
         }
    }


    // Constrain it to move in the z = constant plane only
    // Set the forces in the z direction to zero
    BL_PROFILE_VAR("adv_contstrain z", CONSTRAINZ);

    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = ib_mc.GetParticles(ib_lev).at(index).GetArrayOfStructs();

        long np = markers.size();

        for (int i = 0; i < np; ++i) {

            ParticleType & mark = markers[i];

            // Zero z-force only
            mark.rdata(IBM_realData::forcez) = 0.;

            // fix the first marker in x and y as well. can also use id_1==0
            if (mark.idata(IBM_intData::id_0) == -1) {
                mark.rdata(IBM_realData::forcex) = 0.;
                mark.rdata(IBM_realData::forcey) = 0.;
            }

        }

    }

    BL_PROFILE_VAR_STOP(CONSTRAINZ);


    // Sum predictor forces added to neighbors back to the real markers
    ib_mc.sumNeighbors(IBM_realData::forcex, AMREX_SPACEDIM, 0, 0);



    // Apply Heun's Method time stepping to marker forces
    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = ib_mc.GetParticles(ib_lev).at(index).GetArrayOfStructs();

        long np = markers.size();

        for (int i = 0; i < np; ++i) {

            ParticleType & mark = markers[i];

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                mark.rdata(IBM_realData::forcex + d) = 0.5*mark.rdata(IBM_realData::forcex + d)
                    + 0.5*mark.rdata(IBM_realData::pred_forcex + d);
            }
        }

    }




    //___________________________________________________________________________
    // Spread forces to corrector
    BL_PROFILE_VAR("adv_spread forces to corrector", corrector);

    std::array<MultiFab, AMREX_SPACEDIM> fc_force_corr;
    for (int d=0; d<AMREX_SPACEDIM; ++d){

        fc_force_corr[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
        fc_force_corr[d].setVal(0.);

    }

    // Spread to the `fc_force` multifab
    //ib_mc.fillNeighbors(); // Don't forget to fill neighbor particles

    ib_mc.SpreadMarkers(0, fc_force_corr);
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        fc_force_corr[d].SumBoundary(geom.periodicity());

    BL_PROFILE_VAR_STOP(corrector);


    BL_PROFILE_VAR("adv_fillboundary", fillboundary);

    for (int d=0; d<AMREX_SPACEDIM; d++) {
        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);

        gmres_rhs_u[d].mult(dtinv, 1);
        MultiFab::Add(gmres_rhs_u[d], mfluxdiv_correct[d], 0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], advFluxdivPred[d],   0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], fc_force_corr[d],    0, 0, 1, 0);

        gmres_rhs_u[d].FillBoundary(geom.periodicity());

        MultiFab::Subtract(gmres_rhs_u[d], pg[d], 0, 0, 1, 1);

        MultiFab::Copy(umacNew[d], umac[d], 0, 0, 1, 0);
    }

    BL_PROFILE_VAR_STOP(fillboundary);


    pres.setVal(0.);  // initial guess
    SetPressureBC(pres, geom);

    // call GMRES here
    GMRES(gmres_rhs_u, gmres_rhs_p, umacNew, pres,
          alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
          geom, norm_pre_rhs);

    for (int i=0; i<AMREX_SPACEDIM; i++)
        MultiFab::Copy(umac[i], umacNew[i], 0, 0, 1, 0);
  
    BL_PROFILE_VAR_STOP(advance);
}
