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
    Real spr_k = 10000.0 ; // spring constant

    // initial distance btw markers. TODO: Need to update depending on initial
    // coordinates.
    Real l_db = 0.05;

    // Parameters for calling bending force calculation
    Real driv_k = 10000.0; //bending stiffness
    RealVect driv_u = {0, 0, 1};

    Real driv_period = 0.01;

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


    // TODO: should we constrain it to move in the z = constant plane only?


    //___________________________________________________________________________
    // Move markers according to predictor velocity
    ib_mc.MovePredictor(0, dt);
    // no redistribute needed after updating predictor


    //___________________________________________________________________________
    // Update forces between markers
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

            ParticleType & mark = markers[i];

            // Get previous and next markers connected to current marker (if they exist)
            ParticleType * next_marker = NULL;
            ParticleType * prev_marker = NULL;

            int status = ib_mc.FindConnectedMarkers(markers, mark,
                                                    nbhd_data, nbhd,
                                                    nbhd_index,
                                                    prev_marker, next_marker);


            if (status == -1) Abort("status -1 particle detected in predictor!!! flee for your life!");

            // update spring forces
            if (status == 0) { // has next (p) and prev (m)

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
            }

            // update bending forces for curent, minus/prev, and next/plus
            if (status == 0) { // has next (p) and prev (m)

                // position vectors
                RealVect r, r_m, r_p;
                for(int d=0; d<AMREX_SPACEDIM; ++d) {
                    r[d]   = mark.pos(d) + mark.rdata(IBM_realData::pred_posx + d);
                    r_m[d] = prev_marker->pos(d) + prev_marker->rdata(IBM_realData::pred_posx + d);
                    r_p[d] = next_marker->pos(d) + next_marker->rdata(IBM_realData::pred_posx + d);
                }


                // Set bending forces to zero
                RealVect f_p = RealVect{0., 0., 0.};
                RealVect f   = RealVect{0., 0., 0.};
                RealVect f_m = RealVect{0., 0., 0.};

                // calling the bending force calculation
                Real theta = 0;
                driving_f(f, f_p, f_m, r, r_p, r_m, driv_u, driv_k, theta);

                // updating the force on the minus, current, and plus particles.
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    prev_marker->rdata(IBM_realData::pred_forcex + d) += f_m[d];
                    mark.rdata(IBM_realData::pred_forcex + d)         +=   f[d];
                    next_marker->rdata(IBM_realData::pred_forcex + d) += f_p[d];
                }
            }

            // Increment neighbor list
            int nn      = nbhd[nbhd_index];
            nbhd_index += nn + 1; // +1 <= because the first field contains nn
        }
    }

    // TODO: should we constrain it to move in the z = constant plane only?


    //___________________________________________________________________________
    // Spread forces to predictor
    std::array<MultiFab, AMREX_SPACEDIM> fc_force_pred;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        fc_force_pred[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
        fc_force_pred[d].setVal(0.);
    }

    // Spread predictor forces
    //ib_mc.fillNeighbors(); // Don't forget to fill neighbor particles. This may be redundant

    ib_mc.SpreadPredictor(0, fc_force_pred);
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        fc_force_pred[d].SumBoundary(geom.periodicity());


    for (int d=0; d<AMREX_SPACEDIM; d++) {
        Lumac[d].FillBoundary(geom.periodicity());

        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);

        gmres_rhs_u[d].mult(dtinv, 1);
        MultiFab::Add(gmres_rhs_u[d], mfluxdiv_predict[d], 0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], fc_force_pred[d],    0, 0, 1, 0);
    }

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


    // TODO: should we constrain it to move in the z = constant plane only?


    //___________________________________________________________________________
    // Move markers according to velocity
    ib_mc.MoveMarkers(0, dt);
    ib_mc.Redistribute(); // Don't forget to send particles to the right CPU


    //___________________________________________________________________________
    // Update forces between markers
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

            ParticleType & mark = markers[i];

            // Get previous and next markers connected to current marker (if they exist)
            ParticleType * next_marker = NULL;
            ParticleType * prev_marker = NULL;

            int status = ib_mc.FindConnectedMarkers(markers, mark,
                                                    nbhd_data, nbhd,
                                                    nbhd_index,
                                                    prev_marker, next_marker);

            if (status == -1) Abort("status -1 particle detected in corrector!!! flee for your life!");

            // update spring forces
            if (status == 0) { // has next (p) and prev (m)

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
            }

            // update bending forces for curent, minus/prev, and next/plus
            if (status == 0) { // has next (p) and prev (m)

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

                 // calling the bending force calculation
                 Real theta = 0;
                 driving_f(f, f_p, f_m, r, r_p, r_m, driv_u, driv_k, theta);

                 // updating the force on the minus, current, and plus particles.
                 for (int d=0; d<AMREX_SPACEDIM; ++d) {
                     prev_marker->rdata(IBM_realData::forcex + d) += f_m[d];
                     mark.rdata(IBM_realData::forcex + d)         +=   f[d];
                     next_marker->rdata(IBM_realData::forcex + d) += f_p[d];
                 }
            }

            // Increment neighbor list
            int nn      = nbhd[nbhd_index];
            nbhd_index += nn + 1; // +1 <= because the first field contains nn
         }
    }

    // TODO: should we constrain it to move in the z = constant plane only?


    //___________________________________________________________________________
    // Spread forces to corrector
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

    pres.setVal(0.);  // initial guess
    SetPressureBC(pres, geom);

    // call GMRES here
    GMRES(gmres_rhs_u, gmres_rhs_p, umacNew, pres,
          alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
          geom, norm_pre_rhs);

    for (int i=0; i<AMREX_SPACEDIM; i++)
        MultiFab::Copy(umac[i], umacNew[i], 0, 0, 1, 0);
}
