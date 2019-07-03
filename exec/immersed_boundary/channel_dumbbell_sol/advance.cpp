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
    Real l_db = 0.01;

    // Parameters for calling bending force calculation
    Real bend_k = 10000.0; //bending stiffness
    Real cos_theta0 = 1.0; //initial cos_theta value
    int mark1_id, mark1_cpu;

    RealVect f, f_p, f_m; // bending force vectors for current, plus/next, and minus/previous particles
    //RealVect r, r_p, r_m; // position vectors


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

    //Print() << "This is right before the predictor force particle loop" << std::endl;


    // // To simulate a beam bent by perpendicular flow, set the velocity of the FIRST TWO markers to zero
    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        PairIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = ib_mc.GetParticles(ib_lev).at(index).GetArrayOfStructs();

       long np = markers.size();

      // Print() << " # of particles np = " << np << std::endl;


        // search for the first marker created and set its velocity to zero
        for (int i = 0; i < np; ++i) {
            ParticleType & mark = markers[i];              
           
            /////Set all forces to zero. Get ready for updating //////// 
            mark.rdata(IBM_realData::pred_forcex) = 0.;
            mark.rdata(IBM_realData::pred_forcey) = 0.;
            mark.rdata(IBM_realData::pred_forcez) = 0.;

            //Print() << " set pred_forcex = " << mark.rdata(IBM_realData::pred_forcex) << std::endl;
            //Print() << " set pred_forcey = " << mark.rdata(IBM_realData::pred_forcey) << std::endl;
            //Print() << " set pred_forcez = " << mark.rdata(IBM_realData::pred_forcez) << std::endl;


            if (mark.idata(IBM_intData::id_0) == -1 && mark.idata(IBM_intData::cpu_0 == -1)) {
               mark.rdata(IBM_realData::pred_velx) = 0.;
               mark.rdata(IBM_realData::pred_vely) = 0.;
               mark.rdata(IBM_realData::pred_velz) = 0.;

               mark1_id  = mark.id();  // used below for searching for second particle created
               mark1_cpu = mark.cpu();
            }
        }

        // search for the second marker created and set its velocity to zero.
        for (int i = 0; i < np; ++i) {
            ParticleType & mark = markers[i];
            if (mark.idata(IBM_intData::id_0) == mark1_id && mark.idata(IBM_intData::cpu_0) == mark1_cpu) {
               mark.rdata(IBM_realData::pred_velx) = 0.;
               mark.rdata(IBM_realData::pred_vely) = 0.;
               mark.rdata(IBM_realData::pred_velz) = 0.;
            }
        }
    }


    //___________________________________________________________________________
    // Move markers according to predictor velocity
    ib_mc.MovePredictor(0, dt);


    //
    //
    // TODO: Update forces between predictor markers
    //
    //
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
 
        /////Set all forces to zero. Get ready for updating ////////
        //for (int i = 0; i < np; ++i) {
        //    ParticleType & mark = markers[i];
        //    mark.rdata(IBM_realData::pred_forcex) = 0.;
        //    mark.rdata(IBM_realData::pred_forcey) = 0.;
        //    mark.rdata(IBM_realData::pred_forcez) = 0.;
        //}

        // Set bending forces to zero
        f_p = RealVect{0., 0., 0.};
        f   = RealVect{0., 0., 0.};
        f_m = RealVect{0., 0., 0.};


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
            
            Print() << "force Particle #" << i+1 << "," << "force Neighbor Status = " << status << std::endl;
            Print() << "force current particle neighbor list:"<< nbhd[nbhd_index] << std::endl;
            // Print() << "force current particle ID and CPU:"<<mark.id()<<" and "<<mark.cpu()<<std::endl;
            //Print() << "force current particle ID_0 and CPU_0:"<<mark.idata(IBM_intData::id_0)<<" and "<<mark.idata(IBM_intData::cpu_0)<<std::endl; 

            if (status == 2 || status == 0) { // has prev, only update spring forces for current and prev/minus markers

                Real rx = mark.pos(0) + mark.rdata(IBM_realData::pred_posx)
                        - (prev_marker->pos(0) + prev_marker->rdata(IBM_realData::pred_posx));
                Real ry = mark.pos(1) + mark.rdata(IBM_realData::pred_posy)
                        - (prev_marker->pos(1) + prev_marker->rdata(IBM_realData::pred_posy));
                Real rz = mark.pos(2) + mark.rdata(IBM_realData::pred_posz)
                        - (prev_marker->pos(2) + prev_marker->rdata(IBM_realData::pred_posz));

                Real l2 = rx*rx + ry*ry + rz*rz;
                Real lp = std::sqrt(l2);

                Print() << "This is right before updating predictor force" << std::endl;

                Print() << " pred_forcex = " << mark.rdata(IBM_realData::pred_forcex) << std::endl;
                Print() << " pred_forcey = " << mark.rdata(IBM_realData::pred_forcey) << std::endl;
                Print() << " pred_forcez = " << mark.rdata(IBM_realData::pred_forcez) << std::endl;


                mark.rdata(IBM_realData::pred_forcex) -= spr_k * rx/lp*(lp-l_db);
                mark.rdata(IBM_realData::pred_forcey) -= spr_k * ry/lp*(lp-l_db);
                mark.rdata(IBM_realData::pred_forcez) -= spr_k * rz/lp*(lp-l_db);

                prev_marker->rdata(IBM_realData::pred_forcex) += spr_k * rx/lp*(lp-l_db);
                prev_marker->rdata(IBM_realData::pred_forcey) += spr_k * ry/lp*(lp-l_db);
                prev_marker->rdata(IBM_realData::pred_forcez) += spr_k * rz/lp*(lp-l_db);
  
                Print() << " spring forcex = " << spr_k * rx*(lp-l_db)/lp << std::endl;
                Print() << " spring forcey = " << spr_k * ry*(lp-l_db)/lp << std::endl;
                Print() << " spring forcez = " << spr_k * rz*(lp-l_db)/lp << std::endl;

                Print() << "This is right after updating predictor force" << std::endl;

                Print() << " pred_forcex = " << mark.rdata(IBM_realData::pred_forcex) << std::endl;
                Print() << " pred_forcey = " << mark.rdata(IBM_realData::pred_forcey) << std::endl;
                Print() << " pred_forcez = " << mark.rdata(IBM_realData::pred_forcez) << std::endl;


            } 
	    
            if (status == 0) { // has both prev and next, update bending forces for curent, minus/prev, and next/plus markers
                
                // position vectors
                RealVect r = RealVect{mark.rdata(IBM_realData::pred_posx),
                                      mark.rdata(IBM_realData::pred_posy),
                                      mark.rdata(IBM_realData::pred_posz)};
  
                RealVect r_m = RealVect{prev_marker->rdata(IBM_realData::pred_posx),
                                        prev_marker->rdata(IBM_realData::pred_posy),
                                        prev_marker->rdata(IBM_realData::pred_posz)};

                RealVect r_p = RealVect{next_marker->rdata(IBM_realData::pred_posx),
                                        next_marker->rdata(IBM_realData::pred_posy),
                                        next_marker->rdata(IBM_realData::pred_posz)};

                //calling the bending force calculation
                bending_f(f, f_p, f_m, r, r_p, r_m, bend_k, cos_theta0);

                // updating the force on the minus, current, and plus particles.
                mark.rdata(IBM_realData::pred_forcex) += f[0];
                mark.rdata(IBM_realData::pred_forcey) += f[1];
                mark.rdata(IBM_realData::pred_forcez) += f[2];

                prev_marker->rdata(IBM_realData::pred_forcex) += f_m[0];
                prev_marker->rdata(IBM_realData::pred_forcey) += f_m[1];
                prev_marker->rdata(IBM_realData::pred_forcez) += f_m[2];

                next_marker->rdata(IBM_realData::pred_forcex) += f_p[0];
                next_marker->rdata(IBM_realData::pred_forcey) += f_p[1];
                next_marker->rdata(IBM_realData::pred_forcez) += f_p[2];
           
                Print() << " pred bending force f = " << f << std::endl;
                Print() << " pred bending force f_m  = " << f_m << std::endl;
                Print() << " pred bending force f_p = " << f_p << std::endl;

            } 
            
            // Increment neighbor list
            int nn      = nbhd[nbhd_index];
            nbhd_index += nn + 1; // +1 <= because the first field contains nn
        }
    }


//Abort();

    //___________________________________________________________________________
    // Spread forces to predictor
    std::array<MultiFab, AMREX_SPACEDIM> fc_force_pred;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        fc_force_pred[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
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


    //Print() << "This is right before the predictor force particle loop" << std::endl;


    // // To simulate a beam bent by perpendicular flow, set the velocity of the FIRST TWO markers to zero
    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        PairIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = ib_mc.GetParticles(ib_lev).at(index).GetArrayOfStructs();

        long np = markers.size();

        // search for the first marker created and set its velocity to zero
        for (int i = 0; i < np; ++i) {
            ParticleType & mark = markers[i];

            mark.rdata(IBM_realData::forcex) = 0.;
            mark.rdata(IBM_realData::forcey) = 0.;
            mark.rdata(IBM_realData::forcez) = 0.;

            if (mark.idata(IBM_intData::id_0) == -1 && mark.idata(IBM_intData::cpu_0 == -1)) {
               mark.rdata(IBM_realData::velx) = 0.;
               mark.rdata(IBM_realData::vely) = 0.;
               mark.rdata(IBM_realData::velz) = 0.;

               mark1_id  = mark.id();  // used below for searching for second particle created
               mark1_cpu = mark.cpu();
            }
        }

        // search for the second marker created and set its velocity to zero.
        for (int i = 0; i < np; ++i) {
            ParticleType & mark = markers[i];
            if (mark.idata(IBM_intData::id_0) == mark1_id && mark.idata(IBM_intData::cpu_0) == mark1_cpu) {
               mark.rdata(IBM_realData::velx) = 0.;
               mark.rdata(IBM_realData::vely) = 0.;
               mark.rdata(IBM_realData::velz) = 0.;
            }
        }
    }


    //___________________________________________________________________________
    // Move markers according to velocity
    ib_mc.MoveMarkers(0, dt);
    ib_mc.Redistribute(); // Don't forget to send particles to the right CPU


    //
    //
    // TODO: Update forces between markers
    //
    //
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


        //// Set all forces to zero. Get ready for updating ////////
        //for (int i = 0; i < np; ++i) {
        //    ParticleType & mark = markers[i];
        //    mark.rdata(IBM_realData::forcex) = 0.;
        //    mark.rdata(IBM_realData::forcey) = 0.;
        //    mark.rdata(IBM_realData::forcez) = 0.;
       // }

        // Set bending forces to zero
        f_p = RealVect{0., 0., 0.};
        f   = RealVect{0., 0., 0.};
        f_m = RealVect{0., 0., 0.};

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

            if (status == 2 || status == 0) { // has prev, only update spring forces for current and prev/minus markers

                Real rx = mark.pos(0) - prev_marker->pos(0);
                Real ry = mark.pos(1) - prev_marker->pos(1);
                Real rz = mark.pos(2) - prev_marker->pos(2);

                Real l2 = rx*rx + ry*ry + rz*rz;
                Real lp = std::sqrt(l2);

                mark.rdata(IBM_realData::forcex) -= spr_k * rx/lp*(lp-l_db);
                mark.rdata(IBM_realData::forcey) -= spr_k * ry/lp*(lp-l_db);
                mark.rdata(IBM_realData::forcez) -= spr_k * rz/lp*(lp-l_db);

                prev_marker->rdata(IBM_realData::forcex) += spr_k * rx/lp*(lp-l_db);
                prev_marker->rdata(IBM_realData::forcey) += spr_k * ry/lp*(lp-l_db);
                prev_marker->rdata(IBM_realData::forcez) += spr_k * rz/lp*(lp-l_db);

                Print() << " corr forcex = " << mark.rdata(IBM_realData::forcex) << std::endl;
                Print() << " corr forcey = " << mark.rdata(IBM_realData::forcey) << std::endl;
                Print() << " corr forcez = " << mark.rdata(IBM_realData::forcez) << std::endl;

            } 

	    if (status == 0) { // has both prev and next, update bending forces for curent, minus/prev, and next/plus markers

                // position vectors
                RealVect r   = RealVect{mark.pos(0),mark.pos(1),mark.pos(2)};
                RealVect r_m = RealVect{prev_marker->pos(0),prev_marker->pos(1),prev_marker->pos(2)};
                RealVect r_p = RealVect{next_marker->pos(0),next_marker->pos(1),next_marker->pos(2)};

                //calling the bending force calculation
                bending_f(f, f_p, f_m, r, r_p, r_m, bend_k, cos_theta0);

                // updating the force on the minus, current, and plus markers.
                mark.rdata(IBM_realData::forcex) += f[0];
                mark.rdata(IBM_realData::forcey) += f[1];
                mark.rdata(IBM_realData::forcez) += f[2];

                prev_marker->rdata(IBM_realData::forcex) += f_m[0];
                prev_marker->rdata(IBM_realData::forcey) += f_m[1];
                prev_marker->rdata(IBM_realData::forcez) += f_m[2];

                next_marker->rdata(IBM_realData::forcex) += f_p[0];
                next_marker->rdata(IBM_realData::forcey) += f_p[1];
                next_marker->rdata(IBM_realData::forcez) += f_p[2];

                Print() << " corr bending force f = " << f << std::endl;
                Print() << " corr bending force f_m  = " << f_m << std::endl;
                Print() << " corr bending force f_p = " << f_p << std::endl;

                Print() << " corr TOTAL force fx = " << mark.rdata(IBM_realData::forcex) << std::endl;
                Print() << " corr TOTAL force fy = " << mark.rdata(IBM_realData::forcey) << std::endl;
                Print() << " corr TOTAL force fz = " << mark.rdata(IBM_realData::forcez) << std::endl;
 
                Print() << " corr TOTAL force f_mx = " << prev_marker->rdata(IBM_realData::forcex) << std::endl;
                Print() << " corr TOTAL force f_my = " << prev_marker->rdata(IBM_realData::forcey) << std::endl;
                Print() << " corr TOTAL force f_mz = " << prev_marker->rdata(IBM_realData::forcez) << std::endl;

                Print() << " corr TOTAL force f_px = " << next_marker->rdata(IBM_realData::forcex) << std::endl;
                Print() << " corr TOTAL force f_py = " << next_marker->rdata(IBM_realData::forcey) << std::endl;
                Print() << " corr TOTAL force f_pz = " << next_marker->rdata(IBM_realData::forcez) << std::endl;

            }

            // Increment neighbor list
            int nn      = nbhd[nbhd_index];
            nbhd_index += nn + 1; // +1 <= because the first field contains nn
         }
    }



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
