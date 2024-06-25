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
		    FhdParticleContainer & particles,
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

    // Interpolate particles for the cell body
    particles.ResetMarkers(0);
    particles.InterpolateMarkers(0, umacNew_buffer);
    //particles.InterpolateMarkersGpu(0, dx, umac, RealFaceCoords, check); //from Ions codes where it's done after Advance
    //particles.velNorm(); // need this?


    //// check velocity difference btw two anchor markers on flagellum and two anchor particles on cell body////// 
    Vector<Real> anchor_marker_vel(6); //storing 6 velocity components of 2 anchored markers on the flagellum
    Vector<Real> anchor_particle_vel(6); //storing 6 velocity components of 2 anchored particles on the cell body

    // getting velocities from two anchor markers on flagellum
    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {
        // Get marker data (local to current thread)
        TileIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = ib_mc.GetParticles(ib_lev).at(index).GetArrayOfStructs();

        long np = ib_mc.GetParticles(ib_lev).at(index).numParticles();

        for (int i = 0; i < np; ++i) {

            ParticleType & mark = markers[i];

            if(mark.idata(IBMInt::id_1) == 0) {
                anchor_marker_vel[0] = mark.rdata(IBMReal::velx);
                anchor_marker_vel[1] = mark.rdata(IBMReal::vely);
                anchor_marker_vel[2] = mark.rdata(IBMReal::velz);
            }
            if(mark.idata(IBMInt::id_1) == 1) {
                anchor_marker_vel[3] = mark.rdata(IBMReal::velx);
                anchor_marker_vel[4] = mark.rdata(IBMReal::vely);
                anchor_marker_vel[5] = mark.rdata(IBMReal::velz);
            }                
        }
    }

    // getting velocities from two anchor particles on cell body
    for (FhdParIter pti(particles, ib_lev); pti.isValid(); ++pti) {
        // Get particle data (local to current thread)
        TileIndex index(pti.index(), pti.LocalTileIndex());
        auto & markers = particles.GetParticles(ib_lev).at(index).GetArrayOfStructs();

        long np = particles.GetParticles(ib_lev).at(index).numParticles();

        for (int i = 0; i < np; ++i) {

            auto & mark = markers[i];

            if(mark.idata(FHD_intData::id_global) == 0)  { //anchor particle in the inner layer of the cell body
                anchor_particle_vel[0] = mark.rdata(IBMReal::velx);
                anchor_particle_vel[1] = mark.rdata(IBMReal::vely);
                anchor_particle_vel[2] = mark.rdata(IBMReal::velz);
            }

            if(mark.idata(FHD_intData::id_global) == 1)  { //anchor particle in the outer layer of the cell body
                anchor_particle_vel[3] = mark.rdata(IBMReal::velx);
                anchor_particle_vel[4] = mark.rdata(IBMReal::vely);
                anchor_particle_vel[5] = mark.rdata(IBMReal::velz);
            }
        }
    }

    amrex::Print() << "The velocity of the first anchor marker on flagellum: " << anchor_marker_vel[0] << " i + "
                   << anchor_marker_vel[1] << " j + " << anchor_marker_vel[2] << " k" << std::endl;
    amrex::Print() << "The velocity of the second anchor marker on flagellum: " << anchor_marker_vel[3] << " i + "
                   << anchor_marker_vel[4] << " j + " << anchor_marker_vel[5] << " k" << std::endl;

    amrex::Print() << "The velocity of the first anchor particle on cell body: " << anchor_particle_vel[0] << " i + "
                   << anchor_particle_vel[1] << " j + " << anchor_particle_vel[2] << " k" << std::endl;
    amrex::Print() << "The velocity of the second anchor particle on cell body: " << anchor_particle_vel[3] << " i + "
                   << anchor_particle_vel[4] << " j + " << anchor_particle_vel[5] << " k" << std::endl;

    Print() << "Time to check if anchor particles and anchor markers share the same velocity" << std::endl;
    //exit(0);

    //___________________________________________________________________________
    // Move markers according to velocity: x^(n+1) = x^n + dt/2 J(u^(n+1/2))
    // (constrain it to move in the z = constant plane only)
    constrain_ibm_marker(ib_mc, ib_lev, IBMReal::velz);
    constrain_ibm_marker(particles, ib_lev, IBMReal::velz);

    //if(immbdy::contains_fourier)
    //    anchor_first_marker(ib_mc, ib_lev, IBMReal::velx);
    ib_mc.MoveMarkers(0, dt);

    ib_mc.clearNeighbors(); // Important: clear neighbors before Redistribute
    ib_mc.Redistribute();   // Don't forget to send particles to the right CPU


    // Move particles on cell body
    particles.MoveMarkers(0, dt);
    // particles.MoveIonsCPP(dt, dx, dxp, geom, umac, efield, RealFaceCoords, source, sourceTemp, pparamPlaneList,
    //                           paramPlaneCount, 3);

    // Couple forces between ib_mc and paricles
    // 1. Move flagellum markers.  
    // 2. Move particles
    // 3. Make anchor particles have same position as flagellum anchorpoints
    
    // Real anchor_markers[6]; //storing 2 positions or 6 coordinates of 2 anchored markers
    // Vector<Real> anchor_markers = get_anchor_markers(ib_mc, ib_lev, IBMReal::pred_posx);
    // move_anchor_particles(particles, ib_lev, IBMReal::pred_posx, anchor_markers); //working on it
    //above two can be combined into one function...
    particles.Redistribute();
    
    // Print() << "successfully move anchor particles to anchor markers' positions" << std::endl;
    // exit(0);


    /////////////////////////////////////////////////////////////////
    // Force coupling steps between markers and particles:
    // 1. compute forces in ib_mc and particles seperately
    // 2. sum up forces from respective anchor points
    // 2.a add forces from ib_mc anchor to particle anchor
    // 2.b add forces from particle anchor to ib_mc anchor
    // 3. DON'T Double count (use original / independent forces in 2.a and 2.b)

    // Step 1 (flagellar markers):
    // Update forces between flagellar markers: F^(n+1) = f(x^(n+1)) TODO: expensive =>
    // use infrequently, use updateNeighbors for most steps
    ib_mc.clearNeighbors();
    ib_mc.fillNeighbors(); // Does ghost cells
    ib_mc.buildNeighborList(ib_mc.CheckPair);

    update_ibm_marker(driv_u, driv_amp, time, ib_mc, ib_lev,
                      IBMReal::forcex, false,
                      geom);

    // Constrain it to move in the z = constant plane only
    constrain_ibm_marker(ib_mc, ib_lev, IBMReal::forcez);

    //Step 1 (cell body particles):
    // set velx/y/z and forcex/y/z for each particle on cell body to zero
    particles.ResetMarkers(0);

    // Update bond forces on cell body particles
    particles.computeForcesBondGPU(simParticles);  //check Ion codes and update!!!!!!!!!!!!!!
						   
    
    //Step 2 (This may NOT be necessary as all the forces will be spread to fluid grid):??????????????????????
//    Vector<Real> anchor_marker_forces(6); //storing 6 force components of 2 anchored markers on the flagellum
//    Vector<Real> anchor_particle_forces(6); //storing 6 force components of 2 anchored particles on the cell body


    // getting forces from two anchor markers on flagellum
//    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {
//        // Get marker data (local to current thread)
//        TileIndex index(pti.index(), pti.LocalTileIndex());
//        AoS & markers = ib_mc.GetParticles(ib_lev).at(index).GetArrayOfStructs();

//        long np = ib_mc.GetParticles(ib_lev).at(index).numParticles();

//        for (int i = 0; i < np; ++i) {

//            ParticleType & mark = markers[i];

//            if(mark.idata(IBMInt::id_1) == 0) {
//		anchor_marker_forces[0] = mark.rdata(IBMReal::forcex);
//                anchor_marker_forces[1] = mark.rdata(IBMReal::forcey);
//                anchor_marker_forces[2] = mark.rdata(IBMReal::forcez);
//            }
//            if(mark.idata(IBMInt::id_1) == 1) {
//                anchor_marker_forces[3] = mark.rdata(IBMReal::forcex);
//                anchor_marker_forces[4] = mark.rdata(IBMReal::forcey);
//                anchor_marker_forces[5] = mark.rdata(IBMReal::forcez);
//            }                
//        }
//    }

    // getting forces from two anchor particles on cell body
//    for (FhdParIter pti(particles, ib_lev); pti.isValid(); ++pti) {
        // Get particle data (local to current thread)
//        TileIndex index(pti.index(), pti.LocalTileIndex());
//        auto & markers = particles.GetParticles(ib_lev).at(index).GetArrayOfStructs();

//        long np = particles.GetParticles(ib_lev).at(index).numParticles();

        //Real get_anchor_markers[6]; //storing 6 coordinates of 2 anchored markers

//        for (int i = 0; i < np; ++i) {

//            auto & mark = markers[i];

//            if(mark.idata(FHD_intData::id_global) == 0)  { //anchor particle in the inner layer of the cell body
//                anchor_particle_forces[0] = mark.rdata(IBMReal::forcex);
//                anchor_particle_forces[1] = mark.rdata(IBMReal::forcey);
//                anchor_particle_forces[2] = mark.rdata(IBMReal::forcez);
//            }

//            if(mark.idata(FHD_intData::id_global) == 1)  { //anchor particle in the outer layer of the cell body
//                anchor_particle_forces[3] = mark.rdata(IBMReal::forcex);
//                anchor_particle_forces[4] = mark.rdata(IBMReal::forcey);
//                anchor_particle_forces[5] = mark.rdata(IBMReal::forcez);
//            }
//        }
//    }
    
    /////////////////////////////////////////////////////////
    // Continue? Add above forces between two sets of anchors? 
    // //////////////////////////////////////////////////////


    // Sum predictor forces added to neighbors back to the real markers
    ib_mc.sumNeighbors(IBMReal::forcex, AMREX_SPACEDIM, 0, 0);
    //need this for particles???
    particles.sumNeighbors(IBMReal::forcex, AMREX_SPACEDIM, 0, 0); //check and update!!!!!!!!!

    
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

    // Spread particle forces
    particles.SpreadIonsGPU(dx, dxp, geom, umac, RealFaceCoords, efieldCC, source, sourceTemp); //needs update. copied from Ions code


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
