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
void advance(AmrCoreAdv & amr_core_adv,
             std::array<MultiFab, AMREX_SPACEDIM> & umac,
             std::array<MultiFab, AMREX_SPACEDIM> & umacNew,
             MultiFab & pres, MultiFab & tracer,
             std::array<MultiFab, AMREX_SPACEDIM> & force_ibm,
             IBMarkerMap & ib_forces,
             const std::array<MultiFab, AMREX_SPACEDIM> & mfluxdiv_predict,
             const std::array<MultiFab, AMREX_SPACEDIM> & mfluxdiv_correct,
             std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
             const MultiFab & beta, const MultiFab & gamma,
             const std::array<MultiFab, NUM_EDGE> & beta_ed,
             IBParticleContainer & ib_pc,
             IBCore & ib_core,
             const Geometry geom, const Real & dt, Real time)
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
    // Staggered pressure gradient
    std::array< MultiFab, AMREX_SPACEDIM > pg;

    for (int i=0; i<AMREX_SPACEDIM; i++) {
           gmres_rhs_u[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                 Lumac[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
            advFluxdiv[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
        advFluxdivPred[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                  uMom[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                    pg[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
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
     * Apply non-stochastic boundary conditions                                 *
     *                                                                          *
     ***************************************************************************/

    for (int i=0; i<AMREX_SPACEDIM; i++) {
        umac[i].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(umac[i], i, geom, i);
        MultiFabPhysBCMacVel(umac[i], i, geom, i);
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

    int ibpc_lev = 0; // assume single level for now
    int ib_grow  = 6; // using the 6-point stencil

   // spring coefficient for the immersed boundary points
    Real spring_coefficient = 2e2;

    int nstep=0;

   // Store copies of the concentration surface gradient for predictor
        Vector< std::unique_ptr<MultiFab> > Dc_x0(ibpc_lev+1);
        Vector< std::unique_ptr<MultiFab> > Dc_y0(ibpc_lev+1);
#if (AMREX_SPACEDIM == 3)

        Vector< std::unique_ptr<MultiFab> > Dc_z0(ibpc_lev+1);
#endif
        //IBMarkerContainer ib_marker;
        // Get the locations of the interface and catalyst, as well as the level set and face coordinates for the predicted concentratio advection diffuision (AD) code
        const iMultiFab & iface0 = ib_core.get_TagInterface();
        const iMultiFab & catalyst0 = ib_core.get_TagCatalyst();
        const MultiFab  & LevelSet0=ib_core.get_LevelSet();
        const Vector<std::array<MultiFab, AMREX_SPACEDIM>> & FaceCoords=ib_pc.get_face_coords();
        // Indicates that this is a predictor step
        int corrector=0;

        // Copy previous time step's concentration gradient
         amr_core_adv.con_new_copy(ibpc_lev, Dc_x0, 1);
         amr_core_adv.con_new_copy(ibpc_lev, Dc_y0, 2);
#if (AMREX_SPACEDIM == 3)

         amr_core_adv.con_new_copy(ibpc_lev, Dc_z0, 3);
#endif


    const BoxArray & badpx           = Dc_x0[ibpc_lev]->boxArray();
    const DistributionMapping & dmdpx =Dc_x0[ibpc_lev]->DistributionMap();

    const BoxArray & badpy           = Dc_y0[ibpc_lev]->boxArray();
    const DistributionMapping & dmdpy =Dc_y0[ibpc_lev]->DistributionMap();
#if (AMREX_SPACEDIM == 3)

    const BoxArray & badpz           = Dc_z0[ibpc_lev]->boxArray();
    const DistributionMapping & dmdpz =Dc_z0[ibpc_lev]->DistributionMap();
#endif
    // Put all the componets of the concentration surface gradient into an array
    std::array< MultiFab, AMREX_SPACEDIM > DC_s0;

    DC_s0[0].define(badpx, dmdpx, 1, ib_grow);
    DC_s0[1].define(badpy, dmdpy, 1, ib_grow);

    DC_s0[0].setVal(0.);
    DC_s0[1].setVal(0.);

    DC_s0[0].copy(*Dc_x0[ibpc_lev],0,0,1,0,ib_grow);
    DC_s0[1].copy(*Dc_y0[ibpc_lev],0,0,1,0,ib_grow);

    DC_s0[0].FillBoundary(geom.periodicity());
    DC_s0[1].FillBoundary(geom.periodicity());

#if (AMREX_SPACEDIM == 3)
    DC_s0[2].define(badpz, dmdpz, 1, ib_grow);


    DC_s0[2].setVal(0.);

    DC_s0[2].copy(*Dc_z0[ibpc_lev],0,0,1,0,ib_grow);

    DC_s0[2].FillBoundary(geom.periodicity());
#endif



    //___________________________________________________________________________
    // Collect data on the immersed boundaries interacting with this rank

    Vector<IBP_info> ibp_info = ib_pc.IBParticleInfo(0, true);


    //___________________________________________________________________________
    // Storage data structures for immersed boundary markers

    Vector<ParticleIndex> part_indices(ibp_info.size());

    IBMarkerMap marker_pos;
    IBMarkerMap marker_vel;
    IBMarkerMap marker_pos_0;
    IBMarkerMap marker_delta_0;
    IBMarkerMap marker_force_0;
    IBMarkerMap marker_pos_1;
    IBMarkerMap marker_delta_1;
    IBMarkerMap marker_force_1;

    // Storage for concentration gradient inpterpolated to markers
    IBMarkerMap marker_DCs0;
    IBMarkerMap marker_DCs1;

    //___________________________________________________________________________
    // Build coefficients for IB spreading and interpolation operators

    for (int i=0; i<ibp_info.size(); ++i) {
        part_indices[i] = ibp_info[i].asPairIndex();

        // Pre-allocate particle arrays
        const Vector<RealVect> marker_positions = ib_pc.MarkerPositions(0, part_indices[i]);
        // ... initialized to (0..0) by default constructor
        marker_vel[part_indices[i]].resize(marker_positions.size());
        marker_delta_0[part_indices[i]].resize(marker_positions.size());
        marker_force_0[part_indices[i]].resize(marker_positions.size());
        marker_delta_1[part_indices[i]].resize(marker_positions.size());
        marker_force_1[part_indices[i]].resize(marker_positions.size());

        // The interpolated concentration surface gradients go here
        marker_DCs0[part_indices[i]].resize(marker_positions.size());
        marker_DCs1[part_indices[i]].resize(marker_positions.size());

        // Fill these with initial values
        marker_pos[part_indices[i]]   = marker_positions;
        marker_pos_0[part_indices[i]] = marker_positions;
        marker_pos_1[part_indices[i]] = marker_positions;
    }



    /****************************************************************************
     *                                                                          *
     * Advance immersed boundary markers (for predictor's force)                *
     *                                                                          *
     ***************************************************************************/


    // Explicit force terms used by GMRES to solve for predictor (0) and
    // corrector (1)
    std::array<MultiFab, AMREX_SPACEDIM> force_0;
    std::array<MultiFab, AMREX_SPACEDIM> force_1;

    for (int d=0; d<AMREX_SPACEDIM; d++) {
        force_0[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, umac[d].nGrow());
        force_1[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, umac[d].nGrow());
     }


    //___________________________________________________________________________
    // Predictor step: advect immersed boundary markers
    std::array<MultiFab, AMREX_SPACEDIM> umac_buffer;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        umac_buffer[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
        MultiFab::Copy(umac_buffer[d], umac[d], 0, 0, 1, umac[d].nGrow());
        umac_buffer[d].FillBoundary(geom.periodicity());
        DC_s0[d].FillBoundary(geom.periodicity());

    }

    for (const auto & pindex : part_indices) {
        auto & vel = marker_vel.at(pindex);

        ib_pc.InterpolateMarkers(ibpc_lev, pindex, vel, umac_buffer);
        // interpolate  concentration surface gradient on to markers
        auto & dcs = marker_DCs0.at(pindex);

        ib_pc.InterpolateMarkers(ibpc_lev, pindex, dcs, DC_s0);
   }
    for (const auto & pindex : part_indices) {
        const auto & vel   = marker_vel.at(pindex);
        const auto & pos   = marker_pos.at(pindex);
        const auto & f_0   = ib_forces.at(pindex);
              auto & pos_0 = marker_pos_0.at(pindex);
              auto & del_0 = marker_delta_0.at(pindex);
              auto & force = marker_force_0.at(pindex);
              auto & dcs = marker_DCs0.at(pindex);

        for (int i=0; i<vel.size(); ++i) {
            del_0[i] = dt*vel[i];
            pos_0[i] = pos[i] + del_0[i];
            // " force " is defined by the movement of the immersed boundary particles and the interpolated concentration gradient
            force[i] = f_0[i] - spring_coefficient*del_0[i]+scaling_factor*dcs[i];

        }
    }


    //___________________________________________________________________________
    // Add immersed-boundary forces to predictor's RHS

    for (int d=0; d<AMREX_SPACEDIM; ++d)
        force_0[d].setVal(0);

    for (const auto & pindex : part_indices) {
        const auto & force = marker_force_0.at(pindex);

        ib_pc.SpreadMarkers(ibpc_lev, pindex, force, force_0);


    }

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        force_0[d].FillBoundary(geom.periodicity());

        // MultiFab::Add(force_1[d], force_0[d], 0, 0, 1, 1);
        VisMF::Write(force_0[d], "force_0_" + std::to_string(d));
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
        umac_0[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, umac[d].nGrow());
        umac_1[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, umac[d].nGrow());
    }


    MultiFab p_0(ba, dmap, 1, 1);
    p_0.setVal(0.);

    MultiFab p_1(ba, dmap, 1, 1);
    p_1.setVal(0.);


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
        MultiFabPhysBCMacVel(uMom[d], d, geom, d);
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
        MultiFab::Copy(umac_0[d], umac[d], 0, 0, 1, umac[d].nGrow());
        MultiFab::Copy(umac_1[d], umac[d], 0, 0, 1, umac[d].nGrow());
        MultiFab::Add(force_0[d], force_ibm[d], 0, 0, 1, 1);
        MultiFab::Add(force_1[d], force_ibm[d], 0, 0, 1, 1);
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
          alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
          geom, norm_pre_rhs);

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(umac_0[d], umacNew[d], 0, 0, 1, 1);
        umac_0[d].FillBoundary(geom.periodicity());
        VisMF::Write(umac_0[d], "umac_0_"+std::to_string(d));
    }




    /****************************************************************************
     *                                                                          *
     * ADVANCE (CORRECTOR) STEP (crank-nicolson heun's method: part 2)          *
     *                                                                          *
     ***************************************************************************/

    //___________________________________________________________________________
    // Compute corrector's advective term (using predictor's fluid solution)

    // Store copies of the concentration surface gradient for corrector

        Vector< std::unique_ptr<MultiFab> > Dc_x1(ibpc_lev+1);
        Vector< std::unique_ptr<MultiFab> > Dc_y1(ibpc_lev+1);
#if (AMREX_SPACEDIM == 3)

        Vector< std::unique_ptr<MultiFab> > Dc_z1(ibpc_lev+1);
#endif
        //IBMarkerContainer ib_marker;
        // One step of EvolveChem to solve for predicted concentration
        amr_core_adv.EvolveChem(umac,umac,iface0, iface0, catalyst0, catalyst0, LevelSet0, LevelSet0, ibpc_lev, nstep,dt, time, diffcoeff, FaceCoords,corrector,source_strength);

        // Get the locations of the interface and catalyst, as well as the level set and face coordinates for the corrected concentration advection diffuision (AD) code

        const iMultiFab & catalyst1 = ib_core.get_TagCatalyst();

        const iMultiFab & iface1 = ib_core.get_TagInterface();
        const MultiFab  & LevelSet1=ib_core.get_LevelSet();
        const Vector<std::array<MultiFab, AMREX_SPACEDIM>> & FaceCoords1=ib_pc.get_face_coords();
        // Indicates we are on the corrector step
        corrector=1;

        // copying predicted concentration surface gradients
         amr_core_adv.con_new_copy(ibpc_lev, Dc_x1, 1);
         amr_core_adv.con_new_copy(ibpc_lev, Dc_y1, 2);
#if (AMREX_SPACEDIM == 3)

         amr_core_adv.con_new_copy(ibpc_lev, Dc_z1, 3);
#endif
    // putting predicted concentration surface gradients in an array

    std::array< MultiFab, AMREX_SPACEDIM > DC_s1;
    DC_s1[0].define(badpx, dmdpx, 1,ib_grow );
    DC_s1[1].define(badpy, dmdpy, 1, ib_grow);

    DC_s1[0].setVal(0.);
    DC_s1[1].setVal(0.);
    DC_s1[0].copy(*Dc_x1[ibpc_lev],0,0,1,0,ib_grow);
    DC_s1[1].copy(*Dc_y1[ibpc_lev],0,0,1,0,ib_grow);
    DC_s1[0].FillBoundary(geom.periodicity());
    DC_s1[1].FillBoundary(geom.periodicity());

#if (AMREX_SPACEDIM == 3)
    DC_s1[2].define(badpz, dmdpz, 1, ib_grow);

    DC_s1[2].setVal(0.);

    DC_s1[2].copy(*Dc_z1[ibpc_lev],0,0,1,0,ib_grow);

    DC_s1[2].FillBoundary(geom.periodicity());

#endif


     for (int d=0; d<AMREX_SPACEDIM; d++) {
        umacNew[d].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(umacNew[d], d, geom, d);
        MultiFabPhysBCMacVel(umacNew[d], d, geom, d);

        MultiFab::Copy(uMom[d], umacNew[d], 0, 0, 1, 1);

        // let rho = 1
        uMom[d].mult(1.0, 1);

        // momentum boundary conditions
        uMom[d].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(uMom[d], d, geom, d);
        MultiFabPhysBCMacVel(uMom[d], d, geom, d);
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
    // Corrector step: advect immersed boundary markers

    std::array<MultiFab, AMREX_SPACEDIM> umacNew_buffer;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        umacNew_buffer[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
        MultiFab::Copy(umacNew_buffer[d], umacNew[d], 0, 0, 1, umacNew[d].nGrow());
        umacNew_buffer[d].FillBoundary(geom.periodicity());
        DC_s1[d].FillBoundary(geom.periodicity());
    }

    for (const auto & pindex : part_indices) {
        auto & vel = marker_vel.at(pindex);

        ib_pc.InterpolateMarkers(ibpc_lev, pindex, vel, umacNew_buffer);\

        // interpolating predicted concentration surface gradient on to markers
        auto & dcs = marker_DCs1.at(pindex);

        ib_pc.InterpolateMarkers(ibpc_lev, pindex, dcs, DC_s1);

    }


    for (const auto & pindex : part_indices) {
        const auto & vel   = marker_vel.at(pindex);
        const auto & pos   = marker_pos.at(pindex);
              auto & f_0   = ib_forces.at(pindex);
              auto & pos_1 = marker_pos_1.at(pindex);
              auto & del_1 = marker_delta_1.at(pindex);
              auto & force = marker_force_1.at(pindex);
              auto & dcs = marker_DCs1.at(pindex);

        for (int i=0; i<vel.size(); ++i) {
            del_1[i] = dt*vel[i];
            pos_1[i] = pos[i] + del_1[i];
            // " force " is defined by the movement of the immersed boundary particles and the interpolated concentration gradient

            force[i] = f_0[i] - spring_coefficient*del_1[i]+scaling_factor*dcs[i];
            f_0[i]   = force[i];

        }
    }


    //___________________________________________________________________________
    // Add immersed-boundary forces to predictor's RHS

    for (int d=0; d<AMREX_SPACEDIM; ++d)
        force_1[d].setVal(0);

    for (const auto & pindex : part_indices) {
        const auto & force = marker_force_1.at(pindex);

        ib_pc.SpreadMarkers(ibpc_lev, pindex, force, force_1);

    }

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        force_1[d].FillBoundary(geom.periodicity());
       // force_1[d].setVal(0);
        VisMF::Write(force_1[d], "force_1_" + std::to_string(d));
    }


    //_______________________________________________________________________
    // Set up the RHS for the corrector

    for (int d=0; d<AMREX_SPACEDIM; d++) {
        // explicit part
        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);
        gmres_rhs_u[d].mult(dtinv, 1);
        // RHS 1/2(f_1+f_0) from Heun's method
        MultiFab::Add(force_1[d], force_0[d], 0, 0, 1, 1);
        force_1[d].mult(0.5,1);

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
    gmres_rhs_p.setVal(0.);
    GMRES(gmres_rhs_u, gmres_rhs_p, umacNew, p_1,
          alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
          geom, norm_pre_rhs);
        // One step of EvolveChem to solve for corrected concentration
        amr_core_adv.EvolveChem(umac,umac_0,iface0, iface1, catalyst0, catalyst1, LevelSet0, LevelSet1, ibpc_lev, nstep,dt, time, diffcoeff, FaceCoords,corrector,source_strength);

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        // Output velocity solution
        MultiFab::Copy(umac[d],      umacNew[d],   0, 0, 1, 1);

        // Output immersed-boundary forces
        MultiFab::Copy(force_ibm[d], force_1[d],   0, 0, 1, 1);

        // Output pressure solution
        MultiFab::Copy(pres,         p_1,          0, 0, 1, 1);

    }



}
