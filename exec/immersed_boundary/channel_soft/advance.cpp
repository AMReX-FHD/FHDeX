#include "main_driver.H"
#include "main_driver_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "hydro_functions.H"

#include "common_namespace.H"

#include "gmres_functions.H"
#include "gmres_functions_F.H"

#include "gmres_namespace.H"

#include "ib_functions.H"

#include "IBCore.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;
using namespace common;
using namespace gmres;

// argv contains the name of the inputs file entered at the command line
void advance(AmrCoreAdv & amr_core_adv,
             std::array<MultiFab, AMREX_SPACEDIM> & umac,
             std::array<MultiFab, AMREX_SPACEDIM> & umacNew,
             MultiFab & pres, MultiFab & tracer,
             std::array<MultiFab, AMREX_SPACEDIM> & force_ibm,
             std::array<MultiFab, AMREX_SPACEDIM> & DCs_spread,
             IBMarkerMap & ib_forces,
             const std::array<MultiFab, AMREX_SPACEDIM> & mfluxdiv_predict,
             const std::array<MultiFab, AMREX_SPACEDIM> & mfluxdiv_correct,
             const std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
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



    /****************************************************************************
     *                                                                          *
     * Generate immersed boundary data                                          *
     *                                                                          *
     ***************************************************************************/

    int ibpc_lev = 0; // assume single level for now
    int ib_grow  = 6; // using the 6-point stencil

    Real spring_coefficient = 1e4;

    int nstep=0;

    Print() << " Diff Coeff advance " << diffcoeff <<std::endl;

        Vector< std::unique_ptr<MultiFab> > Dc_x(ibpc_lev+1);
        Vector< std::unique_ptr<MultiFab> > Dc_y(ibpc_lev+1);
        Vector< std::unique_ptr<MultiFab> > Dc_z(ibpc_lev+1);

        //IBMarkerContainer ib_marker;
        // advection diffuision (AD) code
        const iMultiFab & iface = ib_core.get_TagInterface();
        const MultiFab  & LevelSet=ib_core.get_LevelSet();
        const Vector<std::array<MultiFab, AMREX_SPACEDIM>> & FaceCoords=ib_pc.get_face_coords();

        amrex::Print() << "Solving AD Eqn" << std::endl;

        amr_core_adv.EvolveChem(umac, iface, LevelSet, ibpc_lev, nstep,dt/2, time, diffcoeff, FaceCoords);

         amr_core_adv.con_new_copy(ibpc_lev, Dc_x, 1);
         amr_core_adv.con_new_copy(ibpc_lev, Dc_y, 2);
         amr_core_adv.con_new_copy(ibpc_lev, Dc_z, 3);

        amrex::Print() << "After Solving AD Eqn" << std::endl;

        amrex::Print() << "Initializing face centered box array" << std::endl;

    const BoxArray & badpx           = Dc_x[ibpc_lev]->boxArray();
    const DistributionMapping & dmdpx =Dc_x[ibpc_lev]->DistributionMap();

    const BoxArray & badpy           = Dc_y[ibpc_lev]->boxArray();
    const DistributionMapping & dmdpy =Dc_y[ibpc_lev]->DistributionMap();

    const BoxArray & badpz           = Dc_z[ibpc_lev]->boxArray();
    const DistributionMapping & dmdpz =Dc_z[ibpc_lev]->DistributionMap();

    amrex::Print() << "Creating an Array of Multifabs" << std::endl;

    std::array< MultiFab, AMREX_SPACEDIM > DC_s;
    amrex::Print() << "Defining Multifabs with face centered box arrays" << std::endl;
#if (AMREX_SPACEDIM == 2)
    amrex::Print() << "1st element" << std::endl;
  
    DC_s[0].define(badpx, dmdpx, 1, 0);
    amrex::Print() << "2nd element" << std::endl;
    DC_s[1].define(badpy, dmdpy, 1, 0);
    amrex::Print() << "3rd element" << std::endl;

#elif (AMREX_SPACEDIM == 3)
    amrex::Print() << "1st element" << std::endl;
    DC_s[0].define(badpx, dmdpx, 1, 0);
    amrex::Print() << "2nd element" << std::endl;
    DC_s[1].define(badpy, dmdpy, 1, 0);
    amrex::Print() << "3rd element" << std::endl;

    std::cout<< " Box array " << badpz << std::endl;
    std::cout<< " Distribution Map " << dmdpz << std::endl;
    DC_s[2].define(badpz, dmdpz, 1, 0);
    amrex::Print() << "After defining multifabs" << std::endl;

#endif


    DC_s[0].setVal(0.);
    DC_s[1].setVal(0.);
    DC_s[2].setVal(0.);
    amrex::Print() << "Copying gradient into array of multifabs" << std::endl;

    DC_s[0].copy(*Dc_x[ibpc_lev],0,0,1,0,0);
    DC_s[1].copy(*Dc_y[ibpc_lev],0,0,1,0,0);
    DC_s[2].copy(*Dc_z[ibpc_lev],0,0,1,0,0);


    std::array< MultiFab, AMREX_SPACEDIM > DCs_spread0;

#if (AMREX_SPACEDIM == 2)
    amrex::Print() << "1st element" << std::endl;

    DCs_spread0[0].define(badpx, dmdpx, 1, 1);
    amrex::Print() << "2nd element" << std::endl;
    DCs_spread0[1].define(badpy, dmdpy, 1, 1);
    amrex::Print() << "3rd element" << std::endl;

#elif (AMREX_SPACEDIM == 3)
    amrex::Print() << "1st element" << std::endl;
    DCs_spread0[0].define(badpx, dmdpx, 1, 1);
    amrex::Print() << "2nd element" << std::endl;
    DCs_spread0[1].define(badpy, dmdpy, 1, 1);
    amrex::Print() << "3rd element" << std::endl;

    std::cout<< " Box array " << badpz << std::endl;
    std::cout<< " Distribution Map " << dmdpz << std::endl;
    DCs_spread0[2].define(badpz, dmdpz, 1, 1);
    amrex::Print() << "After defining multifabs" << std::endl;

#endif


    DCs_spread0[0].setVal(0.);
    DCs_spread0[1].setVal(0.);
    DCs_spread0[2].setVal(0.);
    amrex::Print() << "Copying gradient into array of multifabs" << std::endl;

    DCs_spread0[0].copy(*Dc_x[ibpc_lev],0,0,1,0,1);
    DCs_spread0[1].copy(*Dc_y[ibpc_lev],0,0,1,0,1);
    DCs_spread0[2].copy(*Dc_z[ibpc_lev],0,0,1,0,1);


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
    amrex::Print() << "Creating data structures for ibmarkers to hold Concentration gradient" << std::endl;

    // Storage for concentration gradient inpterpolated to markers
    IBMarkerMap marker_DCs;
    IBMarkerMap marker_DCs0;

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

        amrex::Print() << "Allocating particle arrays for gradient "<<std::endl;
        marker_DCs[part_indices[i]].resize(marker_positions.size());
        marker_DCs0[part_indices[i]].resize(marker_positions.size());

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
    }

    for (const auto & pindex : part_indices) {
        auto & vel = marker_vel.at(pindex);

        ib_pc.InterpolateMarkers(ibpc_lev, pindex, vel, umac_buffer);
        std::cout<< "pindx "<<std::endl;

        auto & dcs = marker_DCs.at(pindex);
        std::cout << "dcs "<<std::endl;

        ib_pc.InterpolateMarkers(ibpc_lev, pindex, dcs, DC_s);
        //ib_pc.InterpolateMarkers(ibpc_lev, pindex, cy, Dc_y);
        //ib_pc.InterpolateMarkers(ibpc_lev, pindex, cz, Dc_z);

    }

    for (const auto & pindex : part_indices) {
        const auto & vel   = marker_vel.at(pindex);
        const auto & pos   = marker_pos.at(pindex);
        const auto & f_0   = ib_forces.at(pindex);
              auto & pos_0 = marker_pos_0.at(pindex);
              auto & del_0 = marker_delta_0.at(pindex);
              auto & force = marker_force_0.at(pindex);

        for (int i=0; i<vel.size(); ++i) {
            del_0[i] = dt*vel[i];
            pos_0[i] = pos[i] + del_0[i];
            force[i] = f_0[i] - spring_coefficient*del_0[i];

            if (i == 10)
                Print() << "predictor force[" << i << "] = " << force[i] << std::endl;
        }
    }


    //___________________________________________________________________________
    // Add immersed-boundary forces to predictor's RHS

    for (int d=0; d<AMREX_SPACEDIM; ++d)
        force_0[d].setVal(0);

    for (const auto & pindex : part_indices) {
        const auto & force = marker_force_0.at(pindex);

        ib_pc.SpreadMarkers(ibpc_lev, pindex, force, force_0);
        auto & dcs = marker_DCs.at(pindex);
        ib_pc.SpreadMarkers(ibpc_lev, pindex, dcs, DCs_spread0);

    }

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        force_0[d].FillBoundary(geom.periodicity());
        DCs_spread0[d].FillBoundary(geom.periodicity());

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
        MultiFABPhysBCDomainVel(uMom[d], d, geom, d);
        MultiFABPhysBCMacVel(uMom[d], d, geom, d);
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
        MultiFab::Copy(umac_0[d], umac[d], 0, 0, 1, umac[d].nGrow());
        MultiFab::Copy(umac_1[d], umac[d], 0, 0, 1, umac[d].nGrow());
        MultiFab::Add(force_0[d], force_ibm[d], 0, 0, 1, 1);
        MultiFab::Add(force_1[d], force_ibm[d], 0, 0, 1, 1);
    }

    MultiFab::Copy(p_0, pres, 0, 0, 1, 1);
    MultiFab::Copy(p_1, pres, 0, 0, 1, 1);


    //___________________________________________________________________________
    // Set up the RHS for the predictor
    amrex::Real scaling_factor=-0.1;
    
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        // explicit part
        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);
        gmres_rhs_u[d].mult(dtinv, 1);
        int cng=DCs_spread0[d].nGrow();
        DCs_spread0[d].mult(scaling_factor, cng);
        MultiFab::Add(gmres_rhs_u[d], mfluxdiv_predict[d], 0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], force_0[d],          0, 0, 1, 1);
        std::cout<<" Check Add "<<std::endl;
        MultiFab::Add(gmres_rhs_u[d], DCs_spread0[d],       0, 0, 1, 1);
        std::cout<<" Check Add after "<<std::endl;

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


        Vector< std::unique_ptr<MultiFab> > Dc_x0(ibpc_lev+1);
        Vector< std::unique_ptr<MultiFab> > Dc_y0(ibpc_lev+1);
        Vector< std::unique_ptr<MultiFab> > Dc_z0(ibpc_lev+1);

        //IBMarkerContainer ib_marker;
        // advection diffuision (AD) code
        const iMultiFab & iface0 = ib_core.get_TagInterface();
        const MultiFab  & LevelSet0=ib_core.get_LevelSet();
        const Vector<std::array<MultiFab, AMREX_SPACEDIM>> & FaceCoords0=ib_pc.get_face_coords();

        amrex::Print() << "Solving AD Eqn" << std::endl;

        amr_core_adv.EvolveChem(umacNew, iface0, LevelSet0, ibpc_lev, nstep,dt, time, diffcoeff, FaceCoords0);

         amr_core_adv.con_new_copy(ibpc_lev, Dc_x0, 1);
         amr_core_adv.con_new_copy(ibpc_lev, Dc_y0, 2);
         amr_core_adv.con_new_copy(ibpc_lev, Dc_z0, 3);

        amrex::Print() << "After Solving AD Eqn" << std::endl;

        amrex::Print() << "Initializing face centered box array" << std::endl;

    amrex::Print() << "Creating an Array of Multifabs" << std::endl;

    std::array< MultiFab, AMREX_SPACEDIM > DC_s0;
    amrex::Print() << "Defining Multifabs with face centered box arrays" << std::endl;
#if (AMREX_SPACEDIM == 2)
    amrex::Print() << "1st element" << std::endl;
  
    DC_s0[0].define(badpx, dmdpx, 1, 0);
    amrex::Print() << "2nd element" << std::endl;
    DC_s0[1].define(badpy, dmdpy, 1, 0);
    amrex::Print() << "3rd element" << std::endl;

#elif (AMREX_SPACEDIM == 3)
    amrex::Print() << "1st element" << std::endl;
    DC_s0[0].define(badpx, dmdpx, 1, 0);
    amrex::Print() << "2nd element" << std::endl;
    DC_s0[1].define(badpy, dmdpy, 1, 0);
    amrex::Print() << "3rd element" << std::endl;

    std::cout<< " Box array " << badpz << std::endl;
    std::cout<< " Distribution Map " << dmdpz << std::endl;
    DC_s0[2].define(badpz, dmdpz, 1, 0);
    amrex::Print() << "After defining multifabs" << std::endl;

#endif


    DC_s0[0].setVal(0.);
    DC_s0[1].setVal(0.);
    DC_s0[2].setVal(0.);
    amrex::Print() << "Copying gradient into array of multifabs" << std::endl;

    DC_s0[0].copy(*Dc_x[ibpc_lev],0,0,1,0,0);
    DC_s0[1].copy(*Dc_y[ibpc_lev],0,0,1,0,0);
    DC_s0[2].copy(*Dc_z[ibpc_lev],0,0,1,0,0);


    std::array< MultiFab, AMREX_SPACEDIM > DCs_spread1;

#if (AMREX_SPACEDIM == 2)
    amrex::Print() << "1st element" << std::endl;

    DCs_spread1[0].define(badpx, dmdpx, 1, 0);
    amrex::Print() << "2nd element" << std::endl;
    DCs_spread1[1].define(badpy, dmdpy, 1, 0);
    amrex::Print() << "3rd element" << std::endl;

#elif (AMREX_SPACEDIM == 3)
    amrex::Print() << "1st element" << std::endl;
    DCs_spread1[0].define(badpx, dmdpx, 1, 1);
    amrex::Print() << "2nd element" << std::endl;
    DCs_spread1[1].define(badpy, dmdpy, 1, 1);
    amrex::Print() << "3rd element" << std::endl;

    std::cout<< " Box array " << badpz << std::endl;
    std::cout<< " Distribution Map " << dmdpz << std::endl;
    DCs_spread1[2].define(badpz, dmdpz, 1, 1);
    amrex::Print() << "After defining multifabs" << std::endl;

#endif


    DCs_spread1[0].setVal(0.);
    DCs_spread1[1].setVal(0.);
    DCs_spread1[2].setVal(0.);
    amrex::Print() << "Copying gradient into array of multifabs" << std::endl;

    DCs_spread1[0].copy(*Dc_x[ibpc_lev],0,0,1,0,1);
    DCs_spread1[1].copy(*Dc_y[ibpc_lev],0,0,1,0,1);
    DCs_spread1[2].copy(*Dc_z[ibpc_lev],0,0,1,0,1);


     for (int d=0; d<AMREX_SPACEDIM; d++) {
        umacNew[d].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(umacNew[d], d, geom, d);
        MultiFABPhysBCMacVel(umacNew[d], d, geom, d);

        MultiFab::Copy(uMom[d], umacNew[d], 0, 0, 1, 1);

        // let rho = 1
        uMom[d].mult(1.0, 1);

        // momentum boundary conditions
        uMom[d].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(uMom[d], d, geom, d);
        MultiFABPhysBCMacVel(uMom[d], d, geom, d);
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
    // Corrector step: advect immersed boundary markers

    std::array<MultiFab, AMREX_SPACEDIM> umacNew_buffer;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        umacNew_buffer[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
        MultiFab::Copy(umacNew_buffer[d], umacNew[d], 0, 0, 1, umacNew[d].nGrow());
        umacNew_buffer[d].FillBoundary(geom.periodicity());
    }

    for (const auto & pindex : part_indices) {
        auto & vel = marker_vel.at(pindex);

        ib_pc.InterpolateMarkers(ibpc_lev, pindex, vel, umacNew_buffer);
        auto & dcs = marker_DCs0.at(pindex);
        std::cout << "dcs "<<std::endl;

        ib_pc.InterpolateMarkers(ibpc_lev, pindex, dcs, DC_s0);

    }

    for (const auto & pindex : part_indices) {
        const auto & vel   = marker_vel.at(pindex);
        const auto & pos   = marker_pos.at(pindex);
              auto & f_0   = ib_forces.at(pindex);
              auto & pos_1 = marker_pos_1.at(pindex);
              auto & del_1 = marker_delta_1.at(pindex);
              auto & force = marker_force_1.at(pindex);

        for (int i=0; i<vel.size(); ++i) {
            del_1[i] = dt*vel[i];
            pos_1[i] = pos[i] + del_1[i];
            force[i] = f_0[i] - spring_coefficient*del_1[i];
            f_0[i]   = force[i];

            if (i == 10)
                Print() << "corrector force[" << i << "] = " << force[i] << std::endl;
        }
    }


    //___________________________________________________________________________
    // Add immersed-boundary forces to predictor's RHS

    for (int d=0; d<AMREX_SPACEDIM; ++d)
        force_1[d].setVal(0);

    for (const auto & pindex : part_indices) {
        const auto & force = marker_force_1.at(pindex);

        ib_pc.SpreadMarkers(ibpc_lev, pindex, force, force_1);

        auto & dcs = marker_DCs.at(pindex);
        ib_pc.SpreadMarkers(ibpc_lev, pindex, dcs, DCs_spread1);

    }

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        force_1[d].FillBoundary(geom.periodicity());
        DCs_spread1[d].FillBoundary(geom.periodicity());
        VisMF::Write(force_1[d], "force_1_" + std::to_string(d));
    }


    //_______________________________________________________________________
    // Set up the RHS for the corrector

    for (int d=0; d<AMREX_SPACEDIM; d++) {
        // explicit part
        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);
        gmres_rhs_u[d].mult(dtinv, 1);
        int cng=DCs_spread1[d].nGrow();
        DCs_spread1[d].mult(scaling_factor, cng);

        MultiFab::Add(force_1[d], force_0[d], 0, 0, 1, 1);
        force_1[d].mult(0.5,1);
        MultiFab::Add(DCs_spread1[d], DCs_spread0[d], 0, 0, 1, 1);
        DCs_spread1[d].mult(0.5,1);


        MultiFab::Add(gmres_rhs_u[d], mfluxdiv_correct[d], 0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], advFluxdivPred[d],   0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], force_1[d],          0, 0, 1, 1);
        MultiFab::Add(gmres_rhs_u[d], DCs_spread1[d],       0, 0, 1, 1);

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


    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        // Output velocity solution
        MultiFab::Copy(umac[d],      umacNew[d],   0, 0, 1, 1);

        // Output immersed-boundary forces
        MultiFab::Copy(force_ibm[d], force_1[d],   0, 0, 1, 1);
        MultiFab::Copy(DCs_spread[d], DCs_spread1[d],   0, 0, 1, 1);

        // Output pressure solution
        MultiFab::Copy(pres,         p_1,          0, 0, 1, 1);
    }
}
