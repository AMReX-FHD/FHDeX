#include "main_driver.H"
#include "main_driver_F.H"

#include "hydro_functions.H"

#include "StochMomFlux.H"
#include "StructFact.H"


#include "common_functions.H"

#include "gmres_functions.H"



#include <AMReX_ParmParse.H>
#include "IBParticleContainer.H"
#include "IBCore.H"

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

//chemistry
#include <iostream>
#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AmrCoreAdv.H>

using namespace amrex;


//! Defines staggered MultiFab arrays (BoxArrays set according to the
//! nodal_flag_[x,y,z]). Each MultiFab has 1 component, and 1 ghost cell
inline void defineFC(std::array< MultiFab, AMREX_SPACEDIM > & mf_in,
                     const BoxArray & ba, const DistributionMapping & dm,
                     int nghost) {

    for (int i=0; i<AMREX_SPACEDIM; i++)
        mf_in[i].define(convert(ba, nodal_flag_dir[i]), dm, 1, nghost);
}

inline void defineEdge(std::array< MultiFab, AMREX_SPACEDIM > & mf_in,
                     const BoxArray & ba, const DistributionMapping & dm,
                     int nghost) {

    for (int i=0; i<AMREX_SPACEDIM; i++)
        mf_in[i].define(convert(ba, nodal_flag_edge[i]), dm, 1, nghost);
}


//! Sets the value for each component of staggered MultiFab
inline void setVal(std::array< MultiFab, AMREX_SPACEDIM > & mf_in,
                   Real set_val) {

    for (int i=0; i<AMREX_SPACEDIM; i++)
        mf_in[i].setVal(set_val);
}



// argv contains the name of the inputs file entered at the command line
void main_driver(const char * argv) {

    BL_PROFILE_VAR("main_driver()",main_driver);


    /****************************************************************************
     *                                                                          *
     * Initialize simulation                                                    *
     *                                                                          *
     ***************************************************************************/

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();


    //___________________________________________________________________________
    // Load parameters from inputs file, and initialize global parameters
    std::string inputs_file = argv;

    // read in parameters from inputs file into F90 modules NOTE: we use "+1"
    // because of amrex_string_c_to_f expects a null char termination
    read_common_namelist(inputs_file.c_str(), inputs_file.size()+1);

    // copy contents of F90 modules to C++ namespaces NOTE: any changes to
    // global settings in fortran/c++ after this point need to be synchronized
    InitializeCommonNamespace();
    InitializeGmresNamespace();


    //___________________________________________________________________________
    // Set boundary conditions

    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i)
        if (bc_vel_lo[i] <= -1 && bc_vel_hi[i] <= -1)
            is_periodic[i] = 1;

    //___________________________________________________________________________
    // Make BoxArray, DistributionMapping, and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(             0,              0,              0));
        IntVect dom_hi(AMREX_D_DECL(n_cells[0] - 1, n_cells[1] - 1, n_cells[2] - 1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);

        // Break up boxarray "ba" into chunks no larger than "max_grid_size"
        // along a direction note we are converting "Vector<int> max_grid_size"
        // to an IntVect
        ba.maxSize(IntVect(max_grid_size));

        // This defines the physical box, [-1, 1] in each direction
        RealBox real_box({AMREX_D_DECL(prob_lo[0], prob_lo[1], prob_lo[2])},
                         {AMREX_D_DECL(prob_hi[0], prob_hi[1], prob_hi[2])});

        // This defines a Geometry object
        geom.define(domain, & real_box, CoordSys::cartesian, is_periodic.data());
    }

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);



    //___________________________________________________________________________
    // Cell size, and time step
    Real dt         = fixed_dt;
    Real dtinv      = 1.0 / dt;
    const Real * dx = geom.CellSize();


    //___________________________________________________________________________
    // Initialize random number generators
    const int n_rngs = 1;

    // this seems really random :P
    int fhdSeed      = 1;
    int particleSeed = 2;
    int selectorSeed = 3;
    int thetaSeed    = 4;
    int phiSeed      = 5;
    int generalSeed  = 6;

    // each CPU gets a different random seed
    const int proc = ParallelDescriptor::MyProc();
    fhdSeed      += proc;
    particleSeed += proc;
    selectorSeed += proc;
    thetaSeed    += proc;
    phiSeed      += proc;
    generalSeed  += proc;

    // initialize rngs
    rng_initialize( & fhdSeed, & particleSeed, & selectorSeed,
                    & thetaSeed, & phiSeed, & generalSeed);



    /****************************************************************************
     *                                                                          *
     * Initialize physical parameters                                           *
     *                                                                          *
     ***************************************************************************/

    //___________________________________________________________________________
    // Set rho, alpha, beta, gamma:

    // rho is cell-centered
    MultiFab rho(ba, dmap, 1, 1);
    rho.setVal(1.);

    // alpha_fc is face-centered
    Real theta_alpha = 1.;
    std::array< MultiFab, AMREX_SPACEDIM > alpha_fc;
    defineFC(alpha_fc, ba, dmap, 1);
    setVal(alpha_fc, dtinv);

    // beta is cell-centered
    MultiFab beta(ba, dmap, 1, 1);
    beta.setVal(visc_coef);

    // beta is on nodes in 2D, and is on edges in 3D
    std::array< MultiFab, NUM_EDGE > beta_ed;
#if (AMREX_SPACEDIM == 2)
    beta_ed[0].define(convert(ba, nodal_flag), dmap, 1, 1);
    beta_ed[0].setVal(visc_coef);
#elif (AMREX_SPACEDIM == 3)
    defineEdge(beta_ed, ba, dmap, 1);
    setVal(beta_ed, visc_coef);
#endif

    // cell-centered gamma
    MultiFab gamma(ba, dmap, 1, 1);
    gamma.setVal(0.);


    //___________________________________________________________________________
    // Define & initialize eta & temperature MultiFabs

    // eta & temperature
    const Real eta_const  = visc_coef;
    const Real temp_const = T_init[0];      // [units: K]


    // NOTE: eta and temperature live on both cell-centers and edges

    // eta & temperature cell centered
    MultiFab  eta_cc(ba, dmap, 1, 1);
    MultiFab temp_cc(ba, dmap, 1, 1);
    // eta & temperature nodal
    std::array< MultiFab, NUM_EDGE >   eta_ed;
    std::array< MultiFab, NUM_EDGE >  temp_ed;

    // eta_ed and temp_ed are on nodes in 2D, and on edges in 3D
#if (AMREX_SPACEDIM == 2)
    eta_ed[0].define(convert(ba,nodal_flag), dmap, 1, 0);
    temp_ed[0].define(convert(ba,nodal_flag), dmap, 1, 0);

    eta_ed[0].setVal(eta_const);
    temp_ed[0].setVal(temp_const);
#elif (AMREX_SPACEDIM == 3)
    defineEdge(eta_ed, ba, dmap, 1);
    defineEdge(temp_ed, ba, dmap, 1);

    setVal(eta_ed, eta_const);
    setVal(temp_ed, temp_const);
#endif

    // eta_cc and temp_cc are always cell-centered
    eta_cc.setVal(eta_const);
    temp_cc.setVal(temp_const);


    //___________________________________________________________________________
    // Define random fluxes mflux (momentum-flux) divergence, staggered in x,y,z

    // mfluxdiv predictor multifabs
    std::array< MultiFab, AMREX_SPACEDIM >  mfluxdiv_predict;
    defineFC(mfluxdiv_predict, ba, dmap, 1);
    setVal(mfluxdiv_predict, 0.);

    // mfluxdiv corrector multifabs
    std::array< MultiFab, AMREX_SPACEDIM >  mfluxdiv_correct;
    defineFC(mfluxdiv_correct, ba, dmap, 1);
    setVal(mfluxdiv_correct, 0.);

    Vector< amrex::Real > weights;
    // weights = {std::sqrt(0.5), std::sqrt(0.5)};
    weights = {1.0};

    // tracer
    MultiFab tracer(ba, dmap, 1,1);
    tracer.setVal(0.);


    //___________________________________________________________________________
    // Define velocities and pressure

    // pressure for GMRES solve
    MultiFab pres(ba, dmap, 1, 1);
    pres.setVal(0.);  // initial guess

    // staggered velocities
    std::array< MultiFab, AMREX_SPACEDIM > umac;
    defineFC(umac, ba, dmap, 1);

    std::array< MultiFab, AMREX_SPACEDIM > umacNew;
    defineFC(umacNew, ba, dmap, 1);


    //___________________________________________________________________________
    // Define structure factor:

    Vector< std::string > var_names;
    var_names.resize(AMREX_SPACEDIM);
    int cnt = 0;
    std::string x;
    for (int d=0; d<var_names.size(); d++) {
        x = "vel";
        x += (120+d);
        var_names[cnt++] = x;
    }

    MultiFab struct_in_cc;
    struct_in_cc.define(ba, dmap, AMREX_SPACEDIM, 0);

    amrex::Vector< int > s_pairA(AMREX_SPACEDIM);
    amrex::Vector< int > s_pairB(AMREX_SPACEDIM);

    // Select which variable pairs to include in structure factor:
    s_pairA[0] = 0;
    s_pairB[0] = 0;
    //
    s_pairA[1] = 1;
    s_pairB[1] = 1;
    //
#if (AMREX_SPACEDIM == 3)
    s_pairA[2] = 2;
    s_pairB[2] = 2;
#endif

    //StructFact structFact(ba, dmap, var_names);
    // StructFact structFact(ba, dmap, var_names, s_pairA, s_pairB);



    /****************************************************************************
     *                                                                          *
     * Set Initial Conditions                                                   *
     *                                                                          *
     ***************************************************************************/

    //___________________________________________________________________________
    // Initialize velocities (fluid and tracers)

    const RealBox& realDomain = geom.ProbDomain();
    int dm;

    for ( MFIter mfi(beta); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        // initialize velocity
        for (int d=0; d<AMREX_SPACEDIM; ++d)
            init_vel(BL_TO_FORTRAN_BOX(bx),
                     BL_TO_FORTRAN_ANYD(umac[d][mfi]), geom.CellSize(),
                     geom.ProbLo(), geom.ProbHi(), & d,
                     ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));

    	// initialize tracer
        init_s_vel(BL_TO_FORTRAN_BOX(bx),
                   BL_TO_FORTRAN_ANYD(tracer[mfi]),
                   dx, ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));

    }


    //___________________________________________________________________________
    // Ensure that ICs satisfy BCs

    pres.FillBoundary(geom.periodicity());
    MultiFabPhysBC(pres, geom, 0, 1, PRES_BC_COMP);

    for (int i=0; i<AMREX_SPACEDIM; i++) {
        umac[i].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(umac[i], i, geom, i);
        MultiFabPhysBCMacVel(umac[i], i, geom, i);
    }


    //___________________________________________________________________________
    // Add random momentum fluctuations

    // Declare object of StochMomFlux class
    StochMomFlux sMflux (ba, dmap, geom, n_rngs);

    // Add initial equilibrium fluctuations
    sMflux.addMomFluctuations(umac, rho, temp_cc, initial_variance_mom);

    // Project umac onto divergence free field
    MultiFab macphi(ba,dmap, 1, 1);
    MultiFab macrhs(ba,dmap, 1, 1);
    macrhs.setVal(0.);
    MacProj_hydro(umac, rho, geom, true); // from MacProj_hydro.cpp

    // initial guess for new solution
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        MultiFab::Copy(umacNew[d], umac[d], 0, 0, 1, 1);

    int step = 0;
    Real time = 0.;

    // Initialize Chemical fields

    //   AmrCoreAdv amr_core_adv;


    /****************************************************************************
     *                                                                          *
     * Build container for immersed particles                                   *
     *                                                                          *
     ***************************************************************************/

    // TODO: 8 is a magic number for now: it needs to be larger than the
    // particle radius (in cells)
    IBParticleContainer ib_pc(geom, dmap, ba, 10);

    // Initializing colloid position
    Vector<RealVect> ib_pos(1);
    ib_pos[0] = RealVect{AMREX_D_DECL(0.5, 0.5, 0.5)};
    Vector<Real>     ib_r(1);
    ib_r[0]   = 0.10;
    Vector<Real>     ib_rho(1);
    ib_rho[0] = 100.00;

    ib_pc.InitList(0, ib_pos, ib_r, ib_rho);


    // DEBUG: Test interface
    //ib_pc.PrintParticleData(0);

    //__________________________________________________________________________
    // Build IB core

    // NOTE: IBCore (based on AmrCore) needs a parm parse database to contain
    // certain fields => add them here
    {
        ParmParse pp("amr");
        pp.add("max_level", 0);
        pp.addarr("n_cell", n_cells);
    }



    IBCore ib_core;
    ib_core.set_IBParticleContainer(& ib_pc);


    std::array<MultiFab, AMREX_SPACEDIM> force_ibm;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        force_ibm[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        force_ibm[d].setVal(0.);
    }

    //__________________________________________________________________________
    // Build AmrCore and initialize chemical multifabs

    AmrCoreAdv amr_core_adv;

    amr_core_adv.InitData( ba, dmap);


    // Need to have only one level for now
    int lev =0;


    /****************************************************************************
     *                                                                          *
     * Advance Time Steps                                                       *
     *                                                                          *
     ***************************************************************************/

    //___________________________________________________________________________
    // Write out initial state
    if (plot_int > 0) {
        WritePlotFile(step, time, geom, umac, tracer, pres, force_ibm, ib_pc,
                      amr_core_adv, lev);
    }

    Print() << "Write Plot File Success "<<std::endl;

    //___________________________________________________________________________
    // FFT test
    // if (struct_fact_int > 0) {
    //     std::array <MultiFab, AMREX_SPACEDIM> mf_cc;
    //     mf_cc[0].define(ba, dmap, 1, 0);
    //     mf_cc[1].define(ba, dmap, 1, 0);
    //     mf_cc[2].define(ba, dmap, 1, 0);
    //     for ( MFIter mfi(beta); mfi.isValid(); ++mfi ) {
    //         const Box& bx = mfi.validbox();
    //         init_s_vel(BL_TO_FORTRAN_BOX(bx),
    //                    BL_TO_FORTRAN_ANYD(mf_cc[0][mfi]),
    //                    dx, ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));
    //     }
    // }


    ib_pc.FillMarkerPositions(0, 64);


    //___________________________________________________________________________
    // Temporary place to store marker forces _outside_ advanve
    // => this will be moved to the IBParticleContainer class ASAP

    Vector<IBP_info> ibp_info = ib_pc.IBParticleInfo(0, true);

    Vector<ParticleIndex> part_indices(ibp_info.size());

    IBMarkerMap marker_force_0;

    for (int i=0; i<ibp_info.size(); ++i) {
        part_indices[i] = ibp_info[i].asPairIndex();

        // Pre-allocate particle arrays
        const Vector<RealVect> marker_positions = ib_pc.MarkerPositions(0, part_indices[i]);
        // ... initialized to (0..0) by default constructor
        marker_force_0[part_indices[i]].resize(marker_positions.size());
    }

    //---------------------------------------------------------------------------




    for(step = 1; step <= max_step; ++step) {
        Real step_strt_time = ParallelDescriptor::second();

        // if(variance_coef_mom != 0.0) {

        //     //___________________________________________________________________
        //     // Fill stochastic terms

        //     sMflux.fillMomStochastic();

        //     // Compute stochastic force terms (and apply to mfluxdiv_*)
        //     // NOTE: StochMomFlux::StochMomFluxDiv fills ghost cells
        //     sMflux.StochMomFluxDiv(mfluxdiv_predict, 0, eta_cc, eta_ed, temp_cc, temp_ed, weights, dt);
        //     sMflux.StochMomFluxDiv(mfluxdiv_correct, 0, eta_cc, eta_ed, temp_cc, temp_ed, weights, dt);
        // }


        //_______________________________________________________________________
        // Fill Immersed-Boundary interface data

        int lev_ib = 0;
        Real t0_ib = 0;

        ib_core.MakeNewLevelFromScratch(lev_ib, t0_ib, ba, dmap);


        //_______________________________________________________________________
        // Advance umac and chemistry
 
        advance(amr_core_adv,
                umac, umacNew, pres, tracer,
                force_ibm, marker_force_0,
                mfluxdiv_predict, mfluxdiv_correct,
                alpha_fc, beta, gamma, beta_ed,
                ib_pc, ib_core, geom, dt, time);


        // Empty force data
        std::map<ParticleIndex, std::array<Real, AMREX_SPACEDIM>> f_trans;
        f_trans.clear();

        // Print() << "Force data BEFORE Interpolation:" << std::endl;
        // for (const auto & f : f_trans) {
        //     std::cout << f.first.first <<", " << f.first.second << ":" << std::endl;
        //     for (int d=0; d<AMREX_SPACEDIM; ++d)
        //         std::cout << f.second[d] << std::endl;
        // }

        ib_pc.InterpolateParticleForces(0, force_ibm, ib_core, f_trans);

        // Print() << "Force data AFTER Interpolation:" << std::endl;
        // for (const auto & f : f_trans) {
        //     std::cout << f.first.first <<", " << f.first.second << ":" << std::endl;
        //     for (int d=0; d<AMREX_SPACEDIM; ++d)
        //         std::cout << f.second[d] << std::endl;
        // }

        // ib_pc.MoveIBParticles(0, dt, f_trans);


        //_______________________________________________________________________
        // Update structure factor

        if (step > n_steps_skip
                && struct_fact_int > 0
                && (step-n_steps_skip-1)%struct_fact_int == 0
            ) {

            for(int d=0; d<AMREX_SPACEDIM; d++)
                ShiftFaceToCC(umac[d], 0, struct_in_cc, d, 1);
          //  Have to comment this out for now
          //  structFact.FortStructure(struct_in_cc);
        }

        Real step_stop_time = ParallelDescriptor::second() - step_strt_time;
        ParallelDescriptor::ReduceRealMax(step_stop_time);

        amrex::Print() << "Advanced step " << step << " in " << step_stop_time << " seconds\n";

        time = time + dt;

        if (plot_int > 0 && step%plot_int == 0) {
            // write out umac, pres, f_ibm, and particle data to a plotfile
            WritePlotFile(step, time, geom, umac, tracer, pres, force_ibm, ib_pc, amr_core_adv,lev);
        }
    }

    ///////////////////////////////////////////
    if (struct_fact_int > 0) {
        Real dVol = dx[0]*dx[1];
        int tot_n_cells = n_cells[0]*n_cells[1];
        if (AMREX_SPACEDIM == 2) {
            dVol *= cell_depth;
        } else if (AMREX_SPACEDIM == 3) {
            dVol *= dx[2];
            tot_n_cells = n_cells[2]*tot_n_cells;
        }

        // let rho = 1
        Real SFscale = dVol/(k_B*temp_const);
        // Have to comment this out for now
        // Print() << "Hack: structure factor scaling = " << SFscale << std::endl;

        //structFact.Finalize(SFscale);
       // structFact.WritePlotFile(step,time,geom,"plt_SF");
    }

    // Call the timer again and compute the maximum difference between the start
    // time and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;

}
