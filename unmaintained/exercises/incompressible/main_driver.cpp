#include "multispec_test_functions.H"
#include "multispec_test_functions_F.H"

//#include "analysis_functions_F.H"
#include "StochMomFlux.H"
//#include "StructFact.H"


#include "common_functions.H"

#include "gmres_functions.H"

#include "hydro_functions.H"

#include "multispec_functions.H"

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{

    BL_PROFILE_VAR("main_driver()",main_driver);

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::string inputs_file = argv;

    // read in parameters from inputs file into F90 modules
    // we use "+1" because of amrex_string_c_to_f expects a null char termination
    read_common_namelist(inputs_file.c_str(),inputs_file.size()+1);
    read_multispec_namelist(inputs_file.c_str(),inputs_file.size()+1);

    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();
    InitializeMultispecNamespace();
    InitializeGmresNamespace();

    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_vel_lo[i] == -1 && bc_vel_hi[i] == -1) {
            is_periodic[i] = 1;
        }
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
        IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);

        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        // note we are converting "Vector<int> max_grid_size" to an IntVect
        ba.maxSize(IntVect(max_grid_size));

       // This defines the physical box, [-1,1] in each direction.
        RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                         {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    Real dt = fixed_dt;
    Real dtinv = 1.0/dt;
    const Real* dx = geom.CellSize();

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);

    /////////////////////////////////////////
    //Initialise rngs
    /////////////////////////////////////////
    const int n_rngs = 1;

    int fhdSeed = 1;
    int particleSeed = 2;
    int selectorSeed = 3;
    int thetaSeed = 4;
    int phiSeed = 5;
    int generalSeed = 6;

    //Initialise rngs
    rng_initialize(&fhdSeed,&particleSeed,&selectorSeed,&thetaSeed,&phiSeed,&generalSeed);
    /////////////////////////////////////////

    ///////////////////////////////////////////
    // rho, alpha, beta, gamma:
    ///////////////////////////////////////////

    MultiFab rhotot(ba, dmap, 1, 1);
    rhotot.setVal(1.);

    MultiFab rhoxav(ba, dmap, nspecies, 1);
    rhoxav.setVal(1.);

    MultiFab rho(ba, dmap, nspecies, 1);
    rho.setVal(1.);

    // alpha_fc arrays
    Real theta_alpha = 1.;
    std::array< MultiFab, AMREX_SPACEDIM > alpha_fc;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        alpha_fc[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
        alpha_fc[d].setVal(dtinv);
    }

    // beta cell centred
    MultiFab beta(ba, dmap, 1, 1);
    beta.setVal(visc_coef);

    // beta on nodes in 2d
    // beta on edges in 3d
    std::array< MultiFab, NUM_EDGE > beta_ed;
#if (AMREX_SPACEDIM == 2)
    beta_ed[0].define(convert(ba,nodal_flag), dmap, 1, 1);
    beta_ed[0].setVal(visc_coef);
#elif (AMREX_SPACEDIM == 3)
    beta_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
    beta_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
    beta_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
    beta_ed[0].setVal(visc_coef);
    beta_ed[1].setVal(visc_coef);
    beta_ed[2].setVal(visc_coef);
#endif

    // cell-centered gamma
    MultiFab gamma(ba, dmap, 1, 1);
    gamma.setVal(0.);

    ///////////////////////////////////////////

    ///////////////////////////////////////////
    // Define & initalize eta & temperature multifabs
    ///////////////////////////////////////////
    // eta & temperature
    const Real eta_const = visc_coef;
    const Real temp_const = T_init[0];      // [units: K]

    // eta & temperature cell centered
    MultiFab  eta_cc;
    MultiFab temp_cc;
    // eta & temperature nodal
    std::array< MultiFab, NUM_EDGE >   eta_ed;
    std::array< MultiFab, NUM_EDGE >  temp_ed;
    // eta cell-centered
    eta_cc.define(ba, dmap, 1, 1);
    // temperature cell-centered
    temp_cc.define(ba, dmap, 1, 1);
#if (AMREX_SPACEDIM == 2)
    // eta nodal
    eta_ed[0].define(convert(ba,nodal_flag), dmap, 1, 0);
    // temperature nodal
    temp_ed[0].define(convert(ba,nodal_flag), dmap, 1, 0);
#elif (AMREX_SPACEDIM == 3)
    // eta nodal
    eta_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 0);
    eta_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
    eta_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 0);
    // temperature nodal
    temp_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 0);
    temp_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
    temp_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 0);
#endif

    // Initalize eta & temperature multifabs
    // eta cell-centered
    eta_cc.setVal(eta_const);
    // temperature cell-centered
    temp_cc.setVal(temp_const);
#if (AMREX_SPACEDIM == 2)
    // eta nodal
    eta_ed[0].setVal(eta_const);
    // temperature nodal
    temp_ed[0].setVal(temp_const);
#elif (AMREX_SPACEDIM == 3)
    // eta nodal
    eta_ed[0].setVal(eta_const);
    eta_ed[1].setVal(eta_const);
    eta_ed[2].setVal(eta_const);
    // temperature nodal
    temp_ed[0].setVal(temp_const);
    temp_ed[1].setVal(temp_const);
    temp_ed[2].setVal(temp_const);
#endif
    ///////////////////////////////////////////

    ///////////////////////////////////////////
    // Mass and momentum flux divergence MultiFabs:
    ///////////////////////////////////////////

    // mflux divergence, staggered in x,y,z

    // Define mfluxdiv predictor multifabs
    std::array< MultiFab, AMREX_SPACEDIM >  mfluxdiv_predict;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        mfluxdiv_predict[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
        mfluxdiv_predict[d].setVal(0.0);
    }

    // Define mfluxdiv corrector multifabs
    std::array< MultiFab, AMREX_SPACEDIM >  mfluxdiv_correct;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        mfluxdiv_correct[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
        mfluxdiv_correct[d].setVal(0.0);
    }

    ///////////////////////////////////////////

    ///////////////////////////////////////////
    // Stochastic flux divergence class
    ///////////////////////////////////////////

    Vector< amrex::Real > weights;
    // weights = {std::sqrt(0.5), std::sqrt(0.5)};
    weights = {1.0};

    // Declare object of StochMomFlux class
    StochMomFlux sMflux (ba,dmap,geom,n_rngs);

    ///////////////////////////////////////////

    // tracer
    MultiFab tracer(ba,dmap,1,1);

    // pressure for GMRES solve
    MultiFab pres(ba,dmap,1,1);
    pres.setVal(0.);  // initial guess

    // staggered velocities
    std::array< MultiFab, AMREX_SPACEDIM > umac;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umac[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
        umac[d].setVal(0.0);
    }

    std::array< MultiFab, AMREX_SPACEDIM > umacNew;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umacNew[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
        umacNew[d].setVal(0.0);
    }

    ///////////////////////////////////////////
    // structure factor:
    ///////////////////////////////////////////

    /*
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
    s_pairA[1] = 1;
    s_pairB[1] = 1;
#if (AMREX_SPACEDIM == 3)
    s_pairA[2] = 2;
    s_pairB[2] = 2;
#endif

    StructFact structFact(ba,dmap,var_names);
    // StructFact structFact(ba,dmap,var_names,s_pairA,s_pairB);
    */

    ///////////////////////////////////////////

    const RealBox& realDomain = geom.ProbDomain();
    int dm;

    for ( MFIter mfi(beta); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        init_rho_and_umac(BL_TO_FORTRAN_BOX(bx),
                          BL_TO_FORTRAN_FAB(rho[mfi]),
                          BL_TO_FORTRAN_ANYD(umac[0][mfi]),
                          BL_TO_FORTRAN_ANYD(umac[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                          BL_TO_FORTRAN_ANYD(umac[2][mfi]),
#endif
                          dx, geom.ProbLo(), geom.ProbHi(),
                          ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));

        // AMREX_D_TERM(dm=0; init_vel(BL_TO_FORTRAN_BOX(bx),
        //                             BL_TO_FORTRAN_ANYD(umac[0][mfi]), geom.CellSize(),
        //                             geom.ProbLo(), geom.ProbHi() ,&dm,
        //                             ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));,
        //              dm=1; init_vel(BL_TO_FORTRAN_BOX(bx),
        //                             BL_TO_FORTRAN_ANYD(umac[1][mfi]), geom.CellSize(),
        //                             geom.ProbLo(), geom.ProbHi() ,&dm,
        //                             ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));,
        //              dm=2; init_vel(BL_TO_FORTRAN_BOX(bx),
        //                             BL_TO_FORTRAN_ANYD(umac[2][mfi]), geom.CellSize(),
        //                             geom.ProbLo(), geom.ProbHi() ,&dm,
        //                             ZFILL(realDomain.lo()), ZFILL(realDomain.hi())););

        // initialize tracer
        init_s_vel(BL_TO_FORTRAN_BOX(bx),
                   BL_TO_FORTRAN_ANYD(tracer[mfi]),
                   dx, ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));
    }

    // Add initial equilibrium fluctuations
    sMflux.addMomFluctuations(umac, rhotot, temp_cc, initial_variance_mom);

    // Project umac onto divergence free field
    MultiFab macphi(ba,dmap,1,1);
    MultiFab macrhs(ba,dmap,1,1);
    macrhs.setVal(0.0);
    MacProj_hydro(umac,rhotot,geom,true);

    // initial guess for new solution
    for (int d; d<AMREX_SPACEDIM; d++) {
        MultiFab::Copy(umacNew[d], umac[d], 0, 0, 1, 0);
    }

    int step = 0;
    Real time = 0.;

    // write out initial state
    if (plot_int > 0)
    {
        WritePlotFile(step,time,geom,umac,rho,rhoxav,tracer,pres);
    }

    //Time stepping loop
    for(step=1;step<=max_step;++step) {

        Real step_strt_time = ParallelDescriptor::second();

        if(variance_coef_mom != 0.0) {

            // Fill stochastic terms
            sMflux.fillMomStochastic();

            // compute stochastic force terms
            sMflux.StochMomFluxDiv(mfluxdiv_predict,0,eta_cc,eta_ed,temp_cc,temp_ed,
                                   weights,dt);
            sMflux.StochMomFluxDiv(mfluxdiv_correct,0,eta_cc,eta_ed,temp_cc,temp_ed,
                                   weights,dt);

        }

        if(step == 500)
        {
            PrintMF (rhotot);
        }

        // Advance umac
        advance(umac,umacNew,pres,tracer,rho,rhotot,
                mfluxdiv_predict,mfluxdiv_correct,
                alpha_fc,beta,gamma,beta_ed,geom,dt);

        //////////////////////////////////////////////////

        ///////////////////////////////////////////
        // Update structure factor
        ///////////////////////////////////////////
        /*
        if (step > n_steps_skip && struct_fact_int > 0 && (step-n_steps_skip-1)%struct_fact_int == 0) {
            for(int d=0; d<AMREX_SPACEDIM; d++) {
                ShiftFaceToCC(umac[d], 0, struct_in_cc, d, 1);
            }
            structFact.FortStructure(struct_in_cc);
        }
        */

        Real step_stop_time = ParallelDescriptor::second() - step_strt_time;
        ParallelDescriptor::ReduceRealMax(step_stop_time);

        amrex::Print() << "Advanced step " << step << " in " << step_stop_time << " seconds\n";

        time = time + dt;

        MaxSpeed(umac);

        XMeanFab(rho, rhoxav, 0);

        if (plot_int > 0 && step%plot_int == 0) {
            // write out umac & pres to a plotfile
            WritePlotFile(step,time,geom,umac,rho, rhoxav,tracer,pres);
        }
    }

    /*
    if (struct_fact_int > 0) {
        Real dVol = dx[0]*dx[1];
        int tot_n_cells = n_cells[0]*n_cells[1];
        if (AMREX_SPACEDIM == 2) {
            dVol *= cell_depth;
        } else if (AMREX_SPACEDIM == 3) {
            dVol *= dx[2];
            tot_n_cells = n_cells[2]*tot_n_cells;
        }

        // let rhotot = 1
        Real SFscale = dVol/(k_B*temp_const);
        // Print() << "Hack: structure factor scaling = " << SFscale << std::endl;

        structFact.Finalize(SFscale);
        structFact.WritePlotFile(step,time,geom,"plt_SF");
    }
    */

    // Call the timer again and compute the maximum difference between the start time
    // and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;

    // exit(0);

}