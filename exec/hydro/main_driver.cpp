
#include "hydro_test_functions.H"
#include "hydro_test_functions_F.H"

#include "hydro_functions.H"

#include "StochMomFlux.H"

#include "StructFact.H"
#include "TurbForcing.H"

#include "rng_functions_F.H"

#include "common_functions.H"

#include "gmres_functions.H"

#include "common_namespace_declarations.H"

#include "gmres_namespace_declarations.H"

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

    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();
    InitializeGmresNamespace();

    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_vel_lo[i] == -1 && bc_vel_hi[i] == -1) {
            is_periodic[i] = 1;
        }
    }

    // This defines the physical box, [-1,1] in each direction.
    RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                     {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});
    
    IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
    IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
    Box domain(dom_lo, dom_hi);

    Geometry geom(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    
    // BoxArray
    BoxArray ba;

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap;

    Real dt = fixed_dt;
    Real dtinv = 1.0/dt;
    const Real* dx = geom.CellSize();

    /////////////////////////////////////////
    //Initialise rngs
    /////////////////////////////////////////
    const int n_rngs = 1;

    if (restart <= 0) {
        int fhdSeed = 1;
        int particleSeed = 2;
        int selectorSeed = 3;
        int thetaSeed = 4;
        int phiSeed = 5;
        int generalSeed = 6;

        //Initialise rngs
        rng_initialize(&fhdSeed,&particleSeed,&selectorSeed,&thetaSeed,&phiSeed,&generalSeed);

        // initializes the seed for C++ random number calls
        InitRandom(seed+ParallelDescriptor::MyProc());
    }
    /////////////////////////////////////////
    
    int step_start;
    amrex::Real time;

    // object for turbulent forcing
    TurbForcing turbforce;

    // tracer
    MultiFab tracer;
    
    // staggered velocities
    std::array< MultiFab, AMREX_SPACEDIM > umac;
    
    if (restart > 0) {
        ReadCheckPoint(step_start,time,umac,tracer,turbforce,ba,dmap);
    }
    else {

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);

        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        // note we are converting "Vector<int> max_grid_size" to an IntVect
        ba.maxSize(IntVect(max_grid_size));

        dmap.define(ba);
    
	if (turbForcing == 1) {
	  turbforce.define(ba,dmap,turb_a,turb_b);
	}

        const RealBox& realDomain = geom.ProbDomain();
        int dm;
        
        AMREX_D_TERM(umac[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                     umac[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                     umac[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););
    
        tracer.define(ba,dmap,1,1);
        tracer.setVal(0.);
        
        for ( MFIter mfi(tracer); mfi.isValid(); ++mfi ) {
            const Box& bx = mfi.validbox();

            AMREX_D_TERM(dm=0; init_vel(BL_TO_FORTRAN_BOX(bx),
                                        BL_TO_FORTRAN_ANYD(umac[0][mfi]), geom.CellSize(),
                                        geom.ProbLo(), geom.ProbHi() ,&dm,
                                        ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));,
                         dm=1; init_vel(BL_TO_FORTRAN_BOX(bx),
                                        BL_TO_FORTRAN_ANYD(umac[1][mfi]), geom.CellSize(),
                                        geom.ProbLo(), geom.ProbHi() ,&dm,
                                        ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));,
                         dm=2; init_vel(BL_TO_FORTRAN_BOX(bx),
                                        BL_TO_FORTRAN_ANYD(umac[2][mfi]), geom.CellSize(),
                                        geom.ProbLo(), geom.ProbHi() ,&dm,
                                        ZFILL(realDomain.lo()), ZFILL(realDomain.hi())););

    	// initialize tracer
        init_s_vel(BL_TO_FORTRAN_BOX(bx),
    		   BL_TO_FORTRAN_ANYD(tracer[mfi]),
    		   dx, ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));

        }

        // temporary for addMomfluctuations and MacProj_hydro
        MultiFab rho(ba, dmap, 1, 1);
        rho.setVal(1.);
        MultiFab temp_cc(ba, dmap, 1, 1);
        temp_cc.setVal(T_init[0]);
    
        // Add initial equilibrium fluctuations
        if(initial_variance_mom != 0.0) {
            addMomFluctuations(umac, rho, temp_cc, initial_variance_mom,geom);
        }

        // Project umac onto divergence free field
        {
            // macrhs only used once at beginning of simulaton
            // put this in braces so it goes out of scope immediately
            MultiFab macrhs(ba,dmap,1,1);
            macrhs.setVal(0.0);
            MacProj_hydro(umac,rho,geom,true);
        }

        step_start = 1;
        time = 0.;

    }

    if (turbForcing == 1) {
      turbforce.Initialize(geom);
    }

    // pressure for GMRES solve
    MultiFab pres(ba,dmap,1,1);
    pres.setVal(0.);  // initial guess
    
    ///////////////////////////////////////////
    // alpha, beta, gamma:
    ///////////////////////////////////////////

    // alpha_fc arrays
    std::array< MultiFab, AMREX_SPACEDIM > alpha_fc;
    AMREX_D_TERM(alpha_fc[0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 alpha_fc[1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 alpha_fc[2].define(convert(ba,nodal_flag_z), dmap, 1, 0););
    AMREX_D_TERM(alpha_fc[0].setVal(dtinv);,
                 alpha_fc[1].setVal(dtinv);,
                 alpha_fc[2].setVal(dtinv););


    // -1 = inertial backward Euler
    //  0 = inertial
    //  1 = overdamped
    Real factor = (algorithm_type == 0) ? 0.5 : 1.0;
    
    // beta cell centred
    MultiFab beta(ba, dmap, 1, 1);
    beta.setVal(factor*visc_coef); // multiply by factor here

    // beta on nodes in 2d
    // beta on edges in 3d
    std::array< MultiFab, NUM_EDGE > beta_ed;
#if (AMREX_SPACEDIM == 2)
    beta_ed[0].define(convert(ba,nodal_flag), dmap, 1, 0);
#elif (AMREX_SPACEDIM == 3)
    beta_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 0);
    beta_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
    beta_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 0);
#endif
    for (int d=0; d<NUM_EDGE; ++d) {
        beta_ed[d].setVal(factor*visc_coef); // multiply by factor here
    }    

    // cell-centered gamma
    MultiFab gamma(ba, dmap, 1, 1);
    gamma.setVal(0.);

    ///////////////////////////////////////////

    ///////////////////////////////////////////
    // Define & initalize eta & temperature multifabs
    ///////////////////////////////////////////
    // eta & temperature
    const Real eta_const = visc_coef;

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
    temp_cc.setVal(T_init[0]);
#if (AMREX_SPACEDIM == 2)
    // eta nodal
    eta_ed[0].setVal(eta_const);
    // temperature nodal
    temp_ed[0].setVal(T_init[0]);
#elif (AMREX_SPACEDIM == 3)
    // eta nodal
    eta_ed[0].setVal(eta_const);
    eta_ed[1].setVal(eta_const);
    eta_ed[2].setVal(eta_const);
    // temperature nodal
    temp_ed[0].setVal(T_init[0]);
    temp_ed[1].setVal(T_init[0]);
    temp_ed[2].setVal(T_init[0]);
#endif
    ///////////////////////////////////////////

    ///////////////////////////////////////////
    // random fluxes:
    ///////////////////////////////////////////

    // mflux divergence, staggered in x,y,z

    // Define mfluxdiv predictor multifabs
    std::array< MultiFab, AMREX_SPACEDIM >  mfluxdiv_stoch;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
      mfluxdiv_stoch[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
      mfluxdiv_stoch[d].setVal(0.0);
    }

    Vector< amrex::Real > weights;
    // weights = {std::sqrt(0.5), std::sqrt(0.5)};
    weights = {1.0};

    // Declare object of StochMomFlux class
    StochMomFlux sMflux (ba,dmap,geom,n_rngs);

    ///////////////////////////////////////////
    // Initialize structure factor object for analysis
    ///////////////////////////////////////////
    
    // variables are velocities
    int structVars = AMREX_SPACEDIM;
    
    Vector< std::string > var_names;
    var_names.resize(structVars);
    
    int cnt = 0;
    std::string x;

    // velx, vely, velz
    for (int d=0; d<AMREX_SPACEDIM; d++) {
      x = "vel";
      x += (120+d);
      var_names[cnt++] = x;
    }

    MultiFab structFactMF(ba, dmap, structVars, 0);

    // need to use dVol for scaling
    Real dVol = (AMREX_SPACEDIM==2) ? dx[0]*dx[1]*cell_depth : dx[0]*dx[1]*dx[2];
    
    Real dProb = (AMREX_SPACEDIM==2) ? n_cells[0]*n_cells[1] : n_cells[0]*n_cells[1]*n_cells[2];
    dProb = 1./dProb;
    
    Vector<Real> var_scaling(structVars*(structVars+1)/2);
    for (int d=0; d<var_scaling.size(); ++d) {
        var_scaling[d] = 1./dVol;
    }

#if 1
    // option to compute all pairs
    StructFact structFact(ba,dmap,var_names,var_scaling);
#else
    // option to compute only specified pairs
    int nPairs = 2;
    amrex::Vector< int > s_pairA(nPairs);
    amrex::Vector< int > s_pairB(nPairs);

    // Select which variable pairs to include in structure factor:
    s_pairA[0] = 0;
    s_pairB[0] = 0;
    s_pairA[1] = 1;
    s_pairB[1] = 1;
    
    StructFact structFact(ba,dmap,var_names,var_scaling,s_pairA,s_pairB);
#endif
    
    ///////////////////////////////////////////
    // Structure factor object to help compute tubulent energy spectra
    ///////////////////////////////////////////
    
    // option to compute only specified pairs
    amrex::Vector< int > s_pairA(AMREX_SPACEDIM);
    amrex::Vector< int > s_pairB(AMREX_SPACEDIM);

    var_scaling.resize(AMREX_SPACEDIM);
    for (int d=0; d<var_scaling.size(); ++d) {
        var_scaling[d] = 1./dVol;
    }
    
    // Select which variable pairs to include in structure factor:
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        s_pairA[d] = d;
        s_pairB[d] = d;
    }    
    StructFact turbStructFact(ba,dmap,var_names,var_scaling,s_pairA,s_pairB);
    
    ///////////////////////////////////////////

    if (restart < 0) {

        // We do the analysis first so we include the initial condition in the files if n_steps_skip=0
        if (n_steps_skip == 0 && struct_fact_int > 0) {

            // add this snapshot to the average in the structure factor

            // copy velocities into structFactMF
            for(int d=0; d<AMREX_SPACEDIM; d++) {
                ShiftFaceToCC(umac[d], 0, structFactMF, d, 1);
            }
            structFact.FortStructure(structFactMF,geom);
        }

        // write out initial state
        // write out umac, tracer, pres, and divergence to a plotfile
        if (plot_int > 0) {
            WritePlotFile(0,time,geom,umac,tracer,pres);
            if (n_steps_skip == 0 && struct_fact_int > 0) {
                structFact.WritePlotFile(0,0.,geom,"plt_SF");
            }
        }
    }
    
    // FIXME need to fill physical boundary condition ghost cells for tracer

    std::array< MultiFab, AMREX_SPACEDIM > umacTemp;
    AMREX_D_TERM(umacTemp[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 umacTemp[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 umacTemp[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););   

    // storage for grad(U) for energy dissipation calculation
    MultiFab gradU(ba,dmap,AMREX_SPACEDIM,0);
    MultiFab ccTemp(ba,dmap,1,0);
    
    ///////////////////////////////////////////

    //Time stepping loop
    for(int step=step_start;step<=max_step;++step) {

        Real step_strt_time = ParallelDescriptor::second();

	if(variance_coef_mom != 0.0) {

	  // Fill stochastic terms
	  sMflux.fillMomStochastic();

	  // compute stochastic force terms
	  sMflux.StochMomFluxDiv(mfluxdiv_stoch,0,eta_cc,eta_ed,temp_cc,temp_ed,weights,dt);
	}

	// Advance umac
        advance(umac,umacTemp,pres,tracer,mfluxdiv_stoch,
                alpha_fc,beta,gamma,beta_ed,geom,dt,turbforce);

	//////////////////////////////////////////////////

        // add a snapshot to the structure factor
	if (step > n_steps_skip && struct_fact_int > 0 && (step-n_steps_skip)%struct_fact_int == 0) {

            // add this snapshot to the average in the structure factor

            // copy velocities into structFactMF
            for(int d=0; d<AMREX_SPACEDIM; d++) {
                ShiftFaceToCC(umac[d], 0, structFactMF, d, 1);
            }
            structFact.FortStructure(structFactMF,geom);
        }
                
        Real step_stop_time = ParallelDescriptor::second() - step_strt_time;
        ParallelDescriptor::ReduceRealMax(step_stop_time);

        amrex::Print() << "Advanced step " << step << " in " << step_stop_time << " seconds\n";

        time = time + dt;

        if (plot_int > 0 && step%plot_int == 0) {
            // write out umac, tracer, pres, and divergence to a plotfile
            WritePlotFile(step,time,geom,umac,tracer,pres);

            // write out structure factor to plotfile
            if (step > n_steps_skip && struct_fact_int > 0) {
                structFact.WritePlotFile(step,time,geom,"plt_SF");
            }

            // snapshot of instantaneous energy spectra
            if (turbForcing == 1) {

                // copy velocities into structFactMF
                for(int d=0; d<AMREX_SPACEDIM; d++) {
                    ShiftFaceToCC(umac[d], 0, structFactMF, d, 1);
                }
                // reset and compute structure factor
                turbStructFact.FortStructure(structFactMF,geom,1);

                // writing the plotfiles does the shifting and copying into cov_mag
                turbStructFact.WritePlotFile(step,time,geom,"plt_Turb");

                // integrate cov_mag over shells in k and write to file
                turbStructFact.IntegratekShells(step,geom);
            }
        }

        if (chk_int > 0 && step%chk_int == 0) {
            // write out umac and tracer to a checkpoint file
            WriteCheckPoint(step,time,umac,tracer,turbforce);
        }

        // compute kinetic energy integral( (1/2) * rho * U dot U dV)
        Vector<Real> udotu(3);
        Vector<Real> skew(3);
        Vector<Real> kurt(3);
        StagInnerProd(geom,umac,0,umac,0,umacTemp,udotu);
        Print() << "Kinetic energy "
		<< time << " "
                << 0.5*dVol*( udotu[0] + udotu[1] + udotu[2] )
                << std::endl;

        // compute energy dissipation integral(eta grad(U) dot grad(U))
        ComputeCentredGradFC(umac,gradU,geom);
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            CCInnerProd(gradU,d,gradU,d,ccTemp,udotu[d]);
        }
        Print() << "Energy dissipation "
		<< time << " "
                << visc_coef*dProb*( udotu[0] + udotu[1] + udotu[2] )
                << std::endl;
        //      << visc_coef*dVol*( udotu[0] + udotu[1] + udotu[2] )
	//
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            CCMoments(gradU,d,ccTemp,3,skew[d]);
        }

        Print() << "Skewness "
		<< time << " "
                << dProb*skew[0]/(pow(dProb*udotu[0],1.5)) << " "
                << dProb*skew[1]/(pow(dProb*udotu[1],1.5)) << " "
#if (AMREX_SPACEDIM == 3)
                << dProb*skew[2]/(pow(dProb*udotu[2],1.5))
#endif
                << std::endl;

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            CCMoments(gradU,d,ccTemp,4,kurt[d]);
        }

        Print() << "Kurtosis "
		<< time << " "
                << dProb*kurt[0]/(pow(dProb*udotu[0],2.)) << " "
                << dProb*kurt[1]/(pow(dProb*udotu[1],2.)) << " "
#if (AMREX_SPACEDIM == 3)
                << dProb*kurt[2]/(pow(dProb*udotu[2],2.))
#endif
                << std::endl;

        // use gradU as a temporary to store averaged velocities
        AverageFaceToCC(umac,gradU,0);
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            Print() << "Sum of umac in direction " << d << "= "
                    << gradU.sum(d) << std::endl;          
        }
        
        // MultiFab memory usage
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        amrex::Long min_fab_megabytes  = amrex::TotalBytesAllocatedInFabsHWM()/1048576;
        amrex::Long max_fab_megabytes  = min_fab_megabytes;

        ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);

        amrex::Print() << "High-water FAB megabyte spread across MPI nodes: ["
                       << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";

        min_fab_megabytes  = amrex::TotalBytesAllocatedInFabs()/1048576;
        max_fab_megabytes  = min_fab_megabytes;

        ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);

        amrex::Print() << "Curent     FAB megabyte spread across MPI nodes: ["
                       << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";
        
    }

    // Call the timer again and compute the maximum difference between the start time
    // and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;

}
