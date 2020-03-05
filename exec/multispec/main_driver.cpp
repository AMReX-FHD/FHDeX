
#include "multispec_test_functions.H"
#include "multispec_test_functions_F.H"

#include "StochMassFlux.H"
#include "StochMomFlux.H"

//#include "analysis_functions_F.H"
#include "StructFact.H"

#include "common_functions.H"
#include "gmres_functions.H"
#include "multispec_functions.H"

#include "common_namespace_declarations.H"
#include "gmres_namespace_declarations.H"
#include "multispec_namespace_declarations.H"

#include "hydro_functions.H"
#include "rng_functions.H"

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>


// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{
  
    BL_PROFILE_VAR("main_driver()",main_driver);

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    //=============================================================
    // Initialization
    //=============================================================
    
    std::string inputs_file = argv;

    // read in parameters from inputs file into F90 modules
    // we use "+1" because of amrex_string_c_to_f expects a null char termination
    read_common_namelist   (inputs_file.c_str(),inputs_file.size()+1);
    read_multispec_namelist(inputs_file.c_str(),inputs_file.size()+1);
    read_gmres_namelist    (inputs_file.c_str(),inputs_file.size()+1);

    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();
    InitializeMultispecNamespace();
    InitializeGmresNamespace();

    // for reservoirs, make sure the Dirichlet conditions for concentration sum to 1
    //
    //
    //

    // one common seed; not split by process yet like the original code
    int fhdSeed = seed;
    
    // these are unused
    int particleSeed = 2;
    int selectorSeed = 3;
    int thetaSeed = 4;
    int phiSeed = 5;
    int generalSeed = 6;

    //Initialise rngs
    rng_initialize(&fhdSeed,&particleSeed,&selectorSeed,&thetaSeed,&phiSeed,&generalSeed);
    
    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_vel_lo[i] == -1 && bc_vel_hi[i] == -1) {
            is_periodic[i] = 1;
        }
    }

    int ng_s; // ghost cells for density MultiFabs
    if (advection_type == 0) {
        ng_s = 2; // centered advection
    }
    else if (advection_type <= 3) {
        ng_s = 3; // bilinear limited, biliniear unlimited, or unlimited quad bds
    }
    else if (advection_type == 4) {
        ng_s = 4; // limited quad bds
    }

    Real dt;

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    
    if (restart >= 0) {
        Abort("restart not implemented yet");
    }
    else {
        
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

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);
    
    MultiFab rho_old   (ba, dmap, nspecies, ng_s);
    MultiFab rhotot_old(ba, dmap, 1       , ng_s);
    MultiFab pi        (ba, dmap, 1       , 1);

    // staggered velocities
    std::array< MultiFab, AMREX_SPACEDIM > umac;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umac[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
        umac[d].setVal(0.);
    }

    // data structures to help with reservoirs
    // 
    //
    //

    // get grid spacing
    const Real* dx = geom.CellSize();

    // build layouts for staggered multigrid solver and macproject within preconditioner
    //
    //
    //

    if (restart < 0) {
    
        // initialize rho and umac in valid region only
        for ( MFIter mfi(rho_old); mfi.isValid(); ++mfi ) {
            const Box& bx = mfi.validbox();

            init_rho_and_umac(BL_TO_FORTRAN_BOX(bx),
                              BL_TO_FORTRAN_FAB(rho_old[mfi]),
                              BL_TO_FORTRAN_ANYD(umac[0][mfi]),
                              BL_TO_FORTRAN_ANYD(umac[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                              BL_TO_FORTRAN_ANYD(umac[2][mfi]),
#endif
                              dx, geom.ProbLo(), geom.ProbHi());
        }

        // initialize pi, including ghost cells
        pi.setVal(0.);
    }

    // compute rhotot from rho in VALID REGION
    //
    //
    //

    // fill rho and rhotot ghost cells
    //
    //
    //

    // pressure ghost cells
    //
    //
    //

    //=======================================================
    // Build multifabs for all the variables
    //=======================================================

    MultiFab rho_new          (ba, dmap, nspecies, ng_s);
    MultiFab rhotot_new       (ba, dmap, 1       , ng_s);
    MultiFab Temp             (ba, dmap, 1       , ng_s);
    MultiFab diff_mass_fluxdiv(ba, dmap, nspecies, 0);
    MultiFab eta              (ba, dmap, 1       , 1);
    MultiFab kappa            (ba, dmap, 1       , 1);
    
    /////////////////////////////////////////

    // eta and Temp on nodes (2d) or edges (3d)
    std::array< MultiFab, NUM_EDGE > eta_ed;
    std::array< MultiFab, NUM_EDGE > Temp_ed;
#if (AMREX_SPACEDIM == 2)
    eta_ed[0].define (convert(ba,nodal_flag), dmap, 1, 0);
    Temp_ed[0].define(convert(ba,nodal_flag), dmap, 1, 0);
#elif (AMREX_SPACEDIM == 3)
    eta_ed[0].define (convert(ba,nodal_flag_xy), dmap, 1, 0);
    eta_ed[1].define (convert(ba,nodal_flag_xz), dmap, 1, 0);
    eta_ed[2].define (convert(ba,nodal_flag_yz), dmap, 1, 0);
    Temp_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 0);
    Temp_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
    Temp_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 0);
#endif

    MultiFab stoch_mass_fluxdiv(ba,dmap,nspecies,0);
    std::array< MultiFab, AMREX_SPACEDIM > stoch_mass_flux;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
      stoch_mass_flux[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
    }

    // random numbers

    int n_rngs_mass;
    int n_rngs_mom;
    
    if (algorithm_type == 2 || algorithm_type == 5) {
        n_rngs_mass = 2;
        n_rngs_mom = 2;  
    }
    else if (algorithm_type == 6) {
        n_rngs_mass = 2;
        n_rngs_mom = 1;  
    }
    else {
        n_rngs_mass = 1;
        n_rngs_mom = 1;        
    }
    
    // Declare object of StochMassFlux and StochMomFlux class es
    StochMassFlux sMassFlux(ba,dmap,geom,n_rngs_mass);
    StochMomFlux  sMomFlux (ba,dmap,geom,n_rngs_mom);

    // save random state for writing checkpoint
    //
    //
    //

    //=====================================================================
    // Initialize values
    //=====================================================================

    // initial Temp and Temp_ed
    //
    //
    //

    // initialize eta and kappa
    //
    //
    //

    // now that we have eta, we can initialize the inhomogeneous velocity bc's
    // set inhomogeneous velocity bc's to values supplied in inhomogeneous_bc_val
    //
    //
    //

    // velocity boundary conditions
    //
    //
    //

    if (restart < 0) {

        if ((algorithm_type != 2) && (initial_variance_mom != 0.)) {
            // Add initial equilibrium fluctuations
            sMomFlux.addMfluctuations(umac, rhotot_old, Temp, initial_variance_mom);

            // velocity boundary conditions
            //
            //
            //
        }

        dt = fixed_dt;
    }

    ///////////////////////////////////////////
    // structure factor:
    ///////////////////////////////////////////

    // variables are velocity and concentrations
    int structVars = AMREX_SPACEDIM+nspecies;
    
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

    // c1, c2, etc.
    for (int d=0; d<nspecies; d++) {
      x = "c";
      x += (49+d);
      var_names[cnt++] = x;
    }
    
    MultiFab SF(ba, dmap, structVars, 0);

    // need to use dVol for scaling
    Real dVol = dx[0]*dx[1];
    if (AMREX_SPACEDIM == 2) {
	dVol *= cell_depth;
    } else if (AMREX_SPACEDIM == 3) {
	dVol *= dx[2];
    }
    
    Vector< Real > var_scaling(structVars*structVars);
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




                                

    
    
    
    // Project umac onto divergence free field
    MacProj(umac,rhotot_old,geom,true);

    int step = 0;
    Real time = 0.;

    // write out initial state
    if (plot_int > 0) {
	WritePlotFile(step,time,geom,umac,rho_old,pi);
    }

    //Time stepping loop
    for(step=1;step<=max_step;++step) {

        Real step_strt_time = ParallelDescriptor::second();

        if (algorithm_type == 0) {
            /*
            advance_timestep_inertial(umac,pi,rho_old,rhotot_old,
                                      alpha_fc,beta,gamma,beta_ed,geom,dt);
            */
        }
        else {
            Print() << "algorithm_type " << algorithm_type << std::endl;
            Abort("algorithm_type not supported");
        }

	//////////////////////////////////////////////////
	
	///////////////////////////////////////////
	// Update structure factor
	///////////////////////////////////////////
	/*
	if (step > n_steps_skip && struct_fact_int > 0 && (step-n_steps_skip-1)%struct_fact_int == 0) {
	  for(int d=0; d<AMREX_SPACEDIM; d++) {
	    ShiftFaceToCC(umac[d], 0, struct_in_cc, d, 1);
	  }
	  structFact.FortStructure(struct_in_cc,geom);
        }
	*/

        Real step_stop_time = ParallelDescriptor::second() - step_strt_time;
        ParallelDescriptor::ReduceRealMax(step_stop_time);
    
        amrex::Print() << "Advanced step " << step << " in " << step_stop_time << " seconds\n";

        time = time + dt;

        if (plot_int > 0 && step%plot_int == 0) {
    	  WritePlotFile(step,time,geom,umac,rho_new,pi);
        }

        // set old state to new state
        // rho
        // rhotot
        
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
      Real SFscale = dVol/(k_B*Temp_const);
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

}
