#include "multispec_test_functions.H"

#include "StochMomFlux.H"

#include "StructFact.H"

#include "common_functions.H"
#include "gmres_functions.H"
#include "multispec_functions.H"


#include "hydro_functions.H"
#include "rng_functions.H"

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

#include "chrono"

using namespace std::chrono;


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


    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();
    InitializeMultispecNamespace();
    InitializeGmresNamespace();

    if (visc_type > 0 && mixture_type != 0) {
        Abort("positive visc_type and non-zero mixture_type not compatible");
    }

    if (algorithm_type == 6) {
        RhototBCInit();
    }

    // for reservoirs, make sure the Dirichlet conditions for concentration sum to 1
    //
    //
    //

    /////////////////////////////////////////
    //Initialise rngs
    /////////////////////////////////////////
    if (restart < 0) {

        if (seed > 0) {
            // initializes the seed for C++ random number calls
            InitRandom(seed+ParallelDescriptor::MyProc(),
                       ParallelDescriptor::NProcs(),
                       seed+ParallelDescriptor::MyProc());
        } else if (seed == 0) {
            // initializes the seed for C++ random number calls based on the clock
            auto now = time_point_cast<nanoseconds>(system_clock::now());
            int randSeed = now.time_since_epoch().count();
            // broadcast the same root seed to all processors
            ParallelDescriptor::Bcast(&randSeed,1,ParallelDescriptor::IOProcessorNumber());
            InitRandom(randSeed+ParallelDescriptor::MyProc(),
                       ParallelDescriptor::NProcs(),
                       randSeed+ParallelDescriptor::MyProc());
        } else {
            Abort("Must supply non-negative seed");
        }

    }
    /////////////////////////////////////////

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

    int init_step;
    Real time;
    Real dt;

    // make BoxArray and Geometry
    BoxArray ba;

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap;

    MultiFab rho_old;
    MultiFab rhotot_old;
    MultiFab pi;
    std::array< MultiFab, AMREX_SPACEDIM > umac;
    MultiFab Epot;
    std::array< MultiFab, AMREX_SPACEDIM > grad_Epot_old;

    if (restart > 0) {
        ReadCheckPoint(init_step,time,dt,
                       rho_old,rhotot_old,pi,umac,Epot,grad_Epot_old,
                       ba,dmap);
    }
    else {

        init_step = 1;
        time = start_time;
        if (fixed_dt <= 0.) {
            Abort("main_driver.cpp: only fixed_dt > 0 supported");
        }

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);

        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        // note we are converting "Vector<int> max_grid_size" to an IntVect
        ba.maxSize(IntVect(max_grid_size));

        // define DistributionMapping
        dmap.define(ba);

        // staggered velocities
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            umac[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
            umac[d].setVal(0.);

            if (use_charged_fluid) {
                grad_Epot_old[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
            }
        }

        //      added by JBB to allow reset of dt on restart
        if (fixed_dt >0) dt=fixed_dt;

        rho_old.define   (ba, dmap, nspecies, ng_s);
        rhotot_old.define(ba, dmap, 1       , ng_s);
        pi.define        (ba, dmap, 1       , 1);
        if (use_charged_fluid) {
            Epot.define(ba, dmap, 1, 1);
        }
    }

    // moved this to here so can change dt from value in checkpoint
    dt = fixed_dt;

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
        InitRhoUmac(umac,rho_old,geom);

        // initialize pi, including ghost cells
        pi.setVal(0.);
    }

    // compute rhotot from rho in VALID REGION
    ComputeRhotot(rho_old,rhotot_old);

    // fill rho and rhotot ghost cells
    FillRhoRhototGhost(rho_old,rhotot_old,geom);

    // pressure ghost cells
    pi.FillBoundary(geom.periodicity());
    MultiFabPhysBC(pi,geom,0,1,PRES_BC_COMP);

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
        stoch_mass_flux[d].define(convert(ba,nodal_flag_dir[d]), dmap, nspecies, 0);
    }

    std::array< MultiFab, AMREX_SPACEDIM > grad_Epot_new;
    MultiFab charge_old;
    MultiFab charge_new;
    MultiFab permittivity;
    if (use_charged_fluid) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            grad_Epot_new[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
        }
        charge_old.define(ba, dmap, 1, 1);
        charge_new.define(ba, dmap, 1, 1);
        permittivity.define(ba, dmap, 1, 1);

        // set these to zero
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            grad_Epot_old[d].setVal(0.);
            grad_Epot_new[d].setVal(0.);
        }
        charge_old.setVal(0.);
        charge_new.setVal(0.);
        Epot.setVal(0.);

        if (electroneutral) {
            Abort("main_driver.cpp: Electroneutral not implemented");
        }

        // compute total charge
        DotWithZ(rho_old,charge_old);

        // multiply by total volume (all 3 dimensions, even for 2D problems)

        // NOTE: we are using rho = 1 here,
        // so the below is a close approximation to debye lenth
        Real sum_temp = 0.;
        for (int d=0; d<nspecies; ++d) {
            sum_temp += c_init_1[d] * molmass[d] * charge_per_mass[d] * charge_per_mass[d];
        }
        Real debye_len =sqrt(dielectric_const*k_B*T_init[0]/(rho0*sum_temp));
        Print() << "Debye length lambda_D is approx: " << debye_len << std::endl;

        Real total_charge = (AMREX_SPACEDIM == 2) ? charge_old.sum() * dx[0] * dx[1] * cell_depth
                                                  : charge_old.sum() * dx[0] * dx[1] * dx[2];
        Print() << "Initial total charge " << total_charge << std::endl;

        // compute permittivity
        if (dielectric_type == 0) {
            permittivity.setVal(dielectric_const);
        }
        else {
            Abort("main_driver.cpp: dielectric_type != 0 not supported");
        }
    }


    // allocate and build MultiFabs that will contain random numbers
    // by declaring StochMassFlux and StochMomFlux objects
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
    StochMassFlux sMassFlux(ba,dmap,geom,n_rngs_mass);
    StochMomFlux  sMomFlux (ba,dmap,geom,n_rngs_mom);

    // save random state for writing checkpoint
    //
    //
    //

    //=====================================================================
    // Initialize values
    //=====================================================================

    /*
    if (use_charged_fluid) {

    }
    */

    // initial Temp and Temp_ed
    Temp.setVal(T_init[0]); // replace with more general initialization routine
    if (AMREX_SPACEDIM == 2) {
        AverageCCToNode(Temp,Temp_ed[0],0,1,TEMP_BC_COMP,geom);
    }
    else {
        AverageCCToEdge(Temp,Temp_ed,0,1,TEMP_BC_COMP,geom);
    }

    /*
    if (barodiffusion_type > 0) {

    }
    */

    // initialize eta and kappa
    ComputeEta(rho_old, rhotot_old, eta);
    kappa.setVal(0.);
    // replace with more general initialization routine
    //
    //
    if (AMREX_SPACEDIM == 2) {
        AverageCCToNode(eta,eta_ed[0],0,1,SPEC_BC_COMP,geom);
    }
    else {
        AverageCCToEdge(eta,eta_ed,0,1,SPEC_BC_COMP,geom);
    }

    // now that we have eta, we can initialize the inhomogeneous velocity bc's
    // set inhomogeneous velocity bc's to values supplied in inhomogeneous_bc_val
    //
    //
    //

    // velocity boundary conditions
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        // set normal velocity of physical domain boundaries
        MultiFabPhysBCDomainVel(umac[i],geom,i);
        // set transverse velocity behind physical boundaries
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umac[i],geom,i,is_inhomogeneous);
        // fill periodic and interior ghost cells
        umac[i].FillBoundary(geom.periodicity());
        // protect against roundoff issues and sync up
        // faces with the same physical location
        umac[i].OverrideSync(geom.periodicity());
    }

    if (restart < 0) {

        if ((algorithm_type != 2) && (initial_variance_mom != 0.)) {
            // Add initial momentum fluctuations
            addMomFluctuations(umac, rhotot_old, Temp, initial_variance_mom, geom);

            for (int i=0; i<AMREX_SPACEDIM; ++i) {
                // set normal velocity of physical domain boundaries
                MultiFabPhysBCDomainVel(umac[i],geom,i);
                // set transverse velocity behind physical boundaries
                int is_inhomogeneous = 1;
                MultiFabPhysBCMacVel(umac[i],geom,i,is_inhomogeneous);
                // fill periodic and interior ghost cells
                umac[i].FillBoundary(geom.periodicity());
            }
        }
    }

    ///////////////////////////////////////////
    // Initialize structure factor object for analysis
    ///////////////////////////////////////////

    // variables are density, velocity and concentrations
    int structVars = AMREX_SPACEDIM+nspecies+1;

    Vector< std::string > var_names;
    var_names.resize(structVars);

    int cnt = 0;
    std::string x;

    // density
    x = "rho";
    var_names[cnt++] = x;

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

    MultiFab structFactMF(ba, dmap, structVars, 0);

    // need to use dVol for scaling
    Real dVol = dx[0]*dx[1];
    if (AMREX_SPACEDIM == 2) {
        dVol *= cell_depth;
    } else if (AMREX_SPACEDIM == 3) {
        dVol *= dx[2];
    }

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

    /*
      this routine is only called for all inertial simulations (both restart and non-restart)
      it does the following:
      1. fill mass random numbers
      2. computes mass fluxes and flux divergences
      if restarting, the subroutine ends; otherwise
      3. perform an initial projection

      overdamped schemes need to do 1. and 2. within the advance_timestep routine
      in principle, performing an initial projection for overdamped will change
      the reference state for the GMRES solver
      For overdamped the first ever solve cannot have a good reference state
      so in general there is the danger it will be less accurate than subsequent solves
      but I do not see how one can avoid that
      From this perspective it may be useful to keep initial_projection even in overdamped
      because different gmres tolerances may be needed in the first step than in the rest
    */
    if (algorithm_type != 2) {
        InitialProjection(umac,rho_old,rhotot_old,diff_mass_fluxdiv,stoch_mass_fluxdiv,
                          stoch_mass_flux,sMassFlux,Temp,eta,eta_ed,dt,time,geom,
                          charge_old,grad_Epot_old,Epot,permittivity);
    }

    if (restart < 0) {

        // We do the analysis first so we include the initial condition in the files if n_steps_skip=0
        if (n_steps_skip == 0 && struct_fact_int > 0) {

            // add this snapshot to the average in the structure factor

            // copy density into structFactMF
            MultiFab::Copy(structFactMF,rhotot_old,0,0,1,0);

            // copy velocities into structFactMF
            for(int d=0; d<AMREX_SPACEDIM; d++) {
                ShiftFaceToCC(umac[d], 0, structFactMF, d+1, 1);
            }
            // copy concentrations into structFactMF
            MultiFab::Copy(structFactMF,rho_old,0,AMREX_SPACEDIM+1,nspecies,0);
            for(int d=0; d<nspecies; d++) {
                MultiFab::Divide(structFactMF,rhotot_old,0,AMREX_SPACEDIM+d+1,1,0);
            }
            structFact.FortStructure(structFactMF);
        }

        // write initial plotfile and structure factor
        if (plot_int > 0) {
            WritePlotFile(0,0.,geom,umac,rhotot_old,rho_old,pi,charge_old,Epot);
            if (n_steps_skip == 0 && struct_fact_int > 0) {
                structFact.WritePlotFile(0,0.,"plt_SF");
            }
        }

        if (chk_int > 0) {
            // write initial checkpoint
            //
            //
            //
        }

        /*
        if (stats_int > 0) {
            // write initial vertical and horizontal averages (hstat and vstat files)
            //
            //
            //
        }
        */

    }


    // Time stepping loop
    for(int istep=init_step; istep<=max_step; ++istep) {

        Real step_strt_time = ParallelDescriptor::second();

        if (algorithm_type == 0) {
            // inertial
            AdvanceTimestepInertial(umac,rho_old,rho_new,rhotot_old,rhotot_new,
                                    pi,eta,eta_ed,kappa,Temp,Temp_ed,
                                    diff_mass_fluxdiv,stoch_mass_fluxdiv,stoch_mass_flux,
                                    grad_Epot_old,grad_Epot_new,
                                    charge_old,charge_new,Epot,permittivity,
                                    sMassFlux,sMomFlux,
                                    dt,time,istep,geom);
        }
        else if (algorithm_type == 6) {
            // boussinesq
            AdvanceTimestepBousq(umac,rho_old,rho_new,rhotot_old,rhotot_new,
                                 pi,eta,eta_ed,kappa,Temp,Temp_ed,
                                 diff_mass_fluxdiv,stoch_mass_fluxdiv,stoch_mass_flux,
                                 grad_Epot_old,grad_Epot_new,
                                 charge_old,charge_new,Epot,permittivity,
                                 sMassFlux,sMomFlux,
                                 dt,time,istep,geom);
        }
        else {
            Print() << "algorithm_type " << algorithm_type << std::endl;
            Abort("algorithm_type not supported");
        }

        //////////////////////////////////////////////////
        if (istep > n_steps_skip && struct_fact_int > 0 && (istep-n_steps_skip)%struct_fact_int == 0) {

            // add this snapshot to the average in the structure factor

            // copy density into structFactMF
            MultiFab::Copy(structFactMF,rhotot_new,0,0,1,0);

            // copy velocities into structFactMF
            for(int d=0; d<AMREX_SPACEDIM; d++) {
                ShiftFaceToCC(umac[d], 0, structFactMF, d+1, 1);
            }
            // copy concentrations into structFactMF
            MultiFab::Copy(structFactMF,rho_new,0,AMREX_SPACEDIM+1,nspecies,0);
            for(int d=0; d<nspecies; d++) {
                MultiFab::Divide(structFactMF,rhotot_new,0,AMREX_SPACEDIM+d+1,1,0);
            }
            structFact.FortStructure(structFactMF);
        }

        Real step_stop_time = ParallelDescriptor::second() - step_strt_time;
        ParallelDescriptor::ReduceRealMax(step_stop_time);

        amrex::Print() << "Advanced step " << istep << " in " << step_stop_time << " seconds\n";

        time = time + dt;

        // write plotfile at specific intervals
        if (plot_int > 0 && istep%plot_int == 0) {
            WritePlotFile(istep,time,geom,umac,rhotot_new,rho_new,pi,charge_new,Epot);
            if (istep > n_steps_skip && struct_fact_int > 0) {
                structFact.WritePlotFile(istep,time,"plt_SF");
            }
        }

        // write checkpoint at specific intervals
        if (chk_int > 0 && istep%chk_int == 0) {
            WriteCheckPoint(istep,time,dt,rho_new,rhotot_new,pi,umac,Epot,grad_Epot_new);
        }

        // set old state to new state
        MultiFab::Copy(rho_old   ,rho_new   ,0,0,nspecies,ng_s);
        MultiFab::Copy(rhotot_old,rhotot_new,0,0,       1,ng_s);

        if (use_charged_fluid) {
            MultiFab::Copy(charge_old, charge_new, 0, 0, 1, charge_old.nGrow());
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                MultiFab::Copy(grad_Epot_old[d], grad_Epot_new[d], 0, 0, 1, grad_Epot_old[d].nGrow());
            }
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