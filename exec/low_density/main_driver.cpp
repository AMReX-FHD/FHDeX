
#include "rng_functions.H"

#include "common_functions.H"

#include "common_namespace_declarations.H"

#include "compressible_functions.H"

#include "exec_functions.H"

#include "StructFact.H"

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFab.H>


using namespace amrex;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::string inputs_file = argv;

    // read in parameters from inputs file into F90 modules
    // we use "+1" because of amrex_string_c_to_f expects a null char termination
    read_common_namelist(inputs_file.c_str(),inputs_file.size()+1);
    // read_gmres_namelist(inputs_file.c_str(),inputs_file.size()+1);

    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();
    // InitializeGmresNamespace();

  
    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_vel_lo[i] == -1 && bc_vel_hi[i] == -1) {
            is_periodic[i] = 1;
            Print() << "Periodic: " << is_periodic[i] << "\n";
        }
        //is_periodic[i] = 0;
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

    // Print() << "Hack: boxarray = " << ba << "\n";

    Real dt = fixed_dt;
    Real dtinv = 1.0/dt;
    const Real* dx = geom.CellSize();
    const RealBox& realDomain = geom.ProbDomain();

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);

    /////////////////////////////////////////
    //Initialise rngs
    /////////////////////////////////////////
    const int n_rngs = 1;

    const int proc = ParallelDescriptor::MyProc();

    int fhdSeed = 0;
    int particleSeed = 2;
    int selectorSeed = 3;
    int thetaSeed = 4;
    int phiSeed = 5;
    int generalSeed = 0;

    //fhdSeed += 10000*proc;
    particleSeed += 20000*proc;
    selectorSeed += 30000*proc;
    thetaSeed += 40000*proc;
    phiSeed += 50000*proc;
    //generalSeed += 60000*proc;

    //Initialise rngs
    rng_initialize(&fhdSeed,&particleSeed,&selectorSeed,&thetaSeed,&phiSeed,&generalSeed);
    /////////////////////////////////////////


    int eval_dt = 0;

    //conserved quantaties
    MultiFab cu  (ba,dmap,nvars,ngc); 
    MultiFab cup  (ba,dmap,nvars,ngc);  //use some of these for RK3, but prob just do Euler first
    MultiFab cup2  (ba,dmap,nvars,ngc);
    MultiFab cup3  (ba,dmap,nvars,ngc);


    //statistics
    MultiFab cuMeans(ba,dmap,nvars,ngc);
    MultiFab cuVars(ba,dmap,nvars,ngc);

    cuMeans.setVal(0.0);
    cuVars.setVal(0.0);

    //Initialize physical parameters from input vals
    //TODO - set varying density in different regions
    IntVect lo1(AMREX_D_DECL(0,0,0));
    IntVect lo2(AMREX_D_DECL(1,0,0));
    IntVect hi1(AMREX_D_DECL(1,1,1));
    IntVect hi2(AMREX_D_DECL(2,1,1));

    Box b1(lo1,hi1);
    Box b2(lo2,hi2);
    cu.setVal(7,b1,0,1);
    cu.setVal(23,b2,0,1);

    //fluxes - not using for now but keep
    std::array< MultiFab, AMREX_SPACEDIM > flux;
    AMREX_D_TERM(flux[0].define(convert(ba,nodal_flag_x), dmap, nvars, 0);,
                 flux[1].define(convert(ba,nodal_flag_y), dmap, nvars, 0);,
                 flux[2].define(convert(ba,nodal_flag_z), dmap, nvars, 0););

    //stochastic fluxes
    std::array< MultiFab, AMREX_SPACEDIM > stochFlux;
    AMREX_D_TERM(stochFlux[0].define(convert(ba,nodal_flag_x), dmap, nvars, 0);,
                 stochFlux[1].define(convert(ba,nodal_flag_y), dmap, nvars, 0);,
                 stochFlux[2].define(convert(ba,nodal_flag_z), dmap, nvars, 0););

    //zero out stoch fluxes
    AMREX_D_TERM(stochFlux[0].setVal(0.0);,
                 stochFlux[1].setVal(0.0);,
                 stochFlux[2].setVal(0.0););

    Real time = 0;

    int step, statsCount;

    statsCount = 1;

    //Time stepping loop
    for(step=1;step<=max_step;++step)
    {

        eulerStep(cu, flux, stochFlux, geom, dx, dt);
        //RK3step(cu, cup, cup2, cup3, flux, stochFlux, geom, dx, dt);

        if(step == n_steps_skip)
        {
            cuMeans.setVal(0.0);
            cuVars.setVal(0.0);

            statsCount = 1;
        }

	if (step > n_steps_skip) {
	  // evaluateStats(cu, cuMeans, cuVars, prim, primMeans, primVars, spatialCross, eta, etaMean, kappa, kappaMean, delHolder1, delHolder2, delHolder3, delHolder4, delHolder5, delHolder6, statsCount, dx);
	}


        statsCount++;

	amrex::Print() << "Advanced step " << step << "\n";

        if (plot_int > 0 && step > 0 && step%plot_int == 0)
        {

           WritePlotFile(step, time, geom, cu, cuMeans, cuVars);
        }

        time = time + dt;
    }

    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;

    // amrex::Finalize();

}

