
#include "common_functions.H"
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

    InitializeCommonNamespace();

    /////////////////////////////////////////
    // Initialize random number seed on all processors/GPUs
    /////////////////////////////////////////

    int mySeed;
    
    if (seed > 0) {
        // use seed from inputs file
        mySeed = seed;

    } else if (seed == 0) {
        // use a clock-based seed
        auto now = time_point_cast<nanoseconds>(system_clock::now());
        mySeed = now.time_since_epoch().count();
        // broadcast the same root seed to all processors
        ParallelDescriptor::Bcast(&mySeed,1,ParallelDescriptor::IOProcessorNumber());
    } else {
        Abort("Must supply non-negative seed");
    }

    // initialize the seed for C++ random number calls
    InitRandom(mySeed+ParallelDescriptor::MyProc(),
               ParallelDescriptor::NProcs(),
               mySeed+ParallelDescriptor::MyProc());

    /////////////////////////////////////////
    
    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,1);  // set to 1 (periodic) by default
        
    // This defines the physical box, [-1,1] in each direction.
    RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                     {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});
        
    IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
    IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
    Box domain(dom_lo, dom_hi);

    Geometry geom(domain,&real_box,CoordSys::cartesian,is_periodic.data());

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // make BoxArray and Geometry
    BoxArray ba;

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap;    

    // Initialize the boxarray "ba" from the single box "bx"
    ba.define(domain);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    // note we are converting "Vector<int> max_grid_size" to an IntVect
    ba.maxSize(IntVect(max_grid_size));

    // define DistributionMapping
    dmap.define(ba);

    MultiFab height(ba, dmap, 1, 1);


    // copy objects from namelist to temporary variables
    // equilibrium height "h0" is stored in rho0
    // surface tension "gamma" is stored in h_bar
    Real h0 = rho0;
    Real gamma = h_bar;
    
    // Physical time constant for dimensional time
    Real t0 = 3.0*visc_coef*h0/gamma;
    
    Real time = 0.;
    Real dt = 0.4 * (t0/std::pow(h0,4)) * std::pow(dx[0],4) / 16.;
    
    // Time stepping loop
    for(int istep=1; istep<=max_step; ++istep) {









        

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
