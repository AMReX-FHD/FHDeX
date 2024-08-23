
#include "common_functions.H"
#include "chemistry_functions.H"
#include "reactDiff_functions.H"

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

#include "chrono"

using namespace std::chrono;
using namespace amrex;
using namespace common;
using namespace chemistry;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{

    BL_PROFILE_VAR("main_driver()",main_driver);

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::string inputs_file = argv;

    // Initialize variables in namespaces
    InitializeCommonNamespace();
    InitializeChemistryNamespace();
    InitializeReactDiffNamespace();

    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_spec_lo[i] == -1 && bc_spec_hi[i] == -1) {
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

    const Real* dx = geom.CellSize();

    /////////////////////////////////////////
    // Initialize seeds for random number generator
    /////////////////////////////////////////
    if (restart < 0) {

        int mySeed;

        if (seed > 0) {
            // initializes the seed for C++ random number calls with a specified root seed
            mySeed = seed;
        } else if (seed == 0) {
            // initializes the root seed for C++ random number calls based on the clock
            auto now = time_point_cast<nanoseconds>(system_clock::now());
            int mySeed = now.time_since_epoch().count();
            // broadcast the same root seed to all processors
            ParallelDescriptor::Bcast(&mySeed,1,ParallelDescriptor::IOProcessorNumber());
        } else {
            Abort("Must supply non-negative seed");
        }

        // MPI ranks > 0 get a seed inremented by the rank
        InitRandom(mySeed+ParallelDescriptor::MyProc(),
                   ParallelDescriptor::NProcs(),
                   mySeed+ParallelDescriptor::MyProc());

    }
        
    BoxArray ba;
    DistributionMapping dmap;
    
    int step_start;
    amrex::Real time;

    if (restart < 0) {

        step_start = 1;
        time = 0.;
        
        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);

        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        // note we are converting "Vector<int> max_grid_size" to an IntVect
        ba.maxSize(IntVect(max_grid_size));

        dmap.define(ba);

    } else {

        // checkpoint restart
        
    }

    MultiFab n_old(ba,dmap,nspecies,1);
    MultiFab n_new(ba,dmap,nspecies,1);
    
    if (model_file_init) {
        Abort("model_file_init not supported yet");
    } else {
        // Initialize n
        // Init();
        n_old.setVal(0.);
    }

    if (std::abs(initial_variance_mass) > 0.) {
        Abort("initial_variance_mass not supported yet");
        // add_init_n_fluctuations()
    }

    Real dt;
    if (fixed_dt > 0.) {
        dt = fixed_dt;
        Print() << "Setting dt using fixed_dt = " << dt << std::endl;
    } else {
        Real D_Fick_max = 0.;
        for (int i=0; i<nspecies; ++i ) {
            D_Fick_max = std::max(D_Fick_max,D_Fick[i]);
        }
        Real dx_min = dx[0];
        for (int i=1; i<AMREX_SPACEDIM; ++i) {
            dx_min = std::min(dx_min,dx[i]);
        }
        dt = cfl * dx_min / (2. * AMREX_SPACEDIM * D_Fick_max);
        Print() << "Setting dt using explicit diffusion cfl condition = " << dt << std::endl;
    }

    if (inhomogeneous_bc_fix == 1 && temporal_integrator > 0) {
        Abort("comput_n_steady not supported for inhomogeneous_bc_fix == 1 && temporal_integrator > 0 yet");
        // compute_n_steady()
    }

    if (temporal_integrator < 0) { // unsplit schemes
        // Donev: The code will work for a single cell also but may not be the most efficient, so issue warning:
        if (n_cells[0] == 1 && n_cells[1] == 1) {
            Print() << "WARNING in advance_reaction_diffusion: use splitting based schemes (temporal_integrator>=0) for single cell" << std::endl;
        }
        if (nreaction < 1) {
            Print() << "WARNING in advance_reaction_diffusion: use splitting based schemes (temporal_integrator>=0) for diffusion only" << std::endl;
        }
    }

    if (stats_int > 0) {
        Abort("Structure factor not implemented yet");
    }
    
    int istep = (restart < 0) ? 0 : restart;
    WritePlotFile(istep,time,geom,n_old);

    ///////////////////////////////////////////

    // time step loop
    for(int step=step_start;step<=max_step;++step) {
        
        AdvanceTimestep(n_old,n_new,dt,time,geom);

        time += dt;
        MultiFab::Copy(n_new,n_old,0,0,nspecies,1);
        
        if (stats_int > 0 && step%stats_int == 0 && step > n_steps_skip) {
            Abort("fix structure factor snapshot");
        }
        
        if (plot_int > 0 && step%plot_int == 0) {

            WritePlotFile(step,time,geom,n_new);

            if (stats_int > 0 && step > n_steps_skip) {
                Abort("fix structure factor plotfile write");
            }
        }

        if (chk_int > 0 && step%chk_int == 0) {
            Abort("fix checkpoint write");
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
