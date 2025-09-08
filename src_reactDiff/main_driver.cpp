#include "common_functions.H"
#include "chemistry_functions.H"
#include "reactDiff_functions.H"
#include "StructFact.H"

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

    if (nreaction > 0 && use_mole_frac_LMA) {
        if (include_discrete_LMA_correction == 1) {
            Abort("Error: currently use_mole_frac_LMA can be used only with include_discrete_LMA_correction=0");
        }
        if (exclude_solvent_comput_rates != -1) {
            Abort("Error: currently use_mole_frac_LMA can be used only with exclude_solvent_comput_rates=-1");
        }
    }

    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_mass_lo[i] == -1 && bc_mass_hi[i] == -1) {
            is_periodic[i] = 1;
        }
    }

    // This defines the physical box, [-1,1] in each direction.
    RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                     {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

    IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
    IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
    Box domain(dom_lo, dom_hi);

    long cell_count = (AMREX_SPACEDIM==2) ? n_cells[0]*n_cells[1] : n_cells[0]*n_cells[1]*n_cells[2];

    Geometry geom(domain,&real_box,CoordSys::cartesian,is_periodic.data());

    const Real* dx = geom.CellSize();

    BoxArray ba;
    DistributionMapping dmap;

    int step_start;
    amrex::Real time;

    MultiFab n_old;
    MultiFab n_new;

    ///////////////////////////////////////////
    // Initialize structure factor object for analysis
    ///////////////////////////////////////////

    Vector< std::string > var_names;
    var_names.resize(nspecies);

    int cnt = 0;
    std::string x;

    // n0, n1, ...
    for (int d=0; d<nspecies; d++) {
        x = "n";
        x += (49+d);
        var_names[cnt++] = x;
    }

    // need to use dv for scaling
    Real dv = (AMREX_SPACEDIM==2) ? dx[0]*dx[1]*cell_depth : dx[0]*dx[1]*dx[2]*cell_depth;

    // 0 = compute only specified pais listed in s_pairA and s_pairB
    // 1 = compute all possible pairs of variables
    int compute_all_pairs = 1;

    int nPairs = (compute_all_pairs) ? nspecies*(nspecies+1)/2 : 2;

    Vector<Real> var_scaling(nPairs);
    for (int d=0; d<var_scaling.size(); ++d) {
        var_scaling[d] = 1./dv;
    }

    StructFact structFact;

    if (restart < 0) {
        step_start = 1;
        time = 0.;

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);

        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        // note we are converting "Vector<int> max_grid_size" to an IntVect
        ba.maxSize(IntVect(max_grid_size));

        dmap.define(ba);

        n_old.define(ba,dmap,nspecies,1);
        n_new.define(ba,dmap,nspecies,1);

        if (model_file_init) {
            Abort("model_file_init not supported yet");
        } else {
            // Initialize n
            InitN(n_old,geom,time);
        }

        if (std::abs(initial_variance_mass) > 0.) {
            if (integer_populations == 0) {
                Abort("add_init_n_fluctuations not supported yet");
                // add_init_n_fluctuations()
            }
        }

        // structure factor
        if (compute_all_pairs) {
            // option to compute all pairs
            structFact.define(ba,dmap,var_names,var_scaling);
        } else {
            // option to compute only specified pairs
            int nPairs = 2;
            amrex::Vector< int > s_pairA(nPairs);
            amrex::Vector< int > s_pairB(nPairs);

            // Select which variable pairs to include in structure factor:
            s_pairA[0] = 0;
            s_pairB[0] = 0;
            s_pairA[1] = 1;
            s_pairB[1] = 1;

            structFact.define(ba,dmap,var_names,var_scaling,s_pairA,s_pairB);
        }
    } else {
        // checkpoint restart
        Abort("checkpoint read not implemented yet");
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

    int istep = (restart < 0) ? 0 : restart;
    WritePlotFile(istep,time,geom,n_old);

    ///////////////////////////////////////////

    // Create output file for averaged density
    std::ofstream outputFile("averagedDensity.txt");
    outputFile << "time ";

    for (int comp = 0; comp < nspecies; ++comp) {
        outputFile << "comp_" << comp << " ";
    }
    outputFile << std::endl;

    // time step loop
    for(int step=step_start;step<=max_step;++step) {
        // store the current time so we can later compute total run time.
        Real step_strt_time = ParallelDescriptor::second();

        AdvanceTimestep(n_old,n_new,dt,time,geom);

        time += dt;
        MultiFab::Copy(n_old,n_new,0,0,nspecies,1);

        outputFile << std::setprecision(12) << time << " ";
        amrex::Print() << "time = " << time << " n_avg = ";

        // Compute average n for each species, print to file?
        for (int comp = 0; comp < nspecies; ++comp) {
            amrex::Real n_sum = n_old.sum(comp);
            amrex::Real n_avg = n_sum / cell_count;
            amrex::Print() << n_avg << " ";
            outputFile << std::setprecision(15) << n_avg << " ";
        }
        amrex::Print() << std::endl;
        outputFile << std::endl;

        // Call the timer again and compute the maximum difference between the start time
        // and stop time over all processors
        Real step_stop_time = ParallelDescriptor::second() - step_strt_time;
        ParallelDescriptor::ReduceRealMax(step_stop_time);
        amrex::Print() << "Time step " << step << " complted in " << step_stop_time << " seconds\n";

        // add a snapshot to the structure factor
        if (step > n_steps_skip && struct_fact_int > 0 && (step-n_steps_skip)%struct_fact_int == 0) {
            // add this snapshot to the average in the structure factor
            structFact.FortStructure(n_new);
        }

        if (plot_int > 0 && step%plot_int == 0) {
            WritePlotFile(step,time,geom,n_new);

            // write out structure factor to plotfile
            if (step > n_steps_skip && struct_fact_int > 0) {
                structFact.WritePlotFile(step,time,"plt_SF");
            }
        }

        if (chk_int > 0 && step%chk_int == 0) {
            Abort("checkpoint write not implemented yet");
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

    outputFile.close();

    // Call the timer again and compute the maximum difference between the start time
    // and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;
}