#include "spectral_functions.H"
#include "common_functions.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;

void main_driver(const char* argv)
{
    BL_PROFILE_VAR("main_driver()", main_driver);
    amrex::ignore_unused(argv);

    // AMReX has already parsed argv[1], the original hydro inputs file.
    // Populate FHDeX's common namespace before using its geometry and BC data.
    InitializeCommonNamespace();

    Real ts1 = ParallelDescriptor::second();

    ParmParse pp;

    std::string restart_file;
    if (!pp.query("restart_file", restart_file)) {
        Abort("SPECTRAL_FILTER requires restart_file=<hydro checkpoint directory>");
    }

    Real kmin = 0.0;
    // Real kmax = 0.0;
    std::vector<Real> kmax_list;
    if (!pp.query("kmin", kmin)) {
        Abort("SPECTRAL_FILTER requires kmin=<minimum integer wavenumber>");
    }
    // if (!pp.query("kmax", kmax)) {
    //     Abort("SPECTRAL_FILTER requires kmax=<maximum integer wavenumber>");
    // }
    if (!pp.queryarr("kmax_list", kmax_list)) {
        Abort("SPECTRAL_FILTER requires kmax_list=<list of maximum integer wavenumbers>");
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(kmin >= 0.0, "kmin must be non-negative");
    // AMREX_ALWAYS_ASSERT_WITH_MESSAGE(kmax >= kmin, "kmax must be greater than or equal to kmin");
    // Loop through the vector to validate each kmax value
    for (Real kmax : kmax_list) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(kmax >= kmin, "Every kmax must be greater than or equal to kmin");
    }
    int plot_filter = 0;
    pp.query("plot_filter", plot_filter);

    int plot_fourier = 0;
    pp.query("plot_fourier", plot_fourier);

    Vector<int> is_periodic(AMREX_SPACEDIM, 1);
    int const nper = pp.countval("is_periodic");
    if (nper > 0) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            nper == AMREX_SPACEDIM,
            "is_periodic must contain AMREX_SPACEDIM entries");
        pp.getarr("is_periodic", is_periodic, 0, AMREX_SPACEDIM);
    }

    Geometry geom;
    BoxArray ba;
    DistributionMapping dmap;
    MultiFab velocity;
    int step = 0;
    Real time = 0.0;

    SpectralReadCheckPoint(restart_file, is_periodic, geom, ba, dmap, velocity, step, time);

    MultiFab velocity_filter(ba, dmap, 3, 0);

#if (AMREX_SPACEDIM == 3)
    int constexpr num_vv_comps = 6;
#else
    int constexpr num_vv_comps = 3;
#endif

    MultiFab vv_filter(ba, dmap, num_vv_comps, 0);
    

    velocity.FillBoundary(geom.periodicity());

    // Loop over all of the kmax values
    for (Real kmax : kmax_list) {
        velocity_filter.setVal(0.0);
        vv_filter.setVal(0.0);
        
        // Filter the velocity
        SpectralVelDecomp(velocity, velocity_filter, kmin, kmax, geom);
        velocity_filter.FillBoundary(geom.periodicity());
    
        // Filter the outer product of the velocity
        SpectralVelProductDecomp(velocity, vv_filter, kmin, kmax, geom);
        vv_filter.FillBoundary(geom.periodicity());
    
        // velocity.FillBoundary(geom.periodicity());
        // SpectralVelDecomp(velocity, velocity_filter, kmin, kmax, geom);
        // velocity_filter.FillBoundary(geom.periodicity());
    
        if (plot_filter != 0) {
            SpectralWritePlotFile(step, kmin, kmax, geom, velocity, velocity_filter, vv_filter);
        }
        if (plot_fourier != 0) {
            SpectralWriteFourierPlotFile(step, kmin, kmax, geom, velocity_filter, vv_filter);
        }
    }

    Real ts2 = ParallelDescriptor::second() - ts1;
    ParallelDescriptor::ReduceRealMax(ts2, ParallelDescriptor::IOProcessorNumber());
    Print() << "Time (spectral filtering) " << ts2 << " seconds\n";

    Long min_fab_megabytes = TotalBytesAllocatedInFabsHWM() / 1048576;
    Long max_fab_megabytes = min_fab_megabytes;

    ParallelDescriptor::ReduceLongMin(min_fab_megabytes, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceLongMax(max_fab_megabytes, ParallelDescriptor::IOProcessorNumber());

    Print() << "High-water FAB megabyte spread across MPI nodes: ["
            << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";
}
