#include <chrono>

#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_ParallelDescriptor.H>

#include <common_functions.H>
#include <ib_functions.H>
#include <immbdy_namespace.H>

#include <IBMarkerContainer.H>

using namespace std::chrono;
using namespace amrex;

using namespace immbdy;
using namespace ib_flagellum;


// argv contains the name of the inputs file entered at the command line
void main_driver(const char * argv) {

    BL_PROFILE_VAR("main_driver()", main_driver);

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();


    //___________________________________________________________________________
    // Load parameters from inputs file, and initialize global parameters
    std::string inputs_file = argv;

    read_immbdy_namelist(inputs_file.c_str(), inputs_file.size() + 1);

    // copy contents of F90 modules to C++ namespaces NOTE: any changes to
    // global settings in fortran/c++ after this point need to be synchronized
    InitializeCommonNamespace();
    InitializeImmbdyNamespace();
    InitializeIBFlagellumNamespace();


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
    // Initialize immersed boundaries
    // Make sure that the nghost (last argument) is big enough!

    BL_PROFILE_VAR("main_create_markers", create_markers);

    // Find the optimal number of ghost cells for the IBMarkerContainer
    Real min_dx = dx[0];
    for (int d=1; d<AMREX_SPACEDIM; ++d)
	    min_dx = amrex::min(min_dx, dx[d]);

    int ib_nghost = 1;
    Print() << "Initializing IBMarkerContainer with "
            << ib_nghost << " ghost cells" << std::endl;
    // Initialize immersed boundary container
    IBMarkerContainer ib_mc(geom, dmap, ba, ib_nghost);

    for (int i_ib=0; i_ib < n_immbdy; ++i_ib) {

        if (n_marker[i_ib] <= 0) continue;

        int N  = n_marker[i_ib];
        Real L = ib_flagellum::length[i_ib];

        Real l_link = L/(N-1);

        const RealVect & x_0 = offset_0[i_ib];

        Print() << "Initializing flagellum:" << std::endl;
        Print() << "N=      " << N           << std::endl;
        Print() << "L=      " << L           << std::endl;
        Print() << "l_link= " << l_link      << std::endl;
        Print() << "x_0=    " << x_0         << std::endl;


        // using fourier modes => first two nodes reserved as "anchor"
        int N_markers = immbdy::contains_fourier ? N+1 : N;

        Vector<RealVect> marker_positions(N_markers);
        for (int i=0; i<marker_positions.size(); ++i) {
            Real x = x_0[0] + i*l_link;
            // Compute periodic offset. Will work as long as winding number = 1
            Real x_period = x < geom.ProbHi(0) ? x : x - geom.ProbLength(0);

            marker_positions[i] = RealVect{x_period, x_0[1], x_0[2]};
        }

        Vector<Real> marker_radii(N_markers);
        for (int i=0; i<marker_radii.size(); ++i) marker_radii[i] = 4*l_link;

        ib_mc.InitList(0, marker_radii, marker_positions, i_ib);
    }

    ib_mc.fillNeighbors();
    ib_mc.PrintMarkerData(0);
    BL_PROFILE_VAR_STOP(create_markers);

    int rank = ParallelDescriptor::MyProc();
    int size = ParallelDescriptor::NProcs();
    Vector<int> mpi_markers(size);
    mpi_markers[rank] = ib_mc.NextID() - 1;
    ParallelDescriptor::ReduceIntSum(
        Vector<std::reference_wrapper<int>>(
            begin(mpi_markers), end(mpi_markers)
        )
    );

    Print() << "mpi_markers = ";
    for (auto & i:mpi_markers) Print() << i << " ";
    Print() << std::endl;


    // Just for fun, print out the max runtime
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    Print() << "Run time = " << stop_time << std::endl;
}