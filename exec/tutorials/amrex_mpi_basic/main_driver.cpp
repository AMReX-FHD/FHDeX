#include <chrono>

#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_ParallelDescriptor.H>

#include <common_functions.H>

using namespace std::chrono;
using namespace amrex;


// argv contains the name of the inputs file entered at the command line
void main_driver(const char * argv) {

    BL_PROFILE_VAR("main_driver()", main_driver);

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    // copy contents of F90 modules to C++ namespaces NOTE: any changes to
    // global settings in fortran/c++ after this point need to be synchronized
    // InitializeCommonNamespace();

    // Get my current MPI rank and communicator size
    int rank   = ParallelDescriptor::MyProc();
    int n_rank = ParallelDescriptor::NProcs();

    // This vector will hold all ranks
    Vector<int> r_vect(n_rank);

    // Starts out as empty
    Print() << "r_vect = ";
    for (auto & i:r_vect) Print() << i << " ";
    Print() << std::endl;

    // Add local rank and communicate with all other ranks
    r_vect[rank] = rank + 10;

    // Only one element (per rank) added
    std::cout << "[" << rank << "] r_vect = ";
    for (auto & i:r_vect) std::cout << i << " ";
    std::cout << std::endl;

    // MPI Comms
    ParallelDescriptor::ReduceIntSum(
        // MPI communicator needs a list of references to where the data goes
        Vector<std::reference_wrapper<int>>(begin(r_vect), end(r_vect))
    );

    // Now each rank as a list of all ranks
    Print() << "r_vect = ";
    for (auto & i:r_vect) Print() << i << " ";
    Print() << std::endl;

    // Just for fun, print out the max runtime
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    Print() << "Run time = " << stop_time << std::endl;
}