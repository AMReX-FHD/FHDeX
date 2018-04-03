
#include <FHD.H>

using namespace amrex;

int main(int argc, char* argv[])
{

    // in AMReX.cpp
    Initialize(argc,argv);

    // timer for profiling
    BL_PROFILE_VAR("main()", main);

    // wallclock time
    const Real strt_total = ParallelDescriptor::second();

    {
        // declare an FHD object to manage data
        FHD fhd;

        // read in C++/F90 parameters
        // define global C++/F90 variables
        // set up boundary conditions
        // set istep, t_new, t_old
        // allocate MultiFabs and initialize MultiFab data
        fhd.Init();

        // advance solution to final time
        fhd.Evolve();
	
        // wallclock time
        Real end_total = ParallelDescriptor::second() - strt_total;
	
        // print wallclock time
        ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
        Print() << "\nTotal Time: " << end_total << '\n';
    }

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(main);

    // in AMReX.cpp
    Finalize();
}
