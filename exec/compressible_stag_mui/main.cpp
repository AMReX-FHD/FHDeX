#include <AMReX.H>
#include <AMReX_MPMD.H>

#ifdef MUI
// additional header files needed by MUI
#include <mpi.h>
#include <lib_mpi_split.h>
#include <mui.h>
#else
#include <AMReX_MPMD.H>
#endif

// function declaration
void main_driver (const char* argv);

int main (int argc, char* argv[])
{
#ifdef MUI
    MPI_Comm comm = mui::mpi_split_by_app( argc, argv );
#else
    MPI_Comm comm = amrex::MPMD::Initialize(argc, argv);
#endif

#ifdef SLURM
    // expected args: exec %t %o inputs_file ...
    if (argc<4) MPI_Abort(comm,0);
    // remove %t and %o from args list
    argc -= 2;
    for (int i=1;i<argc;i++) argv[i] = argv[i+2];
#else
    // expected args: exec inputs_file ...
    if (argc<2) MPI_Abort(comm,0);
#endif

    amrex::Initialize(argc,argv,true,comm);

    // argv[1] contains the name of the inputs file entered at the command line
    main_driver(argv[1]);

    amrex::Finalize();

#ifdef USE_AMREX_MPMD
    amrex::MPMD::Finalize();
#endif

    return 0;
}

