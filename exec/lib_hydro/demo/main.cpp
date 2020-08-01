#include <AMReX.H>

// function declaration
void main_driver (const char* argv);

int main (int argc, char* argv[])
{
    {
        amrex::Initialize(argc,argv);

        // argv[1] contains the name of the inputs file entered at the command line
        main_driver(argv[1]);

        amrex::Finalize();
    }
    return 0;
}
