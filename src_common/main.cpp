
#include "common.H"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    
    main_driver();

    amrex::Finalize();

    return 0;
}
