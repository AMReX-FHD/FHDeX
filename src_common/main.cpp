
#include "common.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    
    main_driver();

    amrex::Finalize();

    return 0;
}

void main_main ()
{
}
