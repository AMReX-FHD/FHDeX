
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void main_driver ();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    
    main_driver();

    amrex::Finalize();

    return 0;
}
