#include <AMReX.H>
//#include <AMReX_ParmParse.H>

// function declaration
void main_driver (const char* argv);

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

// this specific part has been moved to Flagellum/main_driver.cpp
//    {
//        amrex::ParmParse pp("particles");
//#ifdef AMREX_USE_GPU
//        bool particles_do_tiling = true;
//#else
//        bool particles_do_tiling = false;
//#endif
//        pp.queryAdd("do_tiling", particles_do_tiling);
//    }

    // argv[1] contains the name of the inputs file entered at the command line
    main_driver(argv[1]);

    amrex::Finalize();

    return 0;
}
