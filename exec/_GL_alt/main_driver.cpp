
//#include "hydro_test_functions.H"
//#include "hydro_test_functions_F.H"

//#include "hydro_functions.H"
//#include "hydro_functions_F.H"
//#include "StochMFlux.H"
//#include "StructFact.H"

#include "rng_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "common_namespace.H"
#include "common_namespace_declarations.H"

//#include "compressible_functions.H"
//#include "compressible_functions_F.H"

#include "exec_functions.H"
#include "GL_functions.H"
#include "GL_functions_F.H"

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFab.H>

using namespace amrex;
using namespace common;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::string inputs_file = argv;

    // read in parameters from inputs file into F90 modules
    // we use "+1" because of amrex_string_c_to_f expects a null char termination
    read_common_namelist(inputs_file.c_str(),inputs_file.size()+1);
    read_GL_namelist(inputs_file.c_str(),inputs_file.size()+1);

    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();

  
    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_lo[i] == -1 && bc_hi[i] == -1) {
            is_periodic[i] = 1;
            Print() << "Periodic: " << is_periodic[i] << "\n";
        }
        //is_periodic[i] = 0;
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
        IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);

        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        // note we are converting "Vector<int> max_grid_size" to an IntVect
        ba.maxSize(IntVect(max_grid_size));

       // This defines the physical box, [-1,1] in each direction.
        RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                         {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    Real dt = fixed_dt;

    const Real* dx = geom.CellSize();

    setdt(dx, &dt);

    Real dtinv = 1.0/dt;

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);

    /////////////////////////////////////////
    //Initialise rngs
    /////////////////////////////////////////
    const int n_rngs = 1;

    const int proc = ParallelDescriptor::MyProc();

    int fhdSeed = 0;

    fhdSeed += 10000*proc;

    //Initialise rngs
    rng_init(&fhdSeed);
    /////////////////////////////////////////

    //conserved quantaties
    MultiFab phi  (ba,dmap,nvars,ngc);
    MultiFab phin (ba,dmap,nvars,ngc);

    //possibly for later
    MultiFab rannums(ba,dmap,nvars,ngc);


    Real time = 0;

    int step, statsCount;

    step = 0;

    statsCount = 1;

     Init_Phi(phi,dx);

     WritePlotFile(step, time, geom, phi);

     if( proc == 0 ) {
      
        std::cout << "DT = " << dt << std::endl;
        std::cout << "Mesh spacing = " << dx[0] << " , " << dx[1] << std::endl;
        std::cout << " Domain = [ 0 , " << prob_hi[0] << " ] x [ 0 , " << prob_hi[1] <<" ]" << std::endl;

     }

    amrex::Real integral;

    //Time stepping loop
    for(step=1;step<=max_step;++step)
    {

        RK2step(phi, phin, rannums, geom, dx, dt, integral, step);

//        if(step == n_steps_skip)
//     {
//     on the fly statistics setup?
//      }

//        evaluateStats(phi, dx);

//        statsCount++;

        inc_phi0(&step);

        if(step%500 == 0)
        {    
                amrex::Print() << "Integral " << integral << "\n";
                amrex::Print() << "Advanced step " << step << "\n";
        }

        if (plot_int > 0 && step > 0 && step%plot_int == 0)
        {

           WritePlotFile(step, time, geom, phi);
        }

        time = time + dt;
    }

    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;

//    amrex::Finalize();

}

