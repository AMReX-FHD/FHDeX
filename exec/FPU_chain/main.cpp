#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>

#include "FPU.H"

// for clock-based random seed
#include "chrono"
using namespace std::chrono;

using namespace amrex;

int main (int argc, char* argv[])
{
Initialize(argc,argv);
{

    if (AMREX_SPACEDIM != 2) {
        Abort("Must build with DIM=2");
    }

    // **********************************
    // DECLARE SIMULATION PARAMETERS
    // **********************************

    int n_particles = 8000;
    int n_ensembles = 10;
    int n_steps = 1000000;

    Real a_coef = 1.0;
    Real b_coef = 0.0;
    Real c_coef = 1.0;

    Real beta = 1.0;
    Real pressure = 1.0;

    int seed = 1;

    Real dt = 1.e-3;

    int max_ensembles_per_rank = 1000;

    int diag_int = 1000;
    int plot_int = 10000;

    // ***************************************************************
    // READ INPUT PARAMETER VALUES FROM INPUT FILE AND/OR COMMAND LINE
    // ***************************************************************
    {
        // ParmParse is way of reading inputs from the inputs file
        // pp.get means we require the inputs file to have it
        // pp.query means we optionally need the inputs file to have it - but we must supply a default here
        ParmParse pp;

        pp.query("n_particles",n_particles);
        pp.query("n_ensembles",n_ensembles);
        pp.query("n_steps",n_steps);

        pp.query("a_coef",a_coef);
        pp.query("b_coef",b_coef);
        pp.query("c_coef",c_coef);

        pp.query("beta",beta);
        pp.query("pressure",pressure);

        pp.query("seed",seed);

        pp.query("dt",dt);

        pp.query("max_ensembles_per_rank",max_ensembles_per_rank);

        pp.query("diag_int",diag_int);
        pp.query("plot_int",plot_int);
    }

    int the_seed;
    if (seed > 0) {
        // fixed seed
        the_seed = seed;
    } else if (seed == 0) {
        // initializes the seed for C++ random number calls based on the clock
        // broadcast the same seed to all processors
        auto now = time_point_cast<nanoseconds>(system_clock::now());
        the_seed = now.time_since_epoch().count();
        ParallelDescriptor::Bcast(&the_seed,1,ParallelDescriptor::IOProcessorNumber());
    } else {
        Abort("Must supply non-negative seed");
    }
    InitRandom(the_seed+ParallelDescriptor::MyProc(),
               ParallelDescriptor::NProcs(),
               seed+ParallelDescriptor::MyProc());

    BoxArray ba;
    Geometry geom;

    // define lower and upper indices
    IntVect dom_lo(0,0);
    IntVect dom_hi(n_particles-1, n_ensembles-1);

    // Make a single box that is the entire domain
    Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "domain"
    ba.define(domain);

    IntVect max_grid_size(1024000,max_ensembles_per_rank);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // physical box size in this FPU_chain implementation is not relevant
    RealBox real_box({ 0., 0.,}, { 1., 1.,});

    // periodic in all direction
    Array<int,AMREX_SPACEDIM> is_periodic{1,0};

    // This defines a Geometry object
    geom.define(domain, real_box, CoordSys::cartesian, is_periodic);

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // time = starting time in the simulation
    Real time = 0.0;

    // we only need a ghost cell for the x-direction
    IntVect ng_vect(1,0);

    // components are r and p
    MultiFab state(ba,dm,2,ng_vect);

    // for plotfile
    MultiFab plt_mf(ba,dm,2,0);

    // ******************************
    // SAMPLE TO OBTAIN INITIAL STATE
    // ******************************
    init(state, beta, pressure, a_coef, b_coef, c_coef, 0., 10000, 1.e-3, n_particles, n_ensembles, geom);

    // initial plotfile
    if (plot_int > 0) {
        const std::string& pltfile = amrex::Concatenate("plt",0,7);
        amrex::Print() << "Writing plotfile " << pltfile << std::endl;
        WriteSingleLevelPlotfile(pltfile, state, {"r","p"}, geom, time, 0);
    }

    for (int step=1; step<=n_steps; ++step) {

        time += dt;

        // ****************
        // INTEGRATE A STEP
        // ****************
        FPU_RK4(state,a_coef,b_coef,c_coef,dt,n_particles,n_ensembles,geom);
        amrex::Print() << "Completed step " << step << std::endl;

        // ********
        // PLOTFILE
        // ********
        if (plot_int > 0 && step%plot_int == 0) {
            const std::string& pltfile = amrex::Concatenate("plt",step,7);
            amrex::Print() << "Writing plotfile " << pltfile << std::endl;
            WriteSingleLevelPlotfile(pltfile, state, {"p","r"}, geom, time, step);
        }

        // ****************
        // TEXT DIAGNOSTICS
        // ****************
        //
        //
        //

    }



}
Finalize();
return 0;
}

