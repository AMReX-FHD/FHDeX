#include <AMReX.H>
#include <AMReX_MultiFabUtil.H>
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

    // define lower and upper indices
    IntVect dom_lo(0,0);
    IntVect dom_hi(n_particles-1, n_ensembles-1);

    // Make a single box that is the entire domain
    Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "domain"
    BoxArray ba(domain);

    IntVect max_grid_size(1024000,max_ensembles_per_rank);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // physical box size in this FPU_chain implementation is not relevant
    RealBox real_box({ 0., 0.,}, { 1., 1.,});

    // periodic in all direction
    Array<int,AMREX_SPACEDIM> is_periodic{1,0};

    // This defines a Geometry object
    Geometry geom(domain, real_box, CoordSys::cartesian, is_periodic);

    // time = starting time in the simulation
    Real time = 0.0;

    // we only need a ghost cell for the x-direction
    IntVect ng_vect(1,0);

    // components are r, p, and e
    MultiFab state(ba,dm,3,ng_vect);

    MultiFab S_alphaalpha(ba,dm,6,0);
    S_alphaalpha.setVal(0.);
    int samples = 0;
    Gpu::HostVector<Real> S_alphaalpha_00(n_particles);
    Gpu::HostVector<Real> S_alphaalpha_01(n_particles);
    Gpu::HostVector<Real> S_alphaalpha_02(n_particles);
    Gpu::HostVector<Real> S_alphaalpha_11(n_particles);
    Gpu::HostVector<Real> S_alphaalpha_12(n_particles);
    Gpu::HostVector<Real> S_alphaalpha_22(n_particles);

    // BoxArray to store g_alpha(0,0), will use the same distribution map
    IntVect dom_hi_zero(0, n_ensembles-1);
    Box domain_zero(dom_lo, dom_hi_zero);
    BoxArray ba_zero(domain_zero);
    ba_zero.maxSize(max_grid_size);
    Geometry geom_zero(domain_zero, real_box, CoordSys::cartesian, is_periodic);

    MultiFab g_alpha_zero(ba_zero,dm,3,0);

    // ******************************
    // SAMPLE TO OBTAIN INITIAL STATE
    // ******************************
    init(state, beta, pressure, a_coef, b_coef, c_coef, 0., 10000, 1.e-3, n_particles, n_ensembles, geom);
    compute_energy(state,a_coef,b_coef,c_coef);
    compute_means(state,n_particles,n_ensembles,0);

    g_alpha_zero.ParallelCopy(state, 0, 0, 3);

    // initial plotfile
    if (plot_int > 0) {
        const std::string& pltfile = amrex::Concatenate("plt",0,7);
        amrex::Print() << "Writing plotfile " << pltfile << std::endl;
        WriteSingleLevelPlotfile(pltfile, state, {"r","p","e"}, geom, time, 0);
    }

    for (int step=1; step<=n_steps; ++step) {

        time += dt;

        // ****************
        // INTEGRATE A STEP
        // ****************
        FPU_RK4(state,a_coef,b_coef,c_coef,dt,n_particles,n_ensembles,geom);
        amrex::Print() << "Completed step " << step << std::endl;

        compute_S_alphaalpha(state,g_alpha_zero,S_alphaalpha);

        ++samples;
        S_alphaalpha_00 = sumToLine(S_alphaalpha, 0, 1, domain, 0);
        S_alphaalpha_01 = sumToLine(S_alphaalpha, 1, 1, domain, 0);
        S_alphaalpha_02 = sumToLine(S_alphaalpha, 2, 1, domain, 0);
        S_alphaalpha_11 = sumToLine(S_alphaalpha, 3, 1, domain, 0);
        S_alphaalpha_12 = sumToLine(S_alphaalpha, 4, 1, domain, 0);
        S_alphaalpha_22 = sumToLine(S_alphaalpha, 5, 1, domain, 0);

        if (plot_int > 0 && step%plot_int == 0) {
            const std::string& pltfile = amrex::Concatenate("S_alphaalpha",step,7);
            amrex::Print() << "Writing plotfile " << pltfile << std::endl;
            WriteSingleLevelPlotfile(pltfile, S_alphaalpha, {"00","01","02","11","12","22"}, geom, time, step);

            const std::string S_alphaalphafile = amrex::Concatenate("S_alphaalpha_avg",step,7);
            std::ofstream S_alphaalphaout;
            if (ParallelDescriptor::IOProcessor()) {
                S_alphaalphaout.open(S_alphaalphafile, std::ios::out);
                for (int i=0; i<n_particles; ++i) {
                    S_alphaalpha_00[i] /= (samples*n_ensembles);
                    S_alphaalphaout << " S_alphaalpha i = " << i << " "
                                    << S_alphaalpha_00[i] << " "
                                    << S_alphaalpha_01[i] << " "
                                    << S_alphaalpha_02[i] << " "
                                    << S_alphaalpha_11[i] << " "
                                    << S_alphaalpha_12[i] << " "
                                    << S_alphaalpha_22[i] << "\n";
                }
            }
    }

        // ********
        // PLOTFILE
        // ********
        if (plot_int > 0 && step%plot_int == 0) {
            const std::string& pltfile = amrex::Concatenate("plt",step,7);
            amrex::Print() << "Writing plotfile " << pltfile << std::endl;
            WriteSingleLevelPlotfile(pltfile, state, {"r","p","e"}, geom, time, step);
        }

        // ****************
        // TEXT DIAGNOSTICS
        // ****************
        compute_energy(state,a_coef,b_coef,c_coef);
        if (plot_int > 0 && step%diag_int == 0) {
            compute_means(state,n_particles,n_ensembles,step);
        }

    }



}
Finalize();
return 0;
}

