/*
 * An AMReX-based version of the StochasticHeat python code:
 * Intro FHD paper: https://arxiv.org/abs/2406.12157
 * Original python version: https://github.com/AlejGarcia/IntroFHD
 */

#include "common_functions.H"
#include "rng_functions.H"

#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Random.H>
#include <AMReX_FFT.H>

#include "chrono"

using namespace std::chrono;
using namespace amrex;

int main (int argc, char* argv[])
{

amrex::Initialize(argc,argv);
{

#if (AMREX_SPACEDIM != 2)
    amrex::Abort("Only works with DIM=2 (independent 1D pencils in x)");
#endif
    // **********************************
    // SIMULATION PARAMETERS

    // number of cells in each spatial direction
    int n_particles;
    int n_ensembles;
    int max_ensembles_per_rank; // for parallelization purposes

    // total steps in simulation and time step
    int nsteps;
    Real dt;

    // how often to write a plotfile
    int plot_int;

    // random number seed (positive integer=fixed seed; 0=clock-based seed)
    int seed;

    // size of each finite volume cell - all 3 must be defined regardless of dimensionality
    Real cell_dx;
    Real cell_dy;
    Real cell_dz;

    Real r0;
    Real p0;
    Real e0;
    
    Real A_00;
    Real A_01;
    Real A_02;
    Real A_10;
    Real A_11;
    Real A_12;
    Real A_20;
    Real A_21;
    Real A_22;

    Real D_00;
    Real D_01;
    Real D_02;
    Real D_10;
    Real D_11;
    Real D_12;
    Real D_20;
    Real D_21;
    Real D_22;

    Real B_00;
    Real B_01;
    Real B_02;
    Real B_10;
    Real B_11;
    Real B_12;
    Real B_20;
    Real B_21;
    Real B_22;


    // input parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        // pp.get means we require the inputs file to have it
        // pp.query means we optionally need the inputs file to have it - but we must supply a default here
        ParmParse pp;

        // override defaults set above
        pp.get("n_particles",n_particles);
        pp.get("n_ensembles",n_ensembles);
        pp.get("max_ensembles_per_rank",max_ensembles_per_rank);

        pp.get("nsteps",nsteps);
        pp.get("dt",dt);

        pp.get("plot_int",plot_int);

        pp.get("seed",seed);

        pp.get("cell_dx",cell_dx);
        pp.get("cell_dy",cell_dy);
        pp.get("cell_dz",cell_dz);

        pp.get("r0",r0);
        pp.get("p0",p0);
        pp.get("e0",e0);

        pp.get("A_00",A_00);
        pp.get("A_01",A_01);
        pp.get("A_02",A_02);
        pp.get("A_10",A_10);
        pp.get("A_11",A_11);
        pp.get("A_12",A_12);
        pp.get("A_20",A_20);
        pp.get("A_21",A_21);
        pp.get("A_22",A_22);

        pp.get("D_00",D_00);
        pp.get("D_01",D_01);
        pp.get("D_02",D_02);
        pp.get("D_10",D_10);
        pp.get("D_11",D_11);
        pp.get("D_12",D_12);
        pp.get("D_20",D_20);
        pp.get("D_21",D_21);
        pp.get("D_22",D_22);

        pp.get("B_00",B_00);
        pp.get("B_01",B_01);
        pp.get("B_02",B_02);
        pp.get("B_10",B_10);
        pp.get("B_11",B_11);
        pp.get("B_12",B_12);
        pp.get("B_20",B_20);
        pp.get("B_21",B_21);
        pp.get("B_22",B_22);
    }

    // initialize AMReX random seed (positive integer=fixed seed; 0=clock-based seed)
    if (seed > 0) {
        // initializes the seed for C++ random number calls
        InitRandom(seed+ParallelDescriptor::MyProc(),
                   ParallelDescriptor::NProcs(),
                   seed+ParallelDescriptor::MyProc());
    } else if (seed == 0) {
        // initializes the seed for C++ random number calls based on the clock
        auto now = time_point_cast<nanoseconds>(system_clock::now());
        int randSeed = now.time_since_epoch().count();
        // broadcast the same root seed to all processors
        ParallelDescriptor::Bcast(&randSeed,1,ParallelDescriptor::IOProcessorNumber());
        InitRandom(randSeed+ParallelDescriptor::MyProc(),
                   ParallelDescriptor::NProcs(),
                   randSeed+ParallelDescriptor::MyProc());
    } else {
        Abort("Must supply non-negative seed");
    }

    // **********************************
    // Set up simulation grid, which requires a BoxArray and Distribution Mapping
    // BoxArray ba will contain a list of boxes that cover the domain
    // DistributionMapping dm is a mapping of individual boxes to MPI ranks
    BoxArray ba;

    // AMREX_D_DECL means "do the first X of these, where X is the dimensionality of the simulation"
#if (AMREX_SPACEDIM == 2)
    IntVect dom_lo(            0,             0);
    IntVect dom_hi(n_particles-1, n_ensembles-1);
    IntVect max_grid_size(1024000,max_ensembles_per_rank);
#elif (AMREX_SPACEDIM == 3)
    IntVect dom_lo(            0,             0,             0);
    IntVect dom_hi(n_particles-1, n_particles-1, n_ensembles-1);
    IntVect max_grid_size(1024000,1024000,max_ensembles_per_rank);
#endif

    // Make a single box that is the entire domain
    Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "domain"
    ba.define(domain);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // How Boxes are distributed among MPI processes
    DistributionMapping dm(ba);

    // **********************************
    // Create MultiFab data structures needed for simulation and diagnostics

    MultiFab vars(ba, dm, 3, 1);

    // face-centered MultiFabs for noise and flux
    // dimension is one less than compile dimension (independent pencils or planes)
    Array<MultiFab, AMREX_SPACEDIM-1> noise;
    Array<MultiFab, AMREX_SPACEDIM-1> flux;
    for (int dir=0; dir<AMREX_SPACEDIM-1; dir++)
    {
        // noise[dir] and flux[dir] have three components, zero ghost cells, and is nodal in direction dir
        BoxArray edge_ba = ba;
        edge_ba.surroundingNodes(dir);
        noise[dir].define(edge_ba, dm, 3, 0);
        flux[dir].define(edge_ba, dm, 3, 0);

        noise[dir].setVal(0.);
    }

    // **********************************
    // Set up geometry
    // geom contains information such as the physical domain size,
    //               number of points in the domain, and periodicity

#if (AMREX_SPACEDIM == 2)
    RealBox real_box({                 0.,                  0.},
                     {cell_dx*n_particles, cell_dy*n_ensembles});
#elif (AMREX_SPACEDIM == 3)
    RealBox real_box({0.,                                   0.,                  0.},
                     {cell_dx*n_particles, cell_dy*n_particles, cell_dz*n_ensembles});

#endif

#if (AMREX_SPACEDIM == 2)
    // periodic in x direction
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,0,0)};
#elif (AMREX_SPACEDIM == 3)
    // periodic in x and y directions
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,0)};
#endif

    // This defines a Geometry object
    Geometry geom;
    geom.define(domain, real_box, CoordSys::cartesian, is_periodic);

    // extract dx from the geometry object
    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // volume of grid cell
    Real dV;
#if (AMREX_SPACEDIM == 2)
    dV = dx[0] * dx[1] * cell_dz;
#elif (AMREX_SPACEDIM == 3)
    dV = dx[0] * dx[1] * dx[2];
#endif

    // **********************************
    // INITIALIZE DATA

    // loop over boxes
    for (MFIter mfi(vars); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        const Array4<Real>& vars_fab = vars.array(mfi);

        amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k, amrex::RandomEngine const& engine)
        {
            vars_fab(i,j,k,0) = r0;
            vars_fab(i,j,k,1) = p0;
            vars_fab(i,j,k,2) = e0;
        });
    }

    // time = starting time in the simulation
    Real time = 0.0;

    // Write a plotfile of the initial data if plot_int > 0
    if (plot_int > 0)
    {
        int step = 0;
        const std::string& pltfile = amrex::Concatenate("plt",step,7);
        WriteSingleLevelPlotfile(pltfile, vars, {"r","p","e"}, geom, time, 0);
    }

    // **********************************
    // Time step loop
    for (int step = 1; step <= nsteps; ++step)
    {
        // fill periodic ghost cells
        vars.FillBoundary(geom.periodicity());

        // fill random numbers
        for (int d=0; d<AMREX_SPACEDIM-1; ++d) {
            MultiFabFillRandom(noise[d],0,1.,geom);
            MultiFabFillRandom(noise[d],1,1.,geom);
            MultiFabFillRandom(noise[d],2,1.,geom);
            noise[d].mult( 1./std::sqrt(dt*dV), 0, 3);
        }

        // compute fluxes
        for ( MFIter mfi(vars); mfi.isValid(); ++mfi )
        {
            const Array4<const Real>& vars_fab = vars.array(mfi);

            const Box& xbx = mfi.nodaltilebox(0);
            const Array4<Real>& fluxx = flux[0].array(mfi);
            const Array4<Real>& noisex = noise[0].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Box& ybx = mfi.nodaltilebox(1);
            const Array4<Real>& fluxy = flux[1].array(mfi);
            const Array4<Real>& noisey = noise[1].array(mfi);
#endif
            amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // advection
                // n=0; A_00 (var0) + A_01 (var1) + A_02 (var2)
                // n=1; A_10 (var0) + A_11 (var1) + A_12 (var2)
                // n=2; A_20 (var0) + A_21 (var1) + A_22 (var2)
                fluxx(i,j,k,0) = A_00 * 0.5 * (vars_fab(i,j,k,0) + vars_fab(i-1,j,k,0))
                               + A_01 * 0.5 * (vars_fab(i,j,k,1) + vars_fab(i-1,j,k,1))
                               + A_02 * 0.5 * (vars_fab(i,j,k,2) + vars_fab(i-1,j,k,2));
                fluxx(i,j,k,1) = A_10 * 0.5 * (vars_fab(i,j,k,0) + vars_fab(i-1,j,k,0))
                               + A_11 * 0.5 * (vars_fab(i,j,k,1) + vars_fab(i-1,j,k,1))
                               + A_12 * 0.5 * (vars_fab(i,j,k,2) + vars_fab(i-1,j,k,2));
                fluxx(i,j,k,2) = A_20 * 0.5 * (vars_fab(i,j,k,0) + vars_fab(i-1,j,k,0))
                               + A_21 * 0.5 * (vars_fab(i,j,k,1) + vars_fab(i-1,j,k,1))
                               + A_22 * 0.5 * (vars_fab(i,j,k,2) + vars_fab(i-1,j,k,2));

                // deterministic diffusion
                // n=0; D_00 d/dx(var0) + D_01 d/dx(var1) + D_02 d/dx(var2)
                // n=1; D_10 d/dx(var0) + D_11 d/dx(var1) + D_12 d/dx(var2)
                // n=2; D_20 d/dx(var0) + D_21 d/dx(var1) + D_22 d/dx(var2)
                fluxx(i,j,k,0) += (  D_00 * (vars_fab(i,j,k,0) - vars_fab(i-1,j,k,0))
                                   + D_01 * (vars_fab(i,j,k,1) - vars_fab(i-1,j,k,1))
                                   + D_02 * (vars_fab(i,j,k,2) - vars_fab(i-1,j,k,2)) ) / dx[0];
                fluxx(i,j,k,1) += (  D_10 * (vars_fab(i,j,k,0) - vars_fab(i-1,j,k,0))
                                   + D_11 * (vars_fab(i,j,k,1) - vars_fab(i-1,j,k,1))
                                   + D_12 * (vars_fab(i,j,k,2) - vars_fab(i-1,j,k,2)) ) / dx[0];
                fluxx(i,j,k,2) += (  D_20 * (vars_fab(i,j,k,0) - vars_fab(i-1,j,k,0))
                                   + D_21 * (vars_fab(i,j,k,1) - vars_fab(i-1,j,k,1))
                                   + D_22 * (vars_fab(i,j,k,2) - vars_fab(i-1,j,k,2)) ) / dx[0];

                // stochastic
                // n=0; B_00*noise0 + B_01*noise1 + B_02*noise2
                // n=1; B_10*noise0 + B_11*noise1 + B_12*noise2
                // n=2; B_20*noise0 + B_21*noise1 + B_22*noise2
                fluxx(i,j,k,0) += B_00 * noisex(i,j,k,0) + B_01 * noisex(i,j,k,1) + B_02 * noisex(i,j,k,2);
                fluxx(i,j,k,1) += B_10 * noisex(i,j,k,0) + B_11 * noisex(i,j,k,1) + B_12 * noisex(i,j,k,2);
                fluxx(i,j,k,2) += B_20 * noisex(i,j,k,0) + B_21 * noisex(i,j,k,1) + B_22 * noisex(i,j,k,2);
            });
#if (AMREX_SPACEDIM == 3)
            amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // advection
                // n=0; A_00 (var0) + A_01 (var1) + A_02 (var2)
                // n=1; A_10 (var0) + A_11 (var1) + A_12 (var2)
                // n=2; A_20 (var0) + A_21 (var1) + A_22 (var2)
                fluxy(i,j,k,0) = A_00 * 0.5 * (vars_fab(i,j,k,0) + vars_fab(i,j-1,k,0))
                               + A_01 * 0.5 * (vars_fab(i,j,k,1) + vars_fab(i,j-1,k,1))
                               + A_02 * 0.5 * (vars_fab(i,j,k,2) + vars_fab(i,j-1,k,2));
                fluxy(i,j,k,1) = A_10 * 0.5 * (vars_fab(i,j,k,0) + vars_fab(i,j-1,k,0))
                               + A_11 * 0.5 * (vars_fab(i,j,k,1) + vars_fab(i,j-1,k,1))
                               + A_12 * 0.5 * (vars_fab(i,j,k,2) + vars_fab(i,j-1,k,2));
                fluxy(i,j,k,2) = A_20 * 0.5 * (vars_fab(i,j,k,0) + vars_fab(i,j-1,k,0))
                               + A_21 * 0.5 * (vars_fab(i,j,k,1) + vars_fab(i,j-1,k,1))
                               + A_22 * 0.5 * (vars_fab(i,j,k,2) + vars_fab(i,j-1,k,2));

                // deterministic diffusion
                // n=0; D_00 d/dy(var0) + D_01 d/dy(var1) + D_02 d/dy(var2)
                // n=1; D_10 d/dy(var0) + D_11 d/dy(var1) + D_12 d/dy(var2)
                // n=2; D_20 d/dy(var0) + D_21 d/dy(var1) + D_22 d/dy(var2)
                fluxy(i,j,k,0) += (  D_00 * (vars_fab(i,j,k,0) - vars_fab(i,j-1,k,0))
                                   + D_01 * (vars_fab(i,j,k,1) - vars_fab(i,j-1,k,1))
                                   + D_02 * (vars_fab(i,j,k,2) - vars_fab(i,j-1,k,2)) ) / dx[1];
                fluxy(i,j,k,1) += (  D_10 * (vars_fab(i,j,k,0) - vars_fab(i,j-1,k,0))
                                   + D_11 * (vars_fab(i,j,k,1) - vars_fab(i,j-1,k,1))
                                   + D_12 * (vars_fab(i,j,k,2) - vars_fab(i,j-1,k,2)) ) / dx[1];
                fluxy(i,j,k,2) += (  D_20 * (vars_fab(i,j,k,0) - vars_fab(i,j-1,k,0))
                                   + D_21 * (vars_fab(i,j,k,1) - vars_fab(i,j-1,k,1))
                                   + D_22 * (vars_fab(i,j,k,2) - vars_fab(i,j-1,k,2)) ) / dx[1];

                // stochastic
                // n=0; B_00*noise0 + B_01*noise1 + B_02*noise2
                // n=1; B_10*noise0 + B_11*noise1 + B_12*noise2
                // n=2; B_20*noise0 + B_21*noise1 + B_22*noise2
                fluxy(i,j,k,0) += B_00 * noisey(i,j,k,0) + B_01 * noisey(i,j,k,1) + B_02 * noisey(i,j,k,2);
                fluxy(i,j,k,1) += B_10 * noisey(i,j,k,0) + B_11 * noisey(i,j,k,1) + B_12 * noisey(i,j,k,2);
                fluxy(i,j,k,2) += B_20 * noisey(i,j,k,0) + B_21 * noisey(i,j,k,1) + B_22 * noisey(i,j,k,2);
            });
#endif
        }

        // advance the data by dt
        // new_vars = old_vars + dt * div(flux)
        for ( MFIter mfi(vars); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.validbox();

            const Array4<Real>& vars_fab = vars.array(mfi);

            const Array4<Real>& fluxx = flux[0].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<Real>& fluxy = flux[1].array(mfi);
#endif
            amrex::ParallelFor(bx, 3, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                vars_fab(i,j,k,n) = vars_fab(i,j,k,n) + dt *
                      (fluxx(i+1,j,k,n) - fluxx(i,j,k,n)) / dx[0]
#if (AMREX_SPACEDIM == 3)
                    + (fluxy(i,j+1,k,n) - fluxy(i,j,k,n)) / dx[1]
#endif
                        ;

            });
        }

        // update time
        time = time + dt;

        // Tell the I/O Processor to write out which step we're doing
        amrex::Print() << "Advanced step " << step << "\n";

        // Write a plotfile of the current data (plot_int was defined in the inputs file)
        if (plot_int > 0 && step%plot_int == 0)
        {
            const std::string& pltfile = amrex::Concatenate("plt",step,7);
            WriteSingleLevelPlotfile(pltfile, vars, {"r","p","e"}, geom, time, 0);
        }
    }

}
amrex::Finalize();
return 0;

}
