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
    // **********************************
    // SIMULATION PARAMETERS

    // number of cells in each spatial direction
    int n_cell;

    // total steps in simulation
    int nsteps;

    // how often to write a plotfile
    int plot_int;

    // random number seed (positive integer=fixed seed; 0=clock-based seed)
    int seed;

    Real dt;
    
    Real Length;  // System length (m)
    Real Area;    // System cross-sectional area (m^2)

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

    Real B_11;
    Real B_22;

    
    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        // pp.get means we require the inputs file to have it
        // pp.query means we optionally need the inputs file to have it - but we must supply a default here
        ParmParse pp;

        // override defaults set above
        pp.get("n_cell",n_cell);
        pp.get("nsteps",nsteps);
        pp.get("plot_int",plot_int);
        pp.get("seed",seed);
        pp.get("dt",dt);
        pp.get("Length",Length);
        pp.get("Area",Area);
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
        pp.get("B_11",B_11);
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
    IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));

    // Make a single box that is the entire domain
    Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "domain"
    ba.define(domain);

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // **********************************
    // Create MultiFab data structures needed for simulation and diagonstics

    MultiFab vars(ba, dm, 3, 1);

    // face-centered MultiFabs for noise and flux
    Array<MultiFab, AMREX_SPACEDIM> noise;
    Array<MultiFab, AMREX_SPACEDIM> flux;
    for (int dir=0; dir<AMREX_SPACEDIM; dir++)
    {
        // noise[dir] and flux[dir] have one component, zero ghost cells, and is nodal in direction dir
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

    // This defines the physical box, [0,1] in each direction.
    RealBox real_box({AMREX_D_DECL(    0.,    0.,    0.)},
                     {AMREX_D_DECL(Length,Length,Length)});

    // periodic in all direction
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};

    // This defines a Geometry object
    Geometry geom;
    geom.define(domain, real_box, CoordSys::cartesian, is_periodic);

    // extract dx from the geometry object
    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // volume of grid cell
    Real dV = Area*dx[0];
#if (AMREX_SPACEDIM != 1)
    Abort("Fix dV for multidimensional case");
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
            vars_fab(i,j,k,0) = 0.;
            vars_fab(i,j,k,1) = 0.;
            vars_fab(i,j,k,2) = 0.;
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
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            // no noise for component 0; fill components 1 and 2 only
            MultiFabFillRandom(noise[d],1,1.,geom);
            MultiFabFillRandom(noise[d],2,1.,geom);
            noise[d].mult( 1./std::sqrt(dt*dV), 0, 1);
        }

        // compute fluxes
        for ( MFIter mfi(vars); mfi.isValid(); ++mfi )
        {
            const Array4<const Real>& vars_fab = vars.array(mfi);

            const Box& xbx = mfi.nodaltilebox(0);
            const Array4<Real>& fluxx = flux[0].array(mfi);
            const Array4<Real>& noisex = noise[0].array(mfi);
#if (AMREX_SPACEDIM >= 2)
            const Box& ybx = mfi.nodaltilebox(1);
            const Array4<Real>& fluxy = flux[1].array(mfi);
            const Array4<Real>& noisey = noise[1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Box& zbx = mfi.nodaltilebox(2);
            const Array4<Real>& fluxz = flux[1].array(mfi);
            const Array4<Real>& noisez = noise[1].array(mfi);
#endif
#endif
            amrex::ParallelFor(xbx, 3, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
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
                // n=0; add nothing
                // n=1; B_11 noise1
                // n=1; B_22 noise2
                fluxx(i,j,k,1) += B_11 * noisex(i,j,k,1);
                fluxx(i,j,k,2) += B_22 * noisex(i,j,k,2);
            });
#if (AMREX_SPACEDIM >= 2)
            Abort("Fix y-fluxes for multidimensional case");
            amrex::ParallelFor(ybx, 3, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
            });
#if (AMREX_SPACEDIM == 3)
            Abort("Fix z-fluxes for multidimensional case");
            amrex::ParallelFor(zbx, 3, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
            });
#endif
#endif
        }

        // advance the data by dt
        // new_vars = old_vars + dt * div(flux)
        for ( MFIter mfi(vars); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.validbox();

            const Array4<Real>& vars_fab = vars.array(mfi);

            const Array4<Real>& fluxx = flux[0].array(mfi);
#if (AMREX_SPACEDIM >= 2)
            const Array4<Real>& fluxy = flux[1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<Real>& fluxz = flux[1].array(mfi);
#endif
#endif
            amrex::ParallelFor(bx, 3, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                vars_fab(i,j,k,n) = vars_fab(i,j,k,n) + dt *
                      (fluxx(i+1,j,k,n) - fluxx(i,j,k,n)) / dx[0]
#if (AMREX_SPACEDIM >= 2)
                    + (fluxy(i,j+1,k,n) - fluxy(i,j,k,n)) / dx[1]
#if (AMREX_SPACEDIM >= 3)
                    + (fluxz(i,j,k+1,n) - fluxz(i,j,k,n)) / dx[2]
#endif
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


