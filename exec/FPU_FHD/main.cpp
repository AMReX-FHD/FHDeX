/*
 * An AMReX-based version of the StochasticHeat python code:
 * Intro FHD paper: https://arxiv.org/abs/2406.12157
 * Original python version: https://github.com/AlejGarcia/IntroFHD
 */

#include "common_functions.H"
#include "rng_functions.H"

#include "FPU.H"

#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Random.H>
#include <AMReX_FFT.H>

#include <chrono>
#include <cmath>

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

    // how many steps to skip before defining t=0
    int n_steps_skip;
    
    // how often to write a plotfile
    int plot_int;

    // how often to write out correlation diagnostics
    int diag_int;

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

    Real R_00;
    Real R_01;
    Real R_02;
    Real R_10;
    Real R_11;
    Real R_12;
    Real R_20;
    Real R_21;
    Real R_22;

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

        pp.get("n_steps_skip",n_steps_skip);

        pp.get("plot_int",plot_int);
        pp.get("diag_int",diag_int);

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

        pp.get("R_00",R_00);
        pp.get("R_01",R_01);
        pp.get("R_02",R_02);
        pp.get("R_10",R_10);
        pp.get("R_11",R_11);
        pp.get("R_12",R_12);
        pp.get("R_20",R_20);
        pp.get("R_21",R_21);
        pp.get("R_22",R_22);
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
        amrex::Long randSeed = now.time_since_epoch().count();
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

    MultiFab state(ba, dm, 3, 1);

    // storage for phi = R*state
    MultiFab phi (ba,dm,3,0);
    MultiFab phi0(ba,dm,3,0);

    // diagonal elements of C_alphaalpha
    MultiFab C_alphaalpha(ba,dm,3,0);

    Gpu::HostVector<Real> C_alphaalpha_00(n_particles);
    Gpu::HostVector<Real> C_alphaalpha_11(n_particles);
    Gpu::HostVector<Real> C_alphaalpha_22(n_particles);

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
    // INITIALIZE DATA TO ZERO MEAN UNIT VARIANCE RNG
    MultiFabFillRandom(state,0,1.,geom);
    MultiFabFillRandom(state,1,1.,geom);
    MultiFabFillRandom(state,2,1.,geom);

    // time = starting time in the simulation
    Real time = 0.0;

    // Write a plotfile of the initial data if plot_int > 0
    if (plot_int > 0)
    {
        int step = 0;
        const std::string& pltfile = amrex::Concatenate("plt",step,7);
        WriteSingleLevelPlotfile(pltfile, state, {"r","p","e"}, geom, time, 0);
    }

    // **********************************
    // Time step loop
    for (int step = 1; step <= nsteps; ++step)
    {
        if (step == n_steps_skip) {
            Copy(phi0,state,0,0,3,0);
            ComputePhiFromState(phi0,0.,0.,0.,R_00,R_01,R_02,R_10,R_11,R_12,R_20,R_21,R_22);
        }

        // fill periodic ghost cells
        state.FillBoundary(geom.periodicity());

        // fill random numbers
        for (int d=0; d<AMREX_SPACEDIM-1; ++d) {
            MultiFabFillRandom(noise[d],0,1.,geom);
            MultiFabFillRandom(noise[d],1,1.,geom);
            MultiFabFillRandom(noise[d],2,1.,geom);
            noise[d].mult( 1./std::sqrt(dt*dV), 0, 3);
        }

        // compute fluxes
        for ( MFIter mfi(state); mfi.isValid(); ++mfi )
        {
            const Array4<const Real>& state_fab = state.array(mfi);

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
                fluxx(i,j,k,0) = A_00 * 0.5 * (state_fab(i,j,k,0) + state_fab(i-1,j,k,0))
                               + A_01 * 0.5 * (state_fab(i,j,k,1) + state_fab(i-1,j,k,1))
                               + A_02 * 0.5 * (state_fab(i,j,k,2) + state_fab(i-1,j,k,2));
                fluxx(i,j,k,1) = A_10 * 0.5 * (state_fab(i,j,k,0) + state_fab(i-1,j,k,0))
                               + A_11 * 0.5 * (state_fab(i,j,k,1) + state_fab(i-1,j,k,1))
                               + A_12 * 0.5 * (state_fab(i,j,k,2) + state_fab(i-1,j,k,2));
                fluxx(i,j,k,2) = A_20 * 0.5 * (state_fab(i,j,k,0) + state_fab(i-1,j,k,0))
                               + A_21 * 0.5 * (state_fab(i,j,k,1) + state_fab(i-1,j,k,1))
                               + A_22 * 0.5 * (state_fab(i,j,k,2) + state_fab(i-1,j,k,2));

                // deterministic diffusion
                // n=0; D_00 d/dx(var0) + D_01 d/dx(var1) + D_02 d/dx(var2)
                // n=1; D_10 d/dx(var0) + D_11 d/dx(var1) + D_12 d/dx(var2)
                // n=2; D_20 d/dx(var0) + D_21 d/dx(var1) + D_22 d/dx(var2)
                fluxx(i,j,k,0) += (  D_00 * (state_fab(i,j,k,0) - state_fab(i-1,j,k,0))
                                   + D_01 * (state_fab(i,j,k,1) - state_fab(i-1,j,k,1))
                                   + D_02 * (state_fab(i,j,k,2) - state_fab(i-1,j,k,2)) ) / dx[0];
                fluxx(i,j,k,1) += (  D_10 * (state_fab(i,j,k,0) - state_fab(i-1,j,k,0))
                                   + D_11 * (state_fab(i,j,k,1) - state_fab(i-1,j,k,1))
                                   + D_12 * (state_fab(i,j,k,2) - state_fab(i-1,j,k,2)) ) / dx[0];
                fluxx(i,j,k,2) += (  D_20 * (state_fab(i,j,k,0) - state_fab(i-1,j,k,0))
                                   + D_21 * (state_fab(i,j,k,1) - state_fab(i-1,j,k,1))
                                   + D_22 * (state_fab(i,j,k,2) - state_fab(i-1,j,k,2)) ) / dx[0];

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
                fluxy(i,j,k,0) = A_00 * 0.5 * (state_fab(i,j,k,0) + state_fab(i,j-1,k,0))
                               + A_01 * 0.5 * (state_fab(i,j,k,1) + state_fab(i,j-1,k,1))
                               + A_02 * 0.5 * (state_fab(i,j,k,2) + state_fab(i,j-1,k,2));
                fluxy(i,j,k,1) = A_10 * 0.5 * (state_fab(i,j,k,0) + state_fab(i,j-1,k,0))
                               + A_11 * 0.5 * (state_fab(i,j,k,1) + state_fab(i,j-1,k,1))
                               + A_12 * 0.5 * (state_fab(i,j,k,2) + state_fab(i,j-1,k,2));
                fluxy(i,j,k,2) = A_20 * 0.5 * (state_fab(i,j,k,0) + state_fab(i,j-1,k,0))
                               + A_21 * 0.5 * (state_fab(i,j,k,1) + state_fab(i,j-1,k,1))
                               + A_22 * 0.5 * (state_fab(i,j,k,2) + state_fab(i,j-1,k,2));

                // deterministic diffusion
                // n=0; D_00 d/dy(var0) + D_01 d/dy(var1) + D_02 d/dy(var2)
                // n=1; D_10 d/dy(var0) + D_11 d/dy(var1) + D_12 d/dy(var2)
                // n=2; D_20 d/dy(var0) + D_21 d/dy(var1) + D_22 d/dy(var2)
                fluxy(i,j,k,0) += (  D_00 * (state_fab(i,j,k,0) - state_fab(i,j-1,k,0))
                                   + D_01 * (state_fab(i,j,k,1) - state_fab(i,j-1,k,1))
                                   + D_02 * (state_fab(i,j,k,2) - state_fab(i,j-1,k,2)) ) / dx[1];
                fluxy(i,j,k,1) += (  D_10 * (state_fab(i,j,k,0) - state_fab(i,j-1,k,0))
                                   + D_11 * (state_fab(i,j,k,1) - state_fab(i,j-1,k,1))
                                   + D_12 * (state_fab(i,j,k,2) - state_fab(i,j-1,k,2)) ) / dx[1];
                fluxy(i,j,k,2) += (  D_20 * (state_fab(i,j,k,0) - state_fab(i,j-1,k,0))
                                   + D_21 * (state_fab(i,j,k,1) - state_fab(i,j-1,k,1))
                                   + D_22 * (state_fab(i,j,k,2) - state_fab(i,j-1,k,2)) ) / dx[1];

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
        // new_state = old_state + dt * div(flux)
        for ( MFIter mfi(state); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.validbox();

            const Array4<Real>& state_fab = state.array(mfi);

            const Array4<Real>& fluxx = flux[0].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<Real>& fluxy = flux[1].array(mfi);
#endif
            amrex::ParallelFor(bx, 3, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                state_fab(i,j,k,n) = state_fab(i,j,k,n) + dt *
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
            WriteSingleLevelPlotfile(pltfile, state, {"r","p","e"}, geom, time, 0);

        }

        if (diag_int > 0 && step%diag_int == 0 && step > n_steps_skip) {
            Copy(phi,state,0,0,3,0);
            ComputePhiFromState(phi,0.,0.,0.,R_00,R_01,R_02,R_10,R_11,R_12,R_20,R_21,R_22);

            ComputeCalphaalpha(C_alphaalpha,phi,phi0);
            C_alphaalpha_00 = sumToLine(C_alphaalpha, 0, 1, domain, 0);
            C_alphaalpha_11 = sumToLine(C_alphaalpha, 1, 1, domain, 0);
            C_alphaalpha_22 = sumToLine(C_alphaalpha, 2, 1, domain, 0);

            const std::string C_alphaalphafile = amrex::Concatenate("C_alphaalpha",step,7);
            amrex::Print() << "Writing C_alphaalphafile " << C_alphaalphafile << std::endl;
            std::ofstream C_alphaalphaout;
            if (ParallelDescriptor::IOProcessor()) {
                C_alphaalphaout.open(C_alphaalphafile, std::ios::out);
                for (int i=0; i<n_particles; ++i) {
                        
                    C_alphaalphaout << " C_alphaalpha_00/11/22 = " << i << " "
                                    << C_alphaalpha_00[(i+n_particles/2)%n_particles] << " "
                                    << C_alphaalpha_11[(i+n_particles/2)%n_particles] << " "
                                    << C_alphaalpha_22[(i+n_particles/2)%n_particles] << "\n";
                }
            }
        }

    }

}
amrex::Finalize();
return 0;

}
