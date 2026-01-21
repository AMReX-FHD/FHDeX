/*
 * A simplified single file version of the HeatEquation_EX0_C exmaple.
 * This code is designed to be used with Demo_Tutorial.rst.
 *
 */

#include "common_functions.H"
#include "rng_functions.H"

#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Random.H>

#include "chrono"

using namespace std::chrono;
using namespace amrex;

int main (int argc, char* argv[])
{

amrex::Initialize(argc,argv);
{
    // **********************************
    // SIMULATION PARAMETERS
        
    // number of cells on each side of the domain
    int n_cell = 32;

    // total steps in simulation
    int nsteps = 8000000;

    // how often to write a plotfile
    int plot_int = -1;

    // random number seed (1=fixed seed; 0=clock-based seed)
    int seed = 1;

    // Non-equilibrium flag (0: Thermodynamic equilibrium; 1: Temperature gradient)
    int NONEQ_FLAG = 0;

    // Start from perturbed initial condition (0: No; 1: Yes)
    int PERTURB_FLAG = 1;

    // Thermal fluctuations? (0: No, deterministic; 1: Yes, stochastic)
    int STOCH_FLAG = 1;

    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        // pp.get means we require the inputs file to have it
        // pp.query means we optionally need the inputs file to have it - but we must supply a default here
        ParmParse pp;

        // override defaults set above
        pp.query("n_cell",n_cell);
        pp.query("nsteps",nsteps);
        pp.query("plot_int",plot_int);
        pp.query("seed",seed);
        pp.query("NONEQ_FLAG",NONEQ_FLAG);
        pp.query("PERTURB_FLAG",PERTURB_FLAG);
        pp.query("STOCH_FLAG",STOCH_FLAG);
    }

    if (NONEQ_FLAG == 1) {
        Abort("NONEQ_FLAG=1 requires non-periodic boundaries");
    }
    
    //* Set physical parameters for the system (iron bar)
    Real kB = 1.38e-23;              // Boltzmann constant (J/K)
    Real mAtom = 9.27e-26;           // Mass of iron atom (kg)
    Real rho = 7870.;                // Mass density of iron (kg/m^3)
    Real c_V = 450.;                 // Specific heat capacity of iron (J/(kg K))
    Real ThCond = 70.;               // Thermal conductivity of iron (W/(m K))
    Real Length = 2.0e-8;            // System length (m)
    Real Area = std::pow(2.0e-9,2);  // System cross-sectional area (m^2)

    Real kappa = ThCond / (rho*c_V); // Coefficient in deterministic heat equation

    // Coefficient in stochastic heat equation
    Real alpha = (STOCH_FLAG==1) ? std::sqrt(2.*kB*kappa / (rho*c_V)) : 0.;

    Real stabilityFactor = 0.1;      // Numerical stability if stabilityFactor < 1.

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
    // SIMULATION SETUP

    // make BoxArray and Geometry
    // ba will contain a list of boxes that cover the domain
    // geom contains information such as the physical domain size,
    //               number of points in the domain, and periodicity
    BoxArray ba;
    Geometry geom;

    // AMREX_D_DECL means "do the first X of these, where X is the dimensionality of the simulation"
    IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));

    // Make a single box that is the entire domain
    Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "domain"
    ba.define(domain);

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // Temperature and initial temperature
    MultiFab Temp (ba, dm, 1, 1);
    MultiFab Temp0(ba, dm, 1, 0);

    // diagnostics
    MultiFab sumT(ba,dm,1,0);
    MultiFab sumT2(ba,dm,1,0);
    MultiFab sumTT(ba,dm,1,0);
    MultiFab sumSk(ba,dm,1,0);
    sumT.setVal(0.);
    sumT2.setVal(0.);
    sumTT.setVal(0.);
    sumSk.setVal(0.);

    MultiFab aveT(ba,dm,1,0);
    MultiFab varT(ba,dm,1,0);
    MultiFab corrT(ba,dm,1,0);
    aveT.setVal(0.);
    varT.setVal(0.);
    corrT.setVal(0.);

    MultiFab plotfile(ba,dm,3,0);
    
    int iCorr = n_cell/4;

    // This defines the physical box, [0,1] in each direction.
    RealBox real_box({AMREX_D_DECL(    0.,    0.,    0.)},
                     {AMREX_D_DECL(Length,Length,Length)});

    // periodic in all direction
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};

    // This defines a Geometry object
    geom.define(domain, real_box, CoordSys::cartesian, is_periodic);

    // extract dx from the geometry object
    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // volume of grid cell
    Real dV = Area*dx[0];
#if (AMREX_SPACEDIM != 1)
    Abort("Fix dV for multidimensional case");
#endif

    Real dt = stabilityFactor * dx[0] * dx[0] / (2.*kappa);

    Real Tref = 300.;                // Reference temperature (K)
    Real Tdiff = 400.;               // Temperature difference across the system for NONEQ_FLAG=1

    Real T_Left  = (NONEQ_FLAG==1) ? Tref - Tdiff/2. : Tref;
    Real T_Right = (NONEQ_FLAG==1) ? Tref + Tdiff/2. : Tref;

    // Standard deviation of temperature in a cell at the reference temperature
    Real Tref_SD = std::sqrt(kB*Tref*Tref / (rho*c_V*dV));

    // face-centered MultiFabs for noise and flux
    Array<MultiFab, AMREX_SPACEDIM> noise;
    Array<MultiFab, AMREX_SPACEDIM> flux;
    for (int dir=0; dir<AMREX_SPACEDIM; dir++)
    {
        // noise[dir] and flux[dir] have one component, zero ghost cells, and is nodal in direction dir
        BoxArray edge_ba = ba;
        edge_ba.surroundingNodes(dir);
        noise[dir].define(edge_ba, dm, 1, 0);
        flux[dir].define(edge_ba, dm, 1, 0);
    }

    // time = starting time in the simulation
    Real time = 0.0;

    // **********************************
    // INITIALIZE DATA

    // loop over boxes
    for (MFIter mfi(Temp); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        const Array4<Real>& Temp_fab = Temp.array(mfi);
        const Array4<Real>& Temp0_fab = Temp0.array(mfi);

        amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k, amrex::RandomEngine const& engine)
        {
            Real x = (i+0.5) * dx[0];
            Temp0_fab(i,j,k) = T_Left + (T_Right - T_Left) * x / Length;
            Temp_fab(i,j,k) = Temp0_fab(i,j,k);

            if (PERTURB_FLAG == 1) {
                Temp_fab(i,j,k) += Tref_SD*amrex::RandomNormal(0.,1.,engine);
            }
        });
    }

    // normalize initial noise perturbation
    Real avgT = Temp.sum(0) / n_cell;
    Temp.plus(-(avgT-0.5*(T_Left+T_Right)),0);

    // Write a plotfile of the initial data if plot_int > 0
    if (plot_int > 0)
    {
        int step = 0;
        const std::string& pltfile = amrex::Concatenate("plt",step,7);
        MultiFab::Copy(plotfile,Temp,0,0,1,0);
        MultiFab::Copy(plotfile,aveT,0,1,1,0);
        MultiFab::Copy(plotfile,varT,0,2,1,0);
        WriteSingleLevelPlotfile(pltfile, plotfile, {"Temp","avgT","varT"}, geom, time, 0);
    }

    int Nsamp = 0;
    
    for (int step = 1; step <= nsteps; ++step)
    {        
        // fill periodic ghost cells
        Temp.FillBoundary(geom.periodicity());

        // fill random numbers
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFabFillRandom(noise[d],0.,1.,geom);
            noise[d].mult( 1./std::sqrt(dt*dV), 0, 1);
        }
        
        // compute fluxes
        for ( MFIter mfi(Temp); mfi.isValid(); ++mfi )
        {
            const Array4<const Real>& Temp_fab = Temp.array(mfi);

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
            amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real det = kappa*(Temp_fab(i,j,k) - Temp_fab(i-1,j,k) ) / dx[0];
                Real sto = alpha*0.5*(Temp_fab(i,j,k) + Temp_fab(i-1,j,k))*noisex(i,j,k);
                fluxx(i,j,k) = det + sto;
            });
#if (AMREX_SPACEDIM >= 2)
            amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real det = kappa*(Temp_fab(i,j,k) - Temp_fab(i,j-1,k) ) / dx[1];
                Real sto = alpha*0.5*(Temp_fab(i,j,k) + Temp_fab(i,j-1,k))*noisey(i,j,k)*dx[1];
                fluxy(i,j,k) = det + sto;
            });
#if (AMREX_SPACEDIM == 3)
            amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real det = kappa*(Temp_fab(i,j,k) - Temp_fab(i,j,k-1) ) / dx[2];
                Real sto = alpha*0.5*(Temp_fab(i,j,k) + Temp_fab(i,j,k-1))*noisez(i,j,k)*dx[2];
                fluxz(i,j,k) = det + sto;
            });
#endif
#endif
        }

        // advance the data by dt
        // new_Temp = old_Temp + dt * div(flux)
        for ( MFIter mfi(Temp); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.validbox();

            const Array4<Real>& Temp_fab = Temp.array(mfi);

            const Array4<Real>& sumT_fab = sumT.array(mfi);
            const Array4<Real>& sumT2_fab = sumT2.array(mfi);
            const Array4<Real>& sumTT_fab = sumTT.array(mfi);
            
            const Array4<Real>& fluxx = flux[0].array(mfi);
#if (AMREX_SPACEDIM >= 2)
            const Array4<Real>& fluxy = flux[1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<Real>& fluxz = flux[1].array(mfi);
#endif
#endif
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Temp_fab(i,j,k) = Temp_fab(i,j,k) + dt *
                    (fluxx(i+1,j,k) - fluxx(i,j,k)) / dx[0]
#if (AMREX_SPACEDIM >= 2)
                    + (fluxy(i,j+1,k) - fluxy(i,j,k)) / dx[1]
#if (AMREX_SPACEDIM >= 3)
                    + (fluxz(i,j,k+1) - fluxz(i,j,k)) / dx[2]
#endif
#endif
                        ;

                sumT_fab(i,j,k) += Temp_fab(i,j,k);
                sumT2_fab(i,j,k) += Temp_fab(i,j,k)*Temp_fab(i,j,k);
                sumTT_fab(i,j,k) += Temp_fab(i,j,k)*Temp_fab(iCorr,j,k);

            });
        }

        // increment number of samples
        Nsamp++;

        // diagnostics - average and variance
        for ( MFIter mfi(Temp); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.validbox();

            const Array4<Real>& sumT_fab = sumT.array(mfi);
            const Array4<Real>& sumT2_fab = sumT2.array(mfi);
            const Array4<Real>& sumTT_fab = sumTT.array(mfi);

            const Array4<Real>& aveT_fab = aveT.array(mfi);
            const Array4<Real>& varT_fab = varT.array(mfi);
            const Array4<Real>& corrT_fab = corrT.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                aveT_fab(i,j,k) = sumT_fab(i,j,k) / Nsamp;
                varT_fab(i,j,k) = sumT2_fab(i,j,k) / Nsamp - aveT_fab(i,j,k)*aveT_fab(i,j,k);
            });
        }

        // diagnostics - correlation (need averages first)
        for ( MFIter mfi(Temp); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.validbox();

            const Array4<Real>& sumT_fab = sumT.array(mfi);
            const Array4<Real>& sumT2_fab = sumT.array(mfi);
            const Array4<Real>& sumTT_fab = sumT.array(mfi);

            const Array4<Real>& aveT_fab = aveT.array(mfi);
            const Array4<Real>& varT_fab = varT.array(mfi);
            const Array4<Real>& corrT_fab = corrT.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                corrT_fab(i,j,k) = sumTT_fab(i,j,k) / Nsamp - aveT_fab(i,j,k)*aveT_fab(iCorr,j,k);
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
            MultiFab::Copy(plotfile,Temp,0,0,1,0);
            MultiFab::Copy(plotfile,aveT,0,1,1,0);
            MultiFab::Copy(plotfile,varT,0,2,1,0);
            WriteSingleLevelPlotfile(pltfile, plotfile, {"Temp","avgT","varT"}, geom, time, step);
        }
    }

}
amrex::Finalize();
return 0;

}


