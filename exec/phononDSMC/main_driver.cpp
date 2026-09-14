#include "LocalFunctions.H"
#include "species.H"
#include "paramPlane.H"
#include "StructFact.H"
#include "particle_functions.H"
#include "Checkpoint.H"
#include "chrono"
#include "iostream"
#include "fstream"
#include "DsmcParticleContainer.H"
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

using namespace std::chrono;
using namespace std;

void main_driver(const char* argv)
{
    // timer for total simulation time
    Real strt_time = ParallelDescriptor::second();
    std::string inputs_file = argv;

    InitializeCommonNamespace();

    BoxArray ba;
    IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
    IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
    Box domain(dom_lo, dom_hi);
    DistributionMapping dmap;

    int step = 1;
    Real dt = fixed_dt;
    Real time = 0.;
    int statsCount = 1;

    iMultiFab bCell;

    MultiFab cuInst, cuMeans, cuVars;

    if (seed > 0)
    {
        InitRandom(seed+ParallelDescriptor::MyProc(),
                   ParallelDescriptor::NProcs(),
                   seed+ParallelDescriptor::MyProc());
    }
    else if (seed == 0)
    {
        auto now = time_point_cast<nanoseconds>(system_clock::now());
        int randSeed = now.time_since_epoch().count();
        ParallelDescriptor::Bcast(&randSeed,1,ParallelDescriptor::IOProcessorNumber());
        InitRandom(randSeed+ParallelDescriptor::MyProc(),
                   ParallelDescriptor::NProcs(),
                   randSeed+ParallelDescriptor::MyProc());
    }
    else
    {
        Abort("Must supply non-negative seed");
    }


    if (restart < 0)
    {
        if (seed > 0)
        {
            InitRandom(seed+ParallelDescriptor::MyProc(),
                       ParallelDescriptor::NProcs(),
                       seed+ParallelDescriptor::MyProc());
        }
        else if (seed == 0)
        {
            auto now = time_point_cast<nanoseconds>(system_clock::now());
            int randSeed = now.time_since_epoch().count();
            // broadcast the same root seed to all processors
            ParallelDescriptor::Bcast(&randSeed,1,ParallelDescriptor::IOProcessorNumber());
            InitRandom(randSeed+ParallelDescriptor::MyProc(),
                       ParallelDescriptor::NProcs(),
                       randSeed+ParallelDescriptor::MyProc());
        }
        else
        {
            Abort("Must supply non-negative seed");
        }

        ba.define(domain);
        ba.maxSize(IntVect(max_grid_size));
        dmap.define(ba);
        //////////////////////////////////////
        // Conserved/Primitive Var Setup
        //////////////////////////////////////
        /*
            Conserved Vars:
            0  - rho = (1/V) += m
            1  - Jx  = (1/V) += mu
            2  - Jy  = (1/V) += mv
            3  - Jz  = (1/V) += mw
            4  - K   = (1/V) += m|v|^2
            ... (repeat for each species)
        */

        int ncon  = 5;
        cuInst.define(ba, dmap, ncon, 0); cuInst.setVal(0.);
        cuMeans.define(ba, dmap, ncon, 0); cuMeans.setVal(0.);
        cuVars.define(ba,dmap, ncon, 0); cuVars.setVal(0.);



    }
    else
    {
    //  ReadCheckPoint(step, time, dt, statsCount,
    //      cuInst, cuMeans, cuVars,
    //      primInst, primMeans, primVars,
    //      coVars, spatialCross1D, ncross);
    //  dmap = cuInst.DistributionMap();
    //  ba = cuInst.boxArray();

    }

    Vector<int> is_periodic (AMREX_SPACEDIM,0);
    for (int i=0; i<AMREX_SPACEDIM; ++i)
    {
        if (bc_vel_lo[i] == -1 && bc_vel_hi[i] == -1)
        {
            is_periodic [i] = -1;
        }
    }

    // This defines a Geometry object
    RealBox realDomain({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
        {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

    Geometry geom (domain ,&realDomain,CoordSys::cartesian,is_periodic.data());

    std::ifstream planeFile("paramplanes.dat");
    int fileCount = 0;
    if(planeFile.good())
    {
        planeFile >> fileCount;
    }
    planeFile.close();

    int paramPlaneCount = fileCount;
    paramPlane* paramPlaneList;
    paramPlaneList = new paramPlane[paramPlaneCount];
    BuildParamplanesPhonon(paramPlaneList,paramPlaneCount,realDomain.lo(),realDomain.hi());


    bCell.define(ba, dmap, paramPlaneCount, 1); bCell.setVal(-1);



    // Particle tile size
    Vector<int> ts(BL_SPACEDIM);

    for (int d=0; d<AMREX_SPACEDIM; ++d)
    {
        if (max_particle_tile_size[d] > 0)
        {
            ts[d] = max_particle_tile_size[d];
        }
        else
        {
            ts[d] = max_grid_size[d];
        }
    }

    ParmParse pp ("particles");
    pp.addarr("tile_size", ts);

    int cRange = 0;
    FhdParticleContainer particles(geom, dmap, ba, cRange);

    //particles.InitParticles(dt);

    Real init_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(init_time);
    amrex::Print() << "Initialization time = " << init_time << " seconds " << std::endl;


    Real tbegin, tend;
    int writeStep = 0;

    for (int istep=step; istep<=max_step; ++istep)
    {
        tbegin = ParallelDescriptor::second();

        if ((n_steps_skip > 0 && istep == n_steps_skip) || (n_steps_skip < 0 && istep%n_steps_skip == 0) )
        {
            writeStep = istep;
        }

        particles.SourcePhonons(dt, paramPlaneList, paramPlaneCount);

        particles.MovePhononsCPP(dt, paramPlaneList, paramPlaneCount, writeStep, istep, bCell);

        particles.EvaluateStatsPhonon(cuInst,cuMeans,cuVars,statsCount++,time);

        //////////////////////////////////////
        // PlotFile
        //////////////////////////////////////


        if ((plot_int > 0 && istep%plot_int==0))
        {
            writePlotFile(cuInst,cuMeans,cuVars,particles,geom,time,istep);
        }

    //    if (chk_int > 0 && istep%chk_int == 0 && istep > step)
    //    {
    //        WriteCheckPoint(istep, time, dt, statsCount,
    //            cuInst, cuMeans, cuVars, primInst, primMeans, primVars, coVars,
    //            particles, spatialCross1D, ncross);
    //    }

        if ((n_steps_skip > 0 && istep == n_steps_skip) || (n_steps_skip < 0 && istep%n_steps_skip == 0) ) {
            //reset stats
            statsCount = 1;
            cuMeans.setVal(0.);
            cuVars.setVal(0.);
        }
        tend = ParallelDescriptor::second() - tbegin;
        ParallelDescriptor::ReduceRealMax(tend);
        if(istep%100==0)
        {
            amrex::Print() << "Advanced step " << istep << " in " << tend << " seconds\n";
        }

        time += dt;
    }

    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << " seconds" << std::endl;
}