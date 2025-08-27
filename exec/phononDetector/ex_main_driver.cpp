#include "INS_functions.H"
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
    InitializeGmresNamespace();

    BoxArray ba;
    IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
    IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
    Box domain(dom_lo, dom_hi);
    DistributionMapping dmap;

    MultiFab cuInst, cuMeans, cuVars;
    MultiFab primInst, primMeans, primVars;
    MultiFab coVars;

    int step = 0;
    Real dt = 0;
    int statsCount = 1;
    Real time = 0.;

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

    ba.define(domain);
    ba.maxSize(IntVect(max_grid_size));
    dmap.define(ba);

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

    int paramPlaneCount = 6;
    paramPlane paramPlaneList[paramPlaneCount];
    BuildParamplanes(paramPlaneList,paramPlaneCount,realDomain.lo(),realDomain.hi());

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

    // Output all primitives for structure factor
    int nvarstruct = 6+nspecies*2;
    const Real* dx = geom.CellSize();
    int nstruct = std::ceil((double)nvarstruct*(nvarstruct+1)/2);
    // scale SF results by inverse cell volume
    Vector<Real> var_scaling(nstruct);
    for (int d=0; d<var_scaling.size(); ++d)
    {
        var_scaling[d] = 1./(dx[0]*dx[1]*dx[2]);
    }

    // Collision Cell Vars
    particles.mfselect.define(ba, dmap, nspecies*nspecies, 0);
    particles.mfselect.setVal(0.);
    particles.mfphi.define(ba, dmap, nspecies, 0);
    particles.mfphi.setVal(0.);
    particles.mfvrmax.define(ba, dmap, nspecies*nspecies, 0);
    particles.mfvrmax.setVal(0.);

    particles.InitParticles(dt);

    particles.InitCollisionCells();
    amrex::Print() << "Overwritten dt so Courant number <1: " << dt << "\n";

    Real init_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(init_time);
    amrex::Print() << "Initialization time = " << init_time << " seconds " << std::endl;

    max_step += step;
    n_steps_skip += step;
    Real tbegin, tend;
    for (int istep=step; istep<=max_step; ++istep)
    {
        if(istep%IO_int == 0)
        {
            tbegin = ParallelDescriptor::second();
        }

        particles.CalcSelections(dt);
        particles.CollideParticles(dt);
        particles.Source(dt, paramPlaneList, paramPlaneCount);
        //particles.externalForce(dt);
        particles.MoveParticlesCPP(dt, paramPlaneList, paramPlaneCount);
        //particles.updateTimeStep(geom,dt);

        tend = ParallelDescriptor::second() - tbegin;
        ParallelDescriptor::ReduceRealMax(tend);
        amrex::Print() << "Advanced step " << istep << " in " << tend << " seconds\n";

        time += dt;
    }

    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << " seconds" << std::endl;
}