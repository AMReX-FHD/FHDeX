
#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

// The 3D build is the SPDE (x, y, theta = z); the 2D build is the particle
// code with theta as a particle attribute.
#if (AMREX_SPACEDIM == 3)
#include <AmrCoreFlock.H>
#elif (AMREX_SPACEDIM == 2)
#include <FlockParticles.H>
#endif

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
        BL_PROFILE("main()");
        const auto strt_total = amrex::second();

#if (AMREX_SPACEDIM == 3)
        AmrCoreFlock flock;
        flock.InitData();
        flock.Evolve();
#elif (AMREX_SPACEDIM == 2)
        FlockParticles flock;
        flock.InitData();
        flock.Evolve();
#else
        amrex::Abort("flocking: build with DIM=3 (SPDE) or DIM=2 (particles)");
#endif

        auto end_total = amrex::second() - strt_total;
        ParallelDescriptor::ReduceRealMax(end_total, ParallelDescriptor::IOProcessorNumber());
        amrex::Print() << "\nTotal Time: " << end_total << '\n';
    }

    amrex::Finalize();
}
