#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include <AmrCoreInertial.H>

using namespace amrex;

// Inertial Dean-Kawasaki model in (x, v): the SPDE or the particles, chosen by
// inertial.run_mode, on the same phase-space grid. Build with DIM = 2.
int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
        BL_PROFILE("main()");
        const auto strt_total = amrex::second();

        AmrCoreInertial sim;
        sim.InitData();
        sim.Evolve();

        auto end_total = amrex::second() - strt_total;
        ParallelDescriptor::ReduceRealMax(end_total, ParallelDescriptor::IOProcessorNumber());
        amrex::Print() << "\nTotal Time: " << end_total << '\n';
    }

    amrex::Finalize();
}
