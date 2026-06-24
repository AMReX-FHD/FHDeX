#include "common_functions.H"


#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

#include "StructFact.H"

using namespace amrex;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{

    BL_PROFILE_VAR("main_driver()",main_driver);

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::string inputs_file = argv;

    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();

    MultiFab vel;
    StructFact turbStructFact;
    Geometry geom;

    {
        MultiFab mf;

        plot_init_file += "/Level_0/Cell";

        // read in the plotfile
        VisMF::Read(mf, plot_init_file);

        BoxArray ba = mf.boxArray();
        DistributionMapping dmap(ba);

        // copy velocity into vel
        vel.define(ba,dmap,AMREX_SPACEDIM,0);

        // "prob_type" is the starting component of shifted_vel in the input plotfile
        MultiFab::Copy(vel,mf,prob_type,0,AMREX_SPACEDIM,0);

        // use prob_lo/hi from inputs file
        RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                         {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

        // create a single Box that spans the domain of the (coarse) plotfile
        Box domain = ba.minimalBox();

        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            if (domain.bigEnd(i)+1 != n_cells[i]) {
                Print() << "Mistmatch direction " << i << " "
                        << domain.bigEnd(i)+1 << " " << n_cells[i] << std::endl;
                Abort("n_cells in inputs file does not match plotfile");
            }
        }

        Vector<int> is_periodic(AMREX_SPACEDIM,1);  // set to 1 (periodic) by default

        // geometry object required to build structure factor object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());

        const Real* dx = geom.CellSize();
        Real dVol = (AMREX_SPACEDIM==2) ? dx[0]*dx[1]*cell_depth : dx[0]*dx[1]*dx[2];

        // tell structure factor object we only care about these covariances
        amrex::Vector< int > s_pairA(AMREX_SPACEDIM);
        amrex::Vector< int > s_pairB(AMREX_SPACEDIM);
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            s_pairA[d] = d;
            s_pairB[d] = d;
        }

        // scale solution
        Vector<Real> var_scaling(AMREX_SPACEDIM);
        for (int d=0; d<var_scaling.size(); ++d) {
            var_scaling[d] = 1./dVol;
        }

        Vector< std::string > var_names(AMREX_SPACEDIM);

        int cnt = 0;
        std::string x;

        // velx, vely, velz
        for (int d=0; d<AMREX_SPACEDIM; d++) {
            x = "vel";
            x += (120+d);
            var_names[cnt++] = x;
        }

        turbStructFact.define(ba,dmap,var_names,var_scaling,s_pairA,s_pairB);
    }


    // reset and compute structure factor
    turbStructFact.FortStructure(vel,1);
    turbStructFact.CallFinalize();

    // integrate cov_mag over shells in k and write to file
    turbStructFact.IntegratekShells(0,geom);

    // Call the timer again and compute the maximum difference between the start time
    // and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;

}
