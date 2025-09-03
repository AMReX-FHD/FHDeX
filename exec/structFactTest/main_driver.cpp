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

    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,1);  // set to 1 (periodic) by default

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
        IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);

        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        // note we are converting "Vector<int> max_grid_size" to an IntVect
        ba.maxSize(IntVect(max_grid_size));

       // This defines the physical box, [-1,1] in each direction.
        RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                         {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    const Real* dx = geom.CellSize();

    // we may need to use dVol for scaling
    Real dVol = dx[0]*dx[1];
    int tot_n_cells = n_cells[0]*n_cells[1];
    if (AMREX_SPACEDIM == 2) {
        dVol *= cell_depth;
    } else if (AMREX_SPACEDIM == 3) {
        dVol *= dx[2];
        tot_n_cells = n_cells[2]*tot_n_cells;
    }

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);

    /////////////////////////////////////////

    Vector< std::string > var_names(2);
    var_names[0] = "phi1";
    var_names[1] = "phi2";

    // for the 3 pairs
    Vector< Real > var_scaling(3);
    var_scaling[0] = 1./dVol;
    var_scaling[1] = 1./dVol;
    var_scaling[2] = 1./dVol;

    MultiFab struct_cc(ba, dmap, 2, 0);

    // WRITE INIT ROUTINE
    struct_cc.setVal(0.);

    for (MFIter mfi(struct_cc,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real>& struct_fab = struct_cc.array(mfi);

        amrex::ParallelFor(bx, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (n == 0) {
                struct_fab(i,j,k,n) = i+j+k;
            }
            else if (n == 1) {
                struct_fab(i,j,k,n) = 0.5*sqrt(i+j+k);
            }

            // hooks to specify any function
            /*
            Real x = prob_lo[0] + (i+0.5)*dx[0];
            Real y = prob_lo[0] + (i+0.5)*dx[0];
            Real z = prob_lo[0] + (i+0.5)*dx[0];

            if (n == 0) {
                struct_fab(i,j,k,n) = 1.;
            }
            */


        });

    }


    amrex::Vector< int > s_pairA(3);
    amrex::Vector< int > s_pairB(3);

    // Select which variable pairs to include in structure factor:
    s_pairA[0] = 0;
    s_pairB[0] = 0;
    //
    s_pairA[1] = 0;
    s_pairB[1] = 1;
    //
    s_pairA[2] = 1;
    s_pairB[2] = 1;

    StructFact structFact(ba,dmap,var_names,var_scaling,s_pairA,s_pairB);

    /////////////////////////////////
    // take an FFT and write them out
    MultiFab dft_real(ba, dmap, 2, 0);
    MultiFab dft_imag(ba, dmap, 2, 0);
    structFact.ComputeFFT(struct_cc,dft_real,dft_imag);

    WriteSingleLevelPlotfile("plt_real", dft_real, {"var1", "var2"}, geom, 0., 0);
    WriteSingleLevelPlotfile("plt_imag", dft_imag, {"var1", "var2"}, geom, 0., 0);

    // magnitude
    MultiFab::Multiply(dft_real,dft_real,0,0,2,0);
    MultiFab::Multiply(dft_imag,dft_imag,0,0,2,0);
    MultiFab::Add(dft_real,dft_imag,0,0,2,0);
    SqrtMF(dft_real);
    WriteSingleLevelPlotfile("plt_mag", dft_real, {"var1", "var2"}, geom, 0., 0);
    /////////////////////////////////


    structFact.FortStructure(struct_cc);

    structFact.WritePlotFile(0,0.,"plt_SF");

    // Call the timer again and compute the maximum difference between the start time
    // and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;

}