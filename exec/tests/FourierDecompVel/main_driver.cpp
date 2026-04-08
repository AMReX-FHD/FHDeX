#include "common_functions.H"
#include "StructFact.H"

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_VisMF.H>

// These are for FFTW / cuFFT / rocFFT

#ifdef AMREX_USE_CUDA
#include <cufft.h>
#elif AMREX_USE_HIP
#  if __has_include(<rocfft/rocfft.h>)  // ROCm 5.3+
#    include <rocfft/rocfft.h>
#  else
#    include <rocfft.h>
#  endif
#else
#include <fftw3.h>
#include <fftw3-mpi.h>
#endif

#include <AMReX_GpuComplex.H>

#include <string>

#include <AMReX_MultiFabUtil.H>
#include "AMReX_PlotFileUtil.H"
#include "AMReX_BoxArray.H"
#include <AMReX_Vector.H>
#include <AMReX_MPMD.H>
#include <AMReX_VisMF.H>

#include "rng_functions.H"

#include "chrono"

#define ALIGN 16

using namespace std::chrono;
using namespace amrex;

void main_driver (const char* argv)
{

    // **********************************
    // DECLARE SIMULATION PARAMETERS
    // **********************************

    BL_PROFILE_VAR("main_driver()",main_driver);

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::string inputs_file = argv;

    amrex::AllPrint() << "Compiled with support for maximum species = " << MAX_SPECIES << "\n";

    // copy contents of F90 modules to C++ namespaces
    InitializeCommonNamespace();

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    DistributionMapping dmap;

    // define lower and upper indices of domain
    IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
    IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));

    // Make a single box that is the entire domain
    Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "domain"
    ba.define(domain);

    // create IntVect of max_grid_size
    // IntVect max_grid_size(AMREX_D_DECL(max_grid_size_x,max_grid_size_y,max_grid_size_z));

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(IntVect(max_grid_size));

    // How Boxes are distrubuted among MPI processes
    dmap.define(ba);

    // This defines the physical box size in each direction
    RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                     {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

    // periodic in all direction
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};

    // geometry object for real data
    geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());

    // extract dx from the geometry object
    GpuArray<Real,3> dx = geom.CellSizeArray();
    GpuArray<Real,AMREX_SPACEDIM> reallo;
    GpuArray<Real,AMREX_SPACEDIM> realhi;
    GpuArray<Real,AMREX_SPACEDIM> center;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        reallo[d] = prob_lo[d];
        realhi[d] = prob_hi[d];
        center[d] = ( realhi[d] - reallo[d] ) / 2.;
    }

    // Setup MultiFabs
    std::array< MultiFab, 3 > vel; // staggered velocity field
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        vel[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
        vel[d].setVal(0.0);
    }

    MultiFab structFactMFTurb; // storing CC velocities
    structFactMFTurb.define(ba, dmap, 3, 0);
    structFactMFTurb.setVal(0.0);

    // Fill MultiFab
    Real omega = M_PI/2.0;;
    //for (MFIter mfi(structFactMFTurb); mfi.isValid(); ++mfi) {

    //    const Array4<Real>& vel = structFactMFTurb.array(mfi);
    //    const Box& bx = mfi.tilebox();

    //    amrex::ParallelFor(bx,
    //    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    //    {
    //        if (prob_type == 0) {
    //            if ((i==3) and (j==3) and (k==3)) {
    //                vel(i+1,j,k,0) =  1.0;
    //                vel(i,j,k,0)   = -1.0;
    //                vel(i,j+1,k,1) =  1.0;
    //                vel(i,j,k,1)   = -1.0;
    //                vel(i,j,k+1,2) =  1.0;
    //                vel(i,j,k,2)   = -1.0;
    //            }
    //        }
    //        if (prob_type == 1) {
    //            if ((i==3) and (j==3) and (k==3)) {
    //                vel(i+1,j,k,0)     =  1.0;
    //                vel(i+1,j+1,k,0)   = -1.0;
    //                vel(i+1,j+1,k,1)   =  1.0;
    //                vel(i,j+1,k,1)     = -1.0;
    //            }
    //        }
    //        if (prob_type == 2) {
    //            if ((i==3) and (j==3) and (k==3)) {
    //                vel(i,j,k,0)     =   1.0;
    //                vel(i,j,k,1)     =  -1.0;
    //                vel(i+1,j,k,0)   =   1.0;
    //                vel(i+1,j,k,1)   =   1.0;
    //                vel(i+1,j+1,k,0) =  -1.0;
    //                vel(i+1,j+1,k,1) =   1.0;
    //                vel(i,j+1,k,0)   =  -1.0;
    //                vel(i,j+1,k,1)   =  -1.0;
    //            }
    //        }
    //
    //        Real x = (i+0.5)*dx[0] - center[0];
    //        Real y = (j+0.5)*dx[1] - center[1];
    //        Real z = (k+0.5)*dx[2] - center[2];
    //
    //        if (prob_type == 3) { // pure solenoidal
    //            vel(i,j,k,0) = -1.0*y;
    //            vel(i,j,k,1) =  1.0*x;
    //        }
    //
    //        if (prob_type == 4) { // pure dilatational
    //            vel(i,j,k,0) =  1.0*x;
    //            vel(i,j,k,1) =  1.0*y;
    //        }
    //
    //        if (prob_type == 5) { // pure Laplacian
    //            vel(i,j,k,0) =  1.0*x;
    //            vel(i,j,k,1) = -1.0*y;
    //        }
    //
    //        if (prob_type == 6) { // both
    //            vel(i,j,k,0) =  1.0*x - 1.0*y;
    //            vel(i,j,k,1) =  1.0*x + 1.0*y;
    //        }
    //
    //    });
    //}
    for (MFIter mfi(structFactMFTurb); mfi.isValid(); ++mfi) {

        AMREX_D_TERM(const Array4<Real>& velx = vel[0].array(mfi);,
                     const Array4<Real>& vely = vel[1].array(mfi);,
                     const Array4<Real>& velz = vel[2].array(mfi););

        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        amrex::ParallelFor(tbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real x = prob_lo[0] + (i) * dx[0]     - center[0];
            Real y = prob_lo[1] + (j+0.5) * dx[1] - center[1];
            Real z = prob_lo[2] + (k+0.5) * dx[2] - center[2];
            if (prob_type == 3) { // pure solenoidal
                velx(i,j,k) = -1.0*y;
            }

            if (prob_type == 4) { // pure dilatational
                velx(i,j,k) =  1.0*x;
            }

            if (prob_type == 5) { // pure Laplacian
                velx(i,j,k) =  1.0*x;
            }

            if (prob_type == 6) { // both
                velx(i,j,k) =  1.0*x - 1.0*y;
            }
        });
        amrex::ParallelFor(tby,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real x = prob_lo[0] + (i+0.5) * dx[0] - center[0];
            Real y = prob_lo[1] + (j) * dx[1]     - center[1];
            Real z = prob_lo[2] + (k+0.5) * dx[2] - center[2];
            if (prob_type == 3) { // pure solenoidal
                vely(i,j,k) =  1.0*x;
            }

            if (prob_type == 4) { // pure dilatational
                vely(i,j,k) =  1.0*y;
            }

            if (prob_type == 5) { // pure Laplacian
                vely(i,j,k) = -1.0*y;
            }

            if (prob_type == 6) { // both
                vely(i,j,k) =  1.0*x + 1.0*y;
            }
        });
        amrex::ParallelFor(tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real x = prob_lo[0] + (i+0.5) * dx[0] - center[0];
            Real y = prob_lo[1] + (j+0.5) * dx[1] - center[1];
            Real z = prob_lo[2] + (k) * dx[2]     - center[2];
        });
    }

    // Shift face-MF to CC-MF
    for(int d=0; d<3; d++) {
        ShiftFaceToCC(vel[d], 0, structFactMFTurb, d, 1);
    }

    Real time = 0.;
    int step = 0;
    WriteSingleLevelPlotfile("vel_full", structFactMFTurb, {"u","v","w"}, geom, time, step);

    /*******************************************/
    /********** Start FFT Calculations *********/
    StructFact turbStructFact;
    int structVarsTurb = AMREX_SPACEDIM;
    Vector< std::string > var_names_turb;
    var_names_turb.resize(structVarsTurb);
    int cnt = 0;
    std::string x;
    // velx, vely, velz
    for (int d=0; d<AMREX_SPACEDIM; d++) {
      x = "vel";
      x += (120+d);
      var_names_turb[cnt++] = x;
    }
    // need to use dVol for scaling
    Real dVol = (AMREX_SPACEDIM==2) ? dx[0]*dx[1]*cell_depth : dx[0]*dx[1]*dx[2];
    Real dProb = (AMREX_SPACEDIM==2) ? n_cells[0]*n_cells[1] : n_cells[0]*n_cells[1]*n_cells[2];
    dProb = 1./dProb;
    Vector<Real> var_scaling_turb(structVarsTurb);
    for (int d=0; d<var_scaling_turb.size(); ++d) {
        var_scaling_turb[d] = 1./dVol;
    }
    // option to compute only specified pairs
    amrex::Vector< int > s_pairA_turb(AMREX_SPACEDIM); // vx, vy, vz, rho, P , T
    amrex::Vector< int > s_pairB_turb(AMREX_SPACEDIM); // vx, vy, vz, rho, P , T

    // Select which variable pairs to include in structure factor:
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        s_pairA_turb[d] = d;
        s_pairB_turb[d] = d;
    }
    turbStructFact.define(ba,dmap,var_names_turb,var_scaling_turb,s_pairA_turb,s_pairB_turb);

    // Initialize the boxarray "ba_onegrid" from the single box "domain"
    long npts;
    BoxArray ba_onegrid;
    {
      Box domain = geom.Domain();
      ba_onegrid.define(domain);
      npts = (domain.length(0)*domain.length(1)*domain.length(2));
    }
    Real sqrtnpts = std::sqrt(npts);
    DistributionMapping dmap_onegrid(ba_onegrid);

    MultiFab dft_real, dft_imag;
    dft_real.define(ba_onegrid, dmap_onegrid, 3, 0);
    dft_imag.define(ba_onegrid, dmap_onegrid, 3, 0);

    // Forward FFT
    turbStructFact.ComputeFFT(structFactMFTurb,dft_real,dft_imag,geom,false);

    // Decompose Velocity into Solenoidal and Dilatational Components
    MultiFab dft_decomp_real, dft_decomp_imag;
    dft_decomp_real.define(ba_onegrid, dmap_onegrid, 6, 0);
    dft_decomp_imag.define(ba_onegrid, dmap_onegrid, 6, 0);

    for (MFIter mfi(dft_decomp_real); mfi.isValid(); ++mfi) {

        Box bx = mfi.fabbox();

        Array4<const Real> const& real = dft_real.array(mfi);
        Array4<const Real> const& imag = dft_imag.array(mfi);

        Array4<      Real> const& real_decomp = dft_decomp_real.array(mfi);
        Array4<      Real> const& imag_decomp = dft_decomp_imag.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int nx = bx.length(0);
            int ny = bx.length(1);
            int nz = bx.length(2);

            if (i <= bx.length(0)/2) { // only need to do for the first half of k-space

                // Gradient Operators
                Real GxR = (cos(2.0*M_PI*i/nx)-1.0)/dx[0];
                Real GxC = (sin(2.0*M_PI*i/nx)-0.0)/dx[0];
                Real GyR = (cos(2.0*M_PI*j/ny)-1.0)/dx[1];
                Real GyC = (sin(2.0*M_PI*j/ny)-0.0)/dx[1];
                Real GzR = (cos(2.0*M_PI*k/nz)-1.0)/dx[2];
                Real GzC = (sin(2.0*M_PI*k/nz)-0.0)/dx[2];

                // Inverse Laplacian
                Real Lap = GxR*GxR + GxC*GxC + GyR*GyR + GyC*GyC + GzR*GzR + GzC*GzC;

                // Divergence of vel
                Real divR = real(i,j,k,0)*GxR - imag(i,j,k,0)*GxC +
                            real(i,j,k,1)*GyR - imag(i,j,k,1)*GyC +
                            real(i,j,k,2)*GzR - imag(i,j,k,2)*GzC ;
                Real divC = real(i,j,k,0)*GxC + imag(i,j,k,0)*GxR +
                            real(i,j,k,1)*GyC + imag(i,j,k,1)*GyR +
                            real(i,j,k,2)*GzC + imag(i,j,k,2)*GzR ;

                if (Lap < 1.0e-12) { // zero mode for no bulk motion
                    real_decomp(i,j,k,0) = 0.0;
                    real_decomp(i,j,k,1) = 0.0;
                    real_decomp(i,j,k,2) = 0.0;
                    imag_decomp(i,j,k,0) = 0.0;
                    imag_decomp(i,j,k,1) = 0.0;
                    imag_decomp(i,j,k,2) = 0.0;
                }
                else {
                    // Dilatational velocity
                    real_decomp(i,j,k,0) = (divR*GxR + divC*GxC) / Lap;
                    real_decomp(i,j,k,1) = (divR*GyR + divC*GyC) / Lap;
                    real_decomp(i,j,k,2) = (divR*GzR + divC*GzC) / Lap;
                    imag_decomp(i,j,k,0) = (divC*GxR - divR*GxC) / Lap;
                    imag_decomp(i,j,k,1) = (divC*GyR - divR*GyC) / Lap;
                    imag_decomp(i,j,k,2) = (divC*GzR - divR*GzC) / Lap;

                    // Solenoidal velocity
                    real_decomp(i,j,k,3) = real(i,j,k,0) - real_decomp(i,j,k,0);
                    real_decomp(i,j,k,4) = real(i,j,k,1) - real_decomp(i,j,k,1);
                    real_decomp(i,j,k,5) = real(i,j,k,2) - real_decomp(i,j,k,2);
                    imag_decomp(i,j,k,3) = imag(i,j,k,0) - imag_decomp(i,j,k,0);
                    imag_decomp(i,j,k,4) = imag(i,j,k,1) - imag_decomp(i,j,k,1);
                    imag_decomp(i,j,k,5) = imag(i,j,k,2) - imag_decomp(i,j,k,2);
                }
            }
        });
    }

    // Backward FFT
    MultiFab vel_decomp;
    vel_decomp.define(ba, dmap, 6, 0);
    turbStructFact.InverseFFT(vel_decomp,dft_decomp_real,dft_decomp_imag,geom);

    // Write Output
    WriteSingleLevelPlotfile("vel_decomp", vel_decomp, {"uD","vD","wD","uS","vS","wS"}, geom, time, step);

}
