#include "spectral_functions.H"
#include "AMReX_PlotFileUtil.H"
#include "AMReX_PlotFileDataImpl.H"

#include <sys/stat.h>

#include "chrono"
#include <AMReX_MultiFabUtil.H>
#include "AMReX_PlotFileUtil.H"
#include "AMReX_BoxArray.H"

using namespace std::chrono;

namespace {
    void GotoNextLine (std::istream& is)
    {
        constexpr std::streamsize bl_ignore_max { 100000 };
        is.ignore(bl_ignore_max, '\n');
    }
}

void SpectralReadCheckPoint(amrex::Geometry& geom,
                            const amrex::Box& domain,
                            amrex::MultiFab& prim,
                            std::array<MultiFab, 3>& vel,
                            BoxArray& ba, DistributionMapping& dmap,
                            const amrex::Vector<int> n_cells,
                            const int nprimvars,
                            const amrex::Vector<int> max_grid_size,
                            const amrex::IntVect ngc,
                            const int restart)
{
    // timer for profiling
    BL_PROFILE_VAR("SpectralReadCheckPoint()",SpectralReadCheckPoint);

    // checkpoint file name, e.g., chk0000010
    const std::string& checkpointname = amrex::Concatenate("chk",restart,9);

    amrex::Print() << "Restart from checkpoint " << checkpointname << "\n";

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    std::string line, word;

    // read in old boxarray, and create old distribution map (this is to read in MFabs)
    BoxArray ba_old;
    DistributionMapping dmap_old;

    // initialize new boxarray
    ba.define(domain);
    ba.maxSize(IntVect(max_grid_size));
    dmap.define(ba, ParallelDescriptor::NProcs());

    amrex::Vector<amrex::IntVect> nodal_flag_dir;
    amrex::IntVect                nodal_flag_x;
    amrex::IntVect                nodal_flag_y;
    amrex::IntVect                nodal_flag_z;
    nodal_flag_dir.resize(3);

    for (int i=0; i<3; ++i) {
        nodal_flag_x[i] = int(i==0);
        nodal_flag_y[i] = int(i==1);
        nodal_flag_z[i] = int(i==2);
        AMREX_D_TERM(nodal_flag_dir[0][i] = nodal_flag_x[i];,
                     nodal_flag_dir[1][i] = nodal_flag_y[i];,
                     nodal_flag_dir[2][i] = nodal_flag_z[i];);
    }

    // Header
    {
        std::string File(checkpointname + "/Header");
        Vector<char> fileCharPtr;
        ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
        std::string fileCharPtrString(fileCharPtr.dataPtr());
        std::istringstream is(fileCharPtrString, std::istringstream::in);

        // read in title line
        std::getline(is, line);

        // read in time step number
        int step;
        is >> step;
        GotoNextLine(is);

        // read in time
        Real time;
        is >> time;
        GotoNextLine(is);

        // read in statsCount
        int statsCount;
        is >> statsCount;
        GotoNextLine(is);

        // read in BoxArray (fluid) from Header
        ba_old.readFrom(is);
        GotoNextLine(is);

        // create old distribution mapping
        dmap_old.define(ba_old, ParallelDescriptor::NProcs());

        prim.define(ba,dmap,nprimvars,ngc);
        // velocity and momentum (instantaneous, means, variances)
        for (int d=0; d<3; d++) {
            vel[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ngc);
        }
    }

    // C++ random number engine
    // each MPI process reads in its own file
    int comm_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    int n_ranks;
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

    // read in the MultiFab data
    Read_Copy_MF_Checkpoint(prim,"prim",checkpointname,ba_old,dmap_old,nprimvars,1,ngc);

    Read_Copy_MF_Checkpoint(vel[0],"velx",checkpointname,ba_old,dmap_old,1,1,ngc,0);
    Read_Copy_MF_Checkpoint(vel[1],"vely",checkpointname,ba_old,dmap_old,1,1,ngc,1);
    Read_Copy_MF_Checkpoint(vel[2],"velz",checkpointname,ba_old,dmap_old,1,1,ngc,2);

    // FillBoundaries
    prim.FillBoundary(geom.periodicity());
    vel[0].FillBoundary(geom.periodicity());
    vel[1].FillBoundary(geom.periodicity());
    vel[2].FillBoundary(geom.periodicity());
}

void SpectralVelDecomp(const MultiFab& vel,
                       MultiFab& vel_decomp_filter,
                       const amrex::Real kmin,
                       const amrex::Real kmax,
                       const amrex::Geometry& geom,
                       const amrex::Vector<int> n_cells)
{
    BL_PROFILE_VAR("SpectralVelDecomp()",SpectralVelDecomp);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vel.nComp() == 3,
        "SpectralVelDecomp: must have 3 components of input vel MultiFab");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vel.local_size() == 1,
        "SpectralVelDecomp: Must have one Box per MPI process when using heFFTe");

    const GpuArray<Real, 3> dx = geom.CellSizeArray();

    long npts;
    Box domain = geom.Domain();
    npts = (domain.length(0)*domain.length(1)*domain.length(2));
    Real sqrtnpts = std::sqrt(npts);

    // get box array and distribution map of vel
    DistributionMapping dm = vel.DistributionMap();
    BoxArray ba            = vel.boxArray();

    // since there is 1 MPI rank per box, each MPI rank obtains its local box and the associated boxid
    Box local_box;
    int local_boxid;
    {
        for (int i = 0; i < ba.size(); ++i) {
            Box b = ba[i];
            // each MPI rank has its own local_box Box and local_boxid ID
            if (ParallelDescriptor::MyProc() == dm[i]) {
                local_box = b;
                local_boxid = i;
            }
        }
    }

    // now each MPI rank works on its own box
    // for real->complex fft's, the fft is stored in an (nx/2+1) x ny x nz dataset

    // start by coarsening each box by 2 in the x-direction
    Box c_local_box = amrex::coarsen(local_box, IntVect(AMREX_D_DECL(2,1,1)));

    // if the coarsened box's high-x index is even, we shrink the size in 1 in x
    // this avoids overlap between coarsened boxes
    if (c_local_box.bigEnd(0) * 2 == local_box.bigEnd(0)) {
        c_local_box.setBig(0,c_local_box.bigEnd(0)-1);
    }
    // for any boxes that touch the hi-x domain we
    // increase the size of boxes by 1 in x
    // this makes the overall fft dataset have size (Nx/2+1 x Ny x Nz)
    if (local_box.bigEnd(0) == geom.Domain().bigEnd(0)) {
        c_local_box.growHi(0,1);
    }

    // each MPI rank gets storage for its piece of the fft
    BaseFab<GpuComplex<Real> > spectral_field_Tx(c_local_box, 1, The_Device_Arena()); // totalx
    BaseFab<GpuComplex<Real> > spectral_field_Ty(c_local_box, 1, The_Device_Arena()); // totaly
    BaseFab<GpuComplex<Real> > spectral_field_Tz(c_local_box, 1, The_Device_Arena()); // totalz
    BaseFab<GpuComplex<Real> > spectral_field_Sx(c_local_box, 1, The_Device_Arena()); // solenoidalx
    BaseFab<GpuComplex<Real> > spectral_field_Sy(c_local_box, 1, The_Device_Arena()); // solenoidaly
    BaseFab<GpuComplex<Real> > spectral_field_Sz(c_local_box, 1, The_Device_Arena()); // solenoidalz
    BaseFab<GpuComplex<Real> > spectral_field_Dx(c_local_box, 1, The_Device_Arena()); // dilatationalx
    BaseFab<GpuComplex<Real> > spectral_field_Dy(c_local_box, 1, The_Device_Arena()); // dilatationaly
    BaseFab<GpuComplex<Real> > spectral_field_Dz(c_local_box, 1, The_Device_Arena()); // dilatationalz
    MultiFab vel_single(ba, dm, 1, 0);

    int r2c_direction = 0;

    // ForwardTransform
    // X
    using heffte_complex = typename heffte::fft_output<Real>::type;
    {
#if defined(HEFFTE_CUFFT)
        heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
        heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
        heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
        {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
        {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
        {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
        r2c_direction, ParallelDescriptor::Communicator());

        vel_single.ParallelCopy(vel, 0, 0, 1);
        heffte_complex* spectral_data = (heffte_complex*) spectral_field_Tx.dataPtr();
        fft.forward(vel_single[local_boxid].dataPtr(),spectral_data);
    }
    // Y
    {
#if defined(HEFFTE_CUFFT)
        heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
        heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
        heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
        {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
        {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
        {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
        r2c_direction, ParallelDescriptor::Communicator());

        vel_single.ParallelCopy(vel, 1, 0, 1);
        heffte_complex* spectral_data = (heffte_complex*) spectral_field_Ty.dataPtr();
        fft.forward(vel_single[local_boxid].dataPtr(),spectral_data);
    }
    // Z
    {
#if defined(HEFFTE_CUFFT)
        heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
        heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
        heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
        {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
        {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
        {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
        r2c_direction, ParallelDescriptor::Communicator());

        vel_single.ParallelCopy(vel, 2, 0, 1);
        heffte_complex* spectral_data = (heffte_complex*) spectral_field_Tz.dataPtr();
        fft.forward(vel_single[local_boxid].dataPtr(),spectral_data);
    }

    Gpu::streamSynchronize();

    int nx = n_cells[0];
    int ny = n_cells[1];
    int nz = n_cells[2];

    // Decompose velocity field into solenoidal and dilatational
    Array4< GpuComplex<Real> > spectral_tx = spectral_field_Tx.array();
    Array4< GpuComplex<Real> > spectral_ty = spectral_field_Ty.array();
    Array4< GpuComplex<Real> > spectral_tz = spectral_field_Tz.array();
    Array4< GpuComplex<Real> > spectral_sx = spectral_field_Sx.array();
    Array4< GpuComplex<Real> > spectral_sy = spectral_field_Sy.array();
    Array4< GpuComplex<Real> > spectral_sz = spectral_field_Sz.array();
    Array4< GpuComplex<Real> > spectral_dx = spectral_field_Dx.array();
    Array4< GpuComplex<Real> > spectral_dy = spectral_field_Dy.array();
    Array4< GpuComplex<Real> > spectral_dz = spectral_field_Dz.array();
    ParallelFor(c_local_box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {

        Real GxR = 0.0, GxC = 0.0, GyR = 0.0, GyC = 0.0, GzR = 0.0, GzC = 0.0;

        if (i <= nx/2) {

            // Get the wavevector
            int ki = i;
            int kj = j;
            if (j >= ny/2) kj = ny - j;
            int kk = k;
            if (k >= nz/2) kk = nz - k;

            // Gradient Operators
            GxR = (cos(2.0*M_PI*ki/nx)-1.0)/dx[0];
            GxC = (sin(2.0*M_PI*ki/nx)-0.0)/dx[0];
            GyR = (cos(2.0*M_PI*kj/ny)-1.0)/dx[1];
            GyC = (sin(2.0*M_PI*kj/ny)-0.0)/dx[1];
            GzR = (cos(2.0*M_PI*kk/nz)-1.0)/dx[2];
            GzC = (sin(2.0*M_PI*kk/nz)-0.0)/dx[2];
        }
        else { // conjugate
            amrex::Abort("check the code; i should not go beyond bx.length(0)/2");
        }

        // Get the wavenumber
        int ki = i;
        int kj = j;
        if (j >= ny/2) kj = ny - j;
        int kk = k;
        if (k >= nz/2) kk = nz - k;
        Real knum = (ki*ki + kj*kj + kk*kk);
        knum = std::sqrt(knum);

        // Scale Total velocity FFT components with Filtering
        if ((knum >= kmin) and (knum <= kmax)) {

            spectral_tx(i,j,k) *= (1.0/sqrtnpts);
            spectral_ty(i,j,k) *= (1.0/sqrtnpts);
            spectral_tz(i,j,k) *= (1.0/sqrtnpts);

            // Inverse Laplacian
            Real Lap = GxR*GxR + GxC*GxC + GyR*GyR + GyC*GyC + GzR*GzR + GzC*GzC;

            // Divergence of vel
            Real divR = spectral_tx(i,j,k).real()*GxR - spectral_tx(i,j,k).imag()*GxC +
                       spectral_ty(i,j,k).real()*GyR - spectral_ty(i,j,k).imag()*GyC +
                       spectral_tz(i,j,k).real()*GzR - spectral_tz(i,j,k).imag()*GzC ;
            Real divC = spectral_tx(i,j,k).real()*GxC + spectral_tx(i,j,k).imag()*GxR +
                       spectral_ty(i,j,k).real()*GyC + spectral_ty(i,j,k).imag()*GyR +
                       spectral_tz(i,j,k).real()*GzC + spectral_tz(i,j,k).imag()*GzR ;

            if (Lap < 1.0e-12) { // zero mode for no bulk motion
                spectral_dx(i,j,k) *= 0.0;
                spectral_dy(i,j,k) *= 0.0;
                spectral_dz(i,j,k) *= 0.0;
            }
            else {

                // Dilatational velocity
                GpuComplex<Real> copy_dx((divR*GxR + divC*GxC) / Lap,
                                         (divC*GxR - divR*GxC) / Lap);
                spectral_dx(i,j,k) = copy_dx;

                GpuComplex<Real> copy_dy((divR*GyR + divC*GyC) / Lap,
                                         (divC*GyR - divR*GyC) / Lap);
                spectral_dy(i,j,k) = copy_dy;

                GpuComplex<Real> copy_dz((divR*GzR + divC*GzC) / Lap,
                                         (divC*GzR - divR*GzC) / Lap);
                spectral_dz(i,j,k) = copy_dz;
            }

            // Solenoidal velocity
            spectral_sx(i,j,k) = spectral_tx(i,j,k) - spectral_dx(i,j,k);
            spectral_sy(i,j,k) = spectral_ty(i,j,k) - spectral_dy(i,j,k);
            spectral_sz(i,j,k) = spectral_tz(i,j,k) - spectral_dz(i,j,k);
        }
        else {
            spectral_tx(i,j,k) = 0.0;
            spectral_ty(i,j,k) = 0.0;
            spectral_tz(i,j,k) = 0.0;
            spectral_sx(i,j,k) = 0.0;
            spectral_sy(i,j,k) = 0.0;
            spectral_sz(i,j,k) = 0.0;
            spectral_dx(i,j,k) = 0.0;
            spectral_dy(i,j,k) = 0.0;
            spectral_dz(i,j,k) = 0.0;
        }

    });

    Gpu::streamSynchronize();

    MultiFab vel_decomp_filter_single(ba, dm, 1, 0);
    // inverse Fourier transform filtered total velocity
    {
#if defined(HEFFTE_CUFFT)
        heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
        heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
        heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
        {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
        {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
        {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
        r2c_direction, ParallelDescriptor::Communicator());

        heffte_complex* spectral_data = (heffte_complex*) spectral_field_Tx.dataPtr();
        fft.backward(spectral_data, vel_decomp_filter_single[local_boxid].dataPtr());

        Gpu::streamSynchronize();
        vel_decomp_filter.ParallelCopy(vel_decomp_filter_single, 0, 0, 1);
    }
    {
#if defined(HEFFTE_CUFFT)
        heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
        heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
        heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
        {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
        {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
        {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
        r2c_direction, ParallelDescriptor::Communicator());

        heffte_complex* spectral_data = (heffte_complex*) spectral_field_Ty.dataPtr();
        fft.backward(spectral_data, vel_decomp_filter_single[local_boxid].dataPtr());

        Gpu::streamSynchronize();
        vel_decomp_filter.ParallelCopy(vel_decomp_filter_single, 0, 1, 1);
    }
    {
#if defined(HEFFTE_CUFFT)
        heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
        heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
        heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
        {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
        {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
        {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
        r2c_direction, ParallelDescriptor::Communicator());

        heffte_complex* spectral_data = (heffte_complex*) spectral_field_Tz.dataPtr();
        fft.backward(spectral_data, vel_decomp_filter_single[local_boxid].dataPtr());

        Gpu::streamSynchronize();
        vel_decomp_filter.ParallelCopy(vel_decomp_filter_single, 0, 2, 1);
    }
    // inverse Fourier transform filtered solenoidal and dilatational components
    {
#if defined(HEFFTE_CUFFT)
        heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
        heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
        heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
        {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
        {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
        {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
        r2c_direction, ParallelDescriptor::Communicator());

        heffte_complex* spectral_data = (heffte_complex*) spectral_field_Sx.dataPtr();
        fft.backward(spectral_data, vel_decomp_filter_single[local_boxid].dataPtr());

        Gpu::streamSynchronize();
        vel_decomp_filter.ParallelCopy(vel_decomp_filter_single, 0, 3, 1);
    }
    {
#if defined(HEFFTE_CUFFT)
        heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
        heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
        heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
        {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
        {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
        {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
        r2c_direction, ParallelDescriptor::Communicator());

        heffte_complex* spectral_data = (heffte_complex*) spectral_field_Sy.dataPtr();
        fft.backward(spectral_data, vel_decomp_filter_single[local_boxid].dataPtr());

        Gpu::streamSynchronize();
        vel_decomp_filter.ParallelCopy(vel_decomp_filter_single, 0, 4, 1);
    }
    {
#if defined(HEFFTE_CUFFT)
        heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
        heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
        heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
        {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
        {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
        {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
        r2c_direction, ParallelDescriptor::Communicator());

        heffte_complex* spectral_data = (heffte_complex*) spectral_field_Sz.dataPtr();
        fft.backward(spectral_data, vel_decomp_filter_single[local_boxid].dataPtr());

        Gpu::streamSynchronize();
        vel_decomp_filter.ParallelCopy(vel_decomp_filter_single, 0, 5, 1);
    }
    {
#if defined(HEFFTE_CUFFT)
        heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
        heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
        heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
        {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
        {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
        {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
        r2c_direction, ParallelDescriptor::Communicator());

        heffte_complex* spectral_data = (heffte_complex*) spectral_field_Dx.dataPtr();
        fft.backward(spectral_data, vel_decomp_filter_single[local_boxid].dataPtr());

        Gpu::streamSynchronize();
        vel_decomp_filter.ParallelCopy(vel_decomp_filter_single, 0, 6, 1);
    }
    {
#if defined(HEFFTE_CUFFT)
        heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
        heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
        heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
        {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
        {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
        {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
        r2c_direction, ParallelDescriptor::Communicator());

        heffte_complex* spectral_data = (heffte_complex*) spectral_field_Dy.dataPtr();
        fft.backward(spectral_data, vel_decomp_filter_single[local_boxid].dataPtr());

        Gpu::streamSynchronize();
        vel_decomp_filter.ParallelCopy(vel_decomp_filter_single, 0, 7, 1);
    }
    {
#if defined(HEFFTE_CUFFT)
        heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
        heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
        heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
        {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
        {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
        {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
        r2c_direction, ParallelDescriptor::Communicator());

        heffte_complex* spectral_data = (heffte_complex*) spectral_field_Dz.dataPtr();
        fft.backward(spectral_data, vel_decomp_filter_single[local_boxid].dataPtr());

        Gpu::streamSynchronize();
        vel_decomp_filter.ParallelCopy(vel_decomp_filter_single, 0, 8, 1);
    }


    vel_decomp_filter.mult(1.0/sqrtnpts);

}


void SpectralScalarDecomp(const MultiFab& scalar,
                          MultiFab& scalar_filter,
                          const amrex::Real kmin,
                          const amrex::Real kmax,
                          const amrex::Geometry& geom,
                          const amrex::Vector<int> n_cells)
{
    BL_PROFILE_VAR("SpectralScalarDecomp()",SpectralScalarDecomp);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(scalar.nComp() == 1,
        "SpectralScalarDecomp: must have 1 components of input scalar MultiFab");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(scalar.local_size() == 1,
        "SpectralScalarDecomp: Must have one Box per MPI process when using heFFTe");

    const GpuArray<Real, 3> dx = geom.CellSizeArray();

    long npts;
    Box domain = geom.Domain();
    npts = (domain.length(0)*domain.length(1)*domain.length(2));
    Real sqrtnpts = std::sqrt(npts);

    // get box array and distribution map of vel
    DistributionMapping dm = scalar.DistributionMap();
    BoxArray ba            = scalar.boxArray();

    // since there is 1 MPI rank per box, each MPI rank obtains its local box and the associated boxid
    Box local_box;
    int local_boxid;
    {
        for (int i = 0; i < ba.size(); ++i) {
            Box b = ba[i];
            // each MPI rank has its own local_box Box and local_boxid ID
            if (ParallelDescriptor::MyProc() == dm[i]) {
                local_box = b;
                local_boxid = i;
            }
        }
    }

    // now each MPI rank works on its own box
    // for real->complex fft's, the fft is stored in an (nx/2+1) x ny x nz dataset

    // start by coarsening each box by 2 in the x-direction
    Box c_local_box = amrex::coarsen(local_box, IntVect(AMREX_D_DECL(2,1,1)));

    // if the coarsened box's high-x index is even, we shrink the size in 1 in x
    // this avoids overlap between coarsened boxes
    if (c_local_box.bigEnd(0) * 2 == local_box.bigEnd(0)) {
        c_local_box.setBig(0,c_local_box.bigEnd(0)-1);
    }
    // for any boxes that touch the hi-x domain we
    // increase the size of boxes by 1 in x
    // this makes the overall fft dataset have size (Nx/2+1 x Ny x Nz)
    if (local_box.bigEnd(0) == geom.Domain().bigEnd(0)) {
        c_local_box.growHi(0,1);
    }

    // each MPI rank gets storage for its piece of the fft
    BaseFab<GpuComplex<Real> > spectral_field(c_local_box, 1, The_Device_Arena());
    MultiFab scalar_single(ba, dm, 1, 0);

    int r2c_direction = 0;

    // ForwardTransform
    using heffte_complex = typename heffte::fft_output<Real>::type;
    {
#if defined(HEFFTE_CUFFT)
        heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
        heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
        heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
        {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
        {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
        {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
        r2c_direction, ParallelDescriptor::Communicator());

        scalar_single.ParallelCopy(scalar, 0, 0, 1);
        heffte_complex* spectral_data = (heffte_complex*) spectral_field.dataPtr();
        fft.forward(scalar_single[local_boxid].dataPtr(),spectral_data);
    }

    Gpu::streamSynchronize();

    // filtering
    Array4< GpuComplex<Real> > spectral = spectral_field.array();
    int nx = n_cells[0];
    int ny = n_cells[1];
    int nz = n_cells[2];
    ParallelFor(c_local_box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {

        if (i <= nx/2) {
        }
        else { // conjugate
            amrex::Abort("check the code; i should not go beyond bx.length(0)/2");
        }

        // Get the wavenumber
        int ki = i;
        int kj = j;
        if (j >= ny/2) kj = ny - j;
        int kk = k;
        if (k >= nz/2) kk = nz - k;
        Real knum = (ki*ki + kj*kj + kk*kk);
        knum = std::sqrt(knum);

        // Scale Scalar FFT components with Filtering
        if ((knum >= kmin) and (knum <= kmax)) {
            spectral(i,j,k) *= (1.0/sqrtnpts);
            spectral(i,j,k) *= (1.0/sqrtnpts);
            spectral(i,j,k) *= (1.0/sqrtnpts);
        }
        else {
            spectral(i,j,k) *= 0.0;
            spectral(i,j,k) *= 0.0;
            spectral(i,j,k) *= 0.0;
        }
    });

    Gpu::streamSynchronize();

    MultiFab scalar_filter_single(ba, dm, 1, 0);
    // inverse Fourier transform filtered scalar
    {
#if defined(HEFFTE_CUFFT)
        heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif defined(HEFFTE_ROCFFT)
        heffte::fft3d_r2c<heffte::backend::rocfft> fft
#elif defined(HEFFTE_FFTW)
        heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
        {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
        {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
        {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
        r2c_direction, ParallelDescriptor::Communicator());

        heffte_complex* spectral_data = (heffte_complex*) spectral_field.dataPtr();
        fft.backward(spectral_data, scalar_filter_single[local_boxid].dataPtr());

        Gpu::streamSynchronize();
        scalar_filter.ParallelCopy(scalar_filter_single, 0, 0, 1);
    }

    scalar_filter.mult(1.0/sqrtnpts);

}

void SpectralWritePlotFile(const int step,
                           const amrex::Real& kmin,
                           const amrex::Real& kmax,
                           const amrex::Geometry& geom,
                           const amrex::MultiFab& vel_decomp_in,
                           const amrex::MultiFab& scalar_in,
                           const amrex::MultiFab& vel_total,
                           const amrex::MultiFab& scalar_total)
{

    MultiFab output;

    // Cell-Centered Velocity Gradient Stats (1,2,3 are directions)
    // 0: ux
    // 1: uy
    // 2: uz
    // 3: ux_s
    // 4: uy_s
    // 5: uz_s
    // 6: ux_d
    // 7: uy_d
    // 8: uz_d
    // 9: umag
    // 10: umag_s
    // 11: umag_d
    // 12: scalar
    // 13: divergence = u_1,1 + u_2,2 + u_3,3
    // 14: vorticity w1
    // 15: vorticity w2
    // 16: vorticity w3
    // 17: vorticity mag: sqrt(w1**2 + w2**2 + w3**2)
    // 18: ux_org
    // 19: scalar_org
    output.define(vel_decomp_in.boxArray(), vel_decomp_in.DistributionMap(), 20, 0);
    output.setVal(0.0);

    const GpuArray<Real, 3> dx = geom.CellSizeArray();

    for ( MFIter mfi(output,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.tilebox();

        const Array4<      Real>&             out   = output.array(mfi);

        const Array4<const Real>&  v_decomp         = vel_decomp_in.array(mfi);

        const Array4<const Real>&  sca              = scalar_in.array(mfi);

        const Array4<const Real>&  v_tot            = vel_total.array(mfi);

        const Array4<const Real>&  sca_tot          = scalar_total.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            out(i,j,k,0) = v_decomp(i,j,k,0);
            out(i,j,k,1) = v_decomp(i,j,k,1);
            out(i,j,k,2) = v_decomp(i,j,k,2);
            out(i,j,k,3) = v_decomp(i,j,k,3);
            out(i,j,k,4) = v_decomp(i,j,k,4);
            out(i,j,k,5) = v_decomp(i,j,k,5);
            out(i,j,k,6) = v_decomp(i,j,k,6);
            out(i,j,k,7) = v_decomp(i,j,k,7);
            out(i,j,k,8) = v_decomp(i,j,k,8);

            out(i,j,k,9)  = std::sqrt(out(i,j,k,0)*out(i,j,k,0) + out(i,j,k,1)*out(i,j,k,1) + out(i,j,k,2)*out(i,j,k,2)); // mag
            out(i,j,k,10) = std::sqrt(out(i,j,k,3)*out(i,j,k,3) + out(i,j,k,4)*out(i,j,k,4) + out(i,j,k,5)*out(i,j,k,5)); // mag solednoidal
            out(i,j,k,11) = std::sqrt(out(i,j,k,6)*out(i,j,k,6) + out(i,j,k,7)*out(i,j,k,7) + out(i,j,k,8)*out(i,j,k,8)); // mag solednoidal

            out(i,j,k,12) = sca(i,j,k,0);

            // divergence
            out(i,j,k,13) = 0.5*( (v_decomp(i+1,j,k,0) - v_decomp(i-1,j,k,0))/dx[0] +
                                  (v_decomp(i,j+1,k,1) - v_decomp(i,j-1,k,1))/dx[1] +
                                  (v_decomp(i,j,k+1,2) - v_decomp(i,j,k-1,2))/dx[2] );

            // curl w1 = u_2,1 - u_1,2
            out(i,j,k,14) = 0.5*( (v_decomp(i+1,j,k,1) - v_decomp(i-1,j,k,1))/dx[0] -
                                  (v_decomp(i,j+1,k,0) - v_decomp(i,j-1,k,0))/dx[1] );

            // curl w2 = u_1,3 - u_3,1
            out(i,j,k,15) = 0.5*( (v_decomp(i,j,k+1,0) - v_decomp(i,j,k-1,0))/dx[2] -
                                  (v_decomp(i+1,j,k,2) - v_decomp(i-1,j,k,2))/dx[0] );

            // curl w2 = u_3,2 - u_2,3
            out(i,j,k,16) = 0.5*( (v_decomp(i,j+1,k,2) - v_decomp(i,j-1,k,2))/dx[1] -
                                  (v_decomp(i,j,k+1,1) - v_decomp(i,j,k-1,1))/dx[2] );

            // vorticity magnitude: sqrt(w1*w1 + w2*w2 + w3*w3)
            out(i,j,k,17) = std::sqrt( out(i,j,k,14)*out(i,j,k,14) + out(i,j,k,15)*out(i,j,k,15) + out(i,j,k,16)*out(i,j,k,16) );

            // original velx
            out(i,j,k,18) = v_tot(i,j,k,0);

            // original scalar
            out(i,j,k,19) = sca_tot(i,j,k,0);
        });
    }

    // Write on a plotfile
    std::string plotfilename = amrex::Concatenate("filtered_",step,9);
    std::ostringstream os;
    os << std::setprecision(3) << kmin;
    plotfilename += os.str();;
    plotfilename += "_";
    std::ostringstream oss;
    oss << std::setprecision(3) << kmax;
    plotfilename += oss.str();

    amrex::Vector<std::string> varNames(20);
    varNames[0] = "ux";
    varNames[1] = "uy";
    varNames[2] = "uz";
    varNames[3] = "ux_s";
    varNames[4] = "uy_s";
    varNames[5] = "uz_s";
    varNames[6] = "ux_d";
    varNames[7] = "uy_d";
    varNames[8] = "uz_d";
    varNames[9] = "umag";
    varNames[10] = "umag_s";
    varNames[11] = "umag_d";
    varNames[12] = "rho";
    varNames[13] = "div";
    varNames[14] = "w1";
    varNames[15] = "w2";
    varNames[16] = "w3";
    varNames[17] = "vort";
    varNames[18] = "ux_org";
    varNames[19] = "rho_org";
    WriteSingleLevelPlotfile(plotfilename,output,varNames,geom,0.0,step);
}

void Read_Copy_MF_Checkpoint(amrex::MultiFab& mf, std::string mf_name, const std::string& checkpointname,
                             BoxArray& ba_old, DistributionMapping& dmap_old,
                             int NVARS, int ghost, const amrex::IntVect ngc,
                             int nodal_flag)
{
    // Read into temporary MF from file
    MultiFab mf_temp;
    VisMF::Read(mf_temp,amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", mf_name));

    // Copy temporary MF into the new MF
    if (ghost) {
        mf.ParallelCopy(mf_temp, 0, 0, NVARS, ngc, ngc);
    }
    else {
        mf.ParallelCopy(mf_temp, 0, 0, NVARS, 0, 0);
    }
}

void ShiftFaceToCC(const MultiFab& face_in, int face_comp,
                   MultiFab& cc_in, int cc_comp, int ncomp)
{

    BL_PROFILE_VAR("ShiftFaceToCC()",ShiftFaceToCC);

    if (!face_in.is_nodal(0) && !face_in.is_nodal(1) && !face_in.is_nodal(2)) {
        Abort("ShiftFaceToCC requires a face-centered MultiFab");
    }

    // Loop over boxes (note that mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(cc_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        Array4<Real const> const& face = face_in.array(mfi);

        Array4<Real> const& cc = cc_in.array(mfi);

        amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            cc(i,j,k,cc_comp+n) = face(i,j,k,face_comp+n);
        });
    }
}

void ComputeGrad(const MultiFab & phi_in, std::array<MultiFab, 3> & gphi,
                 int start_incomp, int start_outcomp, int ncomp, int bccomp, const Geometry & geom,
                 int increment)
{
    BL_PROFILE_VAR("ComputeGrad()",ComputeGrad);

    // Physical Domain
    Box dom(geom.Domain());

    const GpuArray<Real, 3> dx = geom.CellSizeArray();

    // if not incrementing, initialize data to zero
    if (increment == 0) {
        for (int dir=0; dir<3; ++dir) {
            gphi[dir].setVal(0.,start_outcomp,ncomp,0);
        }
    }

    // Loop over boxes (note that mfi takes a cell-centered multifab as an argument)
    for ( MFIter mfi(phi_in,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Array4<Real const> & phi = phi_in.array(mfi);

        AMREX_D_TERM(const Array4<Real> & gphix = gphi[0].array(mfi);,
                     const Array4<Real> & gphiy = gphi[1].array(mfi);,
                     const Array4<Real> & gphiz = gphi[2].array(mfi););

        AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                     const Box & bx_y = mfi.nodaltilebox(1);,
                     const Box & bx_z = mfi.nodaltilebox(2););

        amrex::ParallelFor(bx_x, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            gphix(i,j,k,start_outcomp+n) += (phi(i,j,k,start_incomp+n)-phi(i-1,j,k,start_incomp+n)) / dx[0];
        },
                           bx_y, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            gphiy(i,j,k,start_outcomp+n) += (phi(i,j,k,start_incomp+n)-phi(i,j-1,k,start_incomp+n)) / dx[1];
        }
                         , bx_z, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            gphiz(i,j,k,start_outcomp+n) += (phi(i,j,k,start_incomp+n)-phi(i,j,k-1,start_incomp+n)) / dx[2];
        }
        );

    } // end MFIter
}

void FCMoments(const std::array<MultiFab, 3>& m1,
               const amrex::Vector<int>& comps,
               std::array<MultiFab, 3>&  mscr,
               const int& power,
               amrex::Vector<amrex::Real>& prod_val)
{

    BL_PROFILE_VAR("FCMoments()",FCMoments);

    for (int d=0; d<3; ++d) {
        MultiFab::Copy(mscr[d],m1[d],comps[d],0,1,0);
        for(int i=1; i<power; i++){
            MultiFab::Multiply(mscr[d],m1[d],comps[d],0,1,0);
        }
    }
    SumStag(mscr,prod_val);
}

void SumStag(const std::array<MultiFab, 3>& m1,
             amrex::Vector<amrex::Real>& sum)
{
    BL_PROFILE_VAR("SumStag()",SumStag);

    // Initialize to zero
    std::fill(sum.begin(), sum.end(), 0.);

    ReduceOps<ReduceOpSum> reduce_op;

    //////// x-faces

    ReduceData<Real> reduce_datax(reduce_op);
    using ReduceTuple = typename decltype(reduce_datax)::Type;

    for (MFIter mfi(m1[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Box& bx_grid = mfi.validbox();

        auto const& fab = m1[0].array(mfi);

        int xlo = bx_grid.smallEnd(0);
        int xhi = bx_grid.bigEnd(0);

        reduce_op.eval(bx, reduce_datax,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real weight = (i>xlo && i<xhi) ? 1.0 : 0.5;
            return {fab(i,j,k)*weight};
        });
    }

    sum[0] = amrex::get<0>(reduce_datax.value());
    ParallelDescriptor::ReduceRealSum(sum[0]);

    //////// y-faces

    ReduceData<Real> reduce_datay(reduce_op);

    for (MFIter mfi(m1[1],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Box& bx_grid = mfi.validbox();

        auto const& fab = m1[1].array(mfi);

        int ylo = bx_grid.smallEnd(1);
        int yhi = bx_grid.bigEnd(1);

        reduce_op.eval(bx, reduce_datay,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real weight = (j>ylo && j<yhi) ? 1.0 : 0.5;
            return {fab(i,j,k)*weight};
        });
    }

    sum[1] = amrex::get<0>(reduce_datay.value());
    ParallelDescriptor::ReduceRealSum(sum[1]);

    //////// z-faces

    ReduceData<Real> reduce_dataz(reduce_op);

    for (MFIter mfi(m1[2],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Box& bx_grid = mfi.validbox();

        auto const& fab = m1[2].array(mfi);

        int zlo = bx_grid.smallEnd(2);
        int zhi = bx_grid.bigEnd(2);

        reduce_op.eval(bx, reduce_dataz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real weight = (k>zlo && k<zhi) ? 1.0 : 0.5;
            return {fab(i,j,k)*weight};
        });
    }

    sum[2] = amrex::get<0>(reduce_dataz.value());
    ParallelDescriptor::ReduceRealSum(sum[2]);
}

void CCMoments(const amrex::MultiFab& m1,
               const int& comp1,
               amrex::MultiFab& mscr,
               const int& power,
               amrex::Real& prod_val)
{

    BL_PROFILE_VAR("CCMoments()",CCMoments);

    MultiFab::Copy(mscr,m1,comp1,0,1,0);
    for(int i=1; i<power; i++){
        MultiFab::Multiply(mscr,m1,comp1,0,1,0);
    }

    prod_val = 0.;
    SumCC(mscr,0,prod_val,false);
}

void SumCC(const amrex::MultiFab& m1,
           const int& comp,
           amrex::Real& sum,
           const bool& divide_by_ncells)
{
    BL_PROFILE_VAR("SumCC()",SumCC);

    sum = 0.;
    sum = m1.MultiFab::sum(comp, false);

    if (divide_by_ncells == 1) {
        BoxArray ba_temp = m1.boxArray();
        long numpts = ba_temp.numPts();
        sum = sum/(double)(numpts);
    }
}
