#include "common_functions.H"

#include "AMReX_PlotFileUtil.H"

int greatest_common_factor(int,int);
void factor(int,int*,int);

namespace {

// number of cells in the planes transverse to dir
Long NumCellsTransverse (int dir)
{
    Long n = 1;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        if (d != dir) { n *= n_cells[d]; }
    }
    return n;
}

// Sum components incomp..incomp+ncomp-1 of mf_in over the planes transverse to dir.
// Returns a host vector of size n_cells[dir]*ncomp indexed [r*ncomp + n], reduced
// over all MPI ranks. Runs on the device when built for GPU.
// Each bin receives Ny*Nz (or similar) values, so the sums are accumulated in double:
// float atomics lose up to eps*n relative accuracy and are order-dependent.
Vector<double> SumOverPlanes (const MultiFab& mf_in, int dir, int incomp, int ncomp)
{
    const int npts = n_cells[dir];

    Gpu::DeviceVector<double> d_sum(npts*ncomp, 0.0);
    double* psum = d_sum.data();

    for (MFIter mfi(mf_in, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto const& mf = mf_in.const_array(mfi);
        amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            int r = (dir == 0) ? i : ((dir == 1) ? j : k);
            Gpu::Atomic::AddNoRet(&psum[r*ncomp + n], static_cast<double>(mf(i,j,k,incomp+n)));
        });
    }
    Gpu::streamSynchronize();

    Vector<double> h_sum(npts*ncomp);
    Gpu::copy(Gpu::deviceToHost, d_sum.begin(), d_sum.end(), h_sum.begin());

    // sum over all processors
    ParallelDescriptor::ReduceRealSum(h_sum.dataPtr(), npts*ncomp);

    return h_sum;
}

}

void WriteHorizontalAverage(const MultiFab& mf_in, const int& dir, const int& incomp,
                            const int& ncomp, const int& step, const Geometry& geom,
                            const std::string& file_prefix)
{
    // number of points in the averaging direction
    const int npts = n_cells[dir];

    const Real h = geom.CellSize(dir);

    Vector<double> sum = SumOverPlanes(mf_in, dir, incomp, ncomp);

    // divide by the number of cells in each transverse plane
    const double navg_inv = 1.0 / static_cast<double>(NumCellsTransverse(dir));

    if (ParallelDescriptor::IOProcessor()) {
        std::string filename = amrex::Concatenate(file_prefix,step,9);
        std::ofstream outfile;
        outfile.open(filename);

        // write out result; the first column is the physical coordinate
        for (int r=0; r<npts; ++r) {
            outfile << prob_lo[dir] + (r+Real(0.5))*h << " ";
            for (int n=0; n<ncomp; ++n) {
                outfile << Real(sum[r*ncomp + n]*navg_inv) << " ";
            }
            outfile << std::endl;
        }

        outfile.close();
    }
}

void WriteHorizontalAverageToMF(const MultiFab& mf_in, MultiFab& mf_out,
                                const int& dir, const int& incomp,
                                const int& ncomp, int outcomp)
{
    if (outcomp == -1) {
        outcomp = incomp; // default outcomp is incomp unless specified
    }

    // number of points in the averaging direction
    const int npts = n_cells[dir];

    Vector<double> sum = SumOverPlanes(mf_in, dir, incomp, ncomp);

    // divide by the number of cells in each transverse plane (in double), then convert to Real
    const double navg_inv = 1.0 / static_cast<double>(NumCellsTransverse(dir));
    Vector<Real> average(sum.size());
    for (std::size_t m=0; m<sum.size(); ++m) { average[m] = Real(sum[m]*navg_inv); }

    // copy the profile to the device and broadcast it into mf_out
    Gpu::DeviceVector<Real> d_average(npts*ncomp);
    Gpu::copy(Gpu::hostToDevice, average.begin(), average.end(), d_average.begin());
    Real const* pavg = d_average.data();

    for (MFIter mfi(mf_out, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto const& mf = mf_out.array(mfi);
        amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            int r = (dir == 0) ? i : ((dir == 1) ? j : k);
            mf(i,j,k,outcomp+n) = pavg[r*ncomp + n];
        });
    }
    Gpu::streamSynchronize(); // d_average must outlive the kernels
}


void ComputeVerticalAverage(const MultiFab& mf, MultiFab& mf_flat,
                            const int& dir,
                            const int& incomp, const int& ncomp,
                            const int& slablo, const int& slabhi)
{
    BL_PROFILE_VAR("ComputVerticalAverage()",ComputeVerticalAverage);

    // error checking
    if (dir >= AMREX_SPACEDIM) {
        Abort("ComputeVerticalAverage: invalid dir");
    }

    // get a single Box that spans the full domain
    Box domain(mf.boxArray().minimalBox());

    // these are the transverse directions (i.e., NOT the dir direction)
    int dir1=0, dir2=0;
#if (AMREX_SPACEDIM == 2)
    dir1 = 1-dir;
#elif (AMREX_SPACEDIM == 3)
    if (dir == 0) {
        dir1 = 1;
        dir2 = 2;
    } else if (dir == 1) {
        dir1 = 0;
        dir2 = 2;
    } else if (dir == 2) {
        dir1 = 0;
        dir2 = 1;
    }
#endif

    // max_grid_size_pencil will be equal to the number of cells in the domain in the dir direction
    // and uses max_grid_projection to set the non-dir directions
    Vector<int> max_grid_size_pencil(AMREX_SPACEDIM);
    max_grid_size_pencil[dir]  = domain.length(dir);
    max_grid_size_pencil[dir1] = max_grid_projection[0];
#if (AMREX_SPACEDIM == 3)
    max_grid_size_pencil[dir2] = max_grid_projection[1];
#endif

    // this is the inverse of the number of cells in the dir direction we are averaging over
    // by default we average over the entire domain, but one can pass in slab_lo/hi to set bounds
    // sum_domain is the region actually summed over; ReduceToPlane intersects each box with it
    // and collapses the dir index to 0 regardless of its range
    Box sum_domain(domain);
    double ninv;
    if (slablo != -1 && slabhi != 99999) {
        sum_domain.setSmall(dir, slablo);
        sum_domain.setBig(dir, slabhi);
        ninv = 1.0/(slabhi-slablo+1);
    } else {
        ninv = 1.0/(domain.length(dir));
    }

    MultiFab mf_onecomp(mf.boxArray(), mf.DistributionMap(), 1, 0);

    for (int n=0; n<ncomp; ++n) {

        // copy a component of mf into mf_onecomp
        MultiFab::Copy(mf_onecomp,mf,incomp+n,0,1,0);

        // sum up, accumulating in double (column sums of fields with a large mean lose the
        // O(1/sqrt(Nz)) averaged fluctuation to float summation error)
        auto const& ma = mf_onecomp.const_arrays();
        auto fabd = ReduceToPlane<ReduceOpSum,double>(dir, sum_domain, mf_onecomp,
          [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> double
          {
              return static_cast<double>(ma[box_no](i,j,k)); // data at (i,j,k) of Box box_no
          });

        // divide by number of cells in column (in double) to create average, then convert to Real
        FArrayBox fab(fabd.box(), 1);
        {
            auto const& src = fabd.const_array();
            auto const& dst = fab.array();
            const double ninv_d = ninv;
            amrex::ParallelFor(fabd.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                dst(i,j,k) = Real(src(i,j,k)*ninv_d);
            });
            Gpu::streamSynchronize();
        }

        Box dom2d = fab.box();
        Vector<Box> bv(ParallelDescriptor::NProcs(),dom2d);
        BoxArray ba(bv.data(), bv.size());

        Vector<int> pmap(ParallelDescriptor::NProcs());
        std::iota(pmap.begin(), pmap.end(), 0);
        DistributionMapping dm(std::move(pmap));

        MultiFab mftmp(ba, dm, 1, 0, MFInfo().SetAlloc(false));
        mftmp.setFab(ParallelDescriptor::MyProc(),
                     FArrayBox(fab.box(), 1, fab.dataPtr()));

        BoxArray ba2(dom2d);

        ba2.maxSize(IntVect(max_grid_size_pencil));

        if (n==0) {
            mf_flat.define(ba2, DistributionMapping{ba2}, ncomp, 0);
        }

        MultiFab mf_flat_onecomp(ba2, DistributionMapping{ba2}, fab.nComp(), 0);
        mf_flat_onecomp.setVal(0.);
        mf_flat_onecomp.ParallelAdd(mftmp);

        mf_flat.ParallelCopy(mf_flat_onecomp, 0, n, 1);
    }
}

void ExtractSlice(const MultiFab& mf, MultiFab& mf_slice,
                  const int dir, const int slice,
                  const int incomp, const int ncomp)
{
    BL_PROFILE_VAR("ExtractSlice()",ExtractSlice);

    // create BoxArray

    // get lo and hi coordinates of problem domain
    Box domain(mf.boxArray().minimalBox());
    IntVect dom_lo(domain.loVect());
    IntVect dom_hi(domain.hiVect());

    // set lo and hi coordinates in the dir direction to the slice point
    dom_lo[dir] = slice;
    dom_hi[dir] = slice;

    // create a BoxArray with a single box containing the flattened box
    Box domain_slice(dom_lo, dom_hi);
    BoxArray ba_slice(domain_slice);

    // chop up the box based on max_grid_projection
    Vector<int> max_grid_slice(AMREX_SPACEDIM);
#if (AMREX_SPACEDIM == 2)
    max_grid_slice[  dir] = 1;
    max_grid_slice[1-dir] = max_grid_projection[0];
#elif (AMREX_SPACEDIM == 3)
    max_grid_slice[dir] = 1;
    if (dir == 0) {
        max_grid_slice[1] = max_grid_projection[0];
        max_grid_slice[2] = max_grid_projection[1];
    } else if (dir == 1) {
        max_grid_slice[0] = max_grid_projection[0];
        max_grid_slice[2] = max_grid_projection[1];
    } else {
        max_grid_slice[0] = max_grid_projection[0];
        max_grid_slice[1] = max_grid_projection[1];
    }
#endif
    ba_slice.maxSize(IntVect(max_grid_slice));

    // create a new DistributionMapping and define the MultiFab
    DistributionMapping dmap_slice(ba_slice);
    MultiFab mf_slice_tmp(ba_slice,dmap_slice,ncomp,0);

    mf_slice_tmp.ParallelCopy(mf, incomp, 0, ncomp);

    // now copy this into a multifab with index zero in the dir direction rather than slicepoint
    // (structure factor code requires this)
    dom_lo[dir] = 0;
    dom_hi[dir] = 0;

    Box domain_slice2(dom_lo,dom_hi);
    BoxArray ba_slice2(domain_slice2);
    ba_slice2.maxSize(IntVect(max_grid_slice));
    mf_slice.define(ba_slice2,dmap_slice,ncomp,0);

    for ( MFIter mfi(mf_slice_tmp,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.tilebox();

        [[maybe_unused]] const Array4<Real> & slice_local = mf_slice.array(mfi);
        const Array4<Real> & slice_tmp = mf_slice_tmp.array(mfi);

        if (dir == 0) {
            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                slice_local(0,j,k,n) = slice_tmp(i,j,k,n);
            });
        }
        if (dir == 1) {
            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                slice_local(i,0,k,n) = slice_tmp(i,j,k,n);
            });
        }
        if (dir == 2) {
            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                slice_local(i,j,0,n) = slice_tmp(i,j,k,n);
            });
        }
    }
}


void ExtractXPencil(const MultiFab& mf, MultiFab& mf_pencil,
                    const int pencily, const int pencilz,
                    const int incomp, const int ncomp)
{
    BL_PROFILE_VAR("ExtractXPencil()",ExtractXPencil);

    // create BoxArray

    // get lo and hi coordinates of problem domain
    Box domain(mf.boxArray().minimalBox());
    IntVect dom_lo(domain.loVect());
    IntVect dom_hi(domain.hiVect());

    dom_lo[1] = dom_hi[1] = pencily;
    dom_lo[2] = dom_hi[2] = pencilz;

    // create a BoxArray with a single box containing the pencil
    Box domain_pencil(dom_lo, dom_hi);
    BoxArray ba_pencil(domain_pencil);

    // create a new DistributionMapping and define the MultiFab
    DistributionMapping dmap_pencil(ba_pencil);
    MultiFab mf_pencil_tmp(ba_pencil,dmap_pencil,ncomp,0);

    // copy data from full MF into pencil
    mf_pencil_tmp.ParallelCopy(mf, incomp, 0, ncomp);

    // now copy this into a multifab with index zero in x and y
    // (structure factor code requires this)
    dom_lo[1] = dom_hi[1] = 0;
    dom_lo[2] = dom_hi[2] = 0;

    Box domain_pencil2(dom_lo,dom_hi);
    BoxArray ba_pencil2(domain_pencil2);
    mf_pencil.define(ba_pencil2,dmap_pencil,ncomp,0);

    for ( MFIter mfi(mf_pencil_tmp,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.tilebox();

        const Array4<Real> & pencil = mf_pencil.array(mfi);
        const Array4<Real> & pencil_tmp = mf_pencil_tmp.array(mfi);

        amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            pencil(i,0,0,n) = pencil_tmp(i,j,k,n);
        });

    }
}
