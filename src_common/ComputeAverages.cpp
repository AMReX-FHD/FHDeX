#include "common_functions.H"

#include "AMReX_PlotFileUtil.H"

int greatest_common_factor(int,int);
void factor(int,int*,int);

void WriteHorizontalAverage(const MultiFab& mf_in, const int& dir, const int& incomp,
                            const int& ncomp, const int& step, const Geometry& geom,
                            const std::string& file_prefix)
{
    // number of points in the averaging direction
    int npts = n_cells[dir];

    // we use ncomp+1 because th first column is the coorinate
    Vector<Real> average(npts*(ncomp+1),0.);

    Real h = geom.CellSize(dir);

    // dummy variables
    int r=0;
    int comp;

    // no tiling or GPU to easily avoid race conditions
    for (MFIter mfi(mf_in, false); mfi.isValid(); ++mfi) {

        // valid box and the lo/hi coordinates; no ghost cells needed
        const Box& bx = mfi.validbox();
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<const Real> mf = mf_in.array(mfi);

        for (auto n=0; n<ncomp; ++n) {
            comp = incomp+n;
            for (auto k = lo.z; k <= hi.z; ++k) {
                for (auto j = lo.y; j <= hi.y; ++j) {
                    for (auto i = lo.x; i <= hi.x; ++i) {
                        if (dir == 0) {
                            r=i;
                        } else if (dir == 1) {
                            r=j;
                        } else if (dir == 2) {
                            r=k;
                        }
                        // sum up the data
                        // the "+1" is because the first column will store the physical coordinate
                        average[ r*(ncomp+1) + n + 1 ] += mf(i,j,k,comp);
                    }
                }
            }
        }

    } // end MFiter

    // sum over all processors
    ParallelDescriptor::ReduceRealSum(average.dataPtr(),npts*(ncomp+1));

    // divide by the number of cells
    int navg=0;
    if (dir == 0) {
        navg = n_cells[1]*n_cells[2];
    } else if (dir == 1) {
        navg = n_cells[0]*n_cells[2];
    } else if (dir == 2) {
        navg = n_cells[0]*n_cells[1];
    }
    for (r=0; r<npts; ++r) {
        for (auto n=1; n<ncomp+1; ++n) {
            average[r*(ncomp+1) + n] /= navg;
        }
    }

    // compute physical coordinate and store in first column
    for (r=0; r<npts; ++r) {
        average[r*(ncomp+1)] = prob_lo[dir] + (r+0.5)*h;
    }

    if (ParallelDescriptor::IOProcessor()) {
        std::string filename = amrex::Concatenate(file_prefix,step,9);
        std::ofstream outfile;
        outfile.open(filename);

        // write out result
        for (r=0; r<npts; ++r) {
            for (auto n=0; n<ncomp+1; ++n) {
                outfile << average[r*(ncomp+1) + n] << " ";
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
    int npts = n_cells[dir];

    Vector<Real> average(npts*(ncomp),0.);

    // dummy variables
    int r=0;
    int comp;

    // no tiling or GPU to easily avoid race conditions
    for (MFIter mfi(mf_in, false); mfi.isValid(); ++mfi) {

        // valid box and the lo/hi coordinates; no ghost cells needed
        const Box& bx = mfi.validbox();
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<const Real> mf = mf_in.array(mfi);

        for (auto n=0; n<ncomp; ++n) {
            comp = incomp+n;
            for (auto k = lo.z; k <= hi.z; ++k) {
                for (auto j = lo.y; j <= hi.y; ++j) {
                    for (auto i = lo.x; i <= hi.x; ++i) {
                        if (dir == 0) {
                            r=i;
                        } else if (dir == 1) {
                            r=j;
                        } else if (dir == 2) {
                            r=k;
                        }
                        // sum up the data
                        average[ r*(ncomp) + n] += mf(i,j,k,comp);
                    }
                }
            }
        }

    } // end MFiter

    // sum over all processors
    ParallelDescriptor::ReduceRealSum(average.dataPtr(),npts*(ncomp));

    // divide by the number of cells
    int navg=0;
    if (dir == 0) {
        navg = n_cells[1]*n_cells[2];
    } else if (dir == 1) {
        navg = n_cells[0]*n_cells[2];
    } else if (dir == 2) {
        navg = n_cells[0]*n_cells[1];
    }
    for (r=0; r<npts; ++r) {
        for (auto n=0; n<ncomp; ++n) {
            average[r*(ncomp) + n] /= navg;
        }
    }

    // no tiling or GPU to easily avoid race conditions
    for (MFIter mfi(mf_out, false); mfi.isValid(); ++mfi) {

        // valid box and the lo/hi coordinates; no ghost cells needed
        const Box& bx = mfi.validbox();
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<Real> mf = mf_out.array(mfi);

        for (auto n=0; n<ncomp; ++n) {
            comp = outcomp+n;
            for (auto k = lo.z; k <= hi.z; ++k) {
                for (auto j = lo.y; j <= hi.y; ++j) {
                    for (auto i = lo.x; i <= hi.x; ++i) {
                        if (dir == 0) {
                            r=i;
                        } else if (dir == 1) {
                            r=j;
                        } else if (dir == 2) {
                            r=k;
                        }
                        // sum up the data
                        mf(i,j,k,comp) = average[ r*(ncomp) + n];
                    }
                }
            }
        }

    } // end MFiter
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
    Real ninv;
    if (slablo != -1 && slabhi != 99999) {
        ninv = 1./(slabhi-slablo+1);
    } else {
        ninv = 1./(domain.length(dir));
    }

    MultiFab mf_onecomp(mf.boxArray(), mf.DistributionMap(), 1, 0);

    for (int n=0; n<ncomp; ++n) {

        // copy a component of mf into mf_onecomp
        MultiFab::Copy(mf_onecomp,mf,incomp+n,0,1,0);

        // sum up
        auto const& ma = mf_onecomp.const_arrays();
        auto fab = ReduceToPlane<ReduceOpSum,Real>(dir, domain, mf_onecomp,
          [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> Real
          {
              return ma[box_no](i,j,k); // data at (i,j,k) of Box box_no
          });

        Box dom2d = fab.box();
        Vector<Box> bv(ParallelDescriptor::NProcs(),dom2d);
        BoxArray ba(bv.data(), bv.size());

        Vector<int> pmap(ParallelDescriptor::NProcs());
        std::iota(pmap.begin(), pmap.end(), 0);
        DistributionMapping dm(std::move(pmap));

        MultiFab mftmp(ba, dm, 1, 0, MFInfo().SetAlloc(false));
        mftmp.setFab(ParallelDescriptor::MyProc(),
                     FArrayBox(fab.box(), 1, fab.dataPtr()));

        // divide by number of cells in column to create average
        mftmp.mult(ninv);

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

        const Array4<Real> & slice = mf_slice.array(mfi);
        const Array4<Real> & slice_tmp = mf_slice_tmp.array(mfi);

        if (dir == 0) {
            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                slice(0,j,k,n) = slice_tmp(i,j,k,n);
            });
        }
        if (dir == 1) {
            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                slice(i,0,k,n) = slice_tmp(i,j,k,n);
            });
        }
        if (dir == 2) {
            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                slice(i,j,0,n) = slice_tmp(i,j,k,n);
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