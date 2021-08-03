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
    int r;
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
    int navg;
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
                                const int& ncomp)
{
    // number of points in the averaging direction
    int npts = n_cells[dir];

    Vector<Real> average(npts*(ncomp),0.);

    // dummy variables
    int r;
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
    int navg;
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
                mf(i,j,k,comp) = average[ r*(ncomp) + n];
            }
            }
            }
        }

    } // end MFiter
}


void ComputeVerticalAverage(const MultiFab& mf, MultiFab& mf_avg,
			    const Geometry& geom, const int& dir,
			    const int& incomp, const int& ncomp,
                            const int& slablo, const int& slabhi)
{
    BL_PROFILE_VAR("ComputVerticalAverage()",ComputeVerticalAverage);
  
    bool write_data = false;
    std::string plotname;

    MultiFab mf_pencil;

    BoxArray ba_in = mf.boxArray();
    BoxArray ba_pencil;
    BoxArray ba_flat;

    int nbox = ba_in.size();
    Box domain(geom.Domain());
    const DistributionMapping& dmap = mf.DistributionMap();

    Vector<int> max_grid_size_pencil(AMREX_SPACEDIM);
    Vector<int> max_grid_size_flat(AMREX_SPACEDIM);

    int nx[AMREX_SPACEDIM  ], nbx[AMREX_SPACEDIM  ];
    int mx[AMREX_SPACEDIM-1], mbx[AMREX_SPACEDIM-1];

    int nxprod, indlo, indhi, a, b, c;

    // We assume that all grids have the same size hence
    // we have the same nx,ny,nz on all ranks
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        nx[d]  = ba_in[0].size()[d];
        nbx[d] = domain.length(d) / nx[d];
    }
#if (AMREX_SPACEDIM == 2)
    if (nbx[0]*nbx[1] != ba_in.size())
        amrex::Error("ALL GRIDS DO NOT HAVE SAME SIZE");
#elif (AMREX_SPACEDIM == 3)
    if (nbx[0]*nbx[1]*nbx[2] != ba_in.size())
        amrex::Error("ALL GRIDS DO NOT HAVE SAME SIZE");
#endif

#if (AMREX_SPACEDIM == 2)
    if (dir == 0) {
        indlo = 1;
    } else if (dir == 1) {
        indlo = 0;
    } else {
        Abort("ComputeVerticalAverage: invalid dir");
    }
#elif (AMREX_SPACEDIM == 3)
    if (dir == 0) {
        indlo = 1;
        indhi = 2;
    } else if (dir == 1) {
        indlo = 0;
        indhi = 2;
    } else if (dir == 2) {
        indlo = 0;
        indhi = 1;
    } else {
        Abort("ComputeVerticalAverage: invalid dir");
    }
#endif

    IntVect dom_lo(domain.loVect());
    IntVect dom_hi(domain.hiVect());
    dom_hi[dir] = 0;
    Box domain_flat(dom_lo, dom_hi);

#if (AMREX_SPACEDIM == 2)
    mx[0] = max_grid_projection[0];

    mbx[0] = domain.length(indlo)/mx[0];

    Print() << "1D redist: " << mbx[0] << ", grids: " << mx[0] << std::endl;
#elif (AMREX_SPACEDIM == 3)
    mx[0] = max_grid_projection[0];
    mx[1] = max_grid_projection[1];

    mbx[0] = domain.length(indlo)/mx[0];
    mbx[1] = domain.length(indhi)/mx[1];

    Print() << "2D redist: " << mbx[0] << "x" << mbx[1] << ", grids: " << mx[0]
            << "x" << mx[1] << std::endl;
#endif

    max_grid_size_pencil[indlo] = mx[0];
#if (AMREX_SPACEDIM == 3)
    max_grid_size_pencil[indhi] = mx[1];
#endif
    max_grid_size_pencil[dir]   = domain.length(dir);
    ba_pencil.define(domain);
    ba_pencil.maxSize(IntVect(max_grid_size_pencil));
    DistributionMapping dmap_pencil(ba_pencil);
    mf_pencil.define(ba_pencil,dmap_pencil,ncomp,0);

    max_grid_size_flat      = max_grid_size_pencil;
    max_grid_size_flat[dir] = 1;
    ba_flat.define(domain_flat);
    ba_flat.maxSize(IntVect(max_grid_size_flat));
    mf_avg.define(ba_flat,dmap_pencil,ncomp,0);
    mf_avg.setVal(0.);

    // copy/redistrubute to pencils

    mf_pencil.ParallelCopy(mf, incomp, 0, ncomp);

    if (write_data) {
        plotname = "mf_pencil";
        VisMF::Write(mf_pencil,plotname);
    }

    Real ninv;
    if (slablo != -1 && slabhi != 99999) {
        ninv = 1./(slabhi-slablo+1);
    } else {
        ninv = 1./(domain.length(dir));
    }
    
    for ( MFIter mfi(mf_pencil); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<Real> meanfab = mf_avg.array(mfi);
        const Array4<Real> inputfab = mf_pencil.array(mfi);

        if (dir == 0) {
        
            for (auto n = incomp; n<incomp+ncomp; ++n) {
            for (auto k = lo.z; k <= hi.z; ++k) {
            for (auto j = lo.y; j <= hi.y; ++j) {
            for (auto i = lo.x; i <= hi.x; ++i) {
                if ((i >= slablo) and (i <= slabhi)) {
                    meanfab(0,j,k,n) = meanfab(0,j,k,n) + ninv*inputfab(i,j,k,n);
                }
            }
            }
            }
            }
            
        } else if (dir == 1) {
        
            for (auto n = incomp; n<incomp+ncomp; ++n) {
            for (auto k = lo.z; k <= hi.z; ++k) {
            for (auto j = lo.y; j <= hi.y; ++j) {
            for (auto i = lo.x; i <= hi.x; ++i) {
                if ((j >= slablo) and (j <= slabhi)) {
                    meanfab(i,0,k,n) = meanfab(i,0,k,n) + ninv*inputfab(i,j,k,n);
                }
            }
            }
            }
            }

        } else if (dir == 2) {
        
            for (auto n = incomp; n<incomp+ncomp; ++n) {
            for (auto k = lo.z; k <= hi.z; ++k) {
            for (auto j = lo.y; j <= hi.y; ++j) {
            for (auto i = lo.x; i <= hi.x; ++i) {
                if ((k >= slablo) and (k <= slabhi)) {
                    meanfab(i,j,0,n) = meanfab(i,j,0,n) + ninv*inputfab(i,j,k,n);
                }
            }
            }
            }
            }
        }
    }

    if (write_data) {
        if (slablo != -1 && slabhi != 99999) {
            plotname = "mf_avg_slab";
        } else {
            plotname = "mf_avg";
        }
        VisMF::Write(mf_avg,plotname);
    }

}

void ExtractSlice(const MultiFab& mf, MultiFab& mf_slice,
                  const Geometry& geom, const int dir,
                  const int incomp, const int ncomp)
{
    BL_PROFILE_VAR("ExtractSlice()",ExtractSlice);

    Box domain(geom.Domain());

    // create BoxArray
    IntVect dom_lo(domain.loVect());
    IntVect dom_hi(domain.hiVect());
    dom_lo[dir] = slicepoint;
    dom_hi[dir] = slicepoint;

    Box domain_slice(dom_lo, dom_hi);

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

    BoxArray ba_slice(domain_slice);
    ba_slice.maxSize(IntVect(max_grid_slice));

    DistributionMapping dmap_slice(ba_slice);

    mf_slice.define(ba_slice,dmap_slice,ncomp,0);

    mf_slice.ParallelCopy(mf, incomp, 0, ncomp);
}

namespace amrex 
{
    Gpu::HostVector<Real> sumToLine (MultiFab const& mf, int icomp, int ncomp,
                                          Box const& domain, int direction, bool local)
    {
        BL_PROFILE_VAR("sumToLine()",sumToLine);

        int n1d = domain.length(direction) * ncomp;
        Gpu::HostVector<Real> hv(n1d);

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            Gpu::DeviceVector<Real> dv(domain.length(direction), Real(0.0));
            Real* p = dv.data();

            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                Box const& b = mfi.validbox();
                const auto lo = lbound(b);
                const auto len = length(b);
                auto const& fab = mf.const_array(mfi);

                int n2d, n2dx;
                if (direction == 0) {
                    n2d = len.y * len.z;
                    n2dx = len.y;
                } else if (direction == 1) {
                    n2d = len.x * len.z;
                    n2dx = len.x;
                } else {
                    n2d = len.x * len.y;
                    n2dx = len.x;
                }
                int n2dblocks = (n2d+AMREX_GPU_MAX_THREADS-1)/AMREX_GPU_MAX_THREADS;
                int nblocks = n2dblocks * b.length(direction);
#ifdef AMREX_USE_DPCPP
                std::size_t shared_mem_byte = sizeof(Real)*Gpu::Device::warp_size;
                launch(nblocks, AMREX_GPU_MAX_THREADS, shared_mem_byte, Gpu::gpuStream(),
                              [=] AMREX_GPU_DEVICE (Gpu::Handler const& h) noexcept
#else
                launch(nblocks, AMREX_GPU_MAX_THREADS, Gpu::gpuStream(),
                              [=] AMREX_GPU_DEVICE () noexcept
#endif
                {
#ifdef AMREX_USE_DPCPP
                     int i1d = h.blockIdx() / n2dblocks;
                     int i2d = h.threadIdx() + h.blockDim()*(h.blockIdx()-i1d*n2dblocks);
#else
                     int i1d = blockIdx.x / n2dblocks;
                     int i2d = threadIdx.x + blockDim.x*(blockIdx.x-i1d*n2dblocks);
#endif
                     int i2dy = i2d / n2dx;
                     int i2dx = i2d - i2dy*n2dx;
                     int i, j, k, idir;
                     if (direction == 0) {
                         i = i1d  + lo.x;
                         j = i2dx + lo.y;
                         k = i2dy + lo.z;
                         idir = i;
                     } else if (direction == 1) {
                         i = i2dx + lo.x;
                         j = i1d  + lo.y;
                         k = i2dy + lo.z;
                         idir = j;
                     } else {
                         i = i2dx + lo.x;
                         j = i2dy + lo.y;
                         k = i1d  + lo.z;
                         idir = k;
                     }
                     for (int n = 0; n < ncomp; ++n) {
                         Real r = (i2d < n2d) ? fab(i,j,k,n+icomp) : Real(0.0);
#ifdef AMREX_USE_DPCPP
                     Gpu::deviceReduceSum_full(p+n+ncomp*idir, r, h);
#else
                     Gpu::deviceReduceSum_full(p+n+ncomp*idir, r);
#endif
                     }
                });
            }

            Gpu::copyAsync(Gpu::deviceToHost, dv.begin(), dv.end(), hv.begin());
            Gpu::streamSynchronize();
        }
        else
#endif
        {
            for (auto& x : hv) {
                x = Real(0.0);
            }

            Vector<Gpu::HostVector<Real> > other_hv(OpenMP::get_max_threads()-1);

            Vector<Real*> pp(OpenMP::get_max_threads());
            pp[0] = hv.data();
            for (int i = 1; i < OpenMP::get_max_threads(); ++i) {
                other_hv[i-1].resize(n1d, Real(0.0));
                pp[i] = other_hv[i-1].data();
            }

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
            for (MFIter mfi(mf,true); mfi.isValid(); ++mfi) {
                Box const& b = mfi.tilebox();
                auto const& fab = mf.const_array(mfi);
                Real * AMREX_RESTRICT p = pp[OpenMP::get_thread_num()];
                if (direction == 0) {
                    LoopOnCpu(b, ncomp, [&] (int i, int j, int k, int n) noexcept
                    {
                        p[n+ncomp*i] += fab(i,j,k,n+icomp);
                    });
                } else if (direction == 1) {
                    LoopOnCpu(b, ncomp, [&] (int i, int j, int k, int n) noexcept
                    {
                        p[n+ncomp*j] += fab(i,j,k,n+icomp);
                    });
                } else {
                    LoopOnCpu(b, ncomp, [&] (int i, int j, int k, int n) noexcept
                    {
                        p[n+ncomp*k] += fab(i,j,k,n+icomp);
                    });
                }
            }

            if (! other_hv.empty()) {
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
                for (int i = 0; i < n1d; ++i) {
                    for (auto const& v : other_hv) {
                        hv[i] += v[i];
                    }
                }
            }
        }

        if (!local) {
            ParallelAllReduce::Sum(hv.data(), hv.size(), ParallelContext::CommunicatorSub());
        }
        return hv;
    }
}
