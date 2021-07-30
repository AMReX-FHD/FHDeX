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
