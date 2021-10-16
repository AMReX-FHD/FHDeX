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


void ComputeVerticalAverage(const MultiFab& mf, MultiFab& mf_flat,
			    const Geometry& geom, const int& dir,
			    const int& incomp, const int& ncomp,
                            const int& slablo, const int& slabhi)
{
    BL_PROFILE_VAR("ComputVerticalAverage()",ComputeVerticalAverage);

    // error checking
    if (dir >= AMREX_SPACEDIM) {
        Abort("ComputeVerticalAverage: invalid dir");
    }
    
    // debugging
    bool write_data = false;

    // this is a full MultiFab with pencil-shaped boxes
    // we will define mf_flat as a flattened MultiFab that
    // has the same BoxArray but flattened in the dir direction
    // and the same DistributionMapping so
    // we can do the averaging from mf_pencil to mf_flat on a box-by-box basis
    MultiFab mf_pencil;

    // get a single Box that spans the full domain
    Box domain(geom.Domain());

    // these are the transverse directions (i.e., NOT the dir direction)
    int dir1, dir2;
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

    // create the BoxArray for the pencil MultiFab
    BoxArray ba_pencil(domain);
    ba_pencil.maxSize(IntVect(max_grid_size_pencil));

    // create DistributionMapping on the pencil BoxArray
    DistributionMapping dmap_pencil(ba_pencil);

    // build pencil MultiFab
    mf_pencil.define(ba_pencil,dmap_pencil,ncomp,0);

    // copy data from full MultiFab to pencil MultiFab
    mf_pencil.ParallelCopy(mf, incomp, 0, ncomp);

    // create a single flattened box with coordinate index 0 in the dir direction
    IntVect dom_lo(domain.loVect());
    IntVect dom_hi(domain.hiVect());
    if (dom_lo[dir] != 0) {
        Abort("ComputeVerticalAverage requires dom_lo[dir]=0");
    }
    dom_hi[dir] = 0;
    Box domain_flat(dom_lo, dom_hi);
    
    // create the BoxArray for the flattened MultiFab
    BoxArray ba_flat(domain_flat);
    ba_flat.maxSize(IntVect(max_grid_size_pencil));

    // build flattened MultiFab and initialize to zero
    mf_flat.define(ba_flat,dmap_pencil,ncomp,0);
    mf_flat.setVal(0.);

    // this is the inverse of the number of cells in the dir direction we are averaging over
    // by default we average over the entire domain, but one can pass in slab_lo/hi to set bounds
    Real ninv;
    if (slablo != -1 && slabhi != 99999) {
        ninv = 1./(slabhi-slablo+1);
    } else {
        ninv = 1./(domain.length(dir));
    }

    // average pencil data onto the flattened MultiFab
    for ( MFIter mfi(mf_pencil); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<Real> meanfab = mf_flat.array(mfi);
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

    // debugging
    if (write_data) {
        VisMF::Write(mf,"mf_full");
        VisMF::Write(mf_pencil,"mf_pencil");
        VisMF::Write(mf_flat,"mf_flat");
    }

}

void ExtractSlice(const MultiFab& mf, MultiFab& mf_slice,
                  const Geometry& geom, const int dir, const int slice,
                  const int incomp, const int ncomp)
{
    BL_PROFILE_VAR("ExtractSlice()",ExtractSlice);
    
    // create BoxArray

    // get lo and hi coordinates of problem domain
    Box domain(geom.Domain());
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
    mf_slice.define(ba_slice,dmap_slice,ncomp,0);
        
    mf_slice.ParallelCopy(mf, incomp, 0, ncomp);
}
