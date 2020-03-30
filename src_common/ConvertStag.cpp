#include "common_functions.H"

void AverageFaceToCC(const std::array<MultiFab, AMREX_SPACEDIM>& face_in,
                     MultiFab& cc_in, int cc_comp)
{

    BL_PROFILE_VAR("AverageFaceToCC()",AverageFaceToCC);

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(cc_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        AMREX_D_TERM(Array4<Real const> const& facex = face_in[0].array(mfi);,
                     Array4<Real const> const& facey = face_in[1].array(mfi);,
                     Array4<Real const> const& facez = face_in[2].array(mfi););

        Array4<Real> const& cc = cc_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            AMREX_D_TERM(cc(i,j,k,cc_comp  ) = 0.5*(facex(i+1,j,k) + facex(i,j,k));,
                         cc(i,j,k,cc_comp+1) = 0.5*(facey(i,j+1,k) + facey(i,j,k));,
                         cc(i,j,k,cc_comp+2) = 0.5*(facez(i,j,k+1) + facez(i,j,k)););
        });
    }
}

void AverageCCToFace(const MultiFab& cc_in, int cc_comp,
                     std::array<MultiFab, AMREX_SPACEDIM>& face_in, int face_comp,
                     int ncomp)
{

    BL_PROFILE_VAR("AverageCCToFace()",AverageCCToFace);

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(cc_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Array4<Real const> & cc = cc_in.array(mfi);

        AMREX_D_TERM(const Array4<Real> & facex = face_in[0].array(mfi);,
                     const Array4<Real> & facey = face_in[1].array(mfi);,
                     const Array4<Real> & facez = face_in[2].array(mfi););

        AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                     const Box & bx_y = mfi.nodaltilebox(1);,
                     const Box & bx_z = mfi.nodaltilebox(2););

        amrex::ParallelFor(bx_x, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            facex(i,j,k,face_comp+n) = 0.5*(cc(i,j,k,cc_comp+n)+cc(i-1,j,k,cc_comp+n));
        },
                           bx_y, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            facey(i,j,k,face_comp+n) = 0.5*(cc(i,j,k,cc_comp+n)+cc(i,j-1,k,cc_comp+n));
        }
#if (AMREX_SPACEDIM == 3)
        ,
                           bx_z, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            facez(i,j,k,face_comp+n) = 0.5*(cc(i,j,k,cc_comp+n)+cc(i,j,k-1,cc_comp+n));
        }
#endif
        );
    }

}


void ShiftFaceToCC(const MultiFab& face_in, int face_comp,
                   MultiFab& cc_in, int cc_comp, int ncomp)
{

    BL_PROFILE_VAR("ShiftFaceToCC()",ShiftFaceToCC);

    if (!face_in.is_nodal(0) && !face_in.is_nodal(1) && !face_in.is_nodal(2)) {
        Abort("ShiftFaceToCC requires a face-centered MultiFab");
    }

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
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

void AverageCCToNode(const MultiFab& cc_in, MultiFab& node_in, int scomp, int ncomp, int varType,
                     const Geometry& geom)
{

    BL_PROFILE_VAR("AverageCCToNode()",AverageCCToNode);

    // Physical Domain
    Box dom(geom.Domain());
    
    Vector<int> bc_lo(AMREX_SPACEDIM);
    Vector<int> bc_hi(AMREX_SPACEDIM);

    // compute mathematical boundary conditions
    BCPhysToMath(varType,bc_lo,bc_hi);

    int ng = node_in.nGrow();

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(node_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Array4<Real const> & cc = cc_in.array(mfi);
        const Array4<Real> & node = node_in.array(mfi);

        const Box& bx = mfi.growntilebox(ng);
        
        amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            // average to nodes
#if (AMREX_SPACEDIM == 2)
            node(i,j,k,scomp+n) = 0.25*(cc(i,j,k,scomp+n)+cc(i-1,j,k,scomp+n)+cc(i,j-1,k,scomp+n)+cc(i-1,j-1,k,scomp+n));
#elif (AMREX_SPACEDIM == 3)
            node(i,j,k,scomp+n) = 0.125*(cc(i,j,k,scomp+n)+cc(i-1,j,k,scomp+n)+cc(i,j-1,k,scomp+n)+cc(i,j,k-1,scomp+n)
                                 +cc(i-1,j-1,k,scomp+n)+cc(i-1,j,k-1,scomp+n)+cc(i,j-1,k-1,scomp+n)+cc(i-1,j-1,k-1,scomp+n));
#endif
        });

        // boundary conditions
        // note: at physical boundaries,
        // the value in ghost cells represents the value ON the boundary
        if (bc_lo[0] == FOEXTRAP || bc_lo[0] == EXT_DIR) {
            if (bx.smallEnd(0) <= dom.smallEnd(0)) {
                int lo = dom.smallEnd(0);
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (i <= lo) {
#if (AMREX_SPACEDIM == 2)
                        node(i,j,k,scomp+n) = 0.5*(cc(lo-1,j,k,scomp+n) + cc(lo-1,j-1,k,scomp+n));
#elif (AMREX_SPACEDIM == 3)
                        node(i,j,k,scomp+n) = 0.25*(cc(lo-1,j,k,scomp+n) + cc(lo-1,j-1,k,scomp+n) + cc(lo-1,j,k-1,scomp+n) + cc(lo-1,j-1,k-1,scomp+n));
#endif
                    }
                });
            }
        }

        if (bc_hi[0] == FOEXTRAP || bc_hi[0] == EXT_DIR) {
            if (bx.bigEnd(0) >= dom.bigEnd(0)+1) {
                int hi = dom.bigEnd(0);
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (i >= hi+1) {
#if (AMREX_SPACEDIM == 2)
                        node(i,j,k,scomp+n) = 0.5*(cc(hi+1,j,k,scomp+n) + cc(hi+1,j-1,k,scomp+n));
#elif (AMREX_SPACEDIM == 3)
                        node(i,j,k,scomp+n) = 0.25*(cc(hi+1,j,k,scomp+n) + cc(hi+1,j-1,k,scomp+n) + cc(hi+1,j,k-1,scomp+n) + cc(hi+1,j-1,k-1,scomp+n));
#endif
                    }
                });
            }
        }
        
        if (bc_lo[1] == FOEXTRAP || bc_lo[1] == EXT_DIR) {
            if (bx.smallEnd(1) <= dom.smallEnd(1)) {
                int lo = dom.smallEnd(1);
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (j <= lo) {
#if (AMREX_SPACEDIM == 2)
                        node(i,j,k,scomp+n) = 0.5*(cc(i,lo-1,k,scomp+n) + cc(i-1,lo-1,k,scomp+n));
#elif (AMREX_SPACEDIM == 3)
                        node(i,j,k,scomp+n) = 0.25*(cc(i,lo-1,k,scomp+n) + cc(i-1,lo-1,k,scomp+n) + cc(i,lo-1,k-1,scomp+n) + cc(i-1,lo-1,k-1,scomp+n));
#endif
                    }
                });
            }
        }

        if (bc_hi[1] == FOEXTRAP || bc_hi[1] == EXT_DIR) {
            if (bx.bigEnd(1) >= dom.bigEnd(1)+1) {
                int hi = dom.bigEnd(1);
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (j >= hi+1) {
#if (AMREX_SPACEDIM == 2)
                        node(i,j,k,scomp+n) = 0.5*(cc(i,hi+1,k,scomp+n) + cc(i-1,hi+1,k,scomp+n));
#elif (AMREX_SPACEDIM == 3)
                        node(i,j,k,scomp+n) = 0.25*(cc(i,hi+1,k,scomp+n) + cc(i-1,hi+1,k,scomp+n) + cc(i,hi+1,k-1,scomp+n) + cc(i-1,hi+1,j-1,scomp+n));
#endif
                    }
                });
            }
        }

#if (AMREX_SPACEDIM == 3)
        
        if (bc_lo[2] == FOEXTRAP || bc_lo[2] == EXT_DIR) {
            if (bx.smallEnd(2) <= dom.smallEnd(2)) {
                int lo = dom.smallEnd(2);
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (k <= lo) {
                        node(i,j,k,scomp+n) = 0.25*(cc(i,j,lo-1,scomp+n) + cc(i,j-1,lo-1,scomp+n) + cc(i-1,j,lo-1,scomp+n) + cc(i-1,j-1,lo-1,scomp+n));
                    }
                });
            }
        }

        if (bc_hi[2] == FOEXTRAP || bc_hi[2] == EXT_DIR) {
            if (bx.bigEnd(2) >= dom.bigEnd(2)+1) {
                int hi = dom.bigEnd(2);
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (k >= hi+1) {
                        node(i,j,k,scomp+n) = 0.25*(cc(i,j,hi+1,scomp+n) + cc(i,j-1,hi+1,scomp+n) + cc(i-1,j,hi+1,scomp+n) + cc(i-1,j-1,hi+1,scomp+n));
                    }
                });
            }
        }
        
#endif
        
    } // end MFIter
}
