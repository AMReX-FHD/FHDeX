#include "common_functions.H"


// Computes divergence at cell centres from velcocities at cell faces
void ComputeDiv(MultiFab& div,
                const std::array<MultiFab, AMREX_SPACEDIM>& phi_fc,
                int start_incomp, int start_outcomp, int ncomp,
                const Geometry& geom, Real increment)
{
    BL_PROFILE_VAR("ComputeDiv()",ComputeDiv);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(div,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        
        const Box& bx = mfi.tilebox();

        const Array4<Real> & div_fab = div.array(mfi);
        AMREX_D_TERM(Array4<Real const> const& phix_fab = phi_fc[0].array(mfi);,
                     Array4<Real const> const& phiy_fab = phi_fc[1].array(mfi);,
                     Array4<Real const> const& phiz_fab = phi_fc[2].array(mfi););

        if (increment == 0.) {
            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                div_fab(i,j,k,start_outcomp+n) =
                    AMREX_D_TERM(  (phix_fab(i+1,j,k,start_incomp+n) - phix_fab(i,j,k,start_incomp+n)) / dx[0],
                                 + (phiy_fab(i,j+1,k,start_incomp+n) - phiy_fab(i,j,k,start_incomp+n)) / dx[1],
                                 + (phiz_fab(i,j,k+1,start_incomp+n) - phiz_fab(i,j,k,start_incomp+n)) / dx[2]);;
            });
        }
        else
        {
            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                Real temp = 
                    AMREX_D_TERM(  (phix_fab(i+1,j,k,start_incomp+n) - phix_fab(i,j,k,start_incomp+n)) / dx[0],
                                 + (phiy_fab(i,j+1,k,start_incomp+n) - phiy_fab(i,j,k,start_incomp+n)) / dx[1],
                                 + (phiz_fab(i,j,k+1,start_incomp+n) - phiz_fab(i,j,k,start_incomp+n)) / dx[2]);;
                
                div_fab(i,j,k,start_outcomp+n) += increment * temp;
            });
        }
    }
}

// Computes gradient at cell faces of cell centred scalar
void ComputeGrad(const MultiFab & phi_in, std::array<MultiFab, AMREX_SPACEDIM> & gphi,
                 int start_incomp, int start_outcomp, int ncomp, int bccomp, const Geometry & geom,
                 int increment)
{
    BL_PROFILE_VAR("ComputeGrad()",ComputeGrad);
    
    // Physical Domain
    Box dom(geom.Domain());
    
    Vector<int> bc_lo(AMREX_SPACEDIM);
    Vector<int> bc_hi(AMREX_SPACEDIM);

    // compute mathematical boundary conditions
    BCPhysToMath(bccomp,bc_lo,bc_hi);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // Loop over boxes (note that mfi takes a cell-centered multifab as an argument)
    for ( MFIter mfi(phi_in,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Array4<Real const> & phi = phi_in.array(mfi);

        AMREX_D_TERM(const Array4<Real> & gphix = gphi[0].array(mfi);,
                     const Array4<Real> & gphiy = gphi[1].array(mfi);,
                     const Array4<Real> & gphiz = gphi[2].array(mfi););

        AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                     const Box & bx_y = mfi.nodaltilebox(1);,
                     const Box & bx_z = mfi.nodaltilebox(2););

        if (increment == 0) {
        
            amrex::ParallelFor(bx_x, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                gphix(i,j,k,start_outcomp+n) = (phi(i,j,k,start_incomp+n)-phi(i-1,j,k,start_incomp+n)) / dx[0];
            },
                               bx_y, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                gphiy(i,j,k,start_outcomp+n) = (phi(i,j,k,start_incomp+n)-phi(i,j-1,k,start_incomp+n)) / dx[1];
            }
#if (AMREX_SPACEDIM == 3)
                               , bx_z, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                gphiz(i,j,k,start_outcomp+n) = (phi(i,j,k,start_incomp+n)-phi(i,j,k-1,start_incomp+n)) / dx[2];
            }
#endif
            );

            // boundary conditions
            // note: at physical boundaries,
            // alter stencil at boundary since ghost value represents value at boundary
            if (bc_lo[0] == FOEXTRAP || bc_lo[0] == EXT_DIR) {
                if (bx_x.smallEnd(0) <= dom.smallEnd(0)) {
                    int lo = dom.smallEnd(0);
                    amrex::ParallelFor(bx_x, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        if (i == lo) {
                            gphix(i,j,k,start_outcomp+n) = (phi(i,j,k,start_incomp+n)-phi(i-1,j,k,start_incomp+n)) / (0.5*dx[0]);
                        }
                    });
                }
            }
            
            if (bc_hi[0] == FOEXTRAP || bc_hi[0] == EXT_DIR) {
                if (bx_x.bigEnd(0) >= dom.bigEnd(0)+1) {
                    int hi = dom.bigEnd(0)+1;
                    amrex::ParallelFor(bx_x, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        if (i == hi) {
                            gphix(i,j,k,start_outcomp+n) = (phi(i,j,k,start_incomp+n)-phi(i-1,j,k,start_incomp+n)) / (0.5*dx[0]);
                        }
                    });
                }
            }
        
            if (bc_lo[1] == FOEXTRAP || bc_lo[1] == EXT_DIR) {
                if (bx_y.smallEnd(1) <= dom.smallEnd(1)) {
                    int lo = dom.smallEnd(1);
                    amrex::ParallelFor(bx_y, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        if (j == lo) {
                            gphiy(i,j,k,start_outcomp+n) = (phi(i,j,k,start_incomp+n)-phi(i,j-1,k,start_incomp+n)) / (0.5*dx[1]);
                        }
                    });
                }
            }
            
            if (bc_hi[1] == FOEXTRAP || bc_hi[1] == EXT_DIR) {
                if (bx_y.bigEnd(1) >= dom.bigEnd(1)+1) {
                    int hi = dom.bigEnd(1)+1;
                    amrex::ParallelFor(bx_y, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        if (j == hi) {
                            gphiy(i,j,k,start_outcomp+n) = (phi(i,j,k,start_incomp+n)-phi(i,j-1,k,start_incomp+n)) / (0.5*dx[1]);
                        }
                    });
                }
            }

#if (AMREX_SPACEDIM == 3)
            if (bc_lo[2] == FOEXTRAP || bc_lo[2] == EXT_DIR) {
                if (bx_z.smallEnd(2) <= dom.smallEnd(2)) {
                    int lo = dom.smallEnd(2);
                    amrex::ParallelFor(bx_z, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        if (k == lo) {
                            gphiz(i,j,k,start_outcomp+n) = (phi(i,j,k,start_incomp+n)-phi(i,j,k-1,start_incomp+n)) / (0.5*dx[2]);
                        }
                    });
                }
            }
            
            if (bc_hi[2] == FOEXTRAP || bc_hi[2] == EXT_DIR) {
                if (bx_z.bigEnd(2) >= dom.bigEnd(2)+1) {
                    int hi = dom.bigEnd(2)+1;
                    amrex::ParallelFor(bx_z, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        if (k == hi) {
                            gphiz(i,j,k,start_outcomp+n) = (phi(i,j,k,start_incomp+n)-phi(i,j,k-1,start_incomp+n)) / (0.5*dx[2]);
                        }
                    });
                }
            }
#endif
        }
        else { // increment == 1
        
            amrex::ParallelFor(bx_x, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                gphix(i,j,k,start_outcomp+n) += (phi(i,j,k,start_incomp+n)-phi(i-1,j,k,start_incomp+n)) / dx[0];
            },
                               bx_y, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                gphiy(i,j,k,start_outcomp+n) += (phi(i,j,k,start_incomp+n)-phi(i,j-1,k,start_incomp+n)) / dx[1];
            }
#if (AMREX_SPACEDIM == 3)
                               , bx_z, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                gphiz(i,j,k,start_outcomp+n) += (phi(i,j,k,start_incomp+n)-phi(i,j,k-1,start_incomp+n)) / dx[2];
            }
#endif
            );

            // boundary conditions
            // note: at physical boundaries,
            // alter stencil at boundary since ghost value represents value at boundary
            if (bc_lo[0] == FOEXTRAP || bc_lo[0] == EXT_DIR) {
                if (bx_x.smallEnd(0) <= dom.smallEnd(0)) {
                    int lo = dom.smallEnd(0);
                    amrex::ParallelFor(bx_x, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        if (i == lo) {
                            gphix(i,j,k,start_outcomp+n) += (phi(i,j,k,start_incomp+n)-phi(i-1,j,k,start_incomp+n)) / (0.5*dx[0]);
                        }
                    });
                }
            }
            
            if (bc_hi[0] == FOEXTRAP || bc_hi[0] == EXT_DIR) {
                if (bx_x.bigEnd(0) >= dom.bigEnd(0)+1) {
                    int hi = dom.bigEnd(0)+1;
                    amrex::ParallelFor(bx_x, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        if (i == hi) {
                            gphix(i,j,k,start_outcomp+n) += (phi(i,j,k,start_incomp+n)-phi(i-1,j,k,start_incomp+n)) / (0.5*dx[0]);
                        }
                    });
                }
            }
        
            if (bc_lo[1] == FOEXTRAP || bc_lo[1] == EXT_DIR) {
                if (bx_y.smallEnd(1) <= dom.smallEnd(1)) {
                    int lo = dom.smallEnd(1);
                    amrex::ParallelFor(bx_y, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        if (j == lo) {
                            gphiy(i,j,k,start_outcomp+n) += (phi(i,j,k,start_incomp+n)-phi(i,j-1,k,start_incomp+n)) / (0.5*dx[1]);
                        }
                    });
                }
            }
            
            if (bc_hi[1] == FOEXTRAP || bc_hi[1] == EXT_DIR) {
                if (bx_y.bigEnd(1) >= dom.bigEnd(1)+1) {
                    int hi = dom.bigEnd(1)+1;
                    amrex::ParallelFor(bx_y, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        if (j == hi) {
                            gphiy(i,j,k,start_outcomp+n) += (phi(i,j,k,start_incomp+n)-phi(i,j-1,k,start_incomp+n)) / (0.5*dx[1]);
                        }
                    });
                }
            }

#if (AMREX_SPACEDIM == 3)
            if (bc_lo[2] == FOEXTRAP || bc_lo[2] == EXT_DIR) {
                if (bx_z.smallEnd(2) <= dom.smallEnd(2)) {
                    int lo = dom.smallEnd(2);
                    amrex::ParallelFor(bx_z, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        if (k == lo) {
                            gphiz(i,j,k,start_outcomp+n) += (phi(i,j,k,start_incomp+n)-phi(i,j,k-1,start_incomp+n)) / (0.5*dx[2]);
                        }
                    });
                }
            }
            
            if (bc_hi[2] == FOEXTRAP || bc_hi[2] == EXT_DIR) {
                if (bx_z.bigEnd(2) >= dom.bigEnd(2)+1) {
                    int hi = dom.bigEnd(2)+1;
                    amrex::ParallelFor(bx_z, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                    {
                        if (k == hi) {
                            gphiz(i,j,k,start_outcomp+n) += (phi(i,j,k,start_incomp+n)-phi(i,j,k-1,start_incomp+n)) / (0.5*dx[2]);
                        }
                    });
                }
            }
#endif
        } // end increment test
    } // end MFIter
}

// Computes gradient at cell centres from cell centred data - ouputs to a three component mf.
void ComputeCentredGrad(const MultiFab & phi,
                        std::array<MultiFab, AMREX_SPACEDIM> & gphi,
                        const Geometry & geom)
{
    BL_PROFILE_VAR("ComputeCentredGrad()",ComputeCentredGrad);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(phi,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.tilebox();

        Array4<Real const> const& phi_fab = phi.array(mfi);

        AMREX_D_TERM(Array4<Real> const& gphix_fab = gphi[0].array(mfi);,
                     Array4<Real> const& gphiy_fab = gphi[1].array(mfi);,
                     Array4<Real> const& gphiz_fab = gphi[2].array(mfi););

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            AMREX_D_TERM(gphix_fab(i,j,k) = (phi_fab(i+1,j,k) - phi_fab(i-1,j,k) ) / (2.*dx[0]);,
                         gphiy_fab(i,j,k) = (phi_fab(i,j+1,k) - phi_fab(i,j-1,k) ) / (2.*dx[1]);,
                         gphiz_fab(i,j,k) = (phi_fab(i,j,k+1) - phi_fab(i,j,k-1) ) / (2.*dx[2]););
        });
    }
}

// Computes gradient at cell centres from face centred data
// Outputs to 3 different components in a cell-centered MultiFab
void ComputeCentredGradFC(std::array<MultiFab, AMREX_SPACEDIM> & phi,
                          MultiFab & gphi,
                          const Geometry & geom)
{
    BL_PROFILE_VAR("ComputeCentredGradFC()",ComputeCentredGradFC);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(gphi,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.tilebox();

        Array4<Real> const& gphi_fab = gphi.array(mfi);

        AMREX_D_TERM(Array4<Real const> const& phix_fab = phi[0].array(mfi);,
                     Array4<Real const> const& phiy_fab = phi[1].array(mfi);,
                     Array4<Real const> const& phiz_fab = phi[2].array(mfi););

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            AMREX_D_TERM(gphi_fab(i,j,k,0) = (phix_fab(i+1,j,k) - phix_fab(i,j,k) ) / dx[0];,
                         gphi_fab(i,j,k,1) = (phiy_fab(i,j+1,k) - phiy_fab(i,j,k) ) / dx[1];,
                         gphi_fab(i,j,k,2) = (phiz_fab(i,j,k+1) - phiz_fab(i,j,k) ) / dx[2];);
        });
    }
}
