#include "common_functions.H"


// Computes divergence at cell centres from velcocities at cell faces
void ComputeDiv(MultiFab& div,
                const std::array<MultiFab, AMREX_SPACEDIM>& phi_fc,
                int start_incomp, int start_outcomp, int ncomp,
                const Geometry& geom, int increment)
{
    BL_PROFILE_VAR("ComputeDiv()",ComputeDiv);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // if not incrementing, initialize data to zero
    if (increment == 0) {
        div.setVal(0.,start_outcomp,ncomp,0);
    }

    for ( MFIter mfi(div,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.tilebox();

        const Array4<Real> & div_fab = div.array(mfi);
        AMREX_D_TERM(Array4<Real const> const& phix_fab = phi_fc[0].array(mfi);,
                     Array4<Real const> const& phiy_fab = phi_fc[1].array(mfi);,
                     Array4<Real const> const& phiz_fab = phi_fc[2].array(mfi););

        amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            div_fab(i,j,k,start_outcomp+n) +=
                AMREX_D_TERM(  (phix_fab(i+1,j,k,start_incomp+n) - phix_fab(i,j,k,start_incomp+n)) / dx[0],
                             + (phiy_fab(i,j+1,k,start_incomp+n) - phiy_fab(i,j,k,start_incomp+n)) / dx[1],
                             + (phiz_fab(i,j,k+1,start_incomp+n) - phiz_fab(i,j,k,start_incomp+n)) / dx[2]);;
        });
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

    // if not incrementing, initialize data to zero
    if (increment == 0) {
        for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
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
        if (bc_lo[0] == amrex::BCType::foextrap || bc_lo[0] == amrex::BCType::ext_dir) {
            if (bx_x.smallEnd(0) <= dom.smallEnd(0)) {
                int lo = dom.smallEnd(0);
                amrex::ParallelFor(bx_x, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (i == lo) {
                        // subtract off incorrect gradient from above, then add the correct gradient
                        gphix(i,j,k,start_outcomp+n) -= (phi(i,j,k,start_incomp+n)-phi(i-1,j,k,start_incomp+n)) / dx[0];
                        gphix(i,j,k,start_outcomp+n) += (phi(i,j,k,start_incomp+n)-phi(i-1,j,k,start_incomp+n)) / (0.5*dx[0]);
                    }
                });
            }
        }

        if (bc_hi[0] == amrex::BCType::foextrap || bc_hi[0] == amrex::BCType::ext_dir) {
            if (bx_x.bigEnd(0) >= dom.bigEnd(0)+1) {
                int hi = dom.bigEnd(0)+1;
                amrex::ParallelFor(bx_x, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (i == hi) {
                        // subtract off incorrect gradient from above, then add the correct gradient
                        gphix(i,j,k,start_outcomp+n) -= (phi(i,j,k,start_incomp+n)-phi(i-1,j,k,start_incomp+n)) / dx[0];
                        gphix(i,j,k,start_outcomp+n) += (phi(i,j,k,start_incomp+n)-phi(i-1,j,k,start_incomp+n)) / (0.5*dx[0]);
                    }
                });
            }
        }

        if (bc_lo[1] == amrex::BCType::foextrap || bc_lo[1] == amrex::BCType::ext_dir) {
            if (bx_y.smallEnd(1) <= dom.smallEnd(1)) {
                int lo = dom.smallEnd(1);
                amrex::ParallelFor(bx_y, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (j == lo) {
                        // subtract off incorrect gradient from above, then add the correct gradient
                        gphiy(i,j,k,start_outcomp+n) -= (phi(i,j,k,start_incomp+n)-phi(i,j-1,k,start_incomp+n)) / dx[1];
                        gphiy(i,j,k,start_outcomp+n) += (phi(i,j,k,start_incomp+n)-phi(i,j-1,k,start_incomp+n)) / (0.5*dx[1]);
                    }
                });
            }
        }

        if (bc_hi[1] == amrex::BCType::foextrap || bc_hi[1] == amrex::BCType::ext_dir) {
            if (bx_y.bigEnd(1) >= dom.bigEnd(1)+1) {
                int hi = dom.bigEnd(1)+1;
                amrex::ParallelFor(bx_y, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (j == hi) {
                        // subtract off incorrect gradient from above, then add the correct gradient
                        gphiy(i,j,k,start_outcomp+n) -= (phi(i,j,k,start_incomp+n)-phi(i,j-1,k,start_incomp+n)) / dx[1];
                        gphiy(i,j,k,start_outcomp+n) += (phi(i,j,k,start_incomp+n)-phi(i,j-1,k,start_incomp+n)) / (0.5*dx[1]);
                    }
                });
            }
        }

#if (AMREX_SPACEDIM == 3)
        if (bc_lo[2] == amrex::BCType::foextrap || bc_lo[2] == amrex::BCType::ext_dir) {
            if (bx_z.smallEnd(2) <= dom.smallEnd(2)) {
                int lo = dom.smallEnd(2);
                amrex::ParallelFor(bx_z, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (k == lo) {
                        // subtract off incorrect gradient from above, then add the correct gradient
                        gphiz(i,j,k,start_outcomp+n) -= (phi(i,j,k,start_incomp+n)-phi(i,j,k-1,start_incomp+n)) / dx[2];
                        gphiz(i,j,k,start_outcomp+n) += (phi(i,j,k,start_incomp+n)-phi(i,j,k-1,start_incomp+n)) / (0.5*dx[2]);
                    }
                });
            }
        }

        if (bc_hi[2] == amrex::BCType::foextrap || bc_hi[2] == amrex::BCType::ext_dir) {
            if (bx_z.bigEnd(2) >= dom.bigEnd(2)+1) {
                int hi = dom.bigEnd(2)+1;
                amrex::ParallelFor(bx_z, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (k == hi) {
                        // subtract off incorrect gradient from above, then add the correct gradient
                        gphiz(i,j,k,start_outcomp+n) -= (phi(i,j,k,start_incomp+n)-phi(i,j,k-1,start_incomp+n)) / dx[2];
                        gphiz(i,j,k,start_outcomp+n) += (phi(i,j,k,start_incomp+n)-phi(i,j,k-1,start_incomp+n)) / (0.5*dx[2]);
                    }
                });
            }
        }
#endif
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

// Computes gradient at cell centres on a component in a given direction
void ComputeCentredGradCompDir(const MultiFab & phi,
                               MultiFab& gphi,
                               int dir,
                               int incomp,
                               int outcomp,
                               const Geometry & geom)
{
    BL_PROFILE_VAR("ComputeCentredGradCompDir()",ComputeCentredGrad);

    int ioff = (dir == 0) ? 1 : 0;
    int joff = (dir == 1) ? 1 : 0;
    int koff = (dir == 2) ? 1 : 0;

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(phi,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.tilebox();

        Array4<Real const> const& phi_fab = phi.array(mfi);
        Array4<Real> const& gphi_fab = gphi.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            gphi_fab(i,j,k,outcomp) = (phi_fab(i+ioff,j+joff,k+koff,incomp) - phi_fab(i-ioff,j-joff,k-koff,incomp)) / (2.*dx[dir]);
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

// Computes Laplacian at cell centres on a component in a given direction
void ComputeLap(const MultiFab & phi_in,
                MultiFab& Lphi_in,
                int incomp,
                int outcomp,
                int numcomp,
                const Geometry & geom)
{
    BL_PROFILE_VAR("ComputeLap()",ComputeLap);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(phi_in,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.tilebox();

        Array4<Real const> const& phi = phi_in.array(mfi);
        Array4<Real> const& Lphi = Lphi_in.array(mfi);

        amrex::ParallelFor(bx, numcomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Lphi(i,j,k,outcomp+n) =
                  ( phi(i-1,j,k,incomp+n) - 2.*phi(i,j,k,incomp+n) + phi(i+1,j,k,incomp+n) ) / (dx[0]*dx[0])
                + ( phi(i,j-1,k,incomp+n) - 2.*phi(i,j,k,incomp+n) + phi(i,j+1,k,incomp+n) ) / (dx[1]*dx[1])
#if (AMREX_SPACEDIM == 3)
                + ( phi(i,j,k-1,incomp+n) - 2.*phi(i,j,k,incomp+n) + phi(i,j,k+1,incomp+n) ) / (dx[2]*dx[2])
#endif
                ;
        });
    }
}

void ComputeStagLap(std::array<MultiFab, AMREX_SPACEDIM> & phi_in,
                    std::array<MultiFab, AMREX_SPACEDIM> & Lphi_in,
                    const Geometry & geom)
{
    BL_PROFILE_VAR("ComputeStagLap()",ComputeStagLap);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(phi_in[0],TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        AMREX_D_TERM(Box bx_x = mfi.tilebox(nodal_flag_x);,
                     Box bx_y = mfi.tilebox(nodal_flag_y);,
                     Box bx_z = mfi.tilebox(nodal_flag_z););

        AMREX_D_TERM(Array4<Real const> const& phix = phi_in[0].array(mfi);,
                     Array4<Real const> const& phiy = phi_in[1].array(mfi);,
                     Array4<Real const> const& phiz = phi_in[2].array(mfi););

        AMREX_D_TERM(Array4<Real> const& Lphix = Lphi_in[0].array(mfi);,
                     Array4<Real> const& Lphiy = Lphi_in[1].array(mfi);,
                     Array4<Real> const& Lphiz = Lphi_in[2].array(mfi););

        amrex::ParallelFor(bx_x, bx_y,
#if (AMREX_SPACEDIM == 3)
                           bx_z,
#endif
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Lphix(i,j,k) =    (phix(i+1,j,k) - 2.*phix(i,j,k) + phix(i-1,j,k)) / (dx[0]*dx[0])
                                + (phix(i,j+1,k) - 2.*phix(i,j,k) + phix(i,j-1,k)) / (dx[1]*dx[1])
#if (AMREX_SPACEDIM == 3)
                                + (phix(i,j,k+1) - 2.*phix(i,j,k) + phix(i,j,k-1)) / (dx[2]*dx[2])
#endif
                    ;

            },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Lphiy(i,j,k) =    (phiy(i+1,j,k) - 2.*phiy(i,j,k) + phiy(i-1,j,k)) / (dx[0]*dx[0])
                                + (phiy(i,j+1,k) - 2.*phiy(i,j,k) + phiy(i,j-1,k)) / (dx[1]*dx[1])
#if (AMREX_SPACEDIM == 3)
                                + (phiy(i,j,k+1) - 2.*phiy(i,j,k) + phiy(i,j,k-1)) / (dx[2]*dx[2])
#endif
                    ;

            }
#if (AMREX_SPACEDIM == 3)
                           , [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Lphiz(i,j,k) =    (phiz(i+1,j,k) - 2.*phiz(i,j,k) + phiz(i-1,j,k)) / (dx[0]*dx[0])
                                + (phiz(i,j+1,k) - 2.*phiz(i,j,k) + phiz(i,j-1,k)) / (dx[1]*dx[1])
                                + (phiz(i,j,k+1) - 2.*phiz(i,j,k) + phiz(i,j,k-1)) / (dx[2]*dx[2]);
            }
#endif
            );

    }
}

void ComputeCurlFaceToEdge(std::array<MultiFab, AMREX_SPACEDIM> & umac_in,
                           std::array<MultiFab, NUM_EDGE> & curl,
                           const Geometry & geom)
{
    BL_PROFILE_VAR("ComputeCurlFaceToEdge()",ComputeCurlFaceToEdge);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(umac_in[0],TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        AMREX_D_TERM(Array4<Real const> const& umac = umac_in[0].array(mfi);,
                     Array4<Real const> const& vmac = umac_in[1].array(mfi);,
                     Array4<Real const> const& wmac = umac_in[2].array(mfi););

#if (AMREX_SPACEDIM == 3)
        AMREX_D_TERM(Box bx_xy = mfi.tilebox(nodal_flag_xy);,
                     Box bx_xz = mfi.tilebox(nodal_flag_xz);,
                     Box bx_yz = mfi.tilebox(nodal_flag_yz););

        AMREX_D_TERM(Array4<Real> const& curlxy = curl[0].array(mfi);,
                     Array4<Real> const& curlxz = curl[1].array(mfi);,
                     Array4<Real> const& curlyz = curl[2].array(mfi););

        amrex::ParallelFor(bx_xy, bx_xz, bx_yz,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // dv/dx - du/dy
                curlxy(i,j,k) = (vmac(i,j,k) - vmac(i-1,j,k))/dx[0] - (umac(i,j,k)-umac(i,j-1,k))/dx[1];
            },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // du/dz - dw/dx
                curlxz(i,j,k) = (umac(i,j,k)-umac(i,j,k-1))/dx[2] - (wmac(i,j,k)-wmac(i-1,j,k))/dx[0];
            },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // dw/dy - dv/dz
                curlyz(i,j,k) = (wmac(i,j,k)-wmac(i,j-1,k))/dx[1] - (vmac(i,j,k)-vmac(i,j,k-1))/dx[2];
            });
#elif (AMREX_SPACEDIM == 2)
        Box bx_xy = mfi.tilebox(nodal_flag_xy);

        Array4<Real> const& curlxy = curl[0].array(mfi);

        amrex::ParallelFor(bx_xy,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // dv/dx - du/dy
                curlxy(i,j,k) = (vmac(i,j,k) - vmac(i-1,j,k))/dx[0] - (umac(i,j,k)-umac(i,j-1,k))/dx[1];
            });
#endif

    }
}

void ComputeCurlCC(const MultiFab& vel_in,
                   int incomp,
                   MultiFab& curl_in,
                   int outcomp,
                   const Geometry & geom)
{
    BL_PROFILE_VAR("ComputeCurlCC()",ComputeCurlCC);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(vel_in,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.tilebox();

        Array4<Real const> const& vel  = vel_in.array(mfi);
        Array4<Real>       const& curl = curl_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // dw/dy - dv/dz
            curl(i,j,k,outcomp) =
                (vel(i,j+1,k,incomp+2) - vel(i,j-1,k,incomp+2)) / (2.*dx[1]) -
                (vel(i,j,k+1,incomp+1) - vel(i,j,k-1,incomp+1)) / (2.*dx[2]);

            // du/dz - dw/dx
            curl(i,j,k,outcomp+1) =
                (vel(i,j,k+1,incomp+0) - vel(i,j,k-1,incomp+0)) / (2.*dx[2]) -
                (vel(i+1,j,k,incomp+2) - vel(i-1,j,k,incomp+2)) / (2.*dx[0]);

            // dv/dx - du/dy
            curl(i,j,k,outcomp+2) =
                (vel(i+1,j,k,incomp+1) - vel(i-1,j,k,incomp+1)) / (2.*dx[0]) -
                (vel(i,j+1,k,incomp+0) - vel(i,j-1,k,incomp+0)) / (2.*dx[1]);
        });
    }
}

void ComputeDivCC(const MultiFab& vel_in,
                   int incomp,
                   MultiFab& div_in,
                   int outcomp,
                   const Geometry & geom)
{
    BL_PROFILE_VAR("ComputeDivCC()",ComputeCurlCC);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(vel_in,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.tilebox();

        Array4<Real const> const& vel  = vel_in.array(mfi);
        Array4<Real>       const& div  = div_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // dw/dy - dv/dz
            div(i,j,k,outcomp) =
                (vel(i+1,j,k,incomp+0) - vel(i-1,j,k,incomp+0)) / (2.*dx[0]) +
                (vel(i,j+1,k,incomp+1) - vel(i,j-1,k,incomp+1)) / (2.*dx[1]) +
                (vel(i,j,k+1,incomp+2) - vel(i,j,k-1,incomp+2)) / (2.*dx[2]);
        });
    }
}
