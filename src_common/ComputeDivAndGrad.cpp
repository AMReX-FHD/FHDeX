#include "common_functions.H"


// Computes divergence at cell centres from velcocities at cell faces
void ComputeDiv(MultiFab& div,
                const std::array<MultiFab, AMREX_SPACEDIM>& phi_fc,
                int start_incomp, int start_outcomp, int ncomp,
                const Geometry& geom, int increment)
{
    BL_PROFILE_VAR("ComputeDiv()",ComputeDiv);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(div,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        
        const Box& bx = mfi.tilebox();

        const Array4<Real> & div_fab = div.array(mfi);
        AMREX_D_TERM(Array4<Real const> const& phix_fab = phi_fc[0].array(mfi);,
                     Array4<Real const> const& phiy_fab = phi_fc[1].array(mfi);,
                     Array4<Real const> const& phiz_fab = phi_fc[2].array(mfi););

        if (increment == 0) {
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
                div_fab(i,j,k,start_outcomp+n) +=
                    AMREX_D_TERM(  (phix_fab(i+1,j,k,start_incomp+n) - phix_fab(i,j,k,start_incomp+n)) / dx[0],
                                 + (phiy_fab(i,j+1,k,start_incomp+n) - phiy_fab(i,j,k,start_incomp+n)) / dx[1],
                                 + (phiz_fab(i,j,k+1,start_incomp+n) - phiz_fab(i,j,k,start_incomp+n)) / dx[2]);;
            });
        }
    }
}

// Computes gradient at cell faces of cell centred scalar
void ComputeGrad(const MultiFab & phi, std::array<MultiFab, AMREX_SPACEDIM> & gphi,
                 int start_incomp, int start_outcomp, int ncomp, const Geometry & geom)
{
    BL_PROFILE_VAR("ComputeGrad()",ComputeGrad);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(phi,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Array4<Real const> & phi_fab = phi.array(mfi);

        AMREX_D_TERM(const Array4<Real> & gphix_fab = gphi[0].array(mfi);,
                     const Array4<Real> & gphiy_fab = gphi[1].array(mfi);,
                     const Array4<Real> & gphiz_fab = gphi[2].array(mfi););

        AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                     const Box & bx_y = mfi.nodaltilebox(1);,
                     const Box & bx_z = mfi.nodaltilebox(2););

        amrex::ParallelFor(bx_x, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            gphix_fab(i,j,k,start_outcomp+n) = (phi_fab(i,j,k,start_incomp+n)-phi_fab(i-1,j,k,start_incomp+n)) / dx[0];
        },
                           bx_y, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            gphiy_fab(i,j,k,start_outcomp+n) = (phi_fab(i,j,k,start_incomp+n)-phi_fab(i,j-1,k,start_incomp+n)) / dx[1];
        }
#if (AMREX_SPACEDIM == 3)
        ,
                           bx_z, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            gphiz_fab(i,j,k,start_outcomp+n) = (phi_fab(i,j,k,start_incomp+n)-phi_fab(i,j,k-1,start_incomp+n)) / dx[2];
        }
#endif
        );
    }
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
