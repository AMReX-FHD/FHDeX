#include "common_functions.H"
#include "common_functions_F.H"


//Computes divergence at cell centres from velcocities at cell faces
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
            AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, n,
            {
                div_fab(i,j,k,start_outcomp+n) =
                    AMREX_D_TERM(  (phix_fab(i+1,j,k,start_incomp+n) - phix_fab(i,j,k,start_incomp+n)) / dx[0],
                                 + (phiy_fab(i,j+1,k,start_incomp+n) - phiy_fab(i,j,k,start_incomp+n)) / dx[1],
                                 + (phiz_fab(i,j,k+1,start_incomp+n) - phiz_fab(i,j,k,start_incomp+n)) / dx[2]);;
            });
        }
        else
        {
            AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, n,
            {
                div_fab(i,j,k,start_outcomp+n) +=
                    AMREX_D_TERM(  (phix_fab(i+1,j,k,start_incomp+n) - phix_fab(i,j,k,start_incomp+n)) / dx[0],
                                 + (phiy_fab(i,j+1,k,start_incomp+n) - phiy_fab(i,j,k,start_incomp+n)) / dx[1],
                                 + (phiz_fab(i,j,k+1,start_incomp+n) - phiz_fab(i,j,k,start_incomp+n)) / dx[2]);;
            });
        }
    }
}


//Kernel for FC Grad
AMREX_GPU_HOST_DEVICE
inline
void compute_grad (const Box & tbx,
                   AMREX_D_DECL(const Box & xbx,
                                const Box & ybx,
                                const Box & zbx),
                   AMREX_D_DECL(const Array4<Real> & gx,
                                const Array4<Real> & gy,
                                const Array4<Real> & gz),
                   const Array4<Real const> & phi,
                   const GpuArray<Real, AMREX_SPACEDIM>& dx,
                   int start_incomp, int start_outcomp, int ncomp) noexcept
{

    // xbx, ybx, and zbx are the face-centered boxes

    // if running on the host: tlo is the minimal box contains the union of the
    // face-centered grid boxes

    // if running on the gpu: tlo is a box with a single point that comes from
    // the union of the face-centered grid boxes

    const auto tlo = lbound(tbx);
    const auto thi = ubound(tbx);

    // if running on the host, x/y/zlo and x/y/zhi are set to the lower/upper
    // bounds of x/y/zbx

    // if running on the gpu, x/y/zlo and x/y/zhi are set to the single point
    // defined by tlo, unless tlo is outside of the union of the face-centered
    // grid boxes, in which case they are set to values that make sure the loop
    // is not entered

    AMREX_D_TERM(const auto xlo = amrex::elemwiseMax(tlo, lbound(xbx));,
                 const auto ylo = amrex::elemwiseMax(tlo, lbound(ybx));,
                 const auto zlo = amrex::elemwiseMax(tlo, lbound(zbx)););

    AMREX_D_TERM(const auto xhi = amrex::elemwiseMin(thi, ubound(xbx));,
                 const auto yhi = amrex::elemwiseMin(thi, ubound(ybx));,
                 const auto zhi = amrex::elemwiseMin(thi, ubound(zbx)););

    for (int n=0; n<ncomp; ++n) {
        for (int k=xlo.z; k<=xhi.z; ++k) {
            for (int j=xlo.y; j<=xhi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i=xlo.x; i<=xhi.x; ++i) {
                    gx(i, j, k, start_outcomp + n) = (phi(i, j, k, start_incomp + n)
                            - phi(i-1, j, k, start_incomp + n) ) / dx[0];
                }
            }
        }
    }

#if (AMREX_SPACEDIM >= 2)
    for (int n = 0; n < ncomp; ++n) {
        for (int k = ylo.z; k <= yhi.z; ++k) {
            for (int j = ylo.y; j <= yhi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = ylo.x; i <= yhi.x; ++i) {
                    gy(i, j, k, start_outcomp + n) = (phi(i, j, k, start_incomp + n)
                            - phi(i, j-1, k, start_incomp + n) ) / dx[1];
                }
            }
        }
    }
#endif

#if (AMREX_SPACEDIM == 3)
    for (int n = 0; n < ncomp; ++n) {
        for (int k = zlo.z; k <= zhi.z; ++k) {
            for (int j = zlo.y; j <= zhi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = zlo.x; i <= zhi.x; ++i) {
                    gz(i, j, k, start_outcomp + n) = (phi(i, j, k, start_incomp + n)
                            - phi(i, j, k-1, start_incomp + n) ) / dx[2];
                }
            }
        }
    }
#endif
}


//Computes gradient at cell faces of cell centred scalar
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

        const Box& index_bounds = amrex::getIndexBounds(AMREX_D_DECL(bx_x, bx_y, bx_z));

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA(index_bounds, tbx,
        {
            compute_grad(tbx, AMREX_D_DECL(bx_x, bx_y, bx_z),
                         AMREX_D_DECL(gphix_fab, gphiy_fab, gphiz_fab),
                         phi_fab, dx,
                         start_incomp, start_outcomp, ncomp);

        });
    }
}


//Computes gradient at cell centres from cell centred data - ouputs to a three
//component mf.
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

        AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
        {
            AMREX_D_TERM(gphix_fab(i,j,k) = (phi_fab(i+1,j,k) - phi_fab(i-1,j,k) ) / (2.*dx[0]);,
                         gphiy_fab(i,j,k) = (phi_fab(i,j+1,k) - phi_fab(i,j-1,k) ) / (2.*dx[1]);,
                         gphiz_fab(i,j,k) = (phi_fab(i,j,k+1) - phi_fab(i,j,k-1) ) / (2.*dx[2]););
        });
    }
}
