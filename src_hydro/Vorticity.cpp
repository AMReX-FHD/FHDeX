#include "hydro_functions.H"

void MagVort(const std::array< MultiFab, AMREX_SPACEDIM >& umac,
             MultiFab& magvort,
             const Geometry& geom,
             int outcomp) {

    BL_PROFILE_VAR("Vorticity()",Vorticity);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(magvort,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.tilebox();

        const Array4<Real> & vort = magvort.array(mfi);
        AMREX_D_TERM(Array4<Real const> const& u = umac[0].array(mfi);,
                     Array4<Real const> const& v = umac[1].array(mfi);,
                     Array4<Real const> const& w = umac[2].array(mfi););

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
#if (AMREX_SPACEDIM == 2)
            vort(i,j,k,outcomp) = (v(i+1,j,k) + v(i+1,j+1,k) - v(i-1,j,k) - v(i-1,j+1,k)) / (4.*dx[0]) +
                                  (u(i,j+1,k) + u(i+1,j+1,k) - u(i,j-1,k) - u(i+1,j-1,k)) / (4.*dx[1]);
#else
            Real x = (w(i,j+1,k) + w(i,j+1,k+1) - w(i,j-1,k) - w(i,j-1,k+1)) / (4.*dx[1]) +
                     (v(i,j,k+1) + v(i,j+1,k+1) - v(i,j,k-1) - v(i,j+1,k-1)) / (4.*dx[2]);
            Real y = (u(i,j,k+1) + u(i+1,j,k+1) - u(i,j,k-1) - u(i+1,j,k-1)) / (4.*dx[2]) +
                     (w(i+1,j,k) + w(i+1,j,k+1) - w(i-1,j,k) - w(i-1,j,k+1)) / (4.*dx[0]);
            Real z = (v(i+1,j,k) + v(i+1,j+1,k) - v(i-1,j,k) - v(i-1,j+1,k)) / (4.*dx[0]) +
                     (u(i,j+1,k) + u(i+1,j+1,k) - u(i,j-1,k) - u(i+1,j-1,k)) / (4.*dx[1]);

            vort(i,j,k,outcomp) = std::sqrt(x*x + y*y + z*z);
#endif
        });
    }
}
