#include "FPU.H"
#include "common_functions.H"

void calculateFlux(const MultiFab& cu,
                   const MultiFab& cumom,
                   MultiFab& faceflux,
                   MultiFab& cenflux,
                   const MultiFab& stochface,
                   const MultiFab& stochcen,
                   const Geometry& geom)
{
    BL_PROFILE_VAR("calculateFlux()", calculateFlux);

    faceflux.setVal(0.0);
    cenflux.setVal(0.0);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    const Real dxinv = 1.0 / dx[0];
    const Real half_dxinv = 0.5 * dxinv;

    const auto Acoef = FPU::A;
    const auto Dcoef = FPU::D;
    const auto Bcoef = FPU::B;
    const int use_noise = FPU::enable_fluctuations;

    for (MFIter mfi(cu, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& xbx = mfi.nodaltilebox(0);
        const Box& bx  = mfi.tilebox();

        Array4<Real const> const& u  = cu.const_array(mfi);
        Array4<Real const> const& mx = cumom.const_array(mfi);

        Array4<Real> const& ff = faceflux.array(mfi);
        Array4<Real> const& cf = cenflux.array(mfi);

        Array4<Real const> const& sf = stochface.const_array(mfi);
        Array4<Real const> const& sc = stochcen.const_array(mfi);

        amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real u0 = 0.5 * (u(i,j,k,0) + u(i-1,j,k,0));
            const Real u1 = mx(i,j,k,0);
            const Real u2 = 0.5 * (u(i,j,k,2) + u(i-1,j,k,2));

            const Real dudx0 = (u(i,j,k,0) - u(i-1,j,k,0)) * dxinv;
            const Real dudx1 = (mx(i+1,j,k,0) - mx(i-1,j,k,0)) * half_dxinv;
            const Real dudx2 = (u(i,j,k,2) - u(i-1,j,k,2)) * dxinv;

            ff(i,j,k,0) =
                Acoef[0]*u0 + Acoef[1]*u1 + Acoef[2]*u2
              - Dcoef[0]*dudx0 - Dcoef[1]*dudx1 - Dcoef[2]*dudx2;

            ff(i,j,k,1) =
                Acoef[6]*u0 + Acoef[7]*u1 + Acoef[8]*u2
              - Dcoef[6]*dudx0 - Dcoef[7]*dudx1 - Dcoef[8]*dudx2;

            if (use_noise) {
                ff(i,j,k,0) += Bcoef[0] * sf(i,j,k,0);
                ff(i,j,k,1) += Bcoef[2] * sf(i,j,k,1);
            }
        });

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real u0 = u(i,j,k,0);
            const Real u1 = 0.5 * (mx(i,j,k,0) + mx(i+1,j,k,0));
            const Real u2 = u(i,j,k,2);

            const Real dudx0 = (u(i+1,j,k,0) - u(i-1,j,k,0)) * half_dxinv;
            const Real dudx1 = (mx(i+1,j,k,0) - mx(i,j,k,0)) * dxinv;
            const Real dudx2 = (u(i+1,j,k,2) - u(i-1,j,k,2)) * half_dxinv;

            cf(i,j,k,0) =
                Acoef[3]*u0 + Acoef[4]*u1 + Acoef[5]*u2
              - Dcoef[3]*dudx0 - Dcoef[4]*dudx1 - Dcoef[5]*dudx2;

            if (use_noise) {
                cf(i,j,k,0) += Bcoef[1] * sc(i,j,k,0);
            }
        });
    }
}
