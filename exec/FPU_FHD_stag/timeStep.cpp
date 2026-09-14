#include "FPU.H"
#include "common_functions.H"
#include "rng_functions.H"

namespace {

void refreshCellCenteredMomentum(MultiFab& cu,
                                  const MultiFab& cumom,
                                  const Geometry& geom)
{
    for (MFIter mfi(cu, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& u = cu.array(mfi);
        Array4<Real const> const& mx = cumom.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            u(i,j,k,1) = 0.5 * (mx(i,j,k,0) + mx(i+1,j,k,0));
        });
    }
    cu.FillBoundary(geom.periodicity());
}

void fillWeightedNoise(MultiFab& stochface,
                       MultiFab& stochcen,
                       const MultiFab& stochface_A,
                       const MultiFab& stochface_B,
                       const MultiFab& stochcen_A,
                       const MultiFab& stochcen_B,
                       Real swgt1,
                       Real swgt2)
{
    MultiFab::LinComb(stochface,
                      swgt1, stochface_A, 0,
                      swgt2, stochface_B, 0,
                      0, 2, 0);

    MultiFab::LinComb(stochcen,
                      swgt1, stochcen_A, 0,
                      swgt2, stochcen_B, 0,
                      0, 1, 0);
}

void updateCellState(MultiFab& dst,
                     const MultiFab& base,
                     const MultiFab& stage,
                     const MultiFab& faceflux,
                     const Geometry& geom,
                     Real dt,
                     Real base_weight,
                     Real stage_weight,
                     Real rhs_weight)
{
    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    const Real dtdx = dt / dx[0];

    for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& dst_fab = dst.array(mfi);
        Array4<Real const> const& base_fab = base.const_array(mfi);
        Array4<Real const> const& stage_fab = stage.const_array(mfi);
        Array4<Real const> const& xflux = faceflux.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            dst_fab(i,j,k,0) = base_weight * base_fab(i,j,k,0)
                             + stage_weight * (stage_fab(i,j,k,0)
                             - rhs_weight * dtdx * (xflux(i+1,j,k,0) - xflux(i,j,k,0)));

            dst_fab(i,j,k,2) = base_weight * base_fab(i,j,k,2)
                             + stage_weight * (stage_fab(i,j,k,2)
                             - rhs_weight * dtdx * (xflux(i+1,j,k,1) - xflux(i,j,k,1)));
        });
    }
}

void updateFaceMomentum(MultiFab& dstmom,
                        const MultiFab& basemom,
                        const MultiFab& stagemom,
                        const MultiFab& cenflux,
                        const Geometry& geom,
                        Real dt,
                        Real base_weight,
                        Real stage_weight,
                        Real rhs_weight)
{
    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    const Real dtdx = dt / dx[0];

    for (MFIter mfi(dstmom, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& xbx = mfi.nodaltilebox(0);
        Array4<Real> const& dst = dstmom.array(mfi);
        Array4<Real const> const& base = basemom.const_array(mfi);
        Array4<Real const> const& stage = stagemom.const_array(mfi);
        Array4<Real const> const& cenx = cenflux.const_array(mfi);

        amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            dst(i,j,k,0) = base_weight * base(i,j,k,0)
                         + stage_weight * (stage(i,j,k,0)
                         - rhs_weight * dtdx * (cenx(i,j,k,0) - cenx(i-1,j,k,0)));
        });
    }
}

}

void RK3step(MultiFab& cu,
             MultiFab& cumom,
             MultiFab& faceflux,
             MultiFab& cenflux,
             const Geometry& geom,
             const Real dt,
             const int /*step*/)
{
    BL_PROFILE_VAR("RK3step()", RK3step);

    MultiFab cup(cu.boxArray(), cu.DistributionMap(), 3, ngc);
    MultiFab cup2(cu.boxArray(), cu.DistributionMap(), 3, ngc);
    cup.setVal(0.0, 0, 3, ngc);
    cup2.setVal(0.0, 0, 3, ngc);

    MultiFab cupmom(convert(cu.boxArray(), nodal_flag_x), cu.DistributionMap(), 1, ngc);
    MultiFab cup2mom(convert(cu.boxArray(), nodal_flag_x), cu.DistributionMap(), 1, ngc);
    cupmom.setVal(0.0, 0, 1, ngc);
    cup2mom.setVal(0.0, 0, 1, ngc);

    MultiFab stochface(convert(cu.boxArray(), nodal_flag_x), cu.DistributionMap(), 2, 0);
    MultiFab stochcen(cu.boxArray(), cu.DistributionMap(), 1, 0);

    MultiFab stochface_A(stochface.boxArray(), stochface.DistributionMap(), 2, 0);
    MultiFab stochface_B(stochface.boxArray(), stochface.DistributionMap(), 2, 0);
    MultiFab stochcen_A(stochcen.boxArray(), stochcen.DistributionMap(), 1, 0);
    MultiFab stochcen_B(stochcen.boxArray(), stochcen.DistributionMap(), 1, 0);

    stochface_A.setVal(0.0);
    stochface_B.setVal(0.0);
    stochcen_A.setVal(0.0);
    stochcen_B.setVal(0.0);

    if (FPU::enable_fluctuations) {
        const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
#if (AMREX_SPACEDIM == 2)
        const Real dV = dx[0] * dx[1] * FPU::cell_dz;
#else
        const Real dV = dx[0] * dx[1] * dx[2];
#endif
        const Real noise_scale = 1.0 / std::sqrt(dt * dV);

        MultiFabFillRandomNormal(stochface_A, 0, 2, 0.0, 1.0, geom, true, true);
        MultiFabFillRandomNormal(stochface_B, 0, 2, 0.0, 1.0, geom, true, true);
        MultiFabFillRandomNormal(stochcen_A, 0, 1, 0.0, 1.0, geom, true, true);
        MultiFabFillRandomNormal(stochcen_B, 0, 1, 0.0, 1.0, geom, true, true);

        stochface_A.mult(noise_scale, 0, 2, 0);
        stochface_B.mult(noise_scale, 0, 2, 0);
        stochcen_A.mult(noise_scale, 0, 1, 0);
        stochcen_B.mult(noise_scale, 0, 1, 0);
    }

    refreshCellCenteredMomentum(cu, cumom, geom);
    cu.FillBoundary(geom.periodicity());
    cumom.FillBoundary(geom.periodicity());

    fillWeightedNoise(stochface, stochcen,
                      stochface_A, stochface_B,
                      stochcen_A, stochcen_B,
                      1.0,
                      (2.0*std::sqrt(2.0) + std::sqrt(3.0)) / 5.0);

    calculateFlux(cu, cumom, faceflux, cenflux, stochface, stochcen, geom);
    cenflux.FillBoundary(geom.periodicity());
    updateCellState(cup, cu, cu, faceflux, geom, dt, 0.0, 1.0, 1.0);
    updateFaceMomentum(cupmom, cumom, cumom, cenflux, geom, dt, 0.0, 1.0, 1.0);
    cupmom.FillBoundary(geom.periodicity());
    refreshCellCenteredMomentum(cup, cupmom, geom);
    cup.FillBoundary(geom.periodicity());

    fillWeightedNoise(stochface, stochcen,
                      stochface_A, stochface_B,
                      stochcen_A, stochcen_B,
                      1.0,
                      (-4.0*std::sqrt(2.0) + 3.0*std::sqrt(3.0)) / 5.0);

    calculateFlux(cup, cupmom, faceflux, cenflux, stochface, stochcen, geom);
    cenflux.FillBoundary(geom.periodicity());
    updateCellState(cup2, cu, cup, faceflux, geom, dt, 0.75, 0.25, 1.0);
    updateFaceMomentum(cup2mom, cumom, cupmom, cenflux, geom, dt, 0.75, 0.25, 1.0);
    cup2mom.FillBoundary(geom.periodicity());
    refreshCellCenteredMomentum(cup2, cup2mom, geom);
    cup2.FillBoundary(geom.periodicity());

    fillWeightedNoise(stochface, stochcen,
                      stochface_A, stochface_B,
                      stochcen_A, stochcen_B,
                      1.0,
                      (std::sqrt(2.0) - 2.0*std::sqrt(3.0)) / 10.0);

    calculateFlux(cup2, cup2mom, faceflux, cenflux, stochface, stochcen, geom);
    cenflux.FillBoundary(geom.periodicity());
    updateCellState(cu, cu, cup2, faceflux, geom, dt, 1.0/3.0, 2.0/3.0, 1.0);
    updateFaceMomentum(cumom, cumom, cup2mom, cenflux, geom, dt, 1.0/3.0, 2.0/3.0, 1.0);
    cumom.FillBoundary(geom.periodicity());
    refreshCellCenteredMomentum(cu, cumom, geom);
    cu.FillBoundary(geom.periodicity());
}
