#include <AMReX_MultiFabUtil.H>

#include "DensityConv.H"

using namespace amrex;

// In the 3D (SPDE) build the convolution is done on the x-y plane only: phi is
// first integrated over theta (ReduceToPlaneMF), the 2D FFT convolution is done
// on a one-cell-thick domain, and the result is copied into every theta plane.
// This costs one 2D FFT pair instead of a 3D one; since W is constant in theta
// only the theta-mean of phi contributes, so the result is the same.

void
DensityConv::define (Geometry const& geom, BoxArray const& ba, DistributionMapping const& dm,
                     FlockingParams const& fp)
{
    m_geom = geom;
    m_ba = ba;
    m_dm = dm;
    const Box& domain = geom.Domain();

    // the x-y plane the convolution is done on
    Box plane = domain;
#if (AMREX_SPACEDIM == 3)
    plane.setRange(2, 0, 1);
    BoxArray pba(plane);
    IntVect mgs = ba.minimalBox().length();
    for (int i = 0; i < ba.size(); ++i) { mgs.min(ba[i].length()); }
    mgs[2] = 1;
    pba.maxSize(mgs);
    DistributionMapping pdm(pba);
    // the plane layout matching the boxes of ba (each box flattened to z = 0)
    BoxList bl;
    for (int i = 0; i < ba.size(); ++i) { Box b = ba[i]; b.setRange(2, 0, 1); bl.push_back(b); }
    m_flat_ba = BoxArray(std::move(bl));
#else
    BoxArray pba = ba;
    DistributionMapping pdm = dm;
#endif
    m_plane_ba = pba;
    m_plane_dm = pdm;

    m_fwd = std::make_unique<FFT::R2C<Real,FFT::Direction::forward>>(plane);
    m_bwd = std::make_unique<FFT::R2C<Real,FFT::Direction::backward>>(plane);
    auto const& [cba, cdm] = m_fwd->getSpectralDataLayout();
    m_what.define(cba, cdm, 1, 0);
    m_shat.define(cba, cdm, 1, 0);
    m_tmp.define(pba, pdm, 1, 0);

    // W at the cell offsets of the plane, centred at cell 0 with periodic wrap
    // (as init_int_pot in interaction)
    const auto dx = geom.CellSizeArray();
    const Real lx = geom.ProbLength(0);
    const Real ly = geom.ProbLength(1);
    const FlockingParams fp_loc = fp;
    MultiFab wk(pba, pdm, 1, 0);
    for (MFIter mfi(wk); mfi.isValid(); ++mfi) {
        auto const& w = wk.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = i*dx[0];
            Real y = j*dx[1];
            if (x + 0.5*dx[0] > 0.5*lx) { x -= lx; }
            if (y + 0.5*dx[1] > 0.5*ly) { y -= ly; }
            w(i,j,k) = flock_kernel_shape(std::sqrt(x*x + y*y), fp_loc);
        });
    }

    // normalize: sum over the plane of W dx dy = 1
    const Real plane_sum = wk.sum(0)*dx[0]*dx[1];
    if (plane_sum <= 0.0) { Abort("DensityConv: the kernel has no support on the grid; increase flock.kernel_R"); }
    wk.mult(1.0/plane_sum, 0, 1);
    m_norm_ratio = plane_sum*flock_kernel_norm(fp);

    m_fwd->forward(wk, m_what);
}

void
DensityConv::apply (MultiFab const& src, MultiFab& dst)
{
    AMREX_ALWAYS_ASSERT(defined());
    const auto dx = m_geom.CellSizeArray();

#if (AMREX_SPACEDIM == 3)
    // rho(x, y) = sum_k phi(x, y, theta_k) dtheta, on the flattened boxes of src
    const Box& domain = m_geom.Domain();
    auto const& sa = src.const_arrays();
    MultiFab rho_flat = ReduceToPlaneMF<ReduceOpSum>(2, domain, src,
        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
        {
            return sa[b](i,j,k);
        });
    rho_flat.mult(dx[2], 0, 1);
    m_tmp.ParallelCopy(rho_flat, 0, 0, 1);
#else
    MultiFab::Copy(m_tmp, src, 0, 0, 1, 0);
#endif

    m_fwd->forward(m_tmp, m_shat);
    for (MFIter mfi(m_shat); mfi.isValid(); ++mfi) {
        auto const& w = m_what.const_array(mfi);
        auto const& s = m_shat.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            s(i,j,k) *= w(i,j,k);
        });
    }
    m_bwd->backward(m_shat, m_tmp);

    // FFTW is unnormalized: divide by the number of plane points, times dx dy
    const Real npts = static_cast<Real>(m_geom.Domain().length(0))*m_geom.Domain().length(1);
    m_tmp.mult(dx[0]*dx[1]/npts, 0, 1);

#if (AMREX_SPACEDIM == 3)
    // copy rho_tilde into every theta plane of dst
    MultiFab rt_flat(m_flat_ba, m_dm, 1, 0);
    rt_flat.ParallelCopy(m_tmp, 0, 0, 1);
    for (MFIter mfi(dst); mfi.isValid(); ++mfi) {
        auto const& d = dst.array(mfi);
        auto const& r = rt_flat.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            d(i,j,k) = r(i,j,0);
        });
    }
#else
    MultiFab::Copy(dst, m_tmp, 0, 0, 1, 0);
#endif
    dst.FillBoundary(m_geom.periodicity());
}
