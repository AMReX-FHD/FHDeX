#include <AMReX_MultiFabUtil.H>

#include "AlignConv.H"

using namespace amrex;

// In the 3D (SPDE) build the convolution is done on the x-y plane only: the
// polarization p = int phi (cos, sin)(theta) dtheta is formed on the plane
// (ReduceToPlaneMF), the 2D FFT convolution is done on a one-cell-thick
// domain, and Q is copied into every theta plane. Since the kernel is
// H(theta - theta') G = -cos(theta) cos(theta') G - sin(theta) sin(theta') G,
// these two moments are all of phi that the alignment sees.

void
AlignConv::define (Geometry const& geom, BoxArray const& ba, DistributionMapping const& dm,
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
    m_ghat.define(cba, cdm, 1, 0);
    m_shat.define(cba, cdm, 1, 0);
    m_tmp.define(pba, pdm, 2, 0);
    m_one.define(pba, pdm, 1, 0);

    // G at the cell offsets of the plane, centred at cell 0 with periodic wrap
    // (as init_int_pot in interaction)
    const auto dx = geom.CellSizeArray();
    const Real lx = geom.ProbLength(0);
    const Real ly = geom.ProbLength(1);
    const Real rinv = 1.0/fp.align_R;
    MultiFab gk(pba, pdm, 1, 0);
    for (MFIter mfi(gk); mfi.isValid(); ++mfi) {
        auto const& g = gk.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = i*dx[0];
            Real y = j*dx[1];
            if (x + 0.5*dx[0] > 0.5*lx) { x -= lx; }
            if (y + 0.5*dx[1] > 0.5*ly) { y -= ly; }
            g(i,j,k) = flock_align_G(std::sqrt(x*x + y*y)*rinv);
        });
    }

    // rescale so that sum G dx dy is the continuous integral of G
    const Real plane_sum = gk.sum(0)*dx[0]*dx[1];
    if (plane_sum <= 0.0) { Abort("AlignConv: the alignment kernel has no support on the grid; increase flock.align_R"); }
    m_norm_ratio = plane_sum/flock_align_G_integral(fp);
    gk.mult(1.0/m_norm_ratio, 0, 1);

    m_fwd->forward(gk, m_ghat);
}

void
AlignConv::apply (MultiFab const& src, MultiFab& dst)
{
    AMREX_ALWAYS_ASSERT(defined());
    AMREX_ALWAYS_ASSERT(dst.nComp() >= 2);
    const auto dx = m_geom.CellSizeArray();

#if (AMREX_SPACEDIM == 3)
    // p(x, y) = sum_k phi(x, y, theta_k) (cos, sin)(theta_k) dtheta, on the
    // flattened boxes of src
    const Box& domain = m_geom.Domain();
    const Real dth  = dx[2];
    const Real thlo = m_geom.ProbLo(2);
    auto const& sa = src.const_arrays();
    for (int c = 0; c < 2; ++c) {
        MultiFab p_flat = ReduceToPlaneMF<ReduceOpSum>(2, domain, src,
            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
            {
                const Real th = thlo + (k+0.5)*dth;
                return sa[b](i,j,k)*((c == 0) ? std::cos(th) : std::sin(th));
            });
        p_flat.mult(dth, 0, 1);
        m_tmp.ParallelCopy(p_flat, 0, c, 1);
    }
#else
    AMREX_ALWAYS_ASSERT(src.nComp() >= 2);
    MultiFab::Copy(m_tmp, src, 0, 0, 2, 0);
#endif

    // FFTW is unnormalized: divide by the number of plane points, times dx dy
    const Real npts = static_cast<Real>(m_geom.Domain().length(0))*m_geom.Domain().length(1);
    for (int c = 0; c < 2; ++c) {
        MultiFab::Copy(m_one, m_tmp, c, 0, 1, 0);
        m_fwd->forward(m_one, m_shat);
        for (MFIter mfi(m_shat); mfi.isValid(); ++mfi) {
            auto const& g = m_ghat.const_array(mfi);
            auto const& s = m_shat.array(mfi);
            amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                s(i,j,k) *= g(i,j,k);
            });
        }
        m_bwd->backward(m_shat, m_one);
        m_one.mult(dx[0]*dx[1]/npts, 0, 1);
        MultiFab::Copy(m_tmp, m_one, 0, c, 1, 0);
    }

#if (AMREX_SPACEDIM == 3)
    // copy Q into every theta plane of dst
    MultiFab q_flat(m_flat_ba, m_dm, 2, 0);
    q_flat.ParallelCopy(m_tmp, 0, 0, 2);
    for (MFIter mfi(dst); mfi.isValid(); ++mfi) {
        auto const& d = dst.array(mfi);
        auto const& q = q_flat.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            d(i,j,k,0) = q(i,j,0,0);
            d(i,j,k,1) = q(i,j,0,1);
        });
    }
#else
    MultiFab::Copy(dst, m_tmp, 0, 0, 2, 0);
#endif
    dst.FillBoundary(m_geom.periodicity());
}
