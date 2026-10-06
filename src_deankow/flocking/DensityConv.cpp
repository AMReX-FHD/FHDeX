#include "DensityConv.H"

using namespace amrex;

void
DensityConv::define (Geometry const& geom, BoxArray const& ba, DistributionMapping const& dm,
                     FlockingParams const& fp)
{
    m_geom = geom;
    const Box& domain = geom.Domain();
    m_fwd = std::make_unique<FFT::R2C<Real,FFT::Direction::forward>>(domain);
    m_bwd = std::make_unique<FFT::R2C<Real,FFT::Direction::backward>>(domain);

    auto const& [cba, cdm] = m_fwd->getSpectralDataLayout();
    m_what.define(cba, cdm, 1, 0);
    m_shat.define(cba, cdm, 1, 0);
    m_tmp.define(ba, dm, 1, 0);

    // W at the cell offsets, centred at cell 0 with periodic wrap (as
    // init_int_pot in interaction), the same in every theta plane in 3D
    const auto dx  = geom.CellSizeArray();
    const Real lx = geom.ProbLength(0);
    const Real ly = geom.ProbLength(1);
    const FlockingParams fp_loc = fp;
    MultiFab wk(ba, dm, 1, 0);
    for (MFIter mfi(wk); mfi.isValid(); ++mfi) {
        auto const& w = wk.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            amrex::ignore_unused(k);
            Real x = i*dx[0];
            Real y = j*dx[1];
            if (x + 0.5*dx[0] > 0.5*lx) { x -= lx; }
            if (y + 0.5*dx[1] > 0.5*ly) { y -= ly; }
            w(i,j,k) = flock_kernel_shape(std::sqrt(x*x + y*y), fp_loc);
        });
    }

    // normalize: sum over one x-y plane of W dx dy = 1
    Real nplanes = 1.0;
#if (AMREX_SPACEDIM == 3)
    nplanes = domain.length(2);
#endif
    const Real plane_sum = wk.sum(0)/nplanes*dx[0]*dx[1];
    if (plane_sum <= 0.0) { Abort("DensityConv: the kernel has no support on the grid; increase flock.kernel_R"); }
    wk.mult(1.0/plane_sum, 0, 1);
    m_norm_ratio = plane_sum*flock_kernel_norm(fp);

    m_fwd->forward(wk, m_what);
}

void
DensityConv::apply (MultiFab const& src, MultiFab& dst)
{
    AMREX_ALWAYS_ASSERT(defined());
    MultiFab::Copy(m_tmp, src, 0, 0, 1, 0);
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

    // FFTW is unnormalized: divide by the number of points, times the cell volume
    const auto dx = m_geom.CellSizeArray();
    const Real cellvol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
    m_tmp.mult(cellvol/m_geom.Domain().d_numPts(), 0, 1);
    MultiFab::Copy(dst, m_tmp, 0, 0, 1, 0);
    dst.FillBoundary(m_geom.periodicity());
}
