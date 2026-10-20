#ifndef _RunningMean_H_
#define _RunningMean_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>
#include <AMReX_CompensatedSum.H>

// Running averages are updated in increment form, mean += (sample - mean)/steps.
// In single precision the increment drops below half an ulp of the mean once
// steps > sigma/(eps*|<sample>|), so accumulators with a nonzero expected value
// (means of rho, rhoE, rhoYk, Xk, T, theta and second moments with nonzero
// expectation) use Kahan compensation. Zero-expectation accumulators shrink like
// sigma/sqrt(steps) and never swamp, so they use the plain update.
// stepsinv should be computed as Real(1.0/static_cast<double>(steps)).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void runningMean (amrex::Real& mean, amrex::Real sample, amrex::Real stepsinv) noexcept
{
    mean += (sample - mean)*stepsinv;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void runningMeanKahan (amrex::Real& mean, amrex::Real& comp, amrex::Real sample, amrex::Real stepsinv) noexcept
{
    amrex::compensatedAdd(mean, comp, (sample - mean)*stepsinv);
}

#endif
