#ifndef INERTIAL_KERNEL_H_
#define INERTIAL_KERNEL_H_

#include <AMReX_FArrayBox.H>
#include "InertialParams.H"

// Kernels for the inertial SPDE on (x, v), x = direction 0 (periodic) and
// v = direction 1 (reflecting walls at +-Vmax). Fluxes are physical fluxes on
// faces (positive in the +d direction); the update is
//   phi_new = phi_old - dt * sum_d (F_d(+face) - F_d(-face)) / dx_d

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void init_phi (int i, int j, int k, amrex::Array4<amrex::Real> const& phi,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
               InertialParams const& ip)
{
    const amrex::Real x = plo[0] + (i+0.5)*dx[0];
    const amrex::Real v = plo[1] + (j+0.5)*dx[1];
    phi(i,j,k) = inertial_ic(x, v, ip);
}

// Limited slope from the left and right differences a, b.
//   limiter 0: minmod, 1: monotonized central (MC)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real limited_slope (amrex::Real a, amrex::Real b, int limiter) noexcept
{
    if (a*b <= 0.0) { return 0.0; }
    const amrex::Real s = (a > 0.0) ? 1.0 : -1.0;
    const amrex::Real aa = std::abs(a), ab = std::abs(b);
    if (limiter == 1) {
        return s*amrex::min(amrex::min(2.0*aa, 2.0*ab), 0.5*(aa+ab));
    }
    return s*amrex::min(aa, ab);
}

// Upwind advective flux u*phi_face on the face between cells m (minus side)
// and p (plus side), given the cells m2 = m-1 and p2 = p+1 for the slopes
// (as in src_deankow/active).
//   adv_order 0: centred, u (m + p)/2. Skew-symmetric, so it adds no numerical
//                diffusion and keeps the discrete fluctuation-dissipation
//                balance; not positivity preserving, and needs a time
//                integrator whose stability region covers the imaginary axis
//   adv_order 1: first-order upwind
//   adv_order 2: MUSCL, limited linear reconstruction from the upwind cell
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real upwind_flux (amrex::Real u, amrex::Real m2, amrex::Real m, amrex::Real p,
                         amrex::Real p2, int adv_order, int limiter) noexcept
{
    if (adv_order == 0) { return 0.5*u*(m + p); }
    if (adv_order == 2) {
        if (u >= 0.0) { return u*(m + 0.5*limited_slope(m - m2, p - m, limiter)); }
        return u*(p - 0.5*limited_slope(p - m, p2 - p, limiter));
    }
    return (u >= 0.0) ? u*m : u*p;
}

// x face (i-1/2, j): free streaming, u = v_j at the cell-centre velocity
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_x (int i, int j, int k, amrex::Array4<amrex::Real> const& fx,
                     amrex::Array4<amrex::Real const> const& phi,
                     amrex::Real vc, int adv_order, int limiter)
{
    fx(i,j,k) = upwind_flux(vc, phi(i-2,j,k), phi(i-1,j,k), phi(i,j,k), phi(i+1,j,k),
                            adv_order, limiter);
}

// v face (i, j-1/2): -D d_v phi + a phi_face with the drift
// a = -(gamma v_face + U'(x_i)), upwinded. Zero on the wall faces (j = jlo,
// j = jhi + 1): reflecting walls.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_v (int i, int j, int k, amrex::Array4<amrex::Real> const& fv,
                     amrex::Array4<amrex::Real const> const& phi,
                     amrex::Real vface, amrex::Real dU, amrex::Real dvinv,
                     int jlo, int jhi, InertialParams const& ip, int adv_order, int limiter)
{
    if (j <= jlo || j >= jhi + 1) { fv(i,j,k) = 0.0; return; }
    const amrex::Real a = -(ip.gamma*vface + dU);
    fv(i,j,k) = -ip.diff_coeff*(phi(i,j,k) - phi(i,j-1,k))*dvinv
              + upwind_flux(a, phi(i,j-2,k), phi(i,j-1,k), phi(i,j,k), phi(i,j+1,k),
                            adv_order, limiter);
}

// Stochastic v flux sqrt(phi_face) * W, with W ~ N(0, 2 D/(N dV dt)) already in
// sv; phi_face is the average of the two sides, clipped at 0. Zero on the walls.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void scale_stoch_flux_v (int i, int j, int k, amrex::Array4<amrex::Real> const& sv,
                         amrex::Array4<amrex::Real const> const& phi, int jlo, int jhi)
{
    if (j <= jlo || j >= jhi + 1) { sv(i,j,k) = 0.0; return; }
    const amrex::Real pf = 0.5*(phi(i,j,k) + phi(i,j-1,k));
    sv(i,j,k) *= std::sqrt(amrex::max(pf, amrex::Real(0.0)));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void update_phi (int i, int j, int k,
                 amrex::Array4<amrex::Real const> const& phiold,
                 amrex::Array4<amrex::Real      > const& phinew,
                 amrex::Array4<amrex::Real const> const& fx,
                 amrex::Array4<amrex::Real const> const& fv,
                 amrex::Real dt, amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxinv)
{
    phinew(i,j,k) = phiold(i,j,k)
        - dt*( (fx(i+1,j,k) - fx(i,j,k))*dxinv[0]
             + (fv(i,j+1,k) - fv(i,j,k))*dxinv[1] );
}

#endif
