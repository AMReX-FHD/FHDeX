#ifndef FLOCK_KERNEL_H_
#define FLOCK_KERNEL_H_

#include <AMReX_FArrayBox.H>
#include "FlockingParams.H"

// Kernels for the flocking SPDE on (x, y, theta), theta = z. Fluxes are
// physical fluxes on faces (positive in the +d direction); the update is
//   phi_new = phi_old - dt * sum_d (F_d(+face) - F_d(-face)) / dx_d

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void init_phi (int i, int j, int k, amrex::Array4<amrex::Real> const& phi,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
               FlockingParams const& fp)
{
    const amrex::Real x  = plo[0] + (i+0.5)*dx[0];
    const amrex::Real y  = plo[1] + (j+0.5)*dx[1];
    const amrex::Real th = plo[2] + (k+0.5)*dx[2];
    phi(i,j,k) = flocking_ic(x, y, th, fp);
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
// and p (plus side), given the cells m2 = m-1 and p2 = p+1 for the slopes.
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

// x face (i-1/2, j, k): u = v s(rho_tilde_face/rho_bar) cos(theta_k), with the
// smoothed density rho_tilde averaged from the two adjacent cells (s = 1 for
// speed_type 0)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_x (int i, int j, int k, amrex::Array4<amrex::Real> const& fx,
                     amrex::Array4<amrex::Real const> const& phi,
                     amrex::Array4<amrex::Real const> const& rhot, amrex::Real rhobar_inv,
                     amrex::Real theta, FlockingParams const& fp, int adv_order, int limiter)
{
    const amrex::Real u = fp.speed*std::cos(theta)
        * flock_speed_factor(0.5*(rhot(i-1,j,k) + rhot(i,j,k))*rhobar_inv, fp);
    fx(i,j,k) = upwind_flux(u, phi(i-2,j,k), phi(i-1,j,k), phi(i,j,k), phi(i+1,j,k),
                            adv_order, limiter);
}

// y face (i, j-1/2, k): u = v s(rho_tilde_face/rho_bar) sin(theta_k)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_y (int i, int j, int k, amrex::Array4<amrex::Real> const& fy,
                     amrex::Array4<amrex::Real const> const& phi,
                     amrex::Array4<amrex::Real const> const& rhot, amrex::Real rhobar_inv,
                     amrex::Real theta, FlockingParams const& fp, int adv_order, int limiter)
{
    const amrex::Real u = fp.speed*std::sin(theta)
        * flock_speed_factor(0.5*(rhot(i,j-1,k) + rhot(i,j,k))*rhobar_inv, fp);
    fy(i,j,k) = upwind_flux(u, phi(i,j-2,k), phi(i,j-1,k), phi(i,j,k), phi(i,j+1,k),
                            adv_order, limiter);
}

// theta face (i, j, k-1/2): -D d_theta phi + omega*phi_face, with the turning
// rate omega upwinded. omega is the slot for the alignment interaction and is
// 0 for now.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_theta (int i, int j, int k, amrex::Array4<amrex::Real> const& fz,
                         amrex::Array4<amrex::Real const> const& phi,
                         amrex::Real dthinv, FlockingParams const& fp, amrex::Real omega,
                         int adv_order, int limiter)
{
    amrex::Real f = -fp.diff_coeff*(phi(i,j,k) - phi(i,j,k-1))*dthinv;
    if (omega != 0.0) {
        f += upwind_flux(omega, phi(i,j,k-2), phi(i,j,k-1), phi(i,j,k), phi(i,j,k+1),
                         adv_order, limiter);
    }
    fz(i,j,k) = f;
}

// Stochastic theta flux sqrt(phi_face) * W, with W ~ N(0, 2 D/(N dV dt))
// already in sz (filled by MultiFabFillRandom). phi_face is the average of the
// two sides, clipped at 0.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void scale_stoch_flux_theta (int i, int j, int k, amrex::Array4<amrex::Real> const& sz,
                             amrex::Array4<amrex::Real const> const& phi)
{
    const amrex::Real pf = 0.5*(phi(i,j,k) + phi(i,j,k-1));
    sz(i,j,k) *= std::sqrt(amrex::max(pf, amrex::Real(0.0)));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void update_phi (int i, int j, int k,
                 amrex::Array4<amrex::Real const> const& phiold,
                 amrex::Array4<amrex::Real      > const& phinew,
                 amrex::Array4<amrex::Real const> const& fx,
                 amrex::Array4<amrex::Real const> const& fy,
                 amrex::Array4<amrex::Real const> const& fz,
                 amrex::Real dt, amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxinv)
{
    phinew(i,j,k) = phiold(i,j,k)
        - dt*( (fx(i+1,j,k) - fx(i,j,k))*dxinv[0]
             + (fy(i,j+1,k) - fy(i,j,k))*dxinv[1]
             + (fz(i,j,k+1) - fz(i,j,k))*dxinv[2] );
}

#endif
