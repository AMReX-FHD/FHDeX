#ifndef FLOCKING_PARAMS_H_
#define FLOCKING_PARAMS_H_

#include <cmath>
#include <string>

#include <AMReX_Algorithm.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>

// Model and initial-condition parameters shared by the SPDE (3D build: x, y,
// theta = z) and the particle code (2D build: theta is a particle attribute).
//
// Particles:  dX_i = g(theta_i) dt,  dtheta_i = sqrt(2 D) dW_i,
//             g(theta) = v (cos theta, sin theta)
// SPDE:       dphi = D d_theta^2 phi dt - div_x(phi g(theta)) dt
//                    + N^(-1/2) d_theta( sqrt(2 D phi) dB )
// phi is normalized to integrate to 1 over [0,Lx) x [0,Ly) x [0,2 pi).
//
// Read from the "flock." ParmParse namespace, except num_part and diff_coeff,
// which have no prefix as in the other src_deankow codes.

struct FlockingParams
{
    amrex::Real speed      = 1.0;   // v
    amrex::Real diff_coeff = 1.0;   // D, rotational diffusion
    amrex::Real num_part   = -1.0;  // N, must be set

    // initial condition, see flocking_ic
    int         ic_type    = 0;
    amrex::Real ic_amp     = 0.5;
    int         ic_kx      = 1;
    int         ic_ky      = 0;
    int         ic_m       = 1;
    amrex::Real ic_theta0  = 0.0;
    amrex::Real ic_sigma   = 0.1;
    amrex::Real ic_x0      = 0.5;
    amrex::Real ic_y0      = 0.5;

    // domain extent in x, y (theta always spans [0, 2 pi))
    amrex::Real xlo = 0.0, ylo = 0.0, lx = 1.0, ly = 1.0;

    // alignment (flocking) interaction, not implemented yet: must stay 0
    amrex::Real align_K    = 0.0;

    // density-dependent speed v = speed * s(u), u = rho_tilde/rho_bar, with
    // rho_tilde = W * rho the kernel-smoothed spatial density
    //   speed_type 0: s = 1;  1: s = exp(-lambda u);  2: s = max(0, 1 - lambda u)
    int         speed_type   = 0;
    amrex::Real speed_lambda = 1.0;
    //   kernel_type 0: Gaussian exp(-r^2/2R^2) (cut off at 4R for the pair sum)
    //               1: top-hat (r < R);  2: parabolic (1 - r^2/R^2)_+
    //               3: bump (1 - r^2/R^2)_+^2
    int         kernel_type  = 0;
    amrex::Real kernel_R     = 0.05;
    // particles only: 0 particle-mesh (CIC deposit, FFT, CIC interpolation),
    //                 1 direct pair sum over neighbours
    int         density_method = 0;
};

constexpr amrex::Real flock_two_pi = amrex::Real(2.0)*amrex::Real(3.14159265358979323846264338327950288);

// Unnormalized initial density f(x, y, theta) >= 0.
//   ic_type 0: uniform
//   ic_type 1: mode  1 + ic_amp cos(2 pi (ic_kx x/Lx + ic_ky y/Ly)) cos(ic_m (theta - ic_theta0)),
//              needs |ic_amp| <= 1
//   ic_type 2: blob  periodic Gaussian of width ic_sigma at (ic_x0, ic_y0), times
//              1 + ic_amp cos(theta - ic_theta0)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real flocking_ic (amrex::Real x, amrex::Real y, amrex::Real theta,
                         FlockingParams const& fp) noexcept
{
    if (fp.ic_type == 1) {
        const amrex::Real ph = flock_two_pi*(fp.ic_kx*(x - fp.xlo)/fp.lx + fp.ic_ky*(y - fp.ylo)/fp.ly);
        return 1.0 + fp.ic_amp*std::cos(ph)*std::cos(fp.ic_m*(theta - fp.ic_theta0));
    } else if (fp.ic_type == 2) {
        amrex::Real dx = x - fp.ic_x0;
        amrex::Real dy = y - fp.ic_y0;
        dx -= fp.lx*std::floor(dx/fp.lx + 0.5);
        dy -= fp.ly*std::floor(dy/fp.ly + 0.5);
        const amrex::Real s2 = fp.ic_sigma*fp.ic_sigma;
        return std::exp(-0.5*(dx*dx + dy*dy)/s2) * (1.0 + fp.ic_amp*std::cos(theta - fp.ic_theta0));
    }
    return 1.0;
}

// Upper bound of flocking_ic, for rejection sampling.
inline amrex::Real flocking_ic_max (FlockingParams const& fp)
{
    if (fp.ic_type == 1 || fp.ic_type == 2) { return 1.0 + std::abs(fp.ic_amp); }
    return 1.0;
}

// Speed factor s(u), u = rho_tilde / rho_bar
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real flock_speed_factor (amrex::Real u, FlockingParams const& fp) noexcept
{
    if (fp.speed_type == 1) { return std::exp(-fp.speed_lambda*u); }
    if (fp.speed_type == 2) { return amrex::max(amrex::Real(0.0), amrex::Real(1.0) - fp.speed_lambda*u); }
    return 1.0;
}

// Kernel shape W(r) up to normalization, zero beyond flock_kernel_cutoff
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real flock_kernel_shape (amrex::Real r, FlockingParams const& fp) noexcept
{
    const amrex::Real q = r/fp.kernel_R;
    if (fp.kernel_type == 0) { return (q < 4.0) ? std::exp(-0.5*q*q) : 0.0; }
    if (q >= 1.0) { return 0.0; }
    if (fp.kernel_type == 1) { return 1.0; }
    if (fp.kernel_type == 2) { return 1.0 - q*q; }
    return (1.0 - q*q)*(1.0 - q*q);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real flock_kernel_cutoff (FlockingParams const& fp) noexcept
{
    return (fp.kernel_type == 0) ? 4.0*fp.kernel_R : fp.kernel_R;
}

// Continuous normalization: W = norm * shape integrates to 1 over the plane
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real flock_kernel_norm (FlockingParams const& fp) noexcept
{
    const amrex::Real a = 0.5*flock_two_pi*fp.kernel_R*fp.kernel_R;   // pi R^2
    if (fp.kernel_type == 0) { return 1.0/(2.0*a*(1.0 - std::exp(-8.0))); }
    if (fp.kernel_type == 1) { return 1.0/a; }
    if (fp.kernel_type == 2) { return 2.0/a; }
    return 3.0/a;
}

// 2D Fourier transform of the normalized (untruncated) kernel
inline amrex::Real flock_kernel_hat (amrex::Real k, FlockingParams const& fp)
{
    const amrex::Real q = k*fp.kernel_R;
    if (fp.kernel_type == 0) { return std::exp(-0.5*q*q); }
    if (q < 1.e-6) { return 1.0; }
    if (fp.kernel_type == 1) { return 2.0*std::cyl_bessel_j(1.0, q)/q; }
    if (fp.kernel_type == 2) { return 8.0*std::cyl_bessel_j(2.0, q)/(q*q); }
    return 48.0*std::cyl_bessel_j(3.0, q)/(q*q*q);
}

// d ln s / d ln u at the mean density u = 1
inline amrex::Real flock_dlns_dlnu (FlockingParams const& fp)
{
    if (fp.speed_type == 1) { return -fp.speed_lambda; }
    if (fp.speed_type == 2) { return -fp.speed_lambda/(1.0 - fp.speed_lambda); }
    return 0.0;
}

// Linear stability of the uniform state in the diffusive limit (k l_p << 1):
// sigma(k) = -D_eff k^2 (1 + W_hat(k) dln s/dln u), D_eff = v(rho_bar)^2/(2 D)
inline void print_mips_stability (FlockingParams const& fp)
{
    if (fp.speed_type == 0) { return; }
    const amrex::Real g = flock_dlns_dlnu(fp);
    const amrex::Real v1 = fp.speed*flock_speed_factor(1.0, fp);
    const amrex::Real deff = (fp.diff_coeff > 0.0) ? v1*v1/(2.0*fp.diff_coeff) : 0.0;
    const amrex::Real lp = (fp.diff_coeff > 0.0) ? v1/fp.diff_coeff : 0.0;
    const amrex::Real kbox = flock_two_pi/amrex::max(fp.lx, fp.ly);
    amrex::Print() << "Density-dependent speed: speed_type = " << fp.speed_type
                   << " lambda = " << fp.speed_lambda << " kernel_type = " << fp.kernel_type
                   << " kernel_R = " << fp.kernel_R << "\n"
                   << "  at the mean density: v = " << v1 << " dln s/dln u = " << g
                   << " (long-wave instability when < -1)\n"
                   << "  persistence length v/D = " << lp << " (box " << fp.lx << " x " << fp.ly
                   << "), D_eff = v^2/(2D) = " << deff << "\n";
    for (int n = 1; n <= 4; ++n) {
        const amrex::Real k = n*kbox;
        const amrex::Real wh = flock_kernel_hat(k, fp);
        amrex::Print() << "  box mode " << n << ": k = " << k << " W_hat = " << wh
                       << " diffusive-limit growth rate = " << -deff*k*k*(1.0 + wh*g)
                       << ((1.0 + wh*g < 0.0) ? "  -> UNSTABLE" : "  -> stable")
                       << " (k l_p = " << k*lp << ")\n";
    }
}

inline void read_flocking_params (FlockingParams& fp)
{
    {
        amrex::ParmParse pp;
        pp.query("num_part", fp.num_part);
        pp.query("diff_coeff", fp.diff_coeff);
    }
    amrex::ParmParse pp("flock");
    pp.query("speed", fp.speed);
    pp.query("ic_type", fp.ic_type);
    pp.query("ic_amp", fp.ic_amp);
    pp.query("ic_kx", fp.ic_kx);
    pp.query("ic_ky", fp.ic_ky);
    pp.query("ic_m", fp.ic_m);
    pp.query("ic_theta0", fp.ic_theta0);
    pp.query("ic_sigma", fp.ic_sigma);
    pp.query("ic_x0", fp.ic_x0);
    pp.query("ic_y0", fp.ic_y0);
    pp.query("align_K", fp.align_K);
    pp.query("speed_type", fp.speed_type);
    pp.query("speed_lambda", fp.speed_lambda);
    pp.query("kernel_type", fp.kernel_type);
    pp.query("kernel_R", fp.kernel_R);
    pp.query("density_method", fp.density_method);

    if (fp.num_part <= 0.0) { amrex::Abort("num_part must be set to a positive value"); }
    if (fp.diff_coeff < 0.0) { amrex::Abort("diff_coeff must be >= 0"); }
    if (fp.ic_type < 0 || fp.ic_type > 2) { amrex::Abort("flock.ic_type must be 0, 1 or 2"); }
    if (fp.ic_type == 1 && std::abs(fp.ic_amp) > 1.0) { amrex::Abort("flock.ic_amp must be in [-1,1] for ic_type 1"); }
    if (fp.ic_type == 2 && (std::abs(fp.ic_amp) > 1.0 || fp.ic_sigma <= 0.0)) {
        amrex::Abort("flock.ic_type 2 needs |ic_amp| <= 1 and ic_sigma > 0");
    }
    if (fp.align_K != 0.0) { amrex::Abort("flock.align_K: the alignment interaction is not implemented yet"); }
    if (fp.speed_type < 0 || fp.speed_type > 2) { amrex::Abort("flock.speed_type must be 0, 1 or 2"); }
    if (fp.speed_lambda < 0.0) { amrex::Abort("flock.speed_lambda must be >= 0"); }
    if (fp.speed_type == 2 && fp.speed_lambda >= 1.0) {
        amrex::Abort("flock.speed_type 2 needs speed_lambda < 1, so that the speed is positive at the mean density");
    }
    if (fp.kernel_type < 0 || fp.kernel_type > 3) { amrex::Abort("flock.kernel_type must be 0, 1, 2 or 3"); }
    if (fp.kernel_R <= 0.0) { amrex::Abort("flock.kernel_R must be positive"); }
    if (fp.density_method != 0 && fp.density_method != 1) { amrex::Abort("flock.density_method must be 0 or 1"); }

    amrex::Print() << "Flocking: v = " << fp.speed << " D = " << fp.diff_coeff
                   << " N = " << fp.num_part << " ic_type = " << fp.ic_type << "\n";
}

#endif
