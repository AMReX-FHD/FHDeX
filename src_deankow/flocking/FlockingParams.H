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
// Particles:  dX_i = g(theta_i) dt,  dtheta_i = omega(X_i, theta_i) dt + sqrt(2 D) dW_i,
//             g(theta) = v (cos theta, sin theta)
// SPDE:       dphi = D d_theta^2 phi dt - div_x(phi g(theta)) dt
//                    + d_theta( phi d_theta V ) dt + N^(-1/2) d_theta( sqrt(2 D phi) dB )
// phi is normalized to integrate to 1 over [0,Lx) x [0,Ly) x [0,2 pi).
//
// Alignment: V = F * phi (convolution in x, y and theta) with the kernel
//   F(x, theta) = eps H(theta) G(|x|/R),  H(theta) = -cos(theta),
// so that V(x, theta) = -eps (cos(theta) Q_x + sin(theta) Q_y), with
// Q = G_R * p the kernel-smoothed polarization p = int phi (cos, sin)(theta) dtheta,
// and the turning rate is omega = -d_theta V = eps (cos(theta) Q_y - sin(theta) Q_x).
// G is not normalized: G = 1 for s < 0.9, 0 for s > 1.1, and the C1 smooth step
// 75 (1.1 - s)^2 - 250 (1.1 - s)^3 in between (s = r/R).
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

    // alignment interaction F = eps H(theta) G(r/R); eps = 0 switches it off
    amrex::Real align_eps  = 0.0;
    amrex::Real align_R    = 0.1;
    // particles only: Q at the particles by 0 particle-mesh (CIC deposit of the
    // polarization, FFT convolution with G, CIC interpolation) or 1 direct pair sum
    int         align_method = 0;
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

// Alignment kernel profile G(s), s = r/R: 1 for s < 0.9, 0 for s > 1.1, and the
// C1 smooth step 3 t^2 - 2 t^3 with t = 5 (1.1 - s) in between
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real flock_align_G (amrex::Real s) noexcept
{
    if (s <= 0.9) { return 1.0; }
    if (s >= 1.1) { return 0.0; }
    const amrex::Real u = 1.1 - s;
    return 75.0*u*u - 250.0*u*u*u;
}

// Cutoff radius of the alignment kernel
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real flock_align_cutoff (FlockingParams const& fp) noexcept
{
    return 1.1*fp.align_R;
}

// int G(|x|/R) dx dy = 2 pi R^2 int_0^1.1 s G(s) ds = 2 pi R^2 (0.405 + 0.1 - 0.004)
//                    = 1.002 pi R^2
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real flock_align_G_integral (FlockingParams const& fp) noexcept
{
    return 1.002*0.5*flock_two_pi*fp.align_R*fp.align_R;
}

// 2D Fourier transform of G(|x|/R): 2 pi R^2 int_0^1.1 s G(s) J0(k R s) ds
inline amrex::Real flock_align_G_hat (amrex::Real k, FlockingParams const& fp)
{
    const int n = 2200;
    const amrex::Real h = 1.1/n;
    amrex::Real sum = 0.0;
    for (int m = 0; m < n; ++m) {
        const amrex::Real s = (m + 0.5)*h;
        sum += s*flock_align_G(s)*std::cyl_bessel_j(0.0, k*fp.align_R*s);
    }
    return flock_two_pi*fp.align_R*fp.align_R*sum*h;
}

// Linear stability of the uniform isotropic state against polarization: for
// p ~ exp(i k.x), the alignment alone gives the rate eps rho_bar G_hat(k)/2 - D.
// At k = 0 this is exact (the polar transition is at eps_c = 2 D/(rho_bar int G));
// at k > 0 self-propulsion couples p to the density and the rate is an estimate.
inline void print_align_stability (FlockingParams const& fp, amrex::Real rhobar)
{
    const amrex::Real gint = flock_align_G_integral(fp);
    const amrex::Real epsc = 2.0*fp.diff_coeff/(rhobar*gint);
    amrex::Print() << "Alignment: eps = " << fp.align_eps << " R = " << fp.align_R
                   << " int G = " << gint << " threshold eps_c = 2 D/(rho_bar int G) = " << epsc
                   << " (eps/eps_c = " << ((epsc > 0.0) ? fp.align_eps/epsc : 0.0) << ")\n";
    if (fp.align_eps == 0.0) { return; }
    const amrex::Real kbox = flock_two_pi/amrex::max(fp.lx, fp.ly);
    for (int n = 0; n <= 4; ++n) {
        const amrex::Real k = n*kbox;
        const amrex::Real gh = flock_align_G_hat(k, fp);
        const amrex::Real rate = 0.5*fp.align_eps*rhobar*gh - fp.diff_coeff;
        amrex::Print() << "  mode " << n << ": k = " << k << " G_hat = " << gh
                       << " polarization rate eps rho_bar G_hat/2 - D = " << rate
                       << ((rate > 0.0) ? "  -> grows" : "  -> decays")
                       << ((n == 0) ? " (exact)" : " (alignment only)") << "\n";
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
    pp.query("align_eps", fp.align_eps);
    pp.query("align_R", fp.align_R);
    pp.query("align_method", fp.align_method);
    for (const char* key : {"speed_type", "speed_lambda", "kernel_type", "kernel_R",
                            "density_method", "density_check", "align_K"}) {
        if (pp.contains(key)) {
            amrex::Abort(std::string("flock.") + key + " is no longer used: the flocking code has a"
                         " constant speed and the alignment kernel (flock.align_eps, flock.align_R);"
                         " the density-dependent speed is in src_deankow/active");
        }
    }

    if (fp.num_part <= 0.0) { amrex::Abort("num_part must be set to a positive value"); }
    if (fp.diff_coeff < 0.0) { amrex::Abort("diff_coeff must be >= 0"); }
    if (fp.ic_type < 0 || fp.ic_type > 2) { amrex::Abort("flock.ic_type must be 0, 1 or 2"); }
    if (fp.ic_type == 1 && std::abs(fp.ic_amp) > 1.0) { amrex::Abort("flock.ic_amp must be in [-1,1] for ic_type 1"); }
    if (fp.ic_type == 2 && (std::abs(fp.ic_amp) > 1.0 || fp.ic_sigma <= 0.0)) {
        amrex::Abort("flock.ic_type 2 needs |ic_amp| <= 1 and ic_sigma > 0");
    }
    if (fp.align_eps < 0.0) { amrex::Abort("flock.align_eps must be >= 0"); }
    if (fp.align_R <= 0.0) { amrex::Abort("flock.align_R must be positive"); }
    if (fp.align_method != 0 && fp.align_method != 1) { amrex::Abort("flock.align_method must be 0 or 1"); }

    amrex::Print() << "Flocking: v = " << fp.speed << " D = " << fp.diff_coeff
                   << " N = " << fp.num_part << " ic_type = " << fp.ic_type << "\n";
}

#endif
