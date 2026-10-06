#ifndef FLOCKING_PARAMS_H_
#define FLOCKING_PARAMS_H_

#include <cmath>
#include <string>

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

    if (fp.num_part <= 0.0) { amrex::Abort("num_part must be set to a positive value"); }
    if (fp.diff_coeff < 0.0) { amrex::Abort("diff_coeff must be >= 0"); }
    if (fp.ic_type < 0 || fp.ic_type > 2) { amrex::Abort("flock.ic_type must be 0, 1 or 2"); }
    if (fp.ic_type == 1 && std::abs(fp.ic_amp) > 1.0) { amrex::Abort("flock.ic_amp must be in [-1,1] for ic_type 1"); }
    if (fp.ic_type == 2 && (std::abs(fp.ic_amp) > 1.0 || fp.ic_sigma <= 0.0)) {
        amrex::Abort("flock.ic_type 2 needs |ic_amp| <= 1 and ic_sigma > 0");
    }
    if (fp.align_K != 0.0) { amrex::Abort("flock.align_K: the alignment interaction is not implemented yet"); }

    amrex::Print() << "Flocking: v = " << fp.speed << " D = " << fp.diff_coeff
                   << " N = " << fp.num_part << " ic_type = " << fp.ic_type << "\n";
}

#endif
