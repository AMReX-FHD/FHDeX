#ifndef INERTIAL_PARAMS_H_
#define INERTIAL_PARAMS_H_

#include <cmath>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>

// Inertial (underdamped) Brownian particles in one space and one velocity
// dimension (mass 1), and their Dean-Kawasaki SPDE in phase space (x, v):
//
//   Particles:  dX_i = v_i dt,  dv_i = -gamma v_i dt - U'(X_i) dt + sqrt(2 D) dW_i
//   SPDE:       dphi = -d_x(v phi) dt + d_v[(gamma v + U') phi] dt + D d_vv phi dt
//                      + N^(-1/2) d_v( sqrt(2 D phi) dB )
//
// phi is normalized to integrate to 1 over [0, L) x [-Vmax, Vmax]; x is
// periodic, v has reflecting walls at +-Vmax. The equilibrium is
// phi_eq ~ exp(-(v^2/2 + U)/T) with T = D/gamma.
//
// Read from the "inertial." ParmParse namespace, except num_part and
// diff_coeff, which have no prefix as in the other src_deankow codes.

struct InertialParams
{
    amrex::Real gamma      = 1.0;   // friction
    amrex::Real diff_coeff = 1.0;   // D, velocity diffusion
    amrex::Real num_part   = -1.0;  // N, must be set

    // potential: pot_type 0 harmonic U = k/2 (x - x0)^2, minimum image in x
    int         pot_type   = 0;
    amrex::Real pot_k      = 100.0;
    amrex::Real pot_x0     = -1.0;  // default: the centre of the x domain

    // initial condition
    //   ic_type 0: equilibrium phi_eq
    //   ic_type 1: Gaussian at (ic_x0, ic_v0) with widths ic_sx, ic_sv
    //   ic_type 2: uniform
    int         ic_type    = 1;
    amrex::Real ic_x0      = -1.0;  // default: pot_x0 + 0.2
    amrex::Real ic_v0      = 0.0;
    amrex::Real ic_sx      = 0.05;
    amrex::Real ic_sv      = 0.5;

    // domain: x in [xlo, xlo + lx), v in [-vmax, vmax]
    amrex::Real xlo = 0.0, lx = 1.0, vmax = 6.0;

    // 0: SPDE, 1: particles
    int         run_mode   = 0;

    // particle integrator: 0 Euler-Maruyama, 1 BAOAB (kick-drift-OU-drift-kick
    // splitting with an exact Ornstein-Uhlenbeck step for friction and noise)
    int         particle_integrator = 0;

    // particle random numbers: 0 per-thread engine (ParallelForRNG), 1 counter
    // based on (particle id, cpu, step, draw), so a particle gets the same
    // increments on any patch or rank and runs do not depend on the
    // decomposition
    int         particle_rng = 0;

    // hybrid: SPDE on level 0, particles on level 1 where N phi dx dv < n_refine
    int         hybrid = 0;
    amrex::Real n_refine = 10.0;
    int         hybrid_nghost = -1;   // ghost ring of the particle region (default 1 EM, 2 BAOAB)
    amrex::Real rng_seed = 0;         // folded into the counter-based random numbers
    // testing: if tag_box_hi >= tag_box_lo (cell indices x0 v0 x1 v1), the
    // particle region is this box instead of the density criterion
    int tag_box[4] = {0, 0, -1, -1};
    // sampling of v within a cell when particles are created from the SPDE:
    // 0 uniform; 1 (default) tilted by the local equilibrium slope, density
    // ~ exp(-v_c (v - v_low)/T) across the cell (as x uses -U'(x_c)/T);
    // 2 linear profile from the MC-limited v slope of phi itself
    int ghost_v_sampling = 1;
};

// Minimum-image displacement x - x0 on the periodic x domain
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real inertial_dx (amrex::Real x, amrex::Real x0, InertialParams const& ip) noexcept
{
    amrex::Real d = x - x0;
    d -= ip.lx*std::floor(d/ip.lx + 0.5);
    return d;
}

// Potential U(x) and force derivative U'(x)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real inertial_U (amrex::Real x, InertialParams const& ip) noexcept
{
    const amrex::Real d = inertial_dx(x, ip.pot_x0, ip);
    return 0.5*ip.pot_k*d*d;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real inertial_dU (amrex::Real x, InertialParams const& ip) noexcept
{
    return ip.pot_k*inertial_dx(x, ip.pot_x0, ip);
}

// Largest |U'| over the domain (at the periodic seam for the harmonic potential)
inline amrex::Real inertial_max_dU (InertialParams const& ip)
{
    return ip.pot_k*0.5*ip.lx;
}

// Temperature T = D/gamma
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real inertial_T (InertialParams const& ip) noexcept
{
    return ip.diff_coeff/ip.gamma;
}

// Unnormalized equilibrium density exp(-(v^2/2 + U)/T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real inertial_phi_eq (amrex::Real x, amrex::Real v, InertialParams const& ip) noexcept
{
    return std::exp(-(0.5*v*v + inertial_U(x, ip))/inertial_T(ip));
}

// Unnormalized initial density
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real inertial_ic (amrex::Real x, amrex::Real v, InertialParams const& ip) noexcept
{
    if (ip.ic_type == 0) { return inertial_phi_eq(x, v, ip); }
    if (ip.ic_type == 1) {
        const amrex::Real d = inertial_dx(x, ip.ic_x0, ip);
        const amrex::Real w = v - ip.ic_v0;
        return std::exp(-0.5*(d*d/(ip.ic_sx*ip.ic_sx) + w*w/(ip.ic_sv*ip.ic_sv)));
    }
    return 1.0;
}

inline void read_inertial_params (InertialParams& ip)
{
    {
        amrex::ParmParse pp;
        pp.query("num_part", ip.num_part);
        pp.query("diff_coeff", ip.diff_coeff);
    }
    amrex::ParmParse pp("inertial");
    pp.query("gamma", ip.gamma);
    pp.query("pot_type", ip.pot_type);
    pp.query("pot_k", ip.pot_k);
    pp.query("pot_x0", ip.pot_x0);
    pp.query("ic_type", ip.ic_type);
    pp.query("ic_x0", ip.ic_x0);
    pp.query("ic_v0", ip.ic_v0);
    pp.query("ic_sx", ip.ic_sx);
    pp.query("ic_sv", ip.ic_sv);
    pp.query("run_mode", ip.run_mode);
    pp.query("particle_integrator", ip.particle_integrator);
    pp.query("particle_rng", ip.particle_rng);
    pp.query("hybrid", ip.hybrid);
    pp.query("n_refine", ip.n_refine);
    pp.query("hybrid_nghost", ip.hybrid_nghost);
    pp.query("ghost_v_sampling", ip.ghost_v_sampling);
    if (ip.ghost_v_sampling < 0 || ip.ghost_v_sampling > 2) {
        amrex::Abort("inertial.ghost_v_sampling must be 0, 1 or 2");
    }
    {
        amrex::Vector<int> tb;
        if (pp.queryarr("tag_box", tb) && tb.size() == 4) {
            for (int q = 0; q < 4; ++q) { ip.tag_box[q] = tb[q]; }
        }
    }

    if (ip.num_part <= 0.0) { amrex::Abort("num_part must be set to a positive value"); }
    if (ip.diff_coeff <= 0.0) { amrex::Abort("diff_coeff must be positive"); }
    if (ip.gamma <= 0.0) { amrex::Abort("inertial.gamma must be positive"); }
    if (ip.pot_type != 0) { amrex::Abort("inertial.pot_type must be 0 (harmonic)"); }
    if (ip.pot_k < 0.0) { amrex::Abort("inertial.pot_k must be >= 0"); }
    if (ip.ic_type < 0 || ip.ic_type > 2) { amrex::Abort("inertial.ic_type must be 0, 1 or 2"); }
    if (ip.ic_type == 1 && (ip.ic_sx <= 0.0 || ip.ic_sv <= 0.0)) {
        amrex::Abort("inertial.ic_type 1 needs ic_sx > 0 and ic_sv > 0");
    }
    if (ip.run_mode != 0 && ip.run_mode != 1) { amrex::Abort("inertial.run_mode must be 0 (SPDE) or 1 (particles)"); }
    if (ip.hybrid) {
        if (ip.run_mode != 0) { amrex::Abort("inertial.hybrid = 1 needs run_mode = 0 (SPDE on level 0)"); }
        ip.particle_rng = 1;
        if (ip.hybrid_nghost < 0) { ip.hybrid_nghost = (ip.particle_integrator == 1) ? 2 : 1; }
        if (ip.hybrid_nghost < 1) { amrex::Abort("inertial.hybrid_nghost must be >= 1"); }
    }
    if (ip.particle_rng != 0 && ip.particle_rng != 1) { amrex::Abort("inertial.particle_rng must be 0 or 1"); }
    if (ip.particle_integrator != 0 && ip.particle_integrator != 1) {
        amrex::Abort("inertial.particle_integrator must be 0 (Euler-Maruyama) or 1 (BAOAB)");
    }
}

// Fill in the defaults that depend on the domain, and print a summary
inline void finish_inertial_params (InertialParams& ip)
{
    if (ip.pot_x0 < 0.0) { ip.pot_x0 = ip.xlo + 0.5*ip.lx; }
    if (ip.ic_x0 < 0.0) { ip.ic_x0 = ip.pot_x0 + 0.2*ip.lx; }
    const amrex::Real T = inertial_T(ip);
    const amrex::Real sx = (ip.pot_k > 0.0) ? std::sqrt(T/ip.pot_k) : 0.0;
    amrex::Print() << "Inertial DK: " << ((ip.run_mode == 0) ? "SPDE" :
                                          ((ip.particle_integrator == 0) ? "particles (Euler-Maruyama)"
                                                                         : "particles (BAOAB)"))
                   << ", N = " << ip.num_part << ", gamma = " << ip.gamma << ", D = " << ip.diff_coeff
                   << ", T = D/gamma = " << T << "\n"
                   << "  harmonic U = k/2 (x - x0)^2: k = " << ip.pot_k << ", x0 = " << ip.pot_x0
                   << ", omega = " << std::sqrt(ip.pot_k) << "; equilibrium sigma_x = " << sx
                   << " (seam at " << ((sx > 0.0) ? 0.5*ip.lx/sx : 0.0) << " sigma_x), sigma_v = "
                   << std::sqrt(T) << " (wall at " << ip.vmax/std::sqrt(T) << " sigma_v)\n"
                   << "  ic_type = " << ip.ic_type;
    if (ip.ic_type == 1) {
        amrex::Print() << " (Gaussian at x = " << ip.ic_x0 << ", v = " << ip.ic_v0
                       << ", widths " << ip.ic_sx << ", " << ip.ic_sv << ")";
    }
    amrex::Print() << "\n";
    if (ip.hybrid) {
        amrex::Print() << "  hybrid: particles on level 1 where N phi dx dv < " << ip.n_refine
                       << ", ghost ring " << ip.hybrid_nghost << " cell(s), counter-based particle random numbers\n";
    }
}

#endif
