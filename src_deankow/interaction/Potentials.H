#ifndef POTENTIALS_H_
#define POTENTIALS_H_

#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>

// Parameters for the interaction potential (prefix ip_) and the external
// potential (prefix ep_). All are read from the "amr." ParmParse namespace and
// shared by the SPDE (mykernel.H) and the particle (StochasticPC) discretizations.

// Functional form of the interaction potential. To add a new form: add an enum
// value here, a name in read_potential_params, and a case in ip_U and
// ip_dUdr_over_r.
enum IntPotType : int {
    GEMA = 0,  // generalized exponential model: U = ip_eps*exp(-(r/ip_R)^ip_alpha)
    HK   = 1,  // U = ip_eps*ip_R*(1-(r/ip_R)^2)/2 for r < ip_R, 0 otherwise
    MORSE = 2, // generalized Morse: U = -ip_eps_att*exp(-(r-ip_re)/ip_R_att)
               //                        + ip_eps_rep*exp(-(r-ip_re)/ip_R_rep)
    PP   = 3   // piecewise parabolic (Gerber et al., arXiv:2510.17629):
               // U = -ip_eps*ip_R*w(r/ip_R) with, for a = ip_pp_a,
               //   w(s) = ip_pp_alpha*(s^2-a^2)/2 + ip_pp_beta*(a^2-1)/2  s <= a
               //   w(s) = ip_pp_beta*(s^2-1)/2                            a < s <= 1
               //   w(s) = 0                                               s > 1
               // ip_pp_alpha = ip_pp_beta = 1 is hk
};

// Below this separation the Morse pair force is set to zero: dU/dr does not
// vanish at r = 0, so (dU/dr)/r would be dominated by round-off.
constexpr amrex::Real ip_morse_rmin = 1.e-6;

struct PotentialParams
{
    // interaction potential
    int         use_int_pot = 1;
    int         ip_type     = IntPotType::GEMA;
    amrex::Real ip_eps      = 666.0;
    amrex::Real ip_R        = 0.1;
    amrex::Real ip_alpha    = 3.0;
    amrex::Real ip_eps_att  = 1.0;   // morse only
    amrex::Real ip_eps_rep  = 0.5;   // morse only
    amrex::Real ip_R_att    = 0.1;   // morse only
    amrex::Real ip_R_rep    = 0.05;  // morse only
    amrex::Real ip_re       = 0.0;   // morse only
    amrex::Real ip_pp_alpha = 1.0;   // pp only: core stiffness
    amrex::Real ip_pp_beta  = 1.0;   // pp only: tail stiffness
    amrex::Real ip_pp_a     = 0.5;   // pp only: core radius / ip_R
    amrex::Real ip_range    = 0.0;   // particle pair cutoff radius; <= 0 means no cutoff

    // external potential
    int         use_ext_pot = 0;
    amrex::Real ep_alpha    = 0.25;
    amrex::Real ep_beta     = 0.75;
    amrex::Real ep_gamma    = 1.1e-3;
};

// Interaction potential U(r), not truncated by ip_range (the SPDE uses the
// full kernel; ip_range is the particle pair cutoff)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ip_U (amrex::Real r, PotentialParams const& pot) noexcept
{
    switch (pot.ip_type) {
    case IntPotType::MORSE:
        return -pot.ip_eps_att * std::exp(-(r - pot.ip_re)/pot.ip_R_att)
               + pot.ip_eps_rep * std::exp(-(r - pot.ip_re)/pot.ip_R_rep);
    case IntPotType::HK:
    {
        if (r >= pot.ip_R) { return 0.0; }
        const amrex::Real s = r/pot.ip_R;
        return 0.5 * pot.ip_eps * pot.ip_R * (1.0 - s*s);
    }
    case IntPotType::PP:
    {
        if (r >= pot.ip_R) { return 0.0; }
        const amrex::Real s = r/pot.ip_R;
        const amrex::Real a = pot.ip_pp_a;
        const amrex::Real w = (s <= a)
            ? 0.5*pot.ip_pp_alpha*(s*s - a*a) + 0.5*pot.ip_pp_beta*(a*a - 1.0)
            : 0.5*pot.ip_pp_beta*(s*s - 1.0);
        return -pot.ip_eps * pot.ip_R * w;
    }
    case IntPotType::GEMA:
    default:
        return pot.ip_eps * std::exp(-std::pow(r/pot.ip_R, pot.ip_alpha));
    }
}

// (dU/dr)/r, so that the force on particle i from particle j is
// (dU/dr)/r * (x_j - x_i)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ip_dUdr_over_r (amrex::Real r, PotentialParams const& pot) noexcept
{
    if (r <= 0.0) { return 0.0; }
    if (pot.ip_range > 0.0 && r >= pot.ip_range) { return 0.0; }
    switch (pot.ip_type) {
    case IntPotType::MORSE:
    {
        if (r < ip_morse_rmin) { return 0.0; }
        const amrex::Real dUdr =
              pot.ip_eps_att/pot.ip_R_att * std::exp(-(r - pot.ip_re)/pot.ip_R_att)
            - pot.ip_eps_rep/pot.ip_R_rep * std::exp(-(r - pot.ip_re)/pot.ip_R_rep);
        return dUdr / r;
    }
    case IntPotType::HK:
        return (r < pot.ip_R) ? -pot.ip_eps / pot.ip_R : 0.0;
    case IntPotType::PP:
        if (r >= pot.ip_R) { return 0.0; }
        return -pot.ip_eps / pot.ip_R
               * ((r <= pot.ip_pp_a*pot.ip_R) ? pot.ip_pp_alpha : pot.ip_pp_beta);
    case IntPotType::GEMA:
    default:
    {
        const amrex::Real s = std::pow(r/pot.ip_R, pot.ip_alpha);
        return -pot.ip_alpha * pot.ip_eps * std::exp(-s) * s / (r*r);
    }
    }
}

// Magnitude scale of U: max |U| sampled on [0, rmax]. This is |U(0)| for the
// monotone kernels (gema, hk, pp); for morse it also sees the well.
inline amrex::Real ip_Uscale (PotentialParams const& pot, amrex::Real rmax)
{
    constexpr int nsamp = 4096;
    amrex::Real umax = std::abs(ip_U(0.0, pot));
    for (int n = 1; n <= nsamp; ++n) {
        umax = std::max(umax, std::abs(ip_U(rmax*n/nsamp, pot)));
    }
    return umax;
}

inline void read_potential_params (PotentialParams& pot)
{
    amrex::ParmParse pp("amr");

    // Old names are rejected rather than silently ignored
    const char* old_names[3][2] = {{"interaction_strength", "ip_eps"},
                                   {"interaction_scale",    "ip_R"},
                                   {"interaction_range",    "ip_range"}};
    for (auto const& nm : old_names) {
        if (pp.contains(nm[0])) {
            amrex::Abort(std::string("amr.") + nm[0] + " is no longer supported; use amr." + nm[1]);
        }
    }

    pp.query("use_int_pot", pot.use_int_pot);
    std::string ip_type_name = "gema";
    pp.query("ip_type", ip_type_name);
    std::transform(ip_type_name.begin(), ip_type_name.end(), ip_type_name.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (ip_type_name == "gema") {
        pot.ip_type = IntPotType::GEMA;
    } else if (ip_type_name == "hk") {
        pot.ip_type = IntPotType::HK;
    } else if (ip_type_name == "morse") {
        pot.ip_type = IntPotType::MORSE;
    } else if (ip_type_name == "pp") {
        pot.ip_type = IntPotType::PP;
    } else {
        amrex::Abort("Unknown amr.ip_type = " + ip_type_name + "; valid types are: gema, hk, morse, pp");
    }
    pp.query("ip_eps", pot.ip_eps);
    pp.query("ip_R", pot.ip_R);
    pp.query("ip_alpha", pot.ip_alpha);
    pp.query("ip_eps_att", pot.ip_eps_att);
    pp.query("ip_eps_rep", pot.ip_eps_rep);
    pp.query("ip_R_att", pot.ip_R_att);
    pp.query("ip_R_rep", pot.ip_R_rep);
    pp.query("ip_re", pot.ip_re);
    pp.query("ip_pp_alpha", pot.ip_pp_alpha);
    pp.query("ip_pp_beta", pot.ip_pp_beta);
    pp.query("ip_pp_a", pot.ip_pp_a);
    pp.query("ip_range", pot.ip_range);

    pp.query("use_ext_pot", pot.use_ext_pot);
    pp.query("ep_alpha", pot.ep_alpha);
    pp.query("ep_beta", pot.ep_beta);
    pp.query("ep_gamma", pot.ep_gamma);

    if (pot.use_int_pot) {
        if (pot.ip_type == IntPotType::MORSE) {
            if (pot.ip_R_rep <= 0.0) { amrex::Abort("amr.ip_R_rep must be positive"); }
            if (pot.ip_R_att <= pot.ip_R_rep) { amrex::Abort("amr.ip_R_att must be greater than amr.ip_R_rep"); }
            if (pot.ip_eps_rep < 0.0) { amrex::Abort("amr.ip_eps_rep must be non-negative"); }
            if (pot.ip_eps_att <= pot.ip_eps_rep) { amrex::Abort("amr.ip_eps_att must be greater than amr.ip_eps_rep"); }
            if (pot.ip_re < 0.0) { amrex::Abort("amr.ip_re must be non-negative"); }
        } else if (pot.ip_R <= 0.0) {
            amrex::Abort("amr.ip_R must be positive");
        }
        if (pot.ip_type == IntPotType::GEMA && pot.ip_alpha <= 0.0) { amrex::Abort("amr.ip_alpha must be positive"); }
        if (pot.ip_type == IntPotType::PP) {
            if (pot.ip_pp_alpha <= 0.0) { amrex::Abort("amr.ip_pp_alpha must be positive"); }
            if (pot.ip_pp_beta <= 0.0) { amrex::Abort("amr.ip_pp_beta must be positive"); }
            if (pot.ip_pp_a <= 0.0 || pot.ip_pp_a >= 1.0) { amrex::Abort("amr.ip_pp_a must be in (0,1)"); }
        }
    }
    // HK and PP have compact support, so their exact particle cutoff is ip_R
    if (pot.use_int_pot && (pot.ip_type == IntPotType::HK || pot.ip_type == IntPotType::PP)
        && pot.ip_range <= 0.0) {
        pot.ip_range = pot.ip_R;
    }
    if (pot.use_ext_pot && pot.ep_gamma <= 0.0) {
        amrex::Abort("amr.ep_gamma must be positive");
    }

    amrex::Print() << "Interaction potential: use_int_pot = " << pot.use_int_pot
                   << " ip_type = " << ip_type_name;
    if (pot.ip_type == IntPotType::MORSE) {
        amrex::Print() << " ip_eps_att = " << pot.ip_eps_att
                       << " ip_eps_rep = " << pot.ip_eps_rep
                       << " ip_R_att = " << pot.ip_R_att
                       << " ip_R_rep = " << pot.ip_R_rep
                       << " ip_re = " << pot.ip_re;
    } else {
        amrex::Print() << " ip_eps = " << pot.ip_eps
                       << " ip_R = " << pot.ip_R;
    }
    if (pot.ip_type == IntPotType::GEMA) {
        amrex::Print() << " ip_alpha = " << pot.ip_alpha;
    }
    if (pot.ip_type == IntPotType::PP) {
        amrex::Print() << " ip_pp_alpha = " << pot.ip_pp_alpha
                       << " ip_pp_beta = " << pot.ip_pp_beta
                       << " ip_pp_a = " << pot.ip_pp_a;
    }
    amrex::Print() << " ip_range = " << pot.ip_range << "\n";
    amrex::Print() << "External potential: use_ext_pot = " << pot.use_ext_pot
                   << " ep_alpha = " << pot.ep_alpha
                   << " ep_beta = " << pot.ep_beta
                   << " ep_gamma = " << pot.ep_gamma << "\n";

    // Size of the potential dropped by the particle cutoff, relative to the
    // scale of U (|U(0)| for gema, hk and pp)
    if (pot.use_int_pot && pot.ip_type == IntPotType::MORSE) {
        // range at which the attractive tail drops to 1e-5 of the scale, a
        // factor of 10 inside the 1e-4 truncation limit below
        const amrex::Real uscale = ip_Uscale(pot, pot.ip_re + 20.*pot.ip_R_att);
        const amrex::Real r_sugg = pot.ip_re + pot.ip_R_att*std::log(1.e5*pot.ip_eps_att/uscale);
        amrex::Print() << "Morse potential: |U| scale = " << uscale
                       << " suggested amr.ip_range >= " << r_sugg << "\n";
    }
    if (pot.use_int_pot && pot.ip_range > 0.0) {
        const amrex::Real uscale = ip_Uscale(pot, pot.ip_range);
        const amrex::Real ucut = std::abs(ip_U(pot.ip_range, pot));
        const amrex::Real trunc_rel = (uscale > 0.0) ? ucut/uscale : ucut;
        amrex::Print() << "Interaction potential truncation at ip_range: |U(ip_range)| = " << ucut
                       << " |U(ip_range)|/max|U| = " << trunc_rel << "\n";
        if (trunc_rel > 1.e-4) {
            amrex::Abort("amr.ip_range truncates the interaction potential by more than 1e-4 of max|U|; "
                         "increase amr.ip_range");
        }
    }
}

#endif
