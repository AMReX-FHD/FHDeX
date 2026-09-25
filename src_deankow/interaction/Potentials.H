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
    GEMA = 0   // generalized exponential model: U = ip_eps*exp(-(r/ip_R)^ip_alpha)
};

struct PotentialParams
{
    // interaction potential
    int         use_int_pot = 1;
    int         ip_type     = IntPotType::GEMA;
    amrex::Real ip_eps      = 666.0;
    amrex::Real ip_R        = 0.1;
    amrex::Real ip_alpha    = 3.0;
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
    case IntPotType::GEMA:
    default:
    {
        const amrex::Real s = std::pow(r/pot.ip_R, pot.ip_alpha);
        return -pot.ip_alpha * pot.ip_eps * std::exp(-s) * s / (r*r);
    }
    }
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
    } else {
        amrex::Abort("Unknown amr.ip_type = " + ip_type_name + "; valid types are: gema");
    }
    pp.query("ip_eps", pot.ip_eps);
    pp.query("ip_R", pot.ip_R);
    pp.query("ip_alpha", pot.ip_alpha);
    pp.query("ip_range", pot.ip_range);

    pp.query("use_ext_pot", pot.use_ext_pot);
    pp.query("ep_alpha", pot.ep_alpha);
    pp.query("ep_beta", pot.ep_beta);
    pp.query("ep_gamma", pot.ep_gamma);

    if (pot.use_int_pot) {
        if (pot.ip_R <= 0.0) { amrex::Abort("amr.ip_R must be positive"); }
        if (pot.ip_alpha <= 0.0) { amrex::Abort("amr.ip_alpha must be positive"); }
    }
    if (pot.use_ext_pot && pot.ep_gamma <= 0.0) {
        amrex::Abort("amr.ep_gamma must be positive");
    }

    amrex::Print() << "Interaction potential: use_int_pot = " << pot.use_int_pot
                   << " ip_type = " << ip_type_name
                   << " ip_eps = " << pot.ip_eps
                   << " ip_R = " << pot.ip_R
                   << " ip_alpha = " << pot.ip_alpha
                   << " ip_range = " << pot.ip_range << "\n";
    amrex::Print() << "External potential: use_ext_pot = " << pot.use_ext_pot
                   << " ep_alpha = " << pot.ep_alpha
                   << " ep_beta = " << pot.ep_beta
                   << " ep_gamma = " << pot.ep_gamma << "\n";

    // Size of the potential dropped by the particle cutoff, relative to U(0)
    if (pot.use_int_pot && pot.ip_range > 0.0) {
        const amrex::Real u0 = std::abs(ip_U(0.0, pot));
        const amrex::Real ucut = std::abs(ip_U(pot.ip_range, pot));
        const amrex::Real trunc_rel = (u0 > 0.0) ? ucut/u0 : ucut;
        amrex::Print() << "Interaction potential truncation at ip_range: |U(ip_range)| = " << ucut
                       << " |U(ip_range)|/|U(0)| = " << trunc_rel << "\n";
        if (trunc_rel > 1.e-4) {
            amrex::Abort("amr.ip_range truncates the interaction potential by more than 1e-4 of U(0); "
                         "increase amr.ip_range");
        }
    }
}

#endif
