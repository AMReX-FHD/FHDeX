#ifndef POT_PARAMS_H_
#define POT_PARAMS_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>

// External potential shared by the SPDE drift (mykernel.H), the particle drift
// (StochasticPC::AdvectWithRandomWalk) and the ghost-cell sampler
// (StochasticPC::AddParticles), so the three descriptions always agree.
//
//   V(x,y) = (x - alpha)^2 (x - beta)^2 + ky (y - 1/2)^4      (ky = 1/2 by default)
//
// The particles follow dX = -(1/gamma) grad V dt + dW (so D = 1/2), the SPDE
// flux carries phi grad V / gamma, and the equilibrium density is
// proportional to exp(-2 V / gamma). The sampler uses exp(-2 dV / gamma)
// across a cell with the cell-averaged slope dV = (V(x+) - V(x-))/dx.
//
// Inputs (prefix adv): ext_pot (0/1), pot_alpha, pot_beta, pot_gamma, pot_ky.
// The force must be single valued on a periodic face: with is_periodic in x
// the domain has to run from one well to the other (prob_lo = alpha, prob_hi =
// beta, where dV/dx = 0), and with is_periodic in y ky must be 0; otherwise the
// two sides of the periodic face see different fluxes and mass is created.
// ReadParameters aborts on such an inconsistent setup.
// pot_sampler_only = 1 (testing only) applies the tilt in the ghost-cell
// sampler while the dynamics have no potential, which is what the code did
// before the potential became an input; it reproduces that bias on demand.
struct PotParams
{
    int ext_pot = 0;
    int sampler_only = 0;
    amrex::Real alpha = 0.25;
    amrex::Real beta  = 0.75;
    amrex::Real gamma = 1.1e-3;
    amrex::Real ky    = 0.5;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real Vx (amrex::Real x) const noexcept
    {
        return (x - alpha)*(x - alpha)*(x - beta)*(x - beta);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real Vy (amrex::Real y) const noexcept
    {
        return ky*(y - 0.5)*(y - 0.5)*(y - 0.5)*(y - 0.5);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real V (amrex::Real x, amrex::Real y) const noexcept { return Vx(x) + Vy(y); }

    // dV/dx and dV/dy
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real dVx (amrex::Real x) const noexcept
    {
        return 2.*(x - beta)*(x - alpha)*(2.*x - alpha - beta);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real dVy (amrex::Real y) const noexcept
    {
        return 4.*ky*(y - 0.5)*(y - 0.5)*(y - 0.5);
    }
};

// Initial condition (prefix adv): ic_type 0 uniform ic_npc particles per cell;
// 1 Boltzmann ic_npc * exp(-2V/gamma) / <exp(-2V/gamma)> (mean ic_npc per
// cell, norm set by the caller); 2 step: ic_npc outside, ic_npc_in inside the
// index box ic_box = i0 j0 [k0] i1 j1 [k1].
struct ICParams
{
    int type = 0;
    amrex::Real npc = 20.;
    amrex::Real npc_in = 0.;
    amrex::Real norm = 1.;   // <exp(-2V/gamma)> over the domain, for type 1
    int box_lo[3] = {0, 0, 0};
    int box_hi[3] = {-1, -1, -1};
};

#endif
