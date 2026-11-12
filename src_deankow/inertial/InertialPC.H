#ifndef INERTIAL_PC_H_
#define INERTIAL_PC_H_

#include <array>

#include <AMReX_iMultiFab.H>
#include <AMReX_Particles.H>

#include "InertialParams.H"

// Real particle attributes: the start-of-step position, set by Advance and
// used by the hybrid to count crossings of the particle-region boundary
struct InertialRealIdx
{
    enum {
        x_old = 0,
        v_old,
        ncomps
    };
};

// Particles in phase space: pos(0) = x (periodic), pos(1) = v (reflecting
// walls at +-Vmax), on the same (x, v) grid as the SPDE.
//   dX = v dt,  dv = -gamma v dt - U'(X) dt + sqrt(2 D) dW
// In the hybrid, the container's BoxArray is the particle region (level 1)
// followed by its ghost ring, refilled from the SPDE every step.
class InertialPC
    : public amrex::ParticleContainer<InertialRealIdx::ncomps, 0, 0, 0>
{
public:

    using Base = amrex::ParticleContainer<InertialRealIdx::ncomps, 0, 0, 0>;

    InertialPC (const amrex::Geometry& geom, const amrex::DistributionMapping& dmap,
                const amrex::BoxArray& ba)
        : Base(geom, dmap, ba)
        {}

    // Sample num_part particles from inertial_ic: each cell gets about
    // num_part times its share of f at the cell centre (the same quadrature as
    // the SPDE's initial cell values), and x, v are uniform within the cell.
    void InitParticles (InertialParams const& ip);

    // Add particles in the cells where sel != 0, from the density phi on this
    // container's BoxArray: n = floor(N phi dx dv + U) per cell (so 1.3 gives 1
    // with probability 0.7 and 2 with 0.3); x from the density
    // ~ exp(-U'(x_c) (x - x_left)/T) across the cell (inverse CDF, uniform when
    // U' is negligible), v uniform across the cell. Returns the number added.
    amrex::Long SampleCells (amrex::MultiFab const& phi, amrex::iMultiFab const& sel,
                             InertialParams const& ip);

    // Counters of the last Advance: moves truncated to one cell in x and in
    // v, and reflections at the v walls
    struct MoveStats { amrex::Long trunc_x = 0, trunc_v = 0, reflect = 0; };

    // One step, by inertial.particle_integrator:
    //   0 Euler-Maruyama: dx = v dt, dv = (-gamma v - U'(x)) dt + sqrt(2 D dt) Z
    //     with the start-of-step x and v;
    //   1 BAOAB: v -= dt/2 U'(x); x += dt/2 v; v = c v + sqrt(T (1 - c^2)) Z
    //     with c = exp(-gamma dt), T = D/gamma; x += dt/2 v; v -= dt/2 U'(x).
    // Z is per-thread (particle_rng 0) or counter based on (id, cpu, step) (1).
    // The total move is then truncated to at most one cell (dx, dv),
    // velocities beyond +-Vmax are reflected, the start-of-step position is
    // kept in x_old, v_old, and (if redistribute) the particles are
    // redistributed (periodic in x).
    MoveStats Advance (amrex::Real dt, InertialParams const& ip, amrex::Long step = 0,
                       bool redistribute = true);

    // Hybrid bookkeeping after Advance(..., false), with pmask = 1 in the
    // particle region (ghost cells filled periodically) on this BoxArray:
    //   start outside, end inside: -1 particle in xfer at the start cell;
    //   start inside, end outside: +1 particle in xfer at the end cell, and
    //     the particle is removed;
    //   start outside, end outside: removed (an unused ghost particle).
    // xfer (one ghost cell) is in particle counts; the caller sums the ghost
    // cells and scales by 1/(N dx dv). Returns the numbers of entries and exits
    // whose move changed the x cell (0, 2) or only the v cell (1, 3).
    std::array<amrex::Long, 4> DepositCrossings (amrex::iMultiFab const& pmask, amrex::MultiFab& xfer);

    // Remove the particles whose cell has keep == 0 (keep on this BoxArray)
    void RemoveWhere (amrex::iMultiFab const& keep);

    // hist (component 0) = count/(N dx dv) per cell: the particle estimate of
    // the cell average of phi, normalized like the SPDE. N is nnorm if > 0,
    // else the current number of particles.
    void DepositHistogram (amrex::MultiFab& hist, amrex::Real nnorm = -1.0) const;

    // Moments over the particles: <x - x0>, <v>, <(x - x0)^2>, <v^2>, with
    // x - x0 the minimum image from the potential centre
    void Moments (InertialParams const& ip, amrex::Real m[4]) const;
};

#endif
