#ifndef INERTIAL_PC_H_
#define INERTIAL_PC_H_

#include <AMReX_Particles.H>

#include "InertialParams.H"

// Particles in phase space: pos(0) = x (periodic), pos(1) = v (reflecting
// walls at +-Vmax), on the same (x, v) grid as the SPDE.
//   dX = v dt,  dv = -gamma v dt - U'(X) dt + sqrt(2 D) dW
class InertialPC
    : public amrex::ParticleContainer<0, 0, 0, 0>
{
public:

    using Base = amrex::ParticleContainer<0, 0, 0, 0>;

    InertialPC (const amrex::Geometry& geom, const amrex::DistributionMapping& dmap,
                const amrex::BoxArray& ba)
        : Base(geom, dmap, ba)
        {}

    // Sample num_part particles from inertial_ic: each cell gets about
    // num_part times its share of f at the cell centre (the same quadrature as
    // the SPDE's initial cell values), and x, v are uniform within the cell.
    void InitParticles (InertialParams const& ip);

    // Counters of the last Advance: moves truncated to one cell in x and in
    // v, and reflections at the v walls
    struct MoveStats { amrex::Long trunc_x = 0, trunc_v = 0, reflect = 0; };

    // One step, by inertial.particle_integrator:
    //   0 Euler-Maruyama: dx = v dt, dv = (-gamma v - U'(x)) dt + sqrt(2 D dt) Z
    //     with the start-of-step x and v;
    //   1 BAOAB: v -= dt/2 U'(x); x += dt/2 v; v = c v + sqrt(T (1 - c^2)) Z
    //     with c = exp(-gamma dt), T = D/gamma; x += dt/2 v; v -= dt/2 U'(x).
    // The total move is then truncated
    // to at most one cell (dx, dv), velocities beyond +-Vmax are reflected,
    // and the particles are redistributed (periodic in x).
    MoveStats Advance (amrex::Real dt, InertialParams const& ip);

    // hist (component 0) = count/(N dx dv) per cell: the particle estimate of
    // the cell average of phi, normalized like the SPDE
    void DepositHistogram (amrex::MultiFab& hist) const;

    // Moments over the particles: <x - x0>, <v>, <(x - x0)^2>, <v^2>, with
    // x - x0 the minimum image from the potential centre
    void Moments (InertialParams const& ip, amrex::Real m[4]) const;
};

#endif
