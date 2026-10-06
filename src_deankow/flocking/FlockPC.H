#ifndef FLOCK_PC_H_
#define FLOCK_PC_H_

#include <AMReX_NeighborParticles.H>

#include "FlockingParams.H"

// Real particle attributes
struct FlockRealIdx
{
    enum {
        theta = 0,   // direction of travel, in [0, 2 pi)
        rhot,        // smoothed density rho_tilde at the particle (sets the speed)
        ncomps
    };
};

// 2D particles with a direction theta:
//   dX = v (cos theta, sin theta) dt,  dtheta = sqrt(2 D) dW
// A NeighborParticleContainer, as interaction's StochasticPC, so that the
// alignment interaction can later loop over neighbours.
class FlockPC
    : public amrex::NeighborParticleContainer<FlockRealIdx::ncomps, 0, 0, 0>
{
public:

    using Base = amrex::NeighborParticleContainer<FlockRealIdx::ncomps, 0, 0, 0>;

    FlockPC (const amrex::Geometry& geom, const amrex::DistributionMapping& dmap,
             const amrex::BoxArray& ba, int neighbor_cells = 1)
        : Base(geom, dmap, ba, neighbor_cells), m_neighbor_cells(neighbor_cells)
        {}

    // Sample num_part particles from flocking_ic. The expected count in each
    // cell is num_part times the cell's share of the integral of f (midpoint in
    // x, y; nq points in theta); within a cell x, y are uniform and theta is
    // drawn from f(x, y, .) by rejection.
    void InitParticles (FlockingParams const& fp, int nq = 64);

    // One Euler-Maruyama step: position with the start-of-step theta and speed
    // v s(rhot/rho_bar), then theta += sqrt(2 D dt) Z, wrapped into [0, 2 pi).
    void Advance (amrex::Real dt, FlockingParams const& fp, amrex::Real rhobar);

    // Cloud-in-cell deposit of the particles onto rho (one ghost cell), as a
    // density normalized to integrate to 1 over the plane.
    void DepositCIC (amrex::MultiFab& rho) const;

    // Cloud-in-cell interpolation of rt (ghost cells filled) into each
    // particle's rhot attribute.
    void InterpolateCIC (amrex::MultiFab const& rt);

    // rhot_i = (1/N) sum_j W(|x_i - x_j|) over neighbours (and j = i) with the
    // continuous kernel: the grid-free particle model. Needs neighbor_cells to
    // cover the kernel cutoff.
    void ComputeRhoTildePairs (FlockingParams const& fp);

    // Compare each particle's rhot attribute with the CIC interpolation of rt:
    // returns max |diff|, rms diff and rms rhot over all particles.
    void CompareWithCIC (amrex::MultiFab const& rt, amrex::Real& maxdiff,
                         amrex::Real& rmsdiff, amrex::Real& rmsval) const;

    // Fill hist (on this container's grid) with the particle density per theta
    // bin: component k (k < ntheta) is count/(N dx dy dtheta), an estimate of
    // phi(x, y, theta_k) normalized like the SPDE; then the moments
    // rho = sum_k phi_k dtheta, px = sum_k phi_k cos(theta_k) dtheta and py.
    void DepositHistogram (amrex::MultiFab& hist, int ntheta) const;

private:

    int m_neighbor_cells = 1;
};

#endif
