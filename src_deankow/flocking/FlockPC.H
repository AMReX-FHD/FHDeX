#ifndef FLOCK_PC_H_
#define FLOCK_PC_H_

#include <AMReX_NeighborParticles.H>

#include "FlockingParams.H"

// Real particle attributes
struct FlockRealIdx
{
    enum {
        theta = 0,   // direction of travel, in [0, 2 pi)
        qx,          // smoothed polarization Q = G_R * p at the particle
        qy,          //   (sets the turning rate)
        ncomps
    };
};

// 2D particles with a direction theta:
//   dX = v (cos theta, sin theta) dt,
//   dtheta = eps (cos theta Q_y - sin theta Q_x) dt + sqrt(2 D) dW
// A NeighborParticleContainer, as interaction's StochasticPC, so that the
// alignment can loop over neighbours (align_method 1).
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

    // One Euler-Maruyama step: position at speed v with the start-of-step
    // theta, then theta += omega dt + sqrt(2 D dt) Z with the turning rate
    // omega = eps (cos theta Q_y - sin theta Q_x) from the qx, qy attributes,
    // wrapped into [0, 2 pi).
    void Advance (amrex::Real dt, FlockingParams const& fp);

    // Cloud-in-cell deposit of the polarization (cos theta, sin theta) onto p
    // (2 components, one ghost cell), as a density normalized like phi:
    // p = (1/(N dx dy)) sum_i (cos, sin)(theta_i) CIC weight.
    void DepositPolarizationCIC (amrex::MultiFab& p) const;

    // Cloud-in-cell interpolation of q (2 components, ghost cells filled) into
    // each particle's qx, qy attributes.
    void InterpolateCIC (amrex::MultiFab const& q);

    // Q_i = (1/N) sum_j G(|x_i - x_j|/R) (cos, sin)(theta_j) over neighbours
    // (and j = i) with the continuous kernel: the grid-free particle model.
    // Needs neighbor_cells to cover the kernel cutoff 1.1 R.
    void ComputeAlignPairs (FlockingParams const& fp);

    // Compare each particle's qx, qy with the CIC interpolation of q: returns
    // max |diff|, rms diff and rms |Q| over all particles.
    void CompareWithCIC (amrex::MultiFab const& q, amrex::Real& maxdiff,
                         amrex::Real& rmsdiff, amrex::Real& rmsval) const;

    // eps max_i |Q_i|, the bound on the turning rate
    amrex::Real MaxTurningRate (FlockingParams const& fp) const;

    // Fill hist (on this container's grid) with the particle density per theta
    // bin: component k (k < ntheta) is count/(N dx dy dtheta), an estimate of
    // phi(x, y, theta_k) normalized like the SPDE; then the moments
    // rho = sum_k phi_k dtheta, px = sum_k phi_k cos(theta_k) dtheta and py.
    void DepositHistogram (amrex::MultiFab& hist, int ntheta) const;

private:

    int m_neighbor_cells = 1;
};

#endif
