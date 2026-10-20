#ifndef FLOCK_ALIGN_CONV_H_
#define FLOCK_ALIGN_CONV_H_

#include <memory>

#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "FlockingParams.H"

// Periodic convolution of the polarization with the alignment kernel G(|x|/R)
// by FFT, on level 0, following the U * phi convolution of
// src_deankow/interaction:
//
//   Q(x) = sum_x' G(|x - x'|/R) p(x') dx dy,  p = int phi (cos, sin)(theta) dtheta
//
//   2D build: src holds the two components of p on the grid (the
//             particle-mesh polarization), dst gets Q_x, Q_y
//   3D build: src = phi; p is formed on the x-y plane with the weights
//             cos(theta_k) dtheta and sin(theta_k) dtheta, and Q_x, Q_y are
//             copied into every theta plane of dst
//
// G is sampled at the cell offsets (with periodic wrap) and rescaled so that
// sum G dx dy equals the continuous integral 1.002 pi R^2 exactly: the SPDE,
// the particle-mesh and the pair-sum particles then have the same coupling.
class AlignConv
{
public:

    void define (amrex::Geometry const& geom, amrex::BoxArray const& ba,
                 amrex::DistributionMapping const& dm, FlockingParams const& fp);

    bool defined () const { return static_cast<bool>(m_fwd); }

    // dst (2 components, Q_x and Q_y) = G * p(src); ghost cells of dst are
    // filled periodically
    void apply (amrex::MultiFab const& src, amrex::MultiFab& dst);

    // sum of the discrete G samples times dx dy, divided by the continuous
    // integral (1 when the grid resolves G)
    amrex::Real discrete_over_continuous () const { return m_norm_ratio; }

private:

    amrex::Geometry m_geom;
    amrex::BoxArray m_ba;                 // layout of the fields passed to apply
    amrex::DistributionMapping m_dm;
    amrex::BoxArray m_flat_ba;            // 3D: the boxes of m_ba flattened to z = 0
    amrex::BoxArray m_plane_ba;           // layout of the x-y plane the FFT works on
    amrex::DistributionMapping m_plane_dm;
    std::unique_ptr<amrex::FFT::R2C<amrex::Real, amrex::FFT::Direction::forward>>  m_fwd;
    std::unique_ptr<amrex::FFT::R2C<amrex::Real, amrex::FFT::Direction::backward>> m_bwd;
    amrex::MultiFab  m_tmp;               // 2 components on the plane
    amrex::MultiFab  m_one;               // 1 component on the plane, for the FFT
    amrex::cMultiFab m_ghat;
    amrex::cMultiFab m_shat;
    amrex::Real m_norm_ratio = 1.0;
};

#endif
