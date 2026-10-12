#ifndef FLOCK_DENSITY_CONV_H_
#define FLOCK_DENSITY_CONV_H_

#include <memory>

#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "FlockingParams.H"

// Periodic convolution with the sensing kernel W(x, y) by FFT, on level 0,
// following the U * phi convolution of src_deankow/interaction.
//
//   2D build: dst(x) = sum_x' W(x - x') src(x') dx dy, with src the spatial
//             density on the grid (the particle-mesh rho_tilde)
//   3D build: dst(x, theta) = sum_x' W(x - x') rho(x') dx dy = rho_tilde(x), the
//             same in every theta plane, with src = phi and
//             rho = sum_theta phi dtheta: phi is integrated over theta first
//             and the FFT convolution is done on the x-y plane only
//
// W is sampled at the cell offsets (with periodic wrap) and normalized so that
// sum W dx dy = 1 exactly; the mean of rho_tilde is then the mean of rho.
class DensityConv
{
public:

    void define (amrex::Geometry const& geom, amrex::BoxArray const& ba,
                 amrex::DistributionMapping const& dm, FlockingParams const& fp);

    bool defined () const { return static_cast<bool>(m_fwd); }

    // dst (1 component) = W * src; ghost cells of dst are filled periodically
    void apply (amrex::MultiFab const& src, amrex::MultiFab& dst);

    // sum of the discrete W samples times dx dy, divided by the continuous
    // integral of the truncated shape (1 when the grid resolves W)
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
    amrex::MultiFab  m_tmp;
    amrex::cMultiFab m_what;
    amrex::cMultiFab m_shat;
    amrex::Real m_norm_ratio = 1.0;
};

#endif
