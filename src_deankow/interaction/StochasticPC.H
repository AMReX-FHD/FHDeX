#ifndef STOCHASTIC_PC_H_
#define STOCHASTIC_PC_H_

#include <AMReX_NeighborParticles.H>

struct RealIdx
{
    enum {
        xold = 0,
        yold, zold,
        ncomps
    };
};

struct IntIdx
{
    enum {
        ncomps = 0
    };
};

/*

  Helper function that converts global 3D index to *local* id index
  in a given box.

 */
struct FlatIndex {
    amrex::Dim3 lo;
    amrex::Dim3 hi;
    FlatIndex (amrex::Box bx)
        : lo(amrex::lbound(bx)), hi(amrex::ubound(bx)) {}

    AMREX_GPU_HOST_DEVICE
    unsigned int operator() (int i, int j, int k) const {
        int ix = i - lo.x;
        int iy = j - lo.y;
        int iz = k - lo.z;
        int nx = hi.x-lo.x+1;
        int ny = hi.y-lo.y+1;
        int nz = hi.z-lo.z+1;
        unsigned int uix = amrex::min(nx-1,amrex::max(0,ix));
        unsigned int uiy = amrex::min(ny-1,amrex::max(0,iy));
        unsigned int uiz = amrex::min(nz-1,amrex::max(0,iz));
        unsigned int cellid = (uix * ny + uiy) * nz + uiz;
        return cellid;
    }
};

class StochasticPC
    : public amrex::NeighborParticleContainer<RealIdx::ncomps, IntIdx::ncomps, 0, 0>
{

public:

    StochasticPC (amrex::ParGDBBase* gdb, int neighbor_cells = 1)
        : amrex::NeighborParticleContainer<RealIdx::ncomps, IntIdx::ncomps, 0, 0>(gdb, neighbor_cells)
        , m_neighbor_cells(neighbor_cells)
        {}

    StochasticPC (const amrex::Geometry            & geom,
                  const amrex::DistributionMapping & dmap,
                  const amrex::BoxArray            & ba,
                  int neighbor_cells = 1)
        : amrex::NeighborParticleContainer<RealIdx::ncomps, IntIdx::ncomps, 0, 0>(geom, dmap, ba, neighbor_cells)
        , m_neighbor_cells(neighbor_cells)
        {}

    void InitParticles (amrex::MultiFab& phi_fine, amrex::Real num_part);
    void  AddParticles (amrex::MultiFab& phi_fine, const amrex::BoxArray& ba_to_exclude,
                        amrex::Real factor, amrex::Real num_part);

    void ColorParticlesWithPhi (amrex::MultiFab const& phi_fine);

    void AdvectParticles (int lev, amrex::Real dt,
                                  amrex::Real interaction_range,
                                  amrex::Real interaction_strength,
                                  amrex::Real interaction_scale,
                                  amrex::Real num_part,
                                  amrex::Real diff_coeff,
                                  int use_ext_pot);

    void AdvectWithRandomWalk (int lev, amrex::Real dt, amrex::Real diff_coeff);

    void RemoveParticlesNotInBA (const amrex::BoxArray& ba_to_keep);
    void RefluxCrseToFine (const amrex::BoxArray& ba_to_keep, amrex::MultiFab& phi_for_reflux);
    void RefluxFineToCrse (const amrex::BoxArray& ba_to_keep, amrex::MultiFab& phi_fine_for_reflux);

protected:
    amrex::ParticleLocator<amrex::DenseBins<amrex::Box> > m_reflux_particle_locator;
    int m_neighbor_cells = 1;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::IntVect getNewCell (StochasticPC::ParticleType const& p,
                           amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                           amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                           const amrex::Box& domain) noexcept;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::IntVect getOldCell (StochasticPC::ParticleType const& p,
                           amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                           amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                           const amrex::Box& domain) noexcept;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::IntVect periodicCorrectOldCell (const amrex::IntVect& old_pos,
                                       const amrex::IntVect& new_pos,
                                       const amrex::GpuArray<int,AMREX_SPACEDIM>& is_periodic,
                                       const amrex::Box& domain) noexcept;
#endif
