#ifndef _Precon_H_
#define _Precon_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>

#include "MacProj.H"

using namespace amrex;

class Precon {

    MultiFab phi;
    MultiFab mac_rhs;

    std::array< MultiFab, AMREX_SPACEDIM > gradp;

    MacProj macproj;

public:

    Precon();

    void Define(const BoxArray& ba_in,
                const DistributionMapping& dmap_in,
                const Geometry& geom_in);

    void Apply(const std::array<MultiFab, AMREX_SPACEDIM> & b_u,
               const MultiFab & b_p,
               std::array<MultiFab, AMREX_SPACEDIM> & x_u,
               MultiFab & x_p,
               const std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
               const std::array<MultiFab, AMREX_SPACEDIM> & alphainv_fc,
               const MultiFab & beta, const std::array<MultiFab, NUM_EDGE> & beta_ed,
               const MultiFab & gamma,
               const Real & theta_alpha,
               const Geometry & geom,
               StagMGSolver& StagSolver);
};

#endif
