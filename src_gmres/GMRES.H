#ifndef _GMRES_H_
#define _GMRES_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>

#include "common_functions.H"
#include "gmres_functions.H"

#include "StagMGSolver.H"
#include "Precon.H"

using namespace amrex;

class GMRES {

    
    std::array< MultiFab, AMREX_SPACEDIM > r_u;
    std::array< MultiFab, AMREX_SPACEDIM > w_u;
    std::array< MultiFab, AMREX_SPACEDIM > tmp_u;
    std::array< MultiFab, AMREX_SPACEDIM > scr_u;
    std::array< MultiFab, AMREX_SPACEDIM > V_u;
    std::array< MultiFab, AMREX_SPACEDIM > alphainv_fc;

    MultiFab r_p;
    MultiFab w_p;
    MultiFab tmp_p;
    MultiFab scr_p;
    MultiFab V_p;
    
    std::array< MultiFab, AMREX_SPACEDIM > rF_u;
    std::array< MultiFab, AMREX_SPACEDIM > wF_u;
    std::array< MultiFab, AMREX_SPACEDIM > tmpF_u;
    std::array< MultiFab, AMREX_SPACEDIM > scrF_u;
    std::array< MultiFab, AMREX_SPACEDIM > VF_u;
    std::array< MultiFab, AMREX_SPACEDIM > alphainvF_fc;

    MultiFab rF_p;
    MultiFab wF_p;
    MultiFab tmpF_p;
    MultiFab scrF_p;
    MultiFab VF_p;

    StagMGSolver StagSolver;
    Precon Pcon;
    
    int nlevels;

public:

    GMRES (const BoxArray& ba_in,
           const DistributionMapping& dmap_in,
           const Geometry& geom_in);

    GMRES (const BoxArray& ba_in,
           const DistributionMapping& dmap_in,
           const Geometry& geom_in,
           const BoxArray& baF_in,
           const DistributionMapping& dmapF_in,
           const Geometry& geomF_in, int nlev);
           
    void define(const BoxArray& ba_in,
           const DistributionMapping& dmap_in,
           const Geometry& geom_in,
           const BoxArray& baF_in,
           const DistributionMapping& dmapF_in,
           const Geometry& geomF_in, int nlev);           

    void Solve (std::array<MultiFab, AMREX_SPACEDIM> & b_u, MultiFab & b_p,
                std::array<MultiFab, AMREX_SPACEDIM> & x_u, MultiFab & x_p,
                std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
                MultiFab & beta, std::array<MultiFab, NUM_EDGE> & beta_ed,
                MultiFab & gamma,
                Real theta_alpha,
                Geometry & geom,
                Real & norm_pre_rhs, int nlev = 1);
};

#endif
