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

    
/*    std::array< MultiFab, AMREX_SPACEDIM > rF_u;*/
/*    std::array< MultiFab, AMREX_SPACEDIM > wF_u;*/
/*    std::array< MultiFab, AMREX_SPACEDIM > tmpF_u;*/
/*    std::array< MultiFab, AMREX_SPACEDIM > scrF_u;*/
/*    std::array< MultiFab, AMREX_SPACEDIM > VF_u;*/
/*    std::array< MultiFab, AMREX_SPACEDIM > alphainvF_fc;*/

/*    MultiFab rF_p;*/
/*    MultiFab wF_p;*/
/*    MultiFab tmpF_p;*/
/*    MultiFab scrF_p;*/
/*    MultiFab VF_p;*/

    StagMGSolver StagSolver;
    Precon Pcon;
    
    int nlevels;

public:

    Vector<std::array< MultiFab, AMREX_SPACEDIM >> tmp_u;
    Vector<std::array< MultiFab, AMREX_SPACEDIM >> r_u;
    Vector<std::array< MultiFab, AMREX_SPACEDIM >> w_u;

    Vector<std::array< MultiFab, AMREX_SPACEDIM >> scr_u;
    Vector<std::array< MultiFab, AMREX_SPACEDIM >> V_u;
    Vector<std::array< MultiFab, AMREX_SPACEDIM >> alphainv_fc;

    Vector<MultiFab> r_p;
    Vector<MultiFab> w_p;
    Vector<MultiFab> tmp_p;
    Vector<MultiFab> scr_p;
    Vector<MultiFab> V_p;

    GMRES (const BoxArray& ba_in,
           const DistributionMapping& dmap_in,
           const Geometry& geom_in);

    GMRES (const Vector<BoxArray>& ba_in,
              const Vector<DistributionMapping>& dmap_in,
              const Vector<Geometry>& geom_in, int nlev);
           
    void define(const Vector<BoxArray>& ba_in,
              const Vector<DistributionMapping>& dmap_in,
              const Vector<Geometry>& geom_in, int nlev);   
              
                      
    void Solve (std::array<MultiFab, AMREX_SPACEDIM>* & b_u, MultiFab* & b_p,
                   std::array<MultiFab, AMREX_SPACEDIM>* & x_u, MultiFab* & x_p,
                   std::array<MultiFab, AMREX_SPACEDIM>* & alpha_fc,
                   MultiFab* & beta, std::array<MultiFab, NUM_EDGE>* & beta_ed,
                   MultiFab* & gamma,
                   MultiFab* & cc_mask,                   
                   std::array<MultiFab, AMREX_SPACEDIM>* & fc_mask,                   
                   Real theta_alpha,
                   Geometry* & geom,
                   Real & norm_pre_rhs, int nlev);

    void Solve (std::array<MultiFab, AMREX_SPACEDIM> & b_u, MultiFab & b_p,
                   std::array<MultiFab, AMREX_SPACEDIM> & x_u, MultiFab & x_p,
                   std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
                   MultiFab & beta, std::array<MultiFab, NUM_EDGE> & beta_ed,
                   MultiFab & gamma,
                   Real theta_alpha,
                   Geometry & geom,
                   Real & norm_pre_rhs);
                
   void Solve (Vector<std::array<MultiFab, AMREX_SPACEDIM>> & b_u, Vector<MultiFab> & b_p,
                   Vector<std::array<MultiFab, AMREX_SPACEDIM>> & x_u, Vector<MultiFab> & x_p,
                   Vector<std::array<MultiFab, AMREX_SPACEDIM>> & alpha_fc,
                   Vector<MultiFab> & beta, Vector<std::array<MultiFab, NUM_EDGE>> & beta_ed,
                   Vector<MultiFab> & gamma,
                   Vector<MultiFab> & cc_mask,
                   Vector<std::array<MultiFab, AMREX_SPACEDIM>> & fc_mask,                   
                   Real theta_alpha,
                   Vector<Geometry> & geom,
                   Real & norm_pre_rhs);
};

#endif
