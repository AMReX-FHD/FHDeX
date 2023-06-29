#ifndef _Precon_H_
#define _Precon_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>

#include "MacProj.H"

using namespace amrex;

class Precon {

    Vector<MultiFab> phi;
    Vector<MultiFab> mac_rhs;

    Vector<std::array< MultiFab, AMREX_SPACEDIM >> gradp;
    
    int nlevels;
    
    MacProj macproj;

public:

    Precon();
   
    void Define(const BoxArray& ba_in,
                const DistributionMapping& dmap_in,
                const Geometry& geom_in);
                
    void Define(const Vector<BoxArray>& ba_in,
              const Vector<DistributionMapping>& dmap_in,
              const Vector<Geometry>& geom_in);


    void Apply(std::array<MultiFab, AMREX_SPACEDIM> & b_u,
               MultiFab & b_p,
               std::array<MultiFab, AMREX_SPACEDIM> & x_u,
               MultiFab & x_p,
               std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
               std::array<MultiFab, AMREX_SPACEDIM> & alphainv_fc,
               MultiFab & beta, std::array<MultiFab, NUM_EDGE> & beta_ed,
               MultiFab & gamma,
               Real & theta_alpha,
               Geometry & geom,
               StagMGSolver& StagSolver);
               
   void Apply(std::array<MultiFab, AMREX_SPACEDIM>* & b_u,
                   MultiFab*                             & b_p,
                   Vector<std::array<MultiFab, AMREX_SPACEDIM>> & x_u,
                   Vector<MultiFab>                             & x_p,
                   std::array<MultiFab, AMREX_SPACEDIM>* & alpha_fc,
                   Vector<std::array<MultiFab, AMREX_SPACEDIM>> & alphainv_fc,
                   MultiFab*                             & beta, 
                   std::array<MultiFab, NUM_EDGE>*      & beta_ed,
                   MultiFab*                            & gamma,
                   MultiFab*                            & cc_mask,
                   std::array<MultiFab, AMREX_SPACEDIM>* & fc_mask,                                    
                   const Real                           & theta_alpha,
                   Geometry*                            & geom,
                   StagMGSolver& StagSolver);
                   
    void Apply    (Vector<std::array<MultiFab, AMREX_SPACEDIM>> & b_u,
                   Vector<MultiFab>                             & b_p,
                   Vector<std::array<MultiFab, AMREX_SPACEDIM>> & x_u,
                   Vector<MultiFab>                             & x_p,
                   std::array<MultiFab, AMREX_SPACEDIM>* & alpha_fc,
                   Vector<std::array<MultiFab, AMREX_SPACEDIM>> & alphainv_fc,
                   MultiFab*                             & beta, 
                   std::array<MultiFab, NUM_EDGE>*      & beta_ed,
                   MultiFab*                            & gamma,
                   MultiFab*                            & cc_mask,
                   std::array<MultiFab, AMREX_SPACEDIM>* & fc_mask,                   
                   const Real                           & theta_alpha,
                   Geometry*                            & geom,
                   StagMGSolver& StagSolver);
               
   void Apply(std::array<MultiFab, AMREX_SPACEDIM>* & b_u,
                   MultiFab*                             & b_p,
                   std::array<MultiFab, AMREX_SPACEDIM>* & x_u,
                   MultiFab*                             & x_p,
                   std::array<MultiFab, AMREX_SPACEDIM>* & alpha_fc,
                   std::array<MultiFab, AMREX_SPACEDIM>* & alphainv_fc,
                   MultiFab*                             & beta, 
                   std::array<MultiFab, NUM_EDGE>*      & beta_ed,
                   MultiFab*                            & gamma,
                   MultiFab*                            & cc_mask,
                   std::array<MultiFab, AMREX_SPACEDIM>* & fc_mask,                     
                   const Real                           & theta_alpha,
                   Geometry*                            & geom,
                   StagMGSolver& StagSolver);
};

#endif
