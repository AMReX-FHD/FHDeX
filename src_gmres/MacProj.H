#ifndef _MacProj_H_
#define _MacProj_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>

#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLLinOp.H>

#include "common_functions.H"

using namespace amrex;

class MacProj {

    MLABecLaplacian mlabec;
    
public:

    MacProj();
    
    void Define(const BoxArray& ba,
                const DistributionMapping& dmap,
                const Geometry& geom);
                
    void Define(const Vector<BoxArray>& ba,
                const Vector<DistributionMapping>& dmap,
                const Vector<Geometry>& geom, int nlev);

    void Solve(std::array<MultiFab, AMREX_SPACEDIM>& alphainv_fc,
               MultiFab& mac_rhs,
               MultiFab& phi,
               Geometry& geom,
               bool full_solve=false);
               
    void Solve(std::array<MultiFab, AMREX_SPACEDIM>* & alphainv_fc,
                    Vector<MultiFab>& mac_rhs,
                    Vector<MultiFab>& phi,
                    Geometry* & geom,
                    bool full_solve=false);
               
   void Solve(std::array<MultiFab, AMREX_SPACEDIM>* & alphainv_fc,
                MultiFab* & mac_rhs,
                MultiFab* & phi,
                Geometry* & geom,
                bool full_solve=false);
               
    int nlevels;
    
};

#endif
