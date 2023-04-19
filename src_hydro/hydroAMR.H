#ifndef _hydroAMR_H_
#define _hydroAMR_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>

#include "common_functions.H"

using namespace amrex;

class hydroAMR {


public:

    Vector<std::array< MultiFab, AMREX_SPACEDIM > > alpha_fc_patches;
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >   source_patches;
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >   sourceTemp_patches;
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >   sourceRFD_patches;
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >     umac_patches;
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >     umacNew_patches;
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >    umacM_patches;
    Vector<std::array< MultiFab, NUM_EDGE       > >  beta_ed_patches; // nodal in 2D, edge in 3D
    Vector<std::array< MultiFab, NUM_EDGE       > >  eta_ed_patches;
    Vector<std::array< MultiFab, NUM_EDGE       > >  temp_ed_patches;        
    

    Vector<MultiFab>  beta_patches;
    Vector<MultiFab> gamma_patches;
    Vector<MultiFab>  pres_patches;

    Vector<Geometry>  geom_patches;
    Vector<Box>  domain_patches;
    Vector<BoxArray>  ba_patches;
    Vector<DistributionMapping>  dmap_patches;        
   
    int npatches;
    
    std::array< MultiFab, AMREX_SPACEDIM > umac;
    std::array< MultiFab, AMREX_SPACEDIM > umacNew;    
    std::array< MultiFab, AMREX_SPACEDIM > umacM;
    std::array< MultiFab, AMREX_SPACEDIM > alpha_fc;
    std::array< MultiFab, AMREX_SPACEDIM > source;
    std::array< MultiFab, AMREX_SPACEDIM > sourceTemp;
    std::array< MultiFab, AMREX_SPACEDIM > sourceRFD;
    std::array< MultiFab, NUM_EDGE>  beta_ed;
    std::array< MultiFab, NUM_EDGE >  eta_ed;
    std::array< MultiFab, NUM_EDGE > temp_ed;
    
    MultiFab pres;
    MultiFab gamma;
    MultiFab beta;    
    MultiFab eta_cc;
    MultiFab temp_cc;    

    Geometry geom;
    Box domain;
    BoxArray ba;
    DistributionMapping dmap;

    
    hydroAMR(int ang, int * is_periodic, Real dt);
    
    void addPatch(IntVect patch_lo, IntVect patch_hi);
    
};

#endif
