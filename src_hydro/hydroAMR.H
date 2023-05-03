#ifndef _hydroAMR_H_
#define _hydroAMR_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_AmrCore.H>

#include "common_functions.H"

using namespace amrex;

class hydroAMR 
    : public amrex::AmrCore
{
public:

    Vector<std::array< MultiFab, AMREX_SPACEDIM >> umac;
    Vector<std::array< MultiFab, AMREX_SPACEDIM >> umacNew;    
    Vector<std::array< MultiFab, AMREX_SPACEDIM >> umacM;
    Vector<std::array< MultiFab, AMREX_SPACEDIM >> alpha_fc;
    Vector<std::array< MultiFab, AMREX_SPACEDIM >> source;
    Vector<std::array< MultiFab, AMREX_SPACEDIM >> sourceTemp;
    Vector<std::array< MultiFab, AMREX_SPACEDIM >> sourceRFD;
    Vector<std::array< MultiFab, NUM_EDGE>>  beta_ed;
    Vector<std::array< MultiFab, NUM_EDGE >>  eta_ed;
    Vector<std::array< MultiFab, NUM_EDGE >> temp_ed;

    iMultiFab cc_mask;
   
    int npatches;
    int nlevels=2;
    int pg;
    int ng;
    int dt;
    
    bool init = true;
    
//    std::array< iMultiFab, AMREX_SPACEDIM > mask_fc;
//    std::array< iMultiFab, NUM_EDGE > mask_ed;
    
    Vector<MultiFab> pres;
    Vector<MultiFab> gamma;
    Vector<MultiFab> beta;    
    Vector<MultiFab> eta_cc;
    Vector<MultiFab> temp_cc;    

//    Vector<Geometry> geom;
//    Vector<Box> domain;
//    Vector<BoxArray> ba;
//    Vector<DistributionMapping> dmap;

    
    hydroAMR(int ang, int * is_periodic, Real dt_in);

    void initData();
    void updateGrid();
    
    void addPatch(IntVect patch_lo, IntVect patch_hi);
    
    virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping and
    // fill with existing fine and coarse data.
    // overrides the pure virtual function in AmrCore
    virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                              const amrex::DistributionMapping& dm) override;

    // Delete level data
    // overrides the pure virtual function in AmrCore
    virtual void ClearLevel (int lev) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping.
    // Only used during initialization.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                          const amrex::DistributionMapping& dm) override;

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;
    
};

#endif
