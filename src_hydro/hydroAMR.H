#ifndef _hydroAMR_H_
#define _hydroAMR_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_AmrCore.H>

#include "common_functions.H"

#ifndef BCFILL_H
#define BCFILL_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

struct AmrCoreFill
{
    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& /*iv*/, amrex::Array4<amrex::Real> const& /*data*/,
                     const int /*dcomp*/, const int /*numcomp*/,
                     amrex::GeometryData const& /*geom*/, const amrex::Real /*time*/,
                     const amrex::BCRec* /*bcr*/, const int /*bcomp*/,
                     const int /*orig_comp*/) const
        {
            // do something for external Dirichlet (BCType::ext_dir)
        }
};



#endif

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
    Vector<std::array< MultiFab, AMREX_SPACEDIM >> stochMfluxdiv;
    Vector<std::array< MultiFab, NUM_EDGE>>  beta_ed;
    Vector<std::array< MultiFab, NUM_EDGE >>  eta_ed;
    Vector<std::array< MultiFab, NUM_EDGE >> temp_ed;
    
    Vector<std::array< MultiFab, AMREX_SPACEDIM >> gmres_rhs_u;
   

    iMultiFab cc_mask;
   
    int npatches;
    int nlevels=2;
    int pg;
    int ng;
    int dt;
    int nplot = 2;

    amrex::Vector<amrex::BCRec> bcs;
    
    bool init = true;
    
//    std::array< iMultiFab, AMREX_SPACEDIM > mask_fc;
//    std::array< iMultiFab, NUM_EDGE > mask_ed;
    
    Vector<MultiFab> pres;
    Vector<MultiFab> gamma;
    Vector<MultiFab> beta;    
    Vector<MultiFab> eta_cc;
    Vector<MultiFab> temp_cc;
    
    Vector<MultiFab> gmres_rhs_p;
    
    Vector<MultiFab> plot_mf;
    
    Vector<Vector<const MultiFab*>> PlotFileMF () const;

//    Vector<Geometry> geom;
//    Vector<Box> domain;
//    Vector<BoxArray> ba;
//    Vector<DistributionMapping> dmap;

    
    hydroAMR(int ang, int * is_periodic, Real dt_in);

    void initData();
    void updateGrid();
    void FillFine();    
    void addPatch(IntVect patch_lo, IntVect patch_hi);
    void advanceStokes();
    void updateAlpha(int lev, Real dt);
    
    void WritePlotFile(Real time, int step);

    void FillCoarsePatchCell(Vector<MultiFab>& mf);
    void FillCoarsePatchFace(Vector<std::array< MultiFab, AMREX_SPACEDIM >>& mf);
    
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
