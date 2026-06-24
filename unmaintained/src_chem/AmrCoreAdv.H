#ifndef AmrCoreAdv_H_
#define AmrCoreAdv_H_

#include <string>
#include <limits>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_BCRec.H>
#include <AMReX_Vector.H>

using namespace amrex;

class AmrCoreAdv
    : public amrex::AmrCore
{
public:

    ////////////////
    // public member functions

    // constructor - reads in parameters from inputs file
    //             - sizes multilevel arrays and data structures
    //AmrCoreAdv ();
    //virtual ~AmrCoreAdv();

    void Initialize();

    // advance solution (of the advection diffusion equation) in time
    void EvolveChem(
            std::array<MultiFab, AMREX_SPACEDIM> & umac_pre,
            std::array<MultiFab, AMREX_SPACEDIM> & umac,
            const iMultiFab & iface_pre,
            const iMultiFab & iface,
            const iMultiFab & catalyst_pre,
            const iMultiFab & catalyst,
            const MultiFab & LevelSet_pre,
            const MultiFab & LevelSet,
            int lev, int step, Real dt_fluid, Real time, Real diffcoeff,
            const Vector<std::array<MultiFab, AMREX_SPACEDIM>> & face_coords,
            int corrector, Real source_strength);

    // initializes multilevel data
    void InitData ( amrex::BoxArray & ba, amrex::DistributionMapping & dm);

    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
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

    //amrex::Vector<std::unique_ptr<MultiFab>> * con_new_copy(int lev, amrex::Vector< std::unique_ptr<MultiFab>> * conNew);
    void con_new_copy(int lev,  amrex::Vector<std::unique_ptr<MultiFab>> & MF, int indicator);

//    amrex::Vector<std::unique_ptr<MultiFab>> con_return();


private:

    ////////////////
    // private member functions

    // read in some parameters from inputs file
    void ReadParameters();

    // set covered coarse cells to be the average of overlying fine cells
    void AverageDown ();

    // more flexible version of AverageDown() that lets you average down across multiple levels
    void AverageDownTo (int crse_lev);

    // compute the number of cells at a level
    long CountCells (int lev);

    // compute a new multifab by coping in con from valid region and filling ghost cells
    // works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
    void FillPatch (int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp);

    // fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
    void FillCoarsePatch (int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp);

    // utility to copy in data from con_old and/or con_new into another multifab
    void GetData (int lev, amrex::Real time, amrex::Vector<amrex::MultiFab*>& data,
                  amrex::Vector<amrex::Real>& datatime);

    // advance a level by dt
    // includes a recursive call for finer levels
//    void timeStep (amrex::Vector< std::unique_ptr<MultiFab> >* uface, amrex::Vector< std::unique_ptr<amrex::MultiFab> >* vface, amrex::Vector< std::unique_ptr<amrex::MultiFab> >* wface, int lev, amrex::Real time, int iteration);
    void timeStep ( int lev, amrex::Real time, int iteration);

    //
    //give the location of sources and their strength
    //
    amrex::Real sourcelist(amrex::Vector<int> & xloc,
                                   amrex::Vector<int> & yloc,
                                   amrex::Vector<int> & zloc,
                                   amrex::Vector<amrex::Real> & strength) ;
    // advance a single level for a single time step, updates flux registers
//    void Advance (amrex::Vector< std::unique_ptr<MultiFab> >* uface, amrex::Vector< std::unique_ptr<amrex::MultiFab> >* vface, amrex::Vector< std::unique_ptr<amrex::MultiFab> >* wface, int lev, amrex::Real time, amrex::Real dt_lev, int iteration, int ncycle);
    void Advance ( int lev, amrex::Real time, amrex::Real dt_lev, int iteration, int ncycle);

    // a wrapper for EstTimeStep(0
    void ComputeDt ();

    // compute dt from CFL considerations
    Real EstTimeStep (int lev, bool local=false) const;

    // get plotfile name
    std::string PlotFileName (int lev) const;

    // put together an array of multifabs for writing
    amrex::Vector<const amrex::MultiFab*> PlotFileMF () const;

    // set plotfile variables names
    amrex::Vector<std::string> PlotFileVarNames () const;

    // write plotfile to disk
    void WritePlotFile () const;

    // write checkpoint file to disk
    void WriteCheckpointFile () const;

    // read checkpoint file from disk
    void ReadCheckpointFile ();

    // utility to skip to next line in Header
    static void GotoNextLine (std::istream& is);

    ////////////////
    // private data members

    amrex::Vector<int> istep;      // which step?
    amrex::Vector<int> nsubsteps;  // how many substeps on each level?

    // keep track of old time, new time, and time step at each level
    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;

    // array of multifabs to store the solution at each level of refinement
    amrex::Vector<std::unique_ptr<amrex:: MultiFab>> con_pre;
    amrex::Vector<std::unique_ptr<amrex:: MultiFab>> con_new;
    amrex::Vector<std::unique_ptr<amrex:: MultiFab>> con_old;

   // array of multifabs to store the fluid velocity at each level of refinement

   Vector< std::unique_ptr<MultiFab> >  uface;
   Vector< std::unique_ptr<MultiFab> >  vface;
   Vector< std::unique_ptr<MultiFab> >  wface;

   Vector< std::unique_ptr<MultiFab> >  uface_pre;
   Vector< std::unique_ptr<MultiFab> >  vface_pre;
   Vector< std::unique_ptr<MultiFab> >  wface_pre;

   // array of multifabs to store the surface gradient at each level of refinement (face centered and cell centered) and the magnitude of the gradient (cell centerd)

   Vector< std::unique_ptr<MultiFab> >  Dcon_x;
   Vector< std::unique_ptr<MultiFab> >  Dcon_y;
   Vector< std::unique_ptr<MultiFab> >  Dcon_z;

   Vector< std::unique_ptr<MultiFab> >  Dconc_x;
   Vector< std::unique_ptr<MultiFab> >  Dconc_y;
   Vector< std::unique_ptr<MultiFab> >  Dconc_z;

   Vector< std::unique_ptr<MultiFab> >  MagDcon;
   // array of multifabs to store the face centerd positions at each level of refinement

   Vector< std::unique_ptr<MultiFab> >  xface;
   Vector< std::unique_ptr<MultiFab> >  yface;
   Vector< std::unique_ptr<MultiFab> >  zface;
   // multifabs to store the location of the point sources, the location of th einterface, and the levelset
   std::unique_ptr<iMultiFab>           source_loc;
   std::unique_ptr<iMultiFab>           source_loc_pre;

   std::unique_ptr<iMultiFab>           interface_loc;
   std::unique_ptr<iMultiFab>           interface_loc_pre;
   std::unique_ptr<MultiFab>            levset;
   std::unique_ptr<MultiFab>            levset_pre;

    // condition types at the lo/hi walls in each direction
    //BCRec bcs;
   amrex::Vector<BCRec> bcs;
    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] and flux_reg[nlevs_max] are never actually
    // used in the reflux operation
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg;

    ////////////////
    // runtime parameters

    // maximum number of steps and stop time
    int max_step = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

    // if >= 0 we restart from a checkpoint
    std::string restart_chkfile = "";

    // advective cfl number - dt = cfl*dx/umax
    amrex::Real cfl = 0.25;

    // diffusion coefficient
    amrex::Real diffcoeff;
   // source strength
    amrex::Real strength;

    // if we are doing a predictor or corrector step
    int Correct=0;
   // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = 2;

    // hyperbolic refluxing as part of multilevel synchronization
    int do_reflux = 1;

    // plotfile prefix and frequency
    std::string plot_file {"plt"};
    int plot_int = -1;

    // checkpoint prefix and frequency
    std::string chk_file {"chk"};
    int chk_int = -1;
};

#endif
