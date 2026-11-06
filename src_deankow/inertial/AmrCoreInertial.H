#ifndef AMR_CORE_INERTIAL_H_
#define AMR_CORE_INERTIAL_H_

#include <limits>
#include <memory>
#include <string>

#include <AMReX_AmrCore.H>
#include <AMReX_MultiFab.H>

#include "InertialParams.H"
#include "InertialPC.H"

// Inertial Dean-Kawasaki model on the phase-space grid (x, v): x is direction
// 0 (periodic), v is direction 1 with reflecting walls at +-Vmax. One class
// runs either the finite-volume SPDE (inertial.run_mode = 0) or the particles
// (run_mode = 1) on the same grid, so the two can be compared cell for cell;
// the particle histogram is the estimate of the cell averages of phi. Single
// level for now; the AmrCore structure is kept so that a hybrid (particles on
// a fine level, SPDE elsewhere, as in src_deankow/multilevel) can be added.
class AmrCoreInertial
    : public amrex::AmrCore
{
public:

    AmrCoreInertial ();
    ~AmrCoreInertial () override;

    void InitData ();
    void Evolve ();

    // AmrCore interface; only level 0 is used
    void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                  const amrex::DistributionMapping& dm) override;
    void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                 const amrex::DistributionMapping& dm) override;
    void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                      const amrex::DistributionMapping& dm) override;
    void ClearLevel (int lev) override;
    void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

private:

    void ReadParameters ();
    void DefineLevel (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    // one SPDE time step of size dt (AdvancePhi.cpp)
    void AdvancePhi (amrex::Real dt);

    amrex::Real EstTimeStep () const;

    // the field written and analysed: phi (SPDE) or the particle histogram
    amrex::MultiFab const& Density ();

    void PrintDiagnostics ();
    void WritePlotFile ();
    void WriteCheckpointFile () const;
    void ReadCheckpointFile ();

    InertialParams ip;

    // SPDE: phi with 2 ghost cells (MUSCL stencil)
    amrex::MultiFab phi_old;
    amrex::MultiFab phi_new;
    static constexpr int nghost = 2;

    // particles and their histogram on the same grid
    std::unique_ptr<InertialPC> pc;
    amrex::MultiFab hist;
    InertialPC::MoveStats move_total;
    int write_particles = 0;

    amrex::Real t_new = 0.0;
    amrex::Real dt_last = 0.0;
    int step_count = 0;

    // discretization (SPDE)
    int adv_order = 1;        // 0: centred, 1: first-order upwind, 2: MUSCL
    int limiter = 0;          // MUSCL limiter: 0 minmod, 1 MC
    int time_integrator = 0;  // 0: Euler-Maruyama, 1: Heun (deterministic part), 2: stochastic RK3
    int dorand = 1;           // 0: no noise (deterministic mean-field equation)
    amrex::Real cfl = 0.5;
    amrex::Real fixed_dt = -1.0;   // if > 0, the time step for both modes

    // run control
    int max_step = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();
    int seed = 0;
    int diag_int = 1;

    // output
    std::string plot_file {"plt_inertial_"};
    int plot_int = -1;
    amrex::Real plot_dt = -1.0;
    std::string chk_file {"chk_inertial_"};
    int chk_int = -1;
    std::string restart_chkfile;
};

#endif
