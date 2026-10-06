#ifndef AMR_CORE_FLOCK_H_
#define AMR_CORE_FLOCK_H_

#include <string>
#include <limits>

#include <AMReX_AmrCore.H>
#include <AMReX_MultiFab.H>

#include "FlockingParams.H"
#include "DensityConv.H"

// Flocking SPDE on (x, y, theta), with theta = z in [0, 2 pi), periodic in all
// three directions, single level:
//   dphi = D d_theta^2 phi dt - div_x(phi g(theta)) dt + N^(-1/2) d_theta(sqrt(2 D phi) dB)
class AmrCoreFlock
    : public amrex::AmrCore
{
public:

    AmrCoreFlock ();
    ~AmrCoreFlock () override;

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

    static amrex::RealBox MakeRealBox ();
    static amrex::Vector<int> ReadNCell ();

    void ReadParameters ();
    void DefineLevel (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    // one time step of size dt (AdvancePhiAtLevel.cpp)
    void AdvancePhi (amrex::Real dt);

    // rhot = W * phi, the smoothed spatial density that sets the speed; only
    // recomputed when the speed depends on it (speed_type != 0)
    void UpdateRhoTilde (amrex::MultiFab const& phi);

    amrex::Real EstTimeStep () const;

    void PrintDiagnostics () const;
    void WritePlotFile ();
    void WriteCheckpointFile () const;
    void ReadCheckpointFile ();

    FlockingParams fp;

    // phi with 2 ghost cells (MUSCL stencil)
    amrex::MultiFab phi_old;
    amrex::MultiFab phi_new;
    static constexpr int nghost = 2;

    // smoothed density rho_tilde (1 ghost cell), its mean 1/(Lx Ly), and the
    // FFT convolution with the sensing kernel
    amrex::MultiFab rhot;
    amrex::Real rhobar = 1.0;
    DensityConv dconv;

    amrex::Real t_new = 0.0;
    amrex::Real dt_last = 0.0;
    int step_count = 0;

    // discretization
    int adv_order = 1;        // 0: centred, 1: first-order upwind, 2: MUSCL
    int limiter = 0;          // MUSCL limiter: 0 minmod, 1 MC
    int time_integrator = 0;  // 0: Euler-Maruyama, 1: Heun (deterministic part), 2: stochastic RK3
    int dorand = 1;           // 0: no noise (deterministic mean-field equation)
    amrex::Real cfl = 0.5;

    // run control
    int max_step = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();
    int seed = 0;
    int diag_int = 1;

    // output
    std::string plot_file {"plt_spde_"};
    int plot_int = -1;
    amrex::Real plot_dt = -1.0;
    std::string chk_file {"chk_spde_"};
    int chk_int = -1;
    std::string restart_chkfile;
};

#endif
