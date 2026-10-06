#ifndef FLOCK_PARTICLES_H_
#define FLOCK_PARTICLES_H_

#include <limits>
#include <memory>
#include <string>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "FlockingParams.H"
#include "FlockPC.H"

// Driver for the 2D particle version of the flocking model. Reads the same
// x, y geometry and amr.n_cell as the SPDE; the theta histogram has
// n_theta_bins bins (default: the third entry of amr.n_cell), so its plotfile
// can be compared cell for cell with the SPDE.
class FlockParticles
{
public:

    FlockParticles ();

    void InitData ();
    void Evolve ();

private:

    void ReadParameters ();
    void WritePlotFile ();

    FlockingParams fp;

    amrex::Geometry geom;
    amrex::BoxArray grids;
    amrex::DistributionMapping dmap;
    std::unique_ptr<FlockPC> pc;
    amrex::MultiFab hist;

    int n_theta_bins = 32;
    amrex::Real t_new = 0.0;
    int step_count = 0;

    amrex::Real dt_fixed = -1.0;   // time step; default cfl*min(dx, dy)/v
    amrex::Real cfl = 0.5;
    int max_step = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();
    int seed = 0;

    std::string plot_file {"plt_part_"};
    int plot_int = -1;
    amrex::Real plot_dt = -1.0;
    int write_particles = 1;
};

#endif
