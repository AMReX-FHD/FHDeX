#include <chrono>
#include <cmath>
#include <cstdio>

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Random.H>

#include "FlockParticles.H"

using namespace amrex;

FlockParticles::FlockParticles ()
{
    ReadParameters();

    if (seed > 0) {
        InitRandom(seed+ParallelDescriptor::MyProc(), ParallelDescriptor::NProcs(),
                   seed+ParallelDescriptor::MyProc());
    } else if (seed == 0) {
        auto now = std::chrono::time_point_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now());
        int randSeed = now.time_since_epoch().count();
        ParallelDescriptor::Bcast(&randSeed, 1, ParallelDescriptor::IOProcessorNumber());
        InitRandom(randSeed+ParallelDescriptor::MyProc(), ParallelDescriptor::NProcs(),
                   randSeed+ParallelDescriptor::MyProc());
    } else {
        Abort("Must supply non-negative seed");
    }
}

void
FlockParticles::ReadParameters ()
{
    read_flocking_params(fp);

    Vector<int> n_cell;
    Vector<int> max_grid_size(1, 64);
    {
        ParmParse pp("amr");
        pp.getarr("n_cell", n_cell);
        pp.queryarr("max_grid_size", max_grid_size);
        pp.query("plot_file", plot_file);
        pp.query("plot_int", plot_int);
        pp.query("plot_dt", plot_dt);
    }
    if (n_cell.size() < 2) { Abort("amr.n_cell needs at least nx ny"); }
    n_theta_bins = (n_cell.size() > 2) ? n_cell[2] : n_theta_bins;

    Vector<Real> lo(2, 0.0), hi(2, 1.0), v;
    {
        ParmParse pp("geometry");
        if (pp.queryarr("prob_lo", v)) { for (int d = 0; d < 2 && d < v.size(); ++d) { lo[d] = v[d]; } }
        v.clear();
        if (pp.queryarr("prob_hi", v)) { for (int d = 0; d < 2 && d < v.size(); ++d) { hi[d] = v[d]; } }
    }
    {
        ParmParse pp;
        pp.query("n_theta_bins", n_theta_bins);
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);
        pp.query("seed", seed);
        pp.query("particle_dt", dt_fixed);
        pp.query("write_particles", write_particles);
    }
    {
        ParmParse pp("adv");
        pp.query("cfl", cfl);
    }
    if (n_theta_bins < 1) { Abort("n_theta_bins must be >= 1"); }

    Box domain(IntVect(AMREX_D_DECL(0, 0, 0)), IntVect(AMREX_D_DECL(n_cell[0]-1, n_cell[1]-1, 0)));
    RealBox rb(lo.data(), hi.data());
    Array<int,AMREX_SPACEDIM> is_per {AMREX_D_DECL(1, 1, 1)};
    geom.define(domain, rb, CoordSys::cartesian, is_per);

    grids.define(domain);
    grids.maxSize(max_grid_size[0]);
    dmap.define(grids);

    fp.xlo = lo[0];  fp.lx = hi[0] - lo[0];
    fp.ylo = lo[1];  fp.ly = hi[1] - lo[1];

    if (dt_fixed <= 0.0) {
        const Real h = amrex::min(geom.CellSize(0), geom.CellSize(1));
        dt_fixed = (fp.speed > 0.0) ? cfl*h/fp.speed : 0.01;
    }

    amrex::Print() << "Flocking particles: " << n_cell[0] << " x " << n_cell[1]
                   << " cells, " << n_theta_bins << " theta bins, dt = " << dt_fixed << "\n";
}

void
FlockParticles::InitData ()
{
    pc = std::make_unique<FlockPC>(geom, dmap, grids);
    pc->InitParticles(fp);
    hist.define(grids, dmap, n_theta_bins + 3, 0);
    t_new = 0.0;
    step_count = 0;
    if (plot_int > 0 || plot_dt > 0.0) { WritePlotFile(); }
}

void
FlockParticles::Evolve ()
{
    Real cur_time = t_new;
    Real next_plot_time = (plot_dt > 0.0) ? plot_dt : std::numeric_limits<Real>::max();
    int last_plot_step = step_count;

    for (int step = step_count; step < max_step && cur_time < stop_time*(1.0 - 1.e-12); ++step)
    {
        Real dt = dt_fixed;
        dt = amrex::min(dt, stop_time - cur_time);
        dt = amrex::min(dt, next_plot_time - cur_time);

        pc->Advance(dt, fp);

        cur_time += dt;
        t_new = cur_time;
        step_count = step + 1;

        bool plot_now = (plot_int > 0 && step_count % plot_int == 0);
        if (plot_dt > 0.0 && cur_time >= next_plot_time*(1.0 - 1.e-12)) {
            plot_now = true;
            next_plot_time += plot_dt;
        }
        if (plot_now) {
            amrex::Print() << "STEP " << step_count << " TIME = " << cur_time << "\n";
            WritePlotFile();
            last_plot_step = step_count;
        }
    }

    if ((plot_int > 0 || plot_dt > 0.0) && step_count > last_plot_step) { WritePlotFile(); }
}

void
FlockParticles::WritePlotFile ()
{
    pc->DepositHistogram(hist, n_theta_bins);

    Vector<std::string> varnames;
    for (int k = 0; k < n_theta_bins; ++k) {
        char name[32];
        std::snprintf(name, sizeof(name), "phi_t%03d", k);
        varnames.push_back(name);
    }
    varnames.push_back("rho");
    varnames.push_back("px");
    varnames.push_back("py");

    const std::string& plotfilename = amrex::Concatenate(plot_file, step_count, 6);
    amrex::Print() << "Writing plotfile " << plotfilename << "\n";
    WriteSingleLevelPlotfile(plotfilename, hist, varnames, geom, t_new, step_count);
    if (write_particles) {
        pc->WritePlotFile(plotfilename, "particles", {"theta"}, {});
    }
}
