#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Random.H>
#include <AMReX_VisMF.H>

#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <AmrCoreFlock.H>
#include <mykernel.H>

using namespace amrex;

// x, y extent from geometry.prob_lo/prob_hi; theta = z is always [0, 2 pi)
RealBox
AmrCoreFlock::MakeRealBox ()
{
    ParmParse pp("geometry");
    Vector<Real> lo(3, 0.0), hi(3, 1.0);
    hi[2] = flock_two_pi;
    Vector<Real> v;
    if (pp.queryarr("prob_lo", v)) {
        for (int d = 0; d < 2 && d < v.size(); ++d) { lo[d] = v[d]; }
        if (v.size() > 2 && std::abs(v[2]) > 1.e-12) {
            Abort("geometry.prob_lo in z (theta) must be 0: theta always spans [0, 2 pi)");
        }
    }
    v.clear();
    if (pp.queryarr("prob_hi", v)) {
        for (int d = 0; d < 2 && d < v.size(); ++d) { hi[d] = v[d]; }
        if (v.size() > 2 && std::abs(v[2] - flock_two_pi) > 1.e-6) {
            Abort("geometry.prob_hi in z (theta) must be 2 pi: theta always spans [0, 2 pi)");
        }
    }
    Vector<int> per;
    if (pp.queryarr("is_periodic", per)) {
        for (int p : per) {
            if (p != 1) { Abort("the flocking SPDE is periodic in x, y and theta; geometry.is_periodic must be 1 1 1"); }
        }
    }
    return RealBox(lo.data(), hi.data());
}

Vector<int>
AmrCoreFlock::ReadNCell ()
{
    ParmParse pp("amr");
    Vector<int> n_cell;
    pp.getarr("n_cell", n_cell);
    if (n_cell.size() != 3) { Abort("amr.n_cell needs three entries: nx ny ntheta"); }
    return n_cell;
}

AmrCoreFlock::AmrCoreFlock ()
    : AmrCore(MakeRealBox(), 0, ReadNCell(), 0, Vector<IntVect>(), Array<int,AMREX_SPACEDIM>{1, 1, 1})
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

AmrCoreFlock::~AmrCoreFlock () = default;

void
AmrCoreFlock::ReadParameters ()
{
    read_flocking_params(fp);
    const RealBox& rb = Geom(0).ProbDomain();
    fp.xlo = rb.lo(0);  fp.lx = rb.length(0);
    fp.ylo = rb.lo(1);  fp.ly = rb.length(1);
    rhobar = 1.0/(fp.lx*fp.ly);
    print_mips_stability(fp);

    {
        ParmParse pp;
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);
        pp.query("seed", seed);
        pp.query("dorand", dorand);
        pp.query("adv_order", adv_order);
        pp.query("limiter", limiter);
        pp.query("time_integrator", time_integrator);
        pp.query("diag_int", diag_int);
        if (adv_order < 0 || adv_order > 2) { Abort("adv_order must be 0 (centred), 1 (upwind) or 2 (MUSCL)"); }
        if (limiter != 0 && limiter != 1) { Abort("limiter must be 0 (minmod) or 1 (MC)"); }
        if (time_integrator < 0 || time_integrator > 2) {
            Abort("time_integrator must be 0 (Euler-Maruyama), 1 (Heun) or 2 (stochastic RK3)");
        }
        if (adv_order == 0 && time_integrator == 0) {
            Abort("adv_order = 0 (centred advection) is unstable with forward Euler; use time_integrator = 2 (or 1)");
        }
    }
    {
        ParmParse pp("amr");
        pp.query("plot_file", plot_file);
        pp.query("plot_int", plot_int);
        pp.query("plot_dt", plot_dt);
        pp.query("chk_file", chk_file);
        pp.query("chk_int", chk_int);
        pp.query("restart", restart_chkfile);
        if (max_level != 0) { Abort("the flocking SPDE is single level: amr.max_level must be 0"); }
    }
    {
        ParmParse pp("adv");
        pp.query("cfl", cfl);
    }

    amrex::Print() << "Flocking SPDE: " << Geom(0).Domain().length(0) << " x "
                   << Geom(0).Domain().length(1) << " x " << Geom(0).Domain().length(2)
                   << " cells, adv_order = " << adv_order << " limiter = " << limiter
                   << " time_integrator = " << time_integrator << " dorand = " << dorand
                   << " cfl = " << cfl << "\n";
}

void
AmrCoreFlock::DefineLevel (const BoxArray& ba, const DistributionMapping& dm)
{
    phi_old.define(ba, dm, 1, nghost);
    phi_new.define(ba, dm, 1, nghost);
    phi_old.setVal(0.0);
    phi_new.setVal(0.0);

    rhot.define(ba, dm, 1, 1);
    rhot.setVal(rhobar);
    dconv.define(Geom(0), ba, dm, fp);
    amrex::Print() << "Sensing kernel: discrete / continuous normalization = "
                   << dconv.discrete_over_continuous() << "\n";
}

void
AmrCoreFlock::UpdateRhoTilde (MultiFab const& phi)
{
    if (fp.speed_type == 0) { return; }
    dconv.apply(phi, rhot);
}

void
AmrCoreFlock::MakeNewLevelFromScratch (int lev, Real /*time*/, const BoxArray& ba,
                                       const DistributionMapping& dm)
{
    AMREX_ALWAYS_ASSERT(lev == 0);
    DefineLevel(ba, dm);

    const auto dx  = Geom(0).CellSizeArray();
    const auto plo = Geom(0).ProbLoArray();
    const FlockingParams fp_loc = fp;
    for (MFIter mfi(phi_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const& phi = phi_new.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            init_phi(i, j, k, phi, dx, plo, fp_loc);
        });
    }

    // normalize so that phi integrates to 1 over (x, y, theta)
    const Real dv = dx[0]*dx[1]*dx[2];
    const Real mass = phi_new.sum(0) * dv;
    phi_new.mult(1.0/mass, 0, 1);
    phi_new.FillBoundary(Geom(0).periodicity());
}

void
AmrCoreFlock::MakeNewLevelFromCoarse (int, Real, const BoxArray&, const DistributionMapping&)
{
    Abort("AmrCoreFlock is single level");
}

void
AmrCoreFlock::RemakeLevel (int, Real, const BoxArray&, const DistributionMapping&)
{
    Abort("AmrCoreFlock is single level");
}

void
AmrCoreFlock::ClearLevel (int /*lev*/)
{
    phi_old.clear();
    phi_new.clear();
}

void
AmrCoreFlock::ErrorEst (int, TagBoxArray&, Real, int)
{
}

void
AmrCoreFlock::InitData ()
{
    if (restart_chkfile.empty()) {
        InitFromScratch(0.0);
        t_new = 0.0;
        step_count = 0;
        if (plot_int > 0 || plot_dt > 0.0) { WritePlotFile(); }
    } else {
        ReadCheckpointFile();
    }
    PrintDiagnostics();
}

Real
AmrCoreFlock::EstTimeStep () const
{
    // Forward Euler with first-order upwind advection and centred theta
    // diffusion keeps phi >= 0 when dt*(v/dx + v/dy + 2D/dtheta^2) <= 1; use
    // the fraction cfl of that (cfl <= 0.5 for MUSCL).
    const auto dxinv = Geom(0).InvCellSizeArray();
    const Real rate = fp.speed*(dxinv[0] + dxinv[1]) + 2.0*fp.diff_coeff*dxinv[2]*dxinv[2];
    return (rate > 0.0) ? cfl/rate : std::numeric_limits<Real>::max();
}

void
AmrCoreFlock::Evolve ()
{
    Real cur_time = t_new;
    Real next_plot_time = (plot_dt > 0.0) ? plot_dt*(std::floor(cur_time/plot_dt + 1.e-9) + 1.0)
                                          : std::numeric_limits<Real>::max();
    int last_plot_step = step_count;

    for (int step = step_count; step < max_step && cur_time < stop_time*(1.0 - 1.e-12); ++step)
    {
        Real dt = EstTimeStep();
        dt = amrex::min(dt, stop_time - cur_time);
        dt = amrex::min(dt, next_plot_time - cur_time);

        AdvancePhi(dt);

        cur_time += dt;
        t_new = cur_time;
        dt_last = dt;
        step_count = step + 1;

        if (diag_int > 0 && step_count % diag_int == 0) { PrintDiagnostics(); }

        bool plot_now = (plot_int > 0 && step_count % plot_int == 0);
        if (plot_dt > 0.0 && cur_time >= next_plot_time*(1.0 - 1.e-12)) {
            plot_now = true;
            next_plot_time += plot_dt;
        }
        if (plot_now) {
            WritePlotFile();
            last_plot_step = step_count;
        }
        if (chk_int > 0 && step_count % chk_int == 0) { WriteCheckpointFile(); }
    }

    if ((plot_int > 0 || plot_dt > 0.0) && step_count > last_plot_step) { WritePlotFile(); }
}

void
AmrCoreFlock::PrintDiagnostics () const
{
    const auto dx = Geom(0).CellSizeArray();
    const Real dv = dx[0]*dx[1]*dx[2];
    amrex::Print() << "STEP " << step_count << " TIME = " << t_new << " DT = " << dt_last
                   << " mass = " << std::setprecision(15) << phi_new.sum(0)*dv
                   << " min phi = " << phi_new.min(0) << std::setprecision(6) << "\n";
}

void
AmrCoreFlock::WritePlotFile ()
{
    const std::string& plotfilename = amrex::Concatenate(plot_file, step_count, 6);
    amrex::Print() << "Writing plotfile " << plotfilename << "\n";
    MultiFab mf(phi_new.boxArray(), phi_new.DistributionMap(), 3, 0);
    MultiFab::Copy(mf, phi_new, 0, 0, 1, 0);

    // rho(x, y) = sum_k phi(x, y, theta_k) dtheta, stored in every theta plane.
    // A box need not span all of theta, so copy phi to a layout whose boxes
    // are whole theta columns, sum there, and copy back.
    const Box& domain = Geom(0).Domain();
    IntVect col_size = maxGridSize(0);
    col_size[2] = domain.length(2);
    BoxArray cba(domain);
    cba.maxSize(col_size);
    DistributionMapping cdm(cba);
    MultiFab col(cba, cdm, 1, 0);
    col.ParallelCopy(phi_new, 0, 0, 1);

    const Real dth = Geom(0).CellSize(2);
    const int klo = domain.smallEnd(2);
    const int khi = domain.bigEnd(2);
    for (MFIter mfi(col); mfi.isValid(); ++mfi) {
        auto const& c = col.array(mfi);
        Box xy = mfi.validbox();
        xy.setRange(2, klo, 1);
        amrex::ParallelFor(xy, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            Real rho = 0.0;
            for (int k = klo; k <= khi; ++k) { rho += c(i,j,k); }
            rho *= dth;
            for (int k = klo; k <= khi; ++k) { c(i,j,k) = rho; }
        });
    }
    mf.ParallelCopy(col, 0, 1, 1);

    // rho_tilde = W * rho of the current state, the same in every theta plane
    // (written for any speed_type, as a diagnostic)
    {
        MultiFab rt(phi_new.boxArray(), phi_new.DistributionMap(), 1, 1);
        dconv.apply(phi_new, rt);
        MultiFab::Copy(mf, rt, 0, 2, 1, 0);
    }

    WriteSingleLevelPlotfile(plotfilename, mf, {"phi", "rho", "rhot"}, Geom(0), t_new, step_count);
}

void
AmrCoreFlock::WriteCheckpointFile () const
{
    const std::string& checkpointname = amrex::Concatenate(chk_file, step_count);
    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", 1, true);

    if (ParallelDescriptor::IOProcessor()) {
        std::string HeaderFileName(checkpointname + "/Header");
        std::ofstream HeaderFile(HeaderFileName.c_str(), std::ofstream::out | std::ofstream::trunc);
        if (!HeaderFile.good()) { amrex::FileOpenFailed(HeaderFileName); }
        HeaderFile.precision(17);
        HeaderFile << "Checkpoint file for AmrCoreFlock\n";
        HeaderFile << step_count << "\n" << dt_last << "\n" << t_new << "\n";
        boxArray(0).writeOn(HeaderFile);
        HeaderFile << '\n';
    }

    VisMF::Write(phi_new, amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "phi"));
}

void
AmrCoreFlock::ReadCheckpointFile ()
{
    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    std::string File(restart_chkfile + "/Header");
    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line;
    std::getline(is, line);   // title
    is >> step_count >> dt_last >> t_new;
    std::getline(is, line);

    BoxArray ba;
    ba.readFrom(is);
    SetBoxArray(0, ba);
    DistributionMapping dm(ba, ParallelDescriptor::NProcs());
    SetDistributionMap(0, dm);
    DefineLevel(ba, dm);

    VisMF::Read(phi_new, amrex::MultiFabFileFullPrefix(0, restart_chkfile, "Level_", "phi"));
    phi_new.FillBoundary(Geom(0).periodicity());
}
