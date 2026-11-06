#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Random.H>
#include <AMReX_VisMF.H>

#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <AmrCoreInertial.H>
#include <mykernel.H>
#include <myfunc.H>

using namespace amrex;

AmrCoreInertial::AmrCoreInertial ()
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

AmrCoreInertial::~AmrCoreInertial () = default;

void
AmrCoreInertial::ReadParameters ()
{
    read_inertial_params(ip);

    const Geometry& g = Geom(0);
    if (!g.isPeriodic(0) || g.isPeriodic(1)) {
        Abort("geometry.is_periodic must be 1 0: periodic in x, reflecting walls in v");
    }
    const RealBox& rb = g.ProbDomain();
    if (std::abs(rb.lo(1) + rb.hi(1)) > 1.e-12*rb.length(1)) {
        Abort("the v domain must be symmetric: geometry.prob_lo[1] = -geometry.prob_hi[1] = -Vmax");
    }
    ip.xlo  = rb.lo(0);
    ip.lx   = rb.length(0);
    ip.vmax = rb.hi(1);
    finish_inertial_params(ip);

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
        if (ip.run_mode == 0 && adv_order == 0 && time_integrator == 0) {
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
        if (max_level != 0) { Abort("the inertial code is single level for now: amr.max_level must be 0"); }
    }
    {
        ParmParse pp("adv");
        pp.query("cfl", cfl);
    }
    {
        ParmParse pp("inertial");
        pp.query("fixed_dt", fixed_dt);
        pp.query("write_particles", write_particles);
    }
    if (ip.run_mode == 1 && !restart_chkfile.empty()) {
        Abort("restart is only implemented for the SPDE (inertial.run_mode = 0)");
    }

    amrex::Print() << "  grid " << g.Domain().length(0) << " (x) x " << g.Domain().length(1)
                   << " (v), dx = " << g.CellSize(0) << ", dv = " << g.CellSize(1);
    if (ip.run_mode == 0) {
        amrex::Print() << ", adv_order = " << adv_order << " limiter = " << limiter
                       << " time_integrator = " << time_integrator << " dorand = " << dorand;
    }
    amrex::Print() << ", cfl = " << cfl << ", dt = " << EstTimeStep() << "\n";
}

void
AmrCoreInertial::DefineLevel (const BoxArray& ba, const DistributionMapping& dm)
{
    if (ip.run_mode == 0) {
        phi_old.define(ba, dm, 1, nghost);
        phi_new.define(ba, dm, 1, nghost);
        phi_old.setVal(0.0);
        phi_new.setVal(0.0);
    } else {
        hist.define(ba, dm, 1, 0);
        pc = std::make_unique<InertialPC>(Geom(0), dm, ba);
    }
}

void
AmrCoreInertial::MakeNewLevelFromScratch (int lev, Real /*time*/, const BoxArray& ba,
                                          const DistributionMapping& dm)
{
    AMREX_ALWAYS_ASSERT(lev == 0);
    DefineLevel(ba, dm);

    if (ip.run_mode == 1) {
        pc->InitParticles(ip);
        return;
    }

    const auto dx  = Geom(0).CellSizeArray();
    const auto plo = Geom(0).ProbLoArray();
    const InertialParams ip_loc = ip;
    for (MFIter mfi(phi_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const& phi = phi_new.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            init_phi(i, j, k, phi, dx, plo, ip_loc);
        });
    }

    // normalize so that phi integrates to 1 over (x, v)
    const Real mass = phi_new.sum(0)*dx[0]*dx[1];
    phi_new.mult(1.0/mass, 0, 1);
    fill_phi_ghosts(phi_new, Geom(0));
}

void
AmrCoreInertial::MakeNewLevelFromCoarse (int, Real, const BoxArray&, const DistributionMapping&)
{
    Abort("AmrCoreInertial is single level");
}

void
AmrCoreInertial::RemakeLevel (int, Real, const BoxArray&, const DistributionMapping&)
{
    Abort("AmrCoreInertial is single level");
}

void
AmrCoreInertial::ClearLevel (int /*lev*/)
{
    phi_old.clear();
    phi_new.clear();
    hist.clear();
    pc.reset();
}

void
AmrCoreInertial::ErrorEst (int, TagBoxArray&, Real, int)
{
}

void
AmrCoreInertial::InitData ()
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
AmrCoreInertial::EstTimeStep () const
{
    if (fixed_dt > 0.0) { return fixed_dt; }
    // Forward Euler with first-order upwind advection and centred v diffusion
    // keeps phi >= 0 when
    //   dt*(Vmax/dx + (gamma Vmax + max|U'|)/dv + 2D/dv^2) <= 1;
    // use the fraction cfl of that (cfl <= 0.5 for MUSCL).
    const auto dxinv = Geom(0).InvCellSizeArray();
    const Real rate = ip.vmax*dxinv[0] + (ip.gamma*ip.vmax + inertial_max_dU(ip))*dxinv[1]
                    + 2.0*ip.diff_coeff*dxinv[1]*dxinv[1];
    Real dt = cfl/rate;
    if (ip.run_mode == 1) {
        // particles: keep the noise step sqrt(2 D dt) below dv/3, so a move is
        // rarely truncated to one cell
        const Real dv = Geom(0).CellSize(1);
        dt = amrex::min(dt, dv*dv/(18.0*ip.diff_coeff));
    }
    return dt;
}

MultiFab const&
AmrCoreInertial::Density ()
{
    if (ip.run_mode == 0) { return phi_new; }
    pc->DepositHistogram(hist);
    return hist;
}

void
AmrCoreInertial::Evolve ()
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

        if (ip.run_mode == 0) {
            AdvancePhi(dt);
        } else {
            const InertialPC::MoveStats st = pc->Advance(dt, ip);
            move_total.trunc_x += st.trunc_x;
            move_total.trunc_v += st.trunc_v;
            move_total.reflect += st.reflect;
        }

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
        if (chk_int > 0 && step_count % chk_int == 0 && ip.run_mode == 0) { WriteCheckpointFile(); }
    }

    if ((plot_int > 0 || plot_dt > 0.0) && step_count > last_plot_step) { WritePlotFile(); }
}

void
AmrCoreInertial::PrintDiagnostics ()
{
    // mass, <x - x0>, <v>, <(x - x0)^2>, <v^2>
    Real m[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    Real phimin = 0.0;
    if (ip.run_mode == 0) {
        const auto dx  = Geom(0).CellSizeArray();
        const auto plo = Geom(0).ProbLoArray();
        const InertialParams ip_loc = ip;
        ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum> rops;
        ReduceData<Real, Real, Real, Real, Real> rdata(rops);
        using ReduceTuple = typename decltype(rdata)::Type;
        for (MFIter mfi(phi_new); mfi.isValid(); ++mfi) {
            auto const& p = phi_new.const_array(mfi);
            rops.eval(mfi.validbox(), rdata, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                const Real d = inertial_dx(plo[0] + (i+0.5)*dx[0], ip_loc.pot_x0, ip_loc);
                const Real v = plo[1] + (j+0.5)*dx[1];
                const Real w = p(i,j,k)*dx[0]*dx[1];
                return {w, w*d, w*v, w*d*d, w*v*v};
            });
        }
        auto r = rdata.value();
        m[0] = amrex::get<0>(r); m[1] = amrex::get<1>(r); m[2] = amrex::get<2>(r);
        m[3] = amrex::get<3>(r); m[4] = amrex::get<4>(r);
        ParallelDescriptor::ReduceRealSum(m, 5);
        for (int c = 1; c < 5; ++c) { m[c] /= m[0]; }
        phimin = phi_new.min(0);
    } else {
        Real pm[4];
        pc->Moments(ip, pm);
        m[0] = static_cast<Real>(pc->TotalNumberOfParticles());
        for (int c = 0; c < 4; ++c) { m[c+1] = pm[c]; }
    }

    amrex::Print() << "STEP " << step_count << " TIME = " << t_new << " DT = " << dt_last
                   << ((ip.run_mode == 0) ? " mass = " : " particles = ")
                   << std::setprecision(15) << m[0] << std::setprecision(6);
    if (ip.run_mode == 0) { amrex::Print() << " min phi = " << phimin; }
    amrex::Print() << " <x-x0> = " << m[1] << " <v> = " << m[2]
                   << " <(x-x0)^2> = " << m[3] << " <v^2> = " << m[4];
    if (ip.run_mode == 1) {
        amrex::Print() << " truncated x/v = " << move_total.trunc_x << "/" << move_total.trunc_v
                       << " reflected = " << move_total.reflect;
    }
    amrex::Print() << "\n";
}

void
AmrCoreInertial::WritePlotFile ()
{
    const std::string& plotfilename = amrex::Concatenate(plot_file, step_count, 6);
    amrex::Print() << "Writing plotfile " << plotfilename << "\n";
    MultiFab const& phi = Density();
    MultiFab mf(phi.boxArray(), phi.DistributionMap(), 4, 0);
    MultiFab::Copy(mf, phi, 0, 0, 1, 0);

    // moments over v, rho(x) = sum_j phi dv, mom(x) = sum_j v phi dv and
    // ke(x) = sum_j v^2/2 phi dv, stored in every v row. A box need not span
    // all of v, so copy phi to a layout whose boxes are whole v columns.
    const Box& domain = Geom(0).Domain();
    IntVect col_size = maxGridSize(0);
    col_size[1] = domain.length(1);
    BoxArray cba(domain);
    cba.maxSize(col_size);
    DistributionMapping cdm(cba);
    MultiFab col(cba, cdm, 3, 0);
    col.ParallelCopy(phi, 0, 0, 1);

    const Real dv  = Geom(0).CellSize(1);
    const Real vlo = Geom(0).ProbLo(1);
    const int jlo = domain.smallEnd(1);
    const int jhi = domain.bigEnd(1);
    for (MFIter mfi(col); mfi.isValid(); ++mfi) {
        auto const& c = col.array(mfi);
        Box xb = mfi.validbox();
        xb.setRange(1, jlo, 1);
        amrex::ParallelFor(xb, [=] AMREX_GPU_DEVICE (int i, int, int k) noexcept
        {
            Real rho = 0.0, mom = 0.0, ke = 0.0;
            for (int j = jlo; j <= jhi; ++j) {
                const Real v = vlo + (j+0.5)*dv;
                rho += c(i,j,k);
                mom += c(i,j,k)*v;
                ke  += c(i,j,k)*0.5*v*v;
            }
            for (int j = jlo; j <= jhi; ++j) {
                c(i,j,k,0) = rho*dv;
                c(i,j,k,1) = mom*dv;
                c(i,j,k,2) = ke*dv;
            }
        });
    }
    mf.ParallelCopy(col, 0, 1, 3);

    WriteSingleLevelPlotfile(plotfilename, mf, {"phi", "rho", "mom", "ke"}, Geom(0), t_new, step_count);
    if (ip.run_mode == 1 && write_particles) {
        pc->WritePlotFile(plotfilename, "particles");
    }
}

void
AmrCoreInertial::WriteCheckpointFile () const
{
    const std::string& checkpointname = amrex::Concatenate(chk_file, step_count);
    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", 1, true);

    if (ParallelDescriptor::IOProcessor()) {
        std::string HeaderFileName(checkpointname + "/Header");
        std::ofstream HeaderFile(HeaderFileName.c_str(), std::ofstream::out | std::ofstream::trunc);
        if (!HeaderFile.good()) { amrex::FileOpenFailed(HeaderFileName); }
        HeaderFile.precision(17);
        HeaderFile << "Checkpoint file for AmrCoreInertial\n";
        HeaderFile << step_count << "\n" << dt_last << "\n" << t_new << "\n";
        boxArray(0).writeOn(HeaderFile);
        HeaderFile << '\n';
    }

    VisMF::Write(phi_new, amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "phi"));
}

void
AmrCoreInertial::ReadCheckpointFile ()
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
    fill_phi_ghosts(phi_new, Geom(0));
}
