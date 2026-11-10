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
        pp.query("regrid_int", regrid_int);
        if (ip.hybrid) {
            if (max_level != 1) { Abort("the hybrid needs amr.max_level = 1"); }
            if (refRatio(0) != IntVect(1)) { Abort("the hybrid needs amr.ref_ratio = 1"); }
        } else if (max_level != 0) {
            Abort("amr.max_level must be 0 unless inertial.hybrid = 1");
        }
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
    ip.rng_seed = static_cast<Real>(seed);
    if ((ip.run_mode == 1 || ip.hybrid) && !restart_chkfile.empty()) {
        Abort("restart is only implemented for the pure SPDE (run_mode = 0, hybrid = 0)");
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
    if (lev == 1) {
        AMREX_ALWAYS_ASSERT(ip.hybrid);
        SetupFine(ba, dm);
        return;
    }
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
AmrCoreInertial::MakeNewLevelFromCoarse (int lev, Real, const BoxArray& ba, const DistributionMapping& dm)
{
    AMREX_ALWAYS_ASSERT(lev == 1 && ip.hybrid);
    SetupFine(ba, dm);
}

void
AmrCoreInertial::RemakeLevel (int lev, Real, const BoxArray& ba, const DistributionMapping& dm)
{
    AMREX_ALWAYS_ASSERT(lev == 1 && ip.hybrid);
    SetupFine(ba, dm);
}

void
AmrCoreInertial::ClearLevel (int lev)
{
    if (lev == 1) {
        // the particle region disappears: its mass is already in the level-0
        // phi (the last AverageDownParticles), so the particles just go
        if (pc) { pc->clearParticles(); }
        fine_ba = BoxArray();
        BuildMask(mask0, fine_ba);
        phi1.clear();
        return;
    }
    phi_old.clear();
    phi_new.clear();
    hist.clear();
    pc.reset();
}

void
AmrCoreInertial::ErrorEst (int lev, TagBoxArray& tags, Real, int)
{
    // hybrid: tag the cells holding few particles, N phi dx dv < n_refine
    if (!ip.hybrid || lev != 0) { return; }
    const auto dx = Geom(0).CellSizeArray();
    const Real thresh = ip.n_refine/(ip.num_part*dx[0]*dx[1]);
    const bool use_box = ip.tag_box[2] >= ip.tag_box[0] && ip.tag_box[3] >= ip.tag_box[1];
    const int b0 = ip.tag_box[0], b1 = ip.tag_box[1], b2 = ip.tag_box[2], b3 = ip.tag_box[3];
    for (MFIter mfi(phi_new); mfi.isValid(); ++mfi) {
        auto const& p = phi_new.const_array(mfi);
        auto const& t = tags.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const bool tag = use_box ? (i >= b0 && i <= b2 && j >= b1 && j <= b3) : (p(i,j,k) < thresh);
            if (tag) { t(i,j,k) = TagBox::SET; }
        });
    }
}

void
AmrCoreInertial::BuildMask (iMultiFab& mask, const BoxArray& fba) const
{
    mask.define(boxArray(0), DistributionMap(0), 1, 1);
    mask.setVal(0);
    if (!fba.empty()) {
        for (MFIter mfi(mask); mfi.isValid(); ++mfi) {
            for (auto const& is : fba.intersections(mfi.validbox())) {
                mask[mfi].setVal<RunOn::Host>(1, is.second, 0, 1);
            }
        }
    }
    mask.FillBoundary(Geom(0).periodicity());
}

void
AmrCoreInertial::SetupFine (const BoxArray& ba, const DistributionMapping& dm)
{
    const Geometry& g = Geom(0);
    const auto dx = g.CellSizeArray();
    const Real dv = dx[0]*dx[1];
    const int ng = ip.hybrid_nghost;
    const Real mass_before = phi_new.sum(0)*dv;

    iMultiFab mask_new;
    BuildMask(mask_new, ba);

    // particles in cells leaving the region: their mass is already in phi
    if (pc && pc->TotalNumberOfParticles() > 0) {
        iMultiFab keep(pba, pdm, 1, 0);
        keep.setVal(0);
        keep.ParallelCopy(mask_new, 0, 0, 1, IntVect(0), IntVect(0), g.periodicity());
        pc->RemoveWhere(keep);
    }

    // cells new to the region
    iMultiFab newc(boxArray(0), DistributionMap(0), 1, 0);
    const bool had = mask0.ok();
    for (MFIter mfi(newc); mfi.isValid(); ++mfi) {
        auto const& nc = newc.array(mfi);
        auto const& mn = mask_new.const_array(mfi);
        auto const& mo = had ? mask0.const_array(mfi) : mn;
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            nc(i,j,k) = (mn(i,j,k) != 0 && (!had || mo(i,j,k) == 0)) ? 1 : 0;
        });
    }

    fine_ba = ba;
    fine_dm = dm;
    mask0 = std::move(mask_new);

    // particle layout: the region boxes, then the ghost ring (the grown boxes,
    // periodic in x and clipped in v, without overlaps, minus the region)
    const Box& dom = g.Domain();
    const int nx = dom.length(0);
    BoxList bl_grown;
    for (int b = 0; b < ba.size(); ++b) {
        Box gb = amrex::grow(ba[b], ng);
        gb.setSmall(1, std::max(gb.smallEnd(1), dom.smallEnd(1)));
        gb.setBig(1, std::min(gb.bigEnd(1), dom.bigEnd(1)));
        for (int sh : {-nx, 0, nx}) {
            Box sb = gb;
            sb.shift(0, sh);
            sb &= dom;
            if (sb.ok()) { bl_grown.push_back(sb); }
        }
    }
    BoxArray gba(std::move(bl_grown));
    gba.removeOverlap();
    BoxList bl_ghost;
    for (int b = 0; b < gba.size(); ++b) { bl_ghost.join(ba.complementIn(gba[b])); }
    BoxArray ghost_ba(std::move(bl_ghost));
    ghost_ba.maxSize(maxGridSize(0));
    BoxList bl_p(ba);
    bl_p.join(ghost_ba.boxList());
    pba = BoxArray(std::move(bl_p));
    pdm = DistributionMapping(pba);
    n_fine_boxes = ba.size();

    if (!pc) {
        pc = std::make_unique<InertialPC>(g, pdm, pba);
    } else {
        pc->SetParticleBoxArray(0, pba);
        pc->SetParticleDistributionMap(0, pdm);
        pc->Redistribute();
    }

    pmask.define(pba, pdm, 1, ng);
    pmask.setVal(0);
    pmask.ParallelCopy(mask0, 0, 0, 1, IntVect(0), IntVect(ng), g.periodicity());
    gsel.define(pba, pdm, 1, 0);
    for (MFIter mfi(gsel); mfi.isValid(); ++mfi) {
        gsel[mfi].setVal<RunOn::Host>((mfi.index() >= n_fine_boxes) ? 1 : 0);
    }

    // sample particles in the new cells from phi
    iMultiFab nsel(pba, pdm, 1, 0);
    nsel.setVal(0);
    nsel.ParallelCopy(newc, 0, 0, 1, IntVect(0), IntVect(0), g.periodicity());
    MultiFab phip(pba, pdm, 1, 1);
    phip.ParallelCopy(phi_new, 0, 0, 1, IntVect(1), IntVect(1), g.periodicity());
    const Long added = pc->SampleCells(phip, nsel, ip);
    pc->Redistribute();
    regrid_sampled += added;

    phi1.define(ba, dm, 1, 0);
    AverageDownParticles();

    const Real dmass = phi_new.sum(0)*dv - mass_before;
    regrid_dmass += dmass;
    Long ncells = 0;
    for (int b = 0; b < ba.size(); ++b) { ncells += ba[b].numPts(); }
    amrex::Print() << "Particle region: " << ba.size() << " boxes, " << ncells << " of "
                   << dom.numPts() << " cells, ghost ring " << ng << ", " << pc->TotalNumberOfParticles()
                   << " particles (" << added << " sampled in new cells), mass change " << dmass << "\n";
}

void
AmrCoreInertial::AverageDownParticles ()
{
    MultiFab h(pba, pdm, 1, 0);
    pc->DepositHistogram(h, ip.num_part);
    phi1.ParallelCopy(h, 0, 0, 1);
    phi_new.ParallelCopy(phi1, 0, 0, 1);
    fill_phi_ghosts(phi_new, Geom(0));
}

void
AmrCoreInertial::HybridStep (Real dt)
{
    const Geometry& g = Geom(0);
    const auto dx = g.CellSizeArray();
    const bool fine = HasFine();

    // ghost ring: particles sampled from the start-of-step SPDE
    if (fine) {
        MultiFab phip(pba, pdm, 1, 1);
        phip.ParallelCopy(phi_new, 0, 0, 1, IntVect(1), IntVect(1), g.periodicity());
        ghost_sampled += pc->SampleCells(phip, gsel, ip);
    }

    // SPDE on all of level 0, keeping the applied face fluxes
    Array<MultiFab, AMREX_SPACEDIM> fsum;
    if (fine) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            fsum[d].define(amrex::convert(boxArray(0), IntVect::TheDimensionVector(d)), DistributionMap(0), 1, 0);
            fsum[d].setVal(0.0);
        }
    }
    AdvancePhi(dt, fine ? &fsum : nullptr);
    if (!fine) { return; }

    // the particles carry the mass across the region boundary
    undo_cf_fluxes(phi_new, fsum, mask0, g, dt);

    const InertialPC::MoveStats st = pc->Advance(dt, ip, step_count, false);
    move_total.trunc_x += st.trunc_x;
    move_total.trunc_v += st.trunc_v;
    move_total.reflect += st.reflect;

    MultiFab xfer(pba, pdm, 1, 1);
    xfer.setVal(0.0);
    const auto c = pc->DepositCrossings(pmask, xfer);
    for (int q = 0; q < 4; ++q) { crossings[q] += c[q]; }
    xfer.SumBoundary(g.periodicity());
    xfer.mult(1.0/(ip.num_part*dx[0]*dx[1]), 0, 1);
    phi_new.ParallelAdd(xfer, 0, 0, 1, IntVect(0), IntVect(0), g.periodicity());
    pc->Redistribute();

    AverageDownParticles();
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
    if (ip.run_mode == 1 || ip.hybrid) {
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

        if (ip.hybrid) {
            if (regrid_int > 0 && step > 0 && step % regrid_int == 0) { regrid(0, cur_time); }
            HybridStep(dt);
        } else if (ip.run_mode == 0) {
            AdvancePhi(dt);
        } else {
            const InertialPC::MoveStats st = pc->Advance(dt, ip, step);
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
    if (ip.run_mode == 1 || ip.hybrid) {
        amrex::Print() << " truncated x/v = " << move_total.trunc_x << "/" << move_total.trunc_v
                       << " reflected = " << move_total.reflect;
    }
    if (ip.hybrid) {
        Long ncells = 0;
        for (int b = 0; b < fine_ba.size(); ++b) { ncells += fine_ba[b].numPts(); }
        amrex::Print() << " | particles " << (pc ? pc->TotalNumberOfParticles() : 0)
                       << " in " << ncells << " cells, ghost sampled " << ghost_sampled
                       << ", regrid sampled " << regrid_sampled << " (mass change " << regrid_dmass << ")"
                       << ", entries x/v " << crossings[0] << "/" << crossings[1]
                       << ", exits x/v " << crossings[2] << "/" << crossings[3];
    }
    amrex::Print() << "\n";
}

void
AmrCoreInertial::WritePlotFile ()
{
    const std::string& plotfilename = amrex::Concatenate(plot_file, step_count, 6);
    amrex::Print() << "Writing plotfile " << plotfilename << "\n";
    MultiFab const& phi = Density();
    const int ncomp = ip.hybrid ? 5 : 4;
    MultiFab mf(phi.boxArray(), phi.DistributionMap(), ncomp, 0);
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

    if (ip.hybrid) {
        // 1 in the particle region
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            auto const& o = mf.array(mfi);
            auto const& m = mask0.const_array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                o(i,j,k,4) = static_cast<Real>(m(i,j,k));
            });
        }
        WriteSingleLevelPlotfile(plotfilename, mf, {"phi", "rho", "mom", "ke", "particles"}, Geom(0), t_new, step_count);
    } else {
        WriteSingleLevelPlotfile(plotfilename, mf, {"phi", "rho", "mom", "ke"}, Geom(0), t_new, step_count);
    }
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
