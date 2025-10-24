
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <AmrCoreAdv.H>
#include <Kernels.H>
#include <mykernel.H>
#include <myfunc.H>
#include "chrono"

using namespace amrex;
using namespace std::chrono;

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRe boundary condition object
AmrCoreAdv::AmrCoreAdv ()
{

    // periodic boundaries
    //int bc_lo[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    //int bc_hi[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    amrex::Vector<int> bc_lo(AMREX_SPACEDIM,0);
    amrex::Vector<int> bc_hi(AMREX_SPACEDIM,0);

/*
    // walls (Neumann)
    int bc_lo[] = {amrex::BCType::foextrap, amrex::BCType::foextrap, amrex::BCType::foextrap};
    int bc_hi[] = {amrex::BCType::foextrap, amrex::BCType::foextrap, amrex::BCType::foextrap};
*/
    // walls Dirichlet
    //int bc_lo[] = {amrex::BCType::ext_dir, amrex::BCType::ext_dir, amrex::BCType::ext_dir};
    //int bc_hi[] = {amrex::BCType::ext_dir, amrex::BCType::ext_dir, amrex::BCType::ext_dir};

    ReadParameters(bc_lo,bc_hi);

    /////////////////////////////////////////
    //Initialise rngs
    /////////////////////////////////////////
    int restart = -1;

    if (restart < 0) {

        if (seed > 0) {
            // initializes the seed for C++ random number calls
            InitRandom(seed+ParallelDescriptor::MyProc(),
                       ParallelDescriptor::NProcs(),
                       seed+ParallelDescriptor::MyProc());
        } else if (seed == 0) {
            // initializes the seed for C++ random number calls based on the clock
            auto now = time_point_cast<nanoseconds>(system_clock::now());
            int randSeed = now.time_since_epoch().count();
            // broadcast the same root seed to all processors
            ParallelDescriptor::Bcast(&randSeed,1,ParallelDescriptor::IOProcessorNumber());
            InitRandom(randSeed+ParallelDescriptor::MyProc(),
                       ParallelDescriptor::NProcs(),
                       randSeed+ParallelDescriptor::MyProc());
        } else {
            Abort("Must supply non-negative seed");
        }
    }

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    if (do_subcycle) {
        for (int lev = 1; lev <= max_level; ++lev) {
            nsubsteps[lev] = MaxRefRatio(lev-1);
        }
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);

    bcs.resize(1);     // Setup 1-component
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        // lo-side BCs
        if (bc_lo[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_lo[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_lo[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs[0].setLo(idim, bc_lo[idim]);
        }
        else {
            amrex::Abort("Invalid bc_lo");
        }

        // hi-side BCSs
        if (bc_hi[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_hi[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_hi[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs[0].setHi(idim, bc_hi[idim]);
        }
        else {
            amrex::Abort("Invalid bc_hi");
        }
    }

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] is never actually used in the reflux operation
    flux_reg.resize(nlevs_max+1);

    // fillpatcher[lev] is for filling data on level lev using the data on
    // lev-1 and lev.
    fillpatcher.resize(nlevs_max+1);
}

AmrCoreAdv::~AmrCoreAdv () = default;

// advance solution to final time
void
AmrCoreAdv::Evolve ()
{
    Real cur_time = t_new[0];
    int last_plot_file_step = 0;

    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        amrex::Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;

        ComputeDt();

        int lev = 0;
        int iteration = 1;
        timeStepNoSubcycling(cur_time, iteration);

        cur_time += dt[0];

        // sum phi to check conservation
        Real sum_phi_old = phi_old[0].sum();
        Real sum_phi_new = phi_new[0].sum();

        amrex::Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                       << " DT = " << dt[0] << " Sum_old Sum_new Diff (Phi) = "    << std::setw(20) << std::setprecision(12)
                       << std::scientific <<  sum_phi_old << " " << std::setw(2l) << std::setprecision(12)
                       << std::scientific <<  sum_phi_new << " " << std::setw(2l) << std::setprecision(12)
                       << std::scientific << (sum_phi_new - sum_phi_old) << std::endl;

        // sync up time
        for (lev = 0; lev <= finest_level; ++lev) {
            t_new[lev] = cur_time;
        }

        if (plot_int > 0 && (step+1) % plot_int == 0) {
            last_plot_file_step = step+1;
            WritePlotFile();
        }

        if (chk_int > 0 && (step+1) % chk_int == 0) {
            WriteCheckpointFile();
        }

#ifdef AMREX_MEM_PROFILING
        {
            std::ostringstream ss;
            ss << "[STEP " << step+1 << "]";
            MemProfiler::report(ss.str());
        }
#endif

        if (cur_time >= stop_time - 1.e-6*dt[0]) { break; }
    }

    if (plot_int > 0 && istep[0] > last_plot_file_step) {
        WritePlotFile();
    }
}

// initializes multilevel data
void
AmrCoreAdv::InitData ()
{
    if (restart_chkfile.empty()) {
        // start simulation from the beginning
        const Real time = 0.0;

        // Note that InitFromScratch allocates the space for phi at each level,
        //      but only initializes phi at level 0, not at level > 0.
        //      So we can't do average down until we create the particles,
        //      then use the particles to define the level 1 phi
        InitFromScratch(time);

#ifdef AMREX_PARTICLES
        if (max_level > 0) {
            particleData.init_particles((amrex::ParGDBBase*)GetParGDB(), grown_fba, phi_new[0], phi_new[1]);
        }
        else {
            particleData.init_particles((amrex::ParGDBBase*)GetParGDB(), phi_new[0]);
        }
#endif
        AverageDown();
        phi_new[0].FillBoundary();

        MultiFab::Copy(phi_old[0], phi_new[0],0,0,1,0);
        phi_old[0].FillBoundary();

        if (chk_int > 0) {
            WriteCheckpointFile();
        }
    }
    else {
        // restart from a checkpoint
        ReadCheckpointFile();
    }
    if (plot_int > 0) {
        WritePlotFile();
    }
}

void AmrCoreAdv::MakeFBA(const BoxArray& ba)
{
    int lev = 1;
    Box domain(Geom(lev).Domain());
    BoxList valid_bl(ba);
    BoxList com_bl = GetBndryCells(ba,1);
#if (AMREX_SPACEDIM == 2)
    Vector<IntVect> pshifts(9);
#else
    Vector<IntVect> pshifts(27);
#endif

    BoxList com_bl_fixed;

    //
    // Loop over boxes created by GetBndryCells call -- note that if periodic
    // some of these boxes may intersect the valid_bl so we remove those intersections
    // by intersecting with the copmlement of the valid ba
    //
    for (auto& b : com_bl) {
        Box bx(b);

        //
        // First intersect the existing box with the domain and keep that
        // Note that GetBndryCells would not include any cells inside the domain
        // that are part of the original ba
        //
        Box b1 = bx & domain;
        if (!b1.isEmpty()) {
            com_bl_fixed.push_back(b1);
        }

        //
        // Next add the pieces that were outside the domain in a periodic direction
        // Note that GetBndryCells DOES include cells outside the domain
        // that are part of the original ba if shifted periodically
        //
        geom[lev].periodicShift(domain, bx, pshifts);
        for (int n = 0; n < pshifts.size(); n++) {
            Box bx_shift(b);
            bx_shift.shift(pshifts[n]);
            Box b2 = bx_shift & domain;
            if (!b2.isEmpty()) {
                // Now we have to make sure we don't include any intersection of this b2
                // with the valid boxArray
                BoxList bl_comp = complementIn(b2,valid_bl);
                for (auto& b_comp : bl_comp) {
                    Box bx_comp(b_comp);
                    if (!bx_comp.isEmpty()) {
                        com_bl_fixed.push_back(bx_comp);
                    }
                }
            }
        }
    }

    //
    // Remove any duplicated regions in the boundary cells
    //
    com_bl_fixed.simplify();

    //
    // Add the valid boxes
    //
    com_bl_fixed.catenate(valid_bl);
    grown_fba.define(com_bl_fixed);
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data (overrides the pure virtual function in AmrCore)
// regrid  --> RemakeLevel            (if level already existed)
// regrid  --> MakeNewLevelFromCoarse (if adding new level)
void
AmrCoreAdv::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
                                    const DistributionMapping& dm)
{
    const int ncomp = phi_new[lev-1].nComp();
    const int ng = phi_new[lev-1].nGrow();

    amrex::Print() << " CREATE LEVEL " << lev << " " << ba << std::endl;

    phi_new[lev].define(ba, dm, ncomp, ng);
    phi_old[lev].define(ba, dm, ncomp, ng);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg[lev] = std::make_unique<FluxRegister>(ba, dm, refRatio(lev-1), lev, ncomp);
    }

    FillCoarsePatch(lev, time, phi_new[lev], 0, ncomp);
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data (overrides the pure virtual function in AmrCore)
// regrid  --> RemakeLevel            (if level already existed)
// regrid  --> MakeNewLevelFromCoarse (if adding new level)
void
AmrCoreAdv::RemakeLevel (int lev, Real time, const BoxArray& ba,
                         const DistributionMapping& dm)
{
    const int ncomp = phi_new[lev].nComp();
    const int ng = phi_new[lev].nGrow();

    BoxArray old_fine_ba = phi_old[1].boxArray();
    amrex::Print() << " REGRIDDING: NEW GRIDS AT LEVEL " << lev << " " << ba << std::endl;

    if (lev == 1) {
        MakeFBA(ba);
    }

    MultiFab new_state(ba, dm, ncomp, ng);
    MultiFab old_state(ba, dm, ncomp, ng);

    // Must use fillpatch_function
    FillPatch(lev, time, new_state, 0, ncomp, FillPatchType::fillpatch_function);

    std::swap(new_state, phi_new[lev]);
    std::swap(old_state, phi_old[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg[lev] = std::make_unique<FluxRegister>(ba, dm, refRatio(lev-1), lev, ncomp);
    }

#ifdef AMREX_PARTICLES
        if (lev == 1) {
            particleData.regrid_particles(grown_fba, ba, old_fine_ba, phi_new[1]);
        }
#endif
}

// Delete level data
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::ClearLevel (int lev)
{
    phi_new[lev].clear();
    phi_old[lev].clear();
    flux_reg[lev].reset(nullptr);
    fillpatcher[lev].reset(nullptr);
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
// main.cpp --> AmrCoreAdv::InitData --> InitFromScratch --> MakeNewGrids --> MakeNewLevelFromScratch
//                                              restart  --> MakeNewGrids --> MakeNewLevelFromScratch
void AmrCoreAdv::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
                                          const DistributionMapping& dm)
{
    const int ng = 1;

    amrex::Print() << " GRIDS AT LEVEL " << lev << " " << ba << std::endl;

    // ncomp = number of components for each array
    int ncomp;
    if (alg_type == 0) {
       ncomp = 1;
    } else {
       ncomp = 2;
    }

    if (lev == 1) {
        MakeFBA(ba);
    }

    phi_new[lev].define(ba, dm, ncomp, ng);
    phi_old[lev].define(ba, dm, ncomp, ng);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg[lev] = std::make_unique<FluxRegister>(ba, dm, refRatio(lev-1), lev, ncomp);
    }

    const auto problo = Geom(lev).ProbLoArray();
    const auto dx     = Geom(lev).CellSizeArray();

    int Ncomp = phi_new[lev].nComp();

    // External Potential related
    int a_ext_pot = m_ext_pot;
    Real a_alpha  = m_ext_pot_alpha;
    Real a_beta   = m_ext_pot_beta;
    Real a_gamma  = m_ext_pot_gamma;

    if (lev == 0) {
        for (MFIter mfi(phi_new[lev]); mfi.isValid(); ++mfi)
        {
            const Box& vbx = mfi.validbox();
            auto const& phi_arr = phi_new[lev].array(mfi);
            auto npts_scale_local = npts_scale;
            amrex::ParallelFor(vbx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                init_phi(i,j,k,phi_arr,dx,problo,npts_scale_local,Ncomp,
                         a_ext_pot, a_alpha, a_beta, a_gamma);
            });
        }
    } else {
        phi_new[lev].ParallelCopy(phi_new[lev-1], 0, 0, phi_new[lev].nComp());
    }
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::ErrorEst (int lev, TagBoxArray& tags, Real /*time*/, int /*ngrow*/)
{
    static bool first = true;
    static Vector<Real> phierr;

    // only do this during the first call to ErrorEst
    if (first)
    {
        first = false;
        // read in an array of "phierr", which is the tagging threshold
        // in this example, we tag values of "phi" which are greater than phierr
        // for that particular level
        // in subroutine state_error, you could use more elaborate tagging, such
        // as more advanced logical expressions, or gradients, etc.
        ParmParse pp("adv");
        int n = pp.countval("phierr");
        if (n > 0) {
            pp.getarr("phierr", phierr, 0, n);
        }
    }

    if (lev >= phierr.size()) { return; }

//    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const MultiFab& state = phi_new[lev];

    const Real* dx  =  geom[lev].CellSize();
#if (AMREX_SPACEDIM == 2)
    const Real cell_vol = dx[0]*dx[1];
#else
    const Real cell_vol = dx[0]*dx[1]*dx[2];
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif
    {

        for (MFIter mfi(state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx  = mfi.tilebox();
            const auto statefab = state.array(mfi);
            const auto tagfab  = tags.array(mfi);
            Real phierror = phierr[lev]/cell_vol;

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                state_error(i, j, k, tagfab, statefab, phierror, tagval);
            });
        }
    }
}

// read in some parameters from inputs file
void
AmrCoreAdv::ReadParameters ( amrex::Vector<int>& bc_lo, amrex::Vector<int>& bc_hi)
{
    {
        ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);

        npts_scale = 1.0;
        pp.queryAdd("npts_scale", npts_scale);

        alg_type = 0;
        pp.queryAdd("alg_type", alg_type);


        // read in BC; see Src/Base/AMReX_BC_TYPES.H for supported types
        pp.queryarr("bc_lo", bc_lo);
        pp.queryarr("bc_hi", bc_hi);

        seed = 0;
        pp.queryAdd("seed", seed);

        // read in if a direction is ensemble direction
        pp.queryarr("is_ensemble_dir", m_ensemble_dir, 0, AMREX_SPACEDIM);
        // Some asserts for m_ensemble_dir
        for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
            if (m_ensemble_dir[idir]) {
                amrex::Print()<<"max_level: "<<max_level<<"\n";
                AMREX_ALWAYS_ASSERT(max_level == 0);
                AMREX_ALWAYS_ASSERT(Geom(0).CellSize(idir) == Real(1.0));
            }
        }
    }

    {
        ParmParse pp("amr"); // Traditionally, these have prefix, amr.

        pp.query("regrid_int", regrid_int);
        pp.query("plot_file", plot_file);
        pp.query("plot_int", plot_int);
        pp.query("chk_file", chk_file);
        pp.query("chk_int", chk_int);
        pp.query("restart",restart_chkfile);
    }

    {
        ParmParse pp("adv");

        pp.query("cfl", cfl);
        pp.query("do_reflux", do_reflux);
        pp.query("do_subcycle", do_subcycle);
    }

    {
        ParmParse pp("ext_pot");

        pp.query("exists", m_ext_pot);
        if (m_ext_pot) {
            if (AMREX_SPACEDIM != 2) {
                amrex::Abort("External Potential is coded for 2D.\n");
            }
            pp.get("alpha", m_ext_pot_alpha);
            pp.get("beta", m_ext_pot_beta);
            pp.get("gamma", m_ext_pot_gamma);
        }
    }

#ifdef AMREX_PARTICLES
        int a_ensemble_dir_exists = 0;
        for (int edir : m_ensemble_dir) {
            a_ensemble_dir_exists += edir;
        }
        if (a_ensemble_dir_exists) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(alg_type != 0,
                "Ensemble mode with particles requires alg_type != 0, i.e., ncomp = 2");
        }
        particleData.init_particle_params(max_level, a_ensemble_dir_exists);
#endif
}

// set covered coarse cells to be the average of overlying fine cells
void
AmrCoreAdv::AverageDown ()
{
    for (int lev = finest_level; lev >= 1; --lev) {
        phi_new[lev-1].ParallelCopy(phi_new[lev], 0, 0, phi_new[lev].nComp());
    }
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void
AmrCoreAdv::AverageDownTo (int crse_lev)
{
    phi_new[crse_lev].ParallelCopy(phi_new[crse_lev+1], 0, 0, phi_new[crse_lev].nComp());
}

// compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
AmrCoreAdv::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp,
                       FillPatchType fptype)
{
    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetData(0, time, smf, stime);

        if(Gpu::inLaunchRegion())
        {
            GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
            PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > physbc(geom[lev],bcs,gpu_bndry_func);
            amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                        geom[lev], physbc, 0);
        }
        else
        {
            CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
            PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev],bcs,bndry_func);
            amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                        geom[lev], physbc, 0);
        }
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetData(lev-1, time, cmf, ctime);
        GetData(lev  , time, fmf, ftime);

        Interpolater* mapper = &cell_cons_interp;

        if (fptype == FillPatchType::fillpatch_class) {
            if (fillpatcher[lev] == nullptr) {
                fillpatcher[lev] = std::make_unique<FillPatcher<MultiFab>>
                    (boxArray(lev  ), DistributionMap(lev  ), Geom(lev  ),
                     boxArray(lev-1), DistributionMap(lev-1), Geom(lev-1),
                     mf.nGrowVect(), mf.nComp(), mapper);
            }
        }

        if(Gpu::inLaunchRegion())
        {
            GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
            PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > cphysbc(geom[lev-1],bcs,gpu_bndry_func);
            PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > fphysbc(geom[lev],bcs,gpu_bndry_func);

            if (fptype == FillPatchType::fillpatch_class) {
                fillpatcher[lev]->fill(mf, mf.nGrowVect(), time,
                                       cmf, ctime, fmf, ftime, 0, icomp, ncomp,
                                       cphysbc, 0, fphysbc, 0, bcs, 0);
            } else {
                amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                          0, icomp, ncomp, geom[lev-1], geom[lev],
                                          cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                          mapper, bcs, 0);
            }
        }
        else
        {
            CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
            PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bndry_func);
            PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bndry_func);

            if (fptype == FillPatchType::fillpatch_class) {
                fillpatcher[lev]->fill(mf, mf.nGrowVect(), time,
                                       cmf, ctime, fmf, ftime, 0, icomp, ncomp,
                                       cphysbc, 0, fphysbc, 0, bcs, 0);
            } else {
                amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                          0, icomp, ncomp, geom[lev-1], geom[lev],
                                          cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                          mapper, bcs, 0);
            }
        }
    }
}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void
AmrCoreAdv::FillCoarsePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    BL_ASSERT(lev > 0);

    Vector<MultiFab*> cmf;
    Vector<Real> ctime;
    GetData(lev-1, time, cmf, ctime);
    Interpolater* mapper = &cell_cons_interp;

    if (cmf.size() != 1) {
        amrex::Abort("FillCoarsePatch: how did this happen?");
    }

    if(Gpu::inLaunchRegion())
    {
        GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > cphysbc(geom[lev-1],bcs,gpu_bndry_func);
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > fphysbc(geom[lev],bcs,gpu_bndry_func);

        amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
                                     cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                     mapper, bcs, 0);
    }
    else
    {
        CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bndry_func);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bndry_func);

        amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
                                     cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                     mapper, bcs, 0);
    }
}

void
AmrCoreAdv::GetData (int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    if (amrex::almostEqual(time, t_new[lev], 5))
    {
        data.push_back(&phi_new[lev]);
        datatime.push_back(t_new[lev]);
    }
    else if (amrex::almostEqual(time, t_old[lev], 5))
    {
        data.push_back(&phi_old[lev]);
        datatime.push_back(t_old[lev]);
    }
    else
    {
        data.push_back(&phi_old[lev]);
        data.push_back(&phi_new[lev]);
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}

// Advance all the levels with the same dt
void
AmrCoreAdv::timeStepNoSubcycling (Real time, int iteration)
{
    if (max_level > 0 && regrid_int > 0)  // We may need to regrid
    {
        if (istep[0] % regrid_int == 0)
        {
            amrex::Print() << "Regridding at step " << istep[0] << std::endl;
            regrid(0, time);

            AverageDown();

            Real sum_phi_reg_new = phi_new[0].sum();
            Real sum_phi_reg_old = phi_old[0].sum();
            amrex::Print() << " Sum(Phi) new / old / diff / %diff  after regrid = " << std::setw(24) <<  std::setprecision(16) << std::scientific <<
                   sum_phi_reg_new << " " << sum_phi_reg_old << " " << (sum_phi_reg_new-sum_phi_reg_old) << " " <<
                   (sum_phi_reg_new-sum_phi_reg_old)/sum_phi_reg_old << std::endl;
        }
    }

    // Advance phi at level 0 only
    int  lev = 0;
    int nsub = 1;
    std::swap(phi_old[lev], phi_new[lev]);

    if (Verbose()) {
       amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
       amrex::Print() << "ADVANCE with time = " << t_new[lev] << " dt = " << dt[0] << std::endl;
    }

    AdvancePhiAtLevel(lev, time, dt[lev], iteration, nsub);

    if (finest_level > 0) {
        flux_reg[lev+1]->Reflux(phi_new[lev], 1.0, 0, 0, phi_new[lev].nComp(), geom[lev]);
    }

    if (Verbose()) {
        amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
        amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
    }

    int lev_for_particles = finest_level;
#ifdef AMREX_PARTICLES
    // We set this to finest_level rather than 1 so that we can run with max_level = 0;
    //    in that case the advance_particle routine will return without doing anything
    if (finest_level > 0) {
        std::swap(phi_old[lev_for_particles], phi_new[lev_for_particles]);
    }
    const Real* dx  =  geom[lev_for_particles].CellSize();
#if (AMREX_SPACEDIM == 2)
    const Real cell_vol = dx[0]*dx[1];
#else
    const Real cell_vol = dx[0]*dx[1]*dx[2];
#endif
    particleData.advance_particles(lev_for_particles, dt[lev_for_particles], cell_vol,
                                   phi_old[0], phi_new[0], phi_new[lev_for_particles]);
    particleData.Redistribute();
#endif

    // Make sure the coarser levels are consistent with the finer levels
    AverageDown ();

    for (auto& fp : fillpatcher) {
        fp.reset(); // Because the data have changed.
    }

    for (int ilev = 0; ilev <= finest_level; ilev++) {
        ++istep[ilev];
    }

    if (Verbose())
    {
        amrex::Print() << "[Level " << lev_for_particles << " step " << istep[lev_for_particles] << "] ";
        amrex::Print() << "Advanced " << CountCells(lev_for_particles) << " cells" << std::endl;
    }
}

// a wrapper for EstTimeStep
void
AmrCoreAdv::ComputeDt ()
{
    Vector<Real> dt_tmp(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        dt_tmp[lev] = EstTimeStep(lev, t_new[lev]);
    }
    ParallelDescriptor::ReduceRealMin(dt_tmp.data(), int(dt_tmp.size()));

    constexpr Real change_max = 1.1;
    Real dt_0 = dt_tmp[0];
    int n_factor = 1;

    for (int lev = 0; lev <= finest_level; ++lev) {
        dt_tmp[lev] = std::min(dt_tmp[lev], change_max*dt[lev]);
        n_factor *= nsubsteps[lev];
        dt_0 = std::min(dt_0, n_factor*dt_tmp[lev]);
    }

    // Limit dt's by the value of stop_time.
    const Real eps = 1.e-3*dt_0;

    if (t_new[0] + dt_0 > stop_time - eps) {
        dt_0 = stop_time - t_new[0];
    }

    dt[0] = dt_0;

    for (int lev = 1; lev <= finest_level; ++lev) {
        dt[lev] = dt[lev-1] / nsubsteps[lev];
    }
}

// compute dt from CFL considerations
Real
AmrCoreAdv::EstTimeStep (int lev, Real /*time*/)
{
    BL_PROFILE("AmrCoreAdv::EstTimeStep()");

    Real dt_est = std::numeric_limits<Real>::max();

    const Real* dx  =  geom[lev].CellSize();

    Real coeff = AMREX_D_TERM(   2./(dx[0]*dx[0]),
                               + 2./(dx[1]*dx[1]),
                               + 2./(dx[2]*dx[2]) );
    Real est = 1.0 / (2.0*coeff);
    dt_est = amrex::min(dt_est, est);

    dt_est *= cfl;

    return dt_est;
}

// write plotfile to disk
void
AmrCoreAdv::WritePlotFile () const
{
    const std::string& plotfilename = amrex::Concatenate(plot_file, istep[0], 6);

    // Vector of MultiFabs
    Vector<MultiFab> mf(finest_level+1);
    int ncomp_mf = 2; int src_comp = 0;
    for (int lev = 0; lev <= finest_level; ++lev) {
        mf[lev].define(grids[lev], dmap[lev], ncomp_mf, 0);
        MultiFab::Copy(mf[lev],phi_new[lev],src_comp,0,1,0);
        MultiFab::Copy(mf[lev],phi_new[lev],src_comp,1,1,0);

        // Set the fine data in "phi0" to -1 so we can test on that value and plot particles over blank space
        if (lev == 1) {
            mf[lev].setVal(-1.0,1,1,0);
        }
    }


    Vector<std::string> varnames = {"phi", "phi0"};

    amrex::Print() << "Writing plotfile " << plotfilename << "\n";

    if (finest_level == 0)
    {
        int fake_finest_level = 0;
        WriteMultiLevelPlotfile(plotfilename, fake_finest_level+1, GetVecOfConstPtrs(mf), varnames,
                                Geom(), t_new[0], istep, refRatio());
    } else {

        PhysBCFunctNoOp null_bc_for_fill;

        Vector<IntVect>   r2(finest_level);
        Vector<Geometry>  g2(finest_level+1);
        Vector<MultiFab> mf2(finest_level+1);

        mf2[0].define(grids[0], dmap[0], ncomp_mf, 0);

        // Copy level 0 as is
        MultiFab::Copy(mf2[0],mf[0],0,0,ncomp_mf,0);

        // Define a new multi-level array of Geometry's so that we pass the new "domain" at lev > 0
        Array<int,AMREX_SPACEDIM> periodicity =
                     {AMREX_D_DECL(Geom()[0].isPeriodic(0),Geom()[0].isPeriodic(1),Geom()[0].isPeriodic(2))};
        g2[0].define(Geom()[0].Domain(),&(Geom()[0].ProbDomain()),0,periodicity.data());

        r2[0] = IntVect(AMREX_D_DECL(2,2,2));
        for (int lev = 1; lev <= finest_level; ++lev) {
            if (lev > 1) {
                r2[lev-1][0] = r2[lev-2][0] * 2;
                r2[lev-1][1] = r2[lev-2][1] * 2;
#if (AMREX_SPACEDIM > 2)
                r2[lev-1][2] = r2[lev-2][2] * 2;
#endif
            }

            mf2[lev].define(refine(grids[lev],r2[lev-1]), dmap[lev], ncomp_mf, 0);

            // Set the new problem domain
            Box d2(Geom()[lev].Domain());
            d2.refine(r2[lev-1]);

            g2[lev].define(d2,&(Geom()[lev].ProbDomain()),0,periodicity.data());
        }

        amrex::Vector<amrex::BCRec> bcs_temp;
        bcs_temp.resize(2);     // Setup for 2 components in mf
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            bcs_temp[0].setLo(idim, bcs[0].lo()[idim]);
            bcs_temp[1].setLo(idim, bcs[0].lo()[idim]);

            bcs_temp[0].setHi(idim, bcs[0].hi()[idim]);
            bcs_temp[1].setHi(idim, bcs[0].hi()[idim]);
        }

        // Do piecewise interpolation of mf into mf2
        for (int lev = 1; lev <= finest_level; ++lev) {
            Interpolater* mapper_c = &pc_interp;
            InterpFromCoarseLevel(mf2[lev], t_new[lev], mf[lev],
                                  0, 0, ncomp_mf,
                                  geom[lev], g2[lev],
                                  null_bc_for_fill, 0, null_bc_for_fill, 0,
                                  r2[lev-1], mapper_c, bcs_temp, 0);
        }

        // Define an effective ref_ratio which is isotropic to be passed into WriteMultiLevelPlotfile
        Vector<IntVect> rr(finest_level);
        for (int lev = 0; lev < finest_level; ++lev) {
            rr[lev] = IntVect(AMREX_D_DECL(2,2,2));
        }

       WriteMultiLevelPlotfile(plotfilename, finest_level+1,
                                   GetVecOfConstPtrs(mf2), varnames,
                                   g2, t_new[0], istep, rr);

    }

#ifdef AMREX_PARTICLES
   particleData.writePlotFile(plotfilename,phi_new[1]);
#endif
}

void
AmrCoreAdv::WriteCheckpointFile () const
{

    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save (e.g., finest_level, t_new, etc.) and also
    //                     the BoxArrays at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data at each level of refinement

    // checkpoint file name, e.g., chk00010
    const std::string& checkpointname = amrex::Concatenate(chk_file,istep[0]);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    const int nlevels = finest_level+1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
   if (ParallelDescriptor::IOProcessor()) {

       std::string HeaderFileName(checkpointname + "/Header");
       VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
       std::ofstream HeaderFile;
       HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
       HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                               std::ofstream::trunc |
                                               std::ofstream::binary);
       if( ! HeaderFile.good()) {
           amrex::FileOpenFailed(HeaderFileName);
       }

       HeaderFile.precision(17);

       // write out title line
       HeaderFile << "Checkpoint file for AmrCoreAdv\n";

       // write out finest_level
       HeaderFile << finest_level << "\n";

       // write out array of istep
       for (int i = 0; i < istep.size(); ++i) {
           HeaderFile << istep[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of dt
       for (int i = 0; i < dt.size(); ++i) {
           HeaderFile << dt[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of t_new
       for (int i = 0; i < t_new.size(); ++i) {
           HeaderFile << t_new[i] << " ";
       }
       HeaderFile << "\n";

       // write the BoxArray at each level
       for (int lev = 0; lev <= finest_level; ++lev) {
           boxArray(lev).writeOn(HeaderFile);
           HeaderFile << '\n';
       }
   }

   // write the MultiFab data to, e.g., chk00010/Level_0/
   for (int lev = 0; lev <= finest_level; ++lev) {
       VisMF::Write(phi_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "phi"));
   }

#ifdef AMREX_PARTICLES
   particleData.Checkpoint(checkpointname);
#endif

}

namespace {
// utility to skip to next line in Header
void GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}
}

void
AmrCoreAdv::ReadCheckpointFile ()
{
    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    // Header
    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    GotoNextLine(is);

    // read in array of istep
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            istep[i++] = std::stoi(word);
        }
    }

    // read in array of dt
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            dt[i++] = std::stod(word);
        }
    }

    // read in array of t_new
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            t_new[i++] = std::stod(word);
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {

        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

        // set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        // build MultiFab and FluxRegister data
        int ncomp = 1;
        int ng = 0;
        phi_old[lev].define(grids[lev], dmap[lev], ncomp, ng);
        phi_new[lev].define(grids[lev], dmap[lev], ncomp, ng);

        if (lev > 0 && do_reflux) {
            flux_reg[lev] = std::make_unique<FluxRegister>(grids[lev], dmap[lev], refRatio(lev-1), lev, ncomp);
        }
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {
        VisMF::Read(phi_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "phi"));
    }

#ifdef AMREX_PARTICLES
    particleData.Restart((amrex::ParGDBBase*)GetParGDB(),restart_chkfile);
#endif


}
