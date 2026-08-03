
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
#include <init_metric.H>
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

    for (int step = istep[0]; step < max_step + ntherm && cur_time < stop_time; ++step)
    {
        amrex::Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;

        ComputeDt();

        int lev = 0;
        int iteration = 1;
/*
       const auto dx     = geom[lev].CellSizeArray();
       const auto problo = geom[lev].ProbLoArray();

       for (MFIter mfi(gmetric); mfi.isValid(); ++mfi)
        {
            const Box& gbx = mfi.validbox();
            auto const& gmet_arr = gmetric.array(mfi);
            auto const& gsqr_arr = sqrgmetric.array(mfi);
            auto const& detg_arr = detg.array(mfi);
            auto const& newgmet_arr = newgmetric.array(mfi);
            auto const& newgsqr_arr = newsqrgmetric.array(mfi);
            auto const& newdetg_arr = newdetg.array(mfi);
            amrex::ParallelFor(gbx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                new_dk_metric(i,j,k,gmet_arr,gsqr_arr,detg_arr,newgmet_arr,newgsqr_arr,newdetg_arr,dx,problo,cur_time);
            });
        }
*/

        timeStepNoSubcycling(cur_time, iteration);

//  accumulate statistics
        if(nstat > 0 && step > ntherm && (step - ntherm) % nstat == 0){

           statpts +=1.;

           const auto dx     = Geom(lev).CellSizeArray();
           amrex::Real cell_vol = dx[0]*dx[1];

           for (MFIter mfi(phi_new[lev]); mfi.isValid(); ++mfi)
           {
               const Box& vbx = mfi.validbox();
               auto const& phi_arr = phi_new[lev].array(mfi);
               auto const& stats_arr = stats.array(mfi);
               auto const& det_arr = detg.array(mfi);
               amrex::ParallelFor(vbx,
               [=] AMREX_GPU_DEVICE(int i, int j, int k)
               {
                  amrex::Real measure, nparts;
                  if(pure_part == 0) {
                     measure = phi_arr(i,j,k,0);
                     nparts = num_part*phi_arr(i,j,k,0)*cell_vol*det_arr(i,j,k,1);
                  } else {
                     nparts = phi_arr(i,j,k,0);
                     measure = phi_arr(i,j,k)/(cell_vol*det_arr(i,j,k,1)*num_part);
                  }

                  stats_arr(i,j,k,0) += measure;
                  stats_arr(i,j,k,1) += measure*measure;
                  stats_arr(i,j,k,2) += nparts;
                  stats_arr(i,j,k,3) += nparts*nparts;
                 });
           }


        }

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
            AverageDown();
        }
#endif
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

    if (lev == 0) {

       const auto problo = Geom(lev).ProbLoArray();
       const auto probhi = Geom(lev).ProbHiArray();
       const auto dx     = Geom(lev).CellSizeArray();

       stats.define(ba,dm,4,0);
       stats.setVal(0.);

       gmetric.define(ba,dm,4,ng);
       sqrgmetric.define(ba,dm,4,ng);
       detg.define(ba,dm,3,ng);

       newgmetric.define(ba,dm,4,ng);
       newsqrgmetric.define(ba,dm,4,ng);
       newdetg.define(ba,dm,3,ng);

       surf_area = 0.;

       for (MFIter mfi(gmetric); mfi.isValid(); ++mfi)
        {
            const Box& gbx = mfi.validbox();
            auto const& gmet_arr = newgmetric.array(mfi);
            auto const& gsqr_arr = newsqrgmetric.array(mfi);
            auto const& detg_arr = newdetg.array(mfi);
            amrex::ParallelFor(gbx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                init_dk_metric(i,j,k,gmet_arr,gsqr_arr,detg_arr,dx,problo);
                surf_area += dx[0]*dx[1]*detg_arr(i,j,k,1);
            });
        }

        ParallelDescriptor:: ReduceRealSum(surf_area);

        amrex::Print() << "total surface area = " << surf_area <<std::endl;
        gmetric.FillBoundary(Geom(lev).periodicity());
        sqrgmetric.FillBoundary(Geom(lev).periodicity());
        detg.FillBoundary(Geom(lev).periodicity());

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

    amrex::Print() << "pure " << pure_part << " " << max_level << std::endl;

    if (lev == 0) {

        MultiFab volume(ba,dm,1,0);
        volume.setVal(dx[0]*dx[1]);
        MultiFab::Multiply(volume, newdetg, 1, 0, 1, 0);

        for (MFIter mfi(phi_new[lev]); mfi.isValid(); ++mfi)
        {
            const Box& vbx = mfi.validbox();
            auto const& phi_arr = phi_new[lev].array(mfi);
            auto const& det_arr = newdetg.array(mfi);
            amrex::ParallelFor(vbx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                init_phi(i,j,k,phi_arr,det_arr,dx,problo,surf_area,pure_part,Ncomp,ext_pot);
            });
        }

                // amrex::Real increm = phi_arr(i,j,k,0)*dx[0]*dx[1]*det_arr(i,j,k,1);
                // amrex::HostDevice::Atomic::Add(phisum, increm);
             //   phisum += phi_arr(i,j,k,0)*dx[0]*dx[1]*det_arr(i,j,k,1);
            // });

        //
        // Note that when we send in local = true, NO ParallelAllReduce::Sum
        //      is called inside the Dot product -- we will do that before we print
        //
#if 0
        amrex::Real phisum2 = 0.;
        for (MFIter mfi(phi_new[lev]); mfi.isValid(); ++mfi)
        {
        auto const& phi_arr = phi_new[lev].array(mfi);
        for(int i=0; i<32; i++){
        for(int j=0; j<32; i++){
            phisum2+= phi_arr(i,j,0,0);
        }
        }
        }

        amrex::Print() << "integral of phi = " << phisum2 <<std::endl;
#endif

   if(pure_part == 0){
        bool local = false;
        amrex::Real phisum = MultiFab::Dot(phi_new[lev], 0, volume, 0, 1, 0, local);
        amrex::Real partsum = phisum*num_part;

        amrex::Print() << "integral of phi = " << phisum <<std::endl;
        amrex::Print() << "integral of particle = " << partsum <<  " " << partsum/phisum <<std::endl;

        if(std::abs(phisum-1.) > 1.e-10){

        amrex::Print{} << " phi sum in init " << phisum << std::endl;

        for (MFIter mfi(phi_new[lev]); mfi.isValid(); ++mfi)
        {
            const Box& vbx = mfi.validbox();
            auto const& phi_arr = phi_new[lev].array(mfi);
            amrex::ParallelFor(vbx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
//                phi_rescale(i,j,k,phi_arr,phisum);
                  phi_arr(i,j,k,0) /= phisum;
            });
        }

          bool local = false;
          amrex::Real phisum = MultiFab::Dot(phi_new[lev], 0, volume, 0, 1, 0, local);

          amrex::Print() << "integral of phi after rescale = " << phisum <<std::endl;


        }

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

        pp.query("num_part",num_part);
        pp.query("dorand",dorand);

        pp.query("nstat",nstat);
        pp.query("ntherm",ntherm);

        alg_type = 0;
        pp.queryAdd("alg_type", alg_type);

        num_flux = 1;
        pp.query("num_flux",num_flux);

        ext_pot = 0;
        pp.query("ext_pot",ext_pot);


        // read in BC; see Src/Base/AMReX_BC_TYPES.H for supported types
        pp.queryarr("bc_lo", bc_lo);
        pp.queryarr("bc_hi", bc_hi);

        seed = 0;
        pp.queryAdd("seed", seed);
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



#ifdef AMREX_PARTICLES
        particleData.init_particle_params(max_level);
     if(max_level == 1){
        pure_part = 1;
     }
#endif
}

// set covered coarse cells to be the average of overlying fine cells
void AmrCoreAdv::AverageDown ()
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

    AdvancePhiAtLevel(lev, t_new[lev], dt[lev], iteration, nsub);

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
    int ncomp_mf = 11; int src_comp = 0;
    for (int lev = 0; lev <= finest_level; ++lev) {
        mf[lev].define(grids[lev], dmap[lev], ncomp_mf, 1);
        MultiFab::Copy(mf[lev],phi_new[lev],src_comp,0,1,0);
        MultiFab::Copy(mf[lev],phi_new[lev],src_comp,1,1,0);

        // Set the fine data in "phi0" to -1 so we can test on that value and plot particles over blank space
        if (lev == 1) {
            mf[lev].setVal(-1.0,1,1,0);
        }
    }

    int lev = 0;

    Vector<MultiFab> mf_nd(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        BoxArray nodal_ba(grids[lev]); nodal_ba.surroundingNodes();
        mf_nd[lev].define(nodal_ba, dmap[lev], 3, 0);
        mf_nd[lev].setVal(0.0);
    }

    const auto dx     = Geom(lev).CellSizeArray();
    amrex::Real cell_vol = dx[0]*dx[1];

    const auto prob_lo = Geom(lev).ProbLoArray();

    for (MFIter mfi(mf[lev]); mfi.isValid(); ++mfi)
    {
            const Box& vbx = mfi.validbox();
            auto const& mf_arr = mf[lev].array(mfi);
            auto const& phi_arr = phi_new[lev].array(mfi);
            auto const& det_arr = newdetg.array(mfi);
            auto const& stat_arr = stats.array(mfi);
            amrex::ParallelFor(vbx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                amrex::Real measure, nparts;
                if(pure_part == 0) {
                   measure = phi_arr(i,j,k,0);
                   nparts = num_part*phi_arr(i,j,k,0)*cell_vol*det_arr(i,j,k,1);
                } else {
                   nparts = phi_arr(i,j,k,0);
                   measure = phi_arr(i,j,k,0)/(cell_vol*det_arr(i,j,k,1)*num_part);
                }

                mf_arr(i,j,k,2) = det_arr(i,j,k,1);
                mf_arr(i,j,k,3) = nparts;

                amrex::Real twopi = 2.*3.14159265358979323846264338;
                amrex::Real pi = 3.14159265358979323846264338;

                amrex::Real twopisq = twopi*twopi;

                amrex::Real divisor = std::max(statpts,1.);
                mf_arr(i,j,k,4) = stat_arr(i,j,k,0)/divisor;
                mf_arr(i,j,k,5) = stat_arr(i,j,k,1)/divisor - mf_arr(i,j,k,4)*mf_arr(i,j,k,4);
                mf_arr(i,j,k,6) = 1./(cell_vol*det_arr(i,j,k,1)*num_part*surf_area);
                mf_arr(i,j,k,7) = stat_arr(i,j,k,2)/divisor;
                mf_arr(i,j,k,8) = stat_arr(i,j,k,3)/divisor - mf_arr(i,j,k,7)*mf_arr(i,j,k,7);
                mf_arr(i,j,k,9) = cell_vol*det_arr(i,j,k,1)*num_part/surf_area;

#if 0
                amrex::Real amp = 3.;


                amrex::Real xloc = prob_lo[0] + (i+0.5) * dx[0] ;
                amrex::Real yloc = prob_lo[1] + (j+0.5) * dx[1] ;
                amrex::Real sinx,siny;

                sinx = std::sin(xloc);
                siny = std::sin(yloc);

                mf_arr(i,j,k,10) = amp*sinx*siny;

#else
                amrex::Real amp = 4.;


                amrex::Real xloc = prob_lo[0] + (i+0.5) * dx[0] ;
                amrex::Real yloc = prob_lo[1] + (j+0.5) * dx[1] ;
/*
                amrex::Real sinx,siny;
                sinx = std::sin(xloc);
                siny = std::sin(yloc);

                amrex::Real caps, capsy, capr, caprx;

                caps  = 0.;

                capr  = std::pow(std::cos(0.5*(xloc-pi)),4);

                for (int iterm  = 0; iterm < 3; iterm++)
                {
                   amrex::Real shift = pi/3. + 2.*iterm*pi/3.;
                   caps += std::pow(std::cos(0.5*(yloc-shift)),16);
                }


//                mf_arr(i,j,k,10) = amp*sinx*sinx*siny*siny;
                mf_arr(i,j,k,10) = amp*caps*capr;
*/
                amrex::Real tscale = .3;


                 amrex::Real sinx,siny,cosx,cosy;

                 sinx = std::sin(xloc);
                 siny = std::sin(yloc);
                 cosx = std::cos(xloc);
                 cosy = std::cos(yloc);

                 amrex::Real aoft = std::cos(0.5*pi * std::min(t_new[0],tscale) / tscale)*std::cos(0.5*pi * std::min(t_new[0],tscale) / tscale);
                 mf_arr(i,j,k,10) = amp*(aoft*sinx*sinx*siny*siny+(1.-aoft)*cosx*cosx*cosy*cosy);
#endif

#if 0

                amrex::Real sinxpy = std::sin(0.5*(xloc+yloc));
                amrex::Real cosxpy = std::cos(0.5*(xloc+yloc));
                amrex::Real sinxmy = std::sin(0.5*(xloc-yloc));
                amrex::Real cosxmy = std::cos(0.5*(xloc-yloc));
                amrex::Real gamma = 0.2;

                mf_arr(i,j,k,1) = std::exp(-sinxpy*sinxpy*sinxmy*sinxmy/gamma);
#else

                amrex::Real alpha = .3;
                amrex::Real beta = 0.;
                amrex::Real delta = 2.;
                amrex::Real gamma = 1.e-1;

                amrex::Real caps, capsy, capr, caprx;

                caps = 0;
                capr  = (delta-alpha*std::pow(std::sin(xloc),2))*(1.+std::pow(std::cos(xloc/2),6));


             for (int iterm  = 0; iterm < 3; iterm++)
             {
                amrex::Real shift = pi/3. + 2.*iterm*pi/3.;
                caps += std::pow(std::cos(0.5*(yloc-shift)),16);
             }

                mf_arr(i,j,k,1) = std::exp(-(caps+capr)/gamma);
#endif

            });
    }

    mf[lev].FillBoundary(geom[lev].periodicity());

    int nd_lev = 0;
    for (MFIter mfi(mf_nd[nd_lev]); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& mf_arr    = mf[nd_lev].const_array(mfi);
        auto const& mf_nd_arr = mf_nd[nd_lev].array(mfi);
        amrex::Real pi = 3.14159265358979323846264338;
        amrex::Print() << " FILLING NODAL MF ON " << vbx << std::endl;
        amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            mf_nd_arr(i,j,k,0) = mf_arr(i,j,k,2);
            mf_nd_arr(i,j,k,1) = mf_arr(i,j,k,3);
//            mf_nd_arr(i,j,k,2) = mf_arr(i,j,k,10);
            amrex::Real amp = 4.;
            amrex::Real tscale = .3;


            amrex::Real xloc = prob_lo[0] + (i) * dx[0] ;
            amrex::Real yloc = prob_lo[1] + (j) * dx[1] ;
            amrex::Real sinx,siny,cosx,cosy;

            sinx = std::sin(xloc);
            siny = std::sin(yloc);
            cosx = std::cos(xloc);
            cosy = std::cos(yloc);

            amrex::Real aoft = std::cos(0.5*pi * std::min(t_new[0],tscale) / tscale)*std::cos(0.5*pi * std::min(t_new[0],tscale) / tscale);
            mf_nd_arr(i,j,k,2) = amp*(aoft*sinx*sinx*siny*siny+(1.-aoft)*cosx*cosx*cosy*cosy);

/*
            amrex::Real caps, capsy, capr, caprx;

            caps  = 0.;

            capr  = std::pow(std::cos(0.5*(xloc-pi)),4);

            for (int iterm  = 0; iterm < 3; iterm++)
            {
               amrex::Real shift = pi/3. + 2.*iterm*pi/3.;
               caps += std::pow(std::cos(0.5*(yloc-shift)),16);
            }

            mf_nd_arr(i,j,k,2) = amp*caps*capr;
*/
//            mf_nd_arr(i,j,k,2) = amp*sinx*sinx*siny*siny;
        });
    }


    Vector<std::string> varnames = {"phi", "pot_eq","metric","number","mean_mu","var_mu","var_mu_th","mean_N","var_N","var_N_th","surface"};

    amrex::Print() << "Writing plotfile " << plotfilename << "\n";

//    if (finest_level == 0)
//    {
          int fake_finest_level = 0;
          WriteMultiLevelPlotfileWithSurface(plotfilename, fake_finest_level+1,
                                             GetVecOfConstPtrs(mf),
                                             GetVecOfConstPtrs(mf_nd),
                                             varnames,
                                             Geom(), t_new[0], istep, refRatio());
          // WriteMultiLevelPlotfile(plotfilename, fake_finest_level+1, GetVecOfConstPtrs(mf), varnames,
          //                         Geom(), t_new[0], istep, refRatio());
#if 0
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
        bcs_temp.resize(4);     // Setup for 2 components in mf
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            bcs_temp[0].setLo(idim, bcs[0].lo()[idim]);
            bcs_temp[1].setLo(idim, bcs[0].lo()[idim]);
            bcs_temp[2].setLo(idim, bcs[0].lo()[idim]);
            bcs_temp[3].setLo(idim, bcs[0].lo()[idim]);

            bcs_temp[0].setHi(idim, bcs[0].hi()[idim]);
            bcs_temp[1].setHi(idim, bcs[0].hi()[idim]);
            bcs_temp[2].setHi(idim, bcs[0].hi()[idim]);
            bcs_temp[3].setHi(idim, bcs[0].hi()[idim]);
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
#endif

#ifdef AMREX_PARTICLES
   if(max_level == 1){
      particleData.writePlotFile(plotfilename,phi_new[1]);
   }
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

void
AmrCoreAdv::WriteMultiLevelPlotfileWithSurface (const std::string& plotfilename, int nlevels,
                                         const Vector<const MultiFab*>& mf,
                                         const Vector<const MultiFab*>& mf_nd,
                                         const Vector<std::string>& varnames,
                                         const Vector<Geometry>& my_geom,
                                         Real time,
                                         const Vector<int>& level_steps,
                                         const Vector<IntVect>& rr,
                                         const std::string &versionName,
                                         const std::string &levelPrefix,
                                         const std::string &mfPrefix,
                                         const Vector<std::string>& extra_dirs) const
{
    BL_PROFILE("WriteMultiLevelPlotfileWithSurface()");

    AMREX_ALWAYS_ASSERT(nlevels <= mf.size());
    AMREX_ALWAYS_ASSERT(nlevels <= rr.size()+1);
    AMREX_ALWAYS_ASSERT(nlevels <= level_steps.size());
    AMREX_ALWAYS_ASSERT(mf[0]->nComp() == varnames.size());

    bool callBarrier(false);
    PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);
    if (!extra_dirs.empty()) {
        for (const auto& d : extra_dirs) {
            const std::string ed = plotfilename+"/"+d;
            PreBuildDirectorHierarchy(ed, levelPrefix, nlevels, callBarrier);
        }
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::MyProc() == ParallelDescriptor::NProcs()-1) {
        Vector<BoxArray> boxArrays(nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = mf[level]->boxArray();
        }

        auto f = [=]() {
            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
            std::string HeaderFileName(plotfilename + "/Header");
            std::ofstream HeaderFile;
            HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
            HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                    std::ofstream::trunc |
                                                    std::ofstream::binary);
            if( ! HeaderFile.good()) FileOpenFailed(HeaderFileName);
            WriteGenericPlotfileHeaderWithSurface(HeaderFile, nlevels, boxArrays, varnames,
                                                  my_geom, time, level_steps, rr, versionName,
                                                  levelPrefix, mfPrefix);
        };

        if (AsyncOut::UseAsyncOut()) {
            AsyncOut::Submit(std::move(f));
        } else {
            f();
        }
    }

    std::string mf_nodal_prefix = "sfc_nd";
    for (int level = 0; level <= finest_level; ++level)
    {
        if (AsyncOut::UseAsyncOut()) {
            VisMF::AsyncWrite(*mf[level],
                              MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix),
                              true);
            VisMF::AsyncWrite(*mf_nd[level],
                              MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_nodal_prefix),
                              true);
        } else {
            const MultiFab* data;
            std::unique_ptr<MultiFab> mf_tmp;
            if (mf[level]->nGrowVect() != 0) {
                mf_tmp = std::make_unique<MultiFab>(mf[level]->boxArray(),
                                                    mf[level]->DistributionMap(),
                                                    mf[level]->nComp(), 0, MFInfo(),
                                                    mf[level]->Factory());
                MultiFab::Copy(*mf_tmp, *mf[level], 0, 0, mf[level]->nComp(), 0);
                data = mf_tmp.get();
            } else {
                data = mf[level];
            }
            VisMF::Write(*data        , MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
            VisMF::Write(*mf_nd[level], MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_nodal_prefix));
        }
    }
}

void
AmrCoreAdv::WriteGenericPlotfileHeaderWithSurface (std::ostream &HeaderFile,
                                            int nlevels,
                                            const Vector<BoxArray> &bArray,
                                            const Vector<std::string> &varnames,
                                            const Vector<Geometry>& my_geom,
                                            Real my_time,
                                            const Vector<int>& level_steps,
                                            const Vector<IntVect>& my_ref_ratio,
                                            const std::string &versionName,
                                            const std::string &levelPrefix,
                                            const std::string &mfPrefix) const
{
    AMREX_ALWAYS_ASSERT(nlevels <= bArray.size());
    AMREX_ALWAYS_ASSERT(nlevels <= my_ref_ratio.size()+1);
    AMREX_ALWAYS_ASSERT(nlevels <= level_steps.size());

    HeaderFile.precision(17);

    // ---- this is the generic plot file type name
    HeaderFile << versionName << '\n';

    HeaderFile << varnames.size() << '\n';

    for (int ivar = 0; ivar < varnames.size(); ++ivar) {
        HeaderFile << varnames[ivar] << "\n";
    }
    HeaderFile << AMREX_SPACEDIM << '\n';
    HeaderFile << my_time << '\n';
    HeaderFile << finest_level << '\n';
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        HeaderFile << my_geom[0].ProbLo(i) << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        HeaderFile << my_geom[0].ProbHi(i) << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i < finest_level; ++i) {
        HeaderFile << my_ref_ratio[i][0] << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i <= finest_level; ++i) {
        HeaderFile << my_geom[i].Domain() << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i <= finest_level; ++i) {
        HeaderFile << level_steps[i] << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i <= finest_level; ++i) {
        for (int k = 0; k < AMREX_SPACEDIM; ++k) {
            HeaderFile << my_geom[i].CellSize()[k] << ' ';
        }
        HeaderFile << '\n';
    }
    HeaderFile << (int) my_geom[0].Coord() << '\n';
    HeaderFile << "0\n";

   for (int level = 0; level <= finest_level; ++level) {
        HeaderFile << level << ' ' << bArray[level].size() << ' ' << my_time << '\n';
        HeaderFile << level_steps[level] << '\n';

        const IntVect& domain_lo = my_geom[level].Domain().smallEnd();
        for (int i = 0; i < bArray[level].size(); ++i)
        {
            // Need to shift because the RealBox ctor we call takes the
            // physical location of index (0,0,0).  This does not affect
            // the usual cases where the domain index starts with 0.
            const Box& b = shift(bArray[level][i], -domain_lo);
            RealBox loc = RealBox(b, my_geom[level].CellSize(), my_geom[level].ProbLo());
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
            }
        }

        HeaderFile << MultiFabHeaderPath(level, levelPrefix, mfPrefix) << '\n';
    }
    HeaderFile << "1" << "\n";
    HeaderFile << "3" << "\n";
    HeaderFile << "mu_nd" << "\n";
    HeaderFile << "num_nd" << "\n";
    HeaderFile << "surf_nd" << "\n";
    std::string mf_nodal_prefix = "sfc_nd";
    for (int level = 0; level <= finest_level; ++level) {
        HeaderFile << MultiFabHeaderPath(level, levelPrefix, mf_nodal_prefix) << '\n';
    }
}
