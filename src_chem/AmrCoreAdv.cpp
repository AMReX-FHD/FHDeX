
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif
#include <AmrCoreAdv_F.H>
#include <AmrCoreAdv.H>

using namespace amrex;


/*******************************************************************************
 *                                                                             *
 * Initializize problem                                                        *
 *                                                                             * 
 *******************************************************************************/

void AmrCoreAdv::Initialize()
{
//    ReadParameters();

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    phi_new.resize(nlevs_max);

    phi_old.resize(nlevs_max);

    Dphi_x.resize(nlevs_max);
    Dphi_y.resize(nlevs_max);
    Dphi_z.resize(nlevs_max);

    bcs.resize(1);
    
    // periodic boundaries
    int bc_lo[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    int bc_hi[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};

/*
    // walls (Neumann)
    int bc_lo[] = {FOEXTRAP, FOEXTRAP, FOEXTRAP};
    int bc_hi[] = {FOEXTRAP, FOEXTRAP, FOEXTRAP};
*/
//    bcs.resize(1);
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

    // stores fluxes at coarse-fine interface for synchronization this will be
    // sized "nlevs_max+1" NOTE: the flux register associated with
    // flux_reg[lev] is associated with the lev/lev-1 interface (and has grid
    // spacing associated with lev-1) therefore flux_reg[0] is never actually
    // used in the reflux operation
    flux_reg.resize(nlevs_max+1);
}


/*******************************************************************************
 *                                                                             *
 * Advance solution (of the advection diffusion equation) in time for a        *
 * a specified level                                                           *
 *                                                                             *
 *******************************************************************************/

void AmrCoreAdv::EvolveChem(
        std::array<MultiFab, AMREX_SPACEDIM> & umac, 
        const iMultiFab & iface, int lev, int nstep,
        Real dt_fluid)
{
    int time=0;
    dt[lev] = dt_fluid;

    // initialize copies of velocities u_g, v_g, w_g and first derivatives of
    // phi Dphi_x, Dphi_y, Dphi_z
    uface.resize(max_level + 1);
    vface.resize(max_level + 1);
    wface.resize(max_level + 1);


    DistributionMapping phidm = phi_new[lev]->DistributionMap();
    BoxArray phiba            = phi_new[lev]->boxArray();

    BoxArray x_face_ba = phiba;
    BoxArray y_face_ba = phiba;
    BoxArray z_face_ba = phiba;

    x_face_ba.surroundingNodes(0);
    y_face_ba.surroundingNodes(1);
    z_face_ba.surroundingNodes(2);

    for (lev = 0; lev <= finest_level; ++lev) {
        uface[lev].reset(new MultiFab(x_face_ba, phidm, 1, 1));
        vface[lev].reset(new MultiFab(y_face_ba, phidm, 1, 1));
        wface[lev].reset(new MultiFab(z_face_ba, phidm, 1, 1));

        Dphi_x[lev].reset(new MultiFab(phiba, phidm, 1, 0));
        Dphi_y[lev].reset(new MultiFab(phiba, phidm, 1, 0));
        Dphi_z[lev].reset(new MultiFab(phiba, phidm, 1, 0));

        Dphi_x[lev]->setVal(0.);
        Dphi_y[lev]->setVal(0.);
        Dphi_z[lev]->setVal(0.);

       uface[lev]->copy(umac[0], 0, 0, 1, 0, 0);
       vface[lev]->copy(umac[1], 0, 0, 1, 0, 0);
       wface[lev]->copy(umac[2], 0, 0, 1, 0, 0);

       uface[lev]->FillBoundary(geom[lev].periodicity());
       vface[lev]->FillBoundary(geom[lev].periodicity());
       wface[lev]->FillBoundary(geom[lev].periodicity());
    }

    source_loc.reset(new iMultiFab(phiba, phidm, 1, 1));
    source_loc->copy(iface, 0, 0, 1, 0, 0);
    source_loc->FillBoundary(geom[0].periodicity());


    /***************************************************************************
     * Evolve chemical field by integrating time step                          *
     ***************************************************************************/
    lev = 0;
    int iteration = 1;
    timeStep( lev, time, iteration);


    // Commenting out code that writes plot and checkpoint files

    // if (plot_int > 0 && (step+1) % plot_int == 0) {
    //     last_plot_file_step = step+1;
    //     WritePlotFile();
    // }

    // if (chk_int > 0 && (step+1) % chk_int == 0) {
    //     WriteCheckpointFile();
    // }

#ifdef BL_MEM_PROFILING
    {
        std::ostringstream ss;
        ss << "[STEP " << step+1 << "]";
        MemProfiler::report(ss.str());
    }
#endif

    //	if (cur_time >= stop_time - 1.e-6*dt[0]) break;

    //  if (plot_int > 0 && istep[0] > last_plot_file_step) {
    //	    WritePlotFile();
    //  }
}


/*******************************************************************************
 *                                                                             *
 * Initialize multilevel data                                                  *
 *                                                                             *
 *******************************************************************************/

void AmrCoreAdv::InitData ()
{
    Initialize();

    if (restart_chkfile == "") {

        // start simulation from the beginning
        const Real time = 0.0;

        // initialize Levels (done in AmrCore)
        InitFromScratch(time);
        AverageDown();

        // DEBUG: write intial checkpoint
        // if (chk_int > 0) WriteCheckpointFile();

    }
    else {
        // restart from a checkpoint
        ReadCheckpointFile();
    }

   // initialize grad phi
   for (int lev = 0; lev <= finest_level; ++lev) {
       DistributionMapping phidm = phi_new[lev]->DistributionMap();
       BoxArray phiba            = phi_new[lev]->boxArray();

       Dphi_x[lev].reset(new MultiFab(phiba, phidm, 1, 0));
       Dphi_y[lev].reset(new MultiFab(phiba, phidm, 1, 0));
       Dphi_z[lev].reset(new MultiFab(phiba, phidm, 1, 0));

       Dphi_x[lev]->setVal(0.);
       Dphi_y[lev]->setVal(0.);
       Dphi_z[lev]->setVal(0.);
   }

    // DEBUG: write intitial plotfile
    // if (plot_int > 0) WritePlotFile();
}


// Make a new level using provided BoxArray and DistributionMapping and fill
// with interpolated coarse level data. Note: overrides the pure virtual
// function in AmrCore
void AmrCoreAdv::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray & ba,
				         const DistributionMapping & dm)
{

    const int ncomp  = phi_new[lev-1]->nComp();
    const int nghost = phi_new[lev-1]->nGrow();

    phi_new[lev]->define(ba, dm, ncomp, nghost);
    phi_old[lev]->define(ba, dm, ncomp, nghost);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }

    FillCoarsePatch(lev, time, *phi_new[lev], 0, ncomp);
}


// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data. Note: overrides the pure virtual
// function in AmrCore
void AmrCoreAdv::RemakeLevel (int lev, Real time, const BoxArray & ba,
                              const DistributionMapping & dm)
{

    const int ncomp = phi_new[lev]->nComp();
    const int nghost =phi_new[lev]->nGrow();

    MultiFab new_state(ba, dm, ncomp, nghost);
    MultiFab old_state(ba, dm, ncomp, nghost);

    FillPatch(lev, time, new_state, 0, ncomp);

    std::swap(new_state, *phi_new[lev]);
    std::swap(old_state, *phi_old[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }
}

// Delete level data
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::ClearLevel (int lev)
{

    phi_new[lev]->clear();
    phi_old[lev]->clear();
    flux_reg[lev].reset(nullptr);
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void AmrCoreAdv::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
					  const DistributionMapping& dm)
{

    const int ncomp = 1;
    const int nghost = 0;

    phi_new[lev].reset(new MultiFab(ba, dm, ncomp, nghost));
    phi_old[lev].reset(new MultiFab(ba, dm, ncomp, nghost));

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    Real cur_time = t_new[lev];

    MultiFab& state = *phi_new[lev];

    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const int* lo  = box.loVect();
        const int* hi  = box.hiVect();

	initdata(&lev, &cur_time, AMREX_ARLIM_3D(lo), AMREX_ARLIM_3D(hi),
		 BL_TO_FORTRAN_3D(state[mfi]), AMREX_ZFILL(dx),
		 AMREX_ZFILL(prob_lo));
    }

}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
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

    if (lev >= phierr.size()) return;

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    const MultiFab& state = *phi_new[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<int>  itags;

	for (MFIter mfi(state,true); mfi.isValid(); ++mfi)
	{
	    const Box& tilebox  = mfi.tilebox();

            TagBox&     tagfab  = tags[mfi];

	    // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
	    // So we are going to get a temporary integer array.
            // set itags initialle re intent(in) :: dx(3iy to 'untagged' everywhere
            // we define itags over the tilebox region
	    tagfab.get_itags(itags, tilebox);

            // data pointer and index space
	    int*        tptr    = itags.dataPtr();
	    const int*  tlo     = tilebox.loVect();
	    const int*  thi     = tilebox.hiVect();

            // tag cells for refinement
	    state_error(tptr,  AMREX_ARLIM_3D(tlo), AMREX_ARLIM_3D(thi),
			BL_TO_FORTRAN_3D(state[mfi]),
			&tagval, &clearval,
			AMREX_ARLIM_3D(tilebox.loVect()), AMREX_ARLIM_3D(tilebox.hiVect()),
			AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo), &time, &phierr[lev]);
	    //
	    // Now update the tags in the TagBox in the tilebox region
            // to be equal to itags
	    //
	    tagfab.tags_and_untags(itags, tilebox);
	}
    }
}

// read in some parameters from inputs file
void AmrCoreAdv::ReadParameters ()
{
    {
        amrex::Print() << "ReadParameters"<< std::endl;
        ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);
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
        pp.query("diffcoeff", diffcoeff);
        pp.query("do_reflux", do_reflux);
    }
}

// set covered coarse cells to be the average of overlying fine cells
void AmrCoreAdv::AverageDown ()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
        amrex::average_down(
                * phi_new[lev+1], * phi_new[lev],
                     geom[lev+1],   geom[lev],
                0, phi_new[lev]->nComp(), refRatio(lev)
            );
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void AmrCoreAdv::AverageDownTo (int crse_lev)
{
    amrex::average_down(
            * phi_new[crse_lev+1], * phi_new[crse_lev],
                 geom[crse_lev+1],   geom[crse_lev],
            0, phi_new[crse_lev]->nComp(), refRatio(crse_lev)
        );
}

// compute the number of cells at a level
long AmrCoreAdv::CountCells (int lev)
{
    const int N = grids[lev].size();

    long cnt = 0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:cnt)
#endif
    for (int i = 0; i < N; ++i)
        cnt += grids[lev][i].numPts();

    return cnt;
}

// compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
AmrCoreAdv::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    if (lev == 0)
    {
	Vector<MultiFab*> smf;
	Vector<Real> stime;
	GetData(0, time, smf, stime);

//	PhysBCFunct physbc(geom[lev],bcs,BndryFunctBase(phifill));
//	amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
//				     geom[lev], physbc);
        BndryFuncArray bfunc(phifill);
        PhysBCFunct<BndryFuncArray> physbc(geom[lev], bcs, bfunc);
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                    geom[lev], physbc, 0);



  
   }
    else
    {
	Vector<MultiFab*> cmf, fmf;
	Vector<Real> ctime, ftime;
	GetData(lev-1, time, cmf, ctime);
	GetData(lev  , time, fmf, ftime);

//        PhysBCFunct cphysbc(geom[lev-1],bcs,BndryFunctBase(phifill));
//        PhysBCFunct fphysbc(geom[lev  ],bcs,BndryFunctBase(phifill));
        BndryFuncArray bfunc(phifill);
        PhysBCFunct<BndryFuncArray> cphysbc(geom[lev-1],bcs,bfunc);
        PhysBCFunct<BndryFuncArray> fphysbc(geom[lev  ],bcs,bfunc);

	Interpolater* mapper = &cell_cons_interp;

	amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
				   0, icomp, ncomp, geom[lev-1], geom[lev],
				   cphysbc,0, fphysbc, 0, refRatio(lev-1),
				   mapper, bcs,0);
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

    if (cmf.size() != 1) {
	amrex::Abort("FillCoarsePatch: how did this happen?");
    }
    BndryFuncArray bfunc(phifill);
    PhysBCFunct<BndryFuncArray> cphysbc(geom[lev-1],bcs,bfunc);
    PhysBCFunct<BndryFuncArray> fphysbc(geom[lev  ],bcs,bfunc);

//    PhysBCFunct cphysbc(geom[lev-1],bcs,BndryFunctBase(phifill));
//    PhysBCFunct fphysbc(geom[lev  ],bcs,BndryFunctBase(phifill));

    Interpolater* mapper = &cell_cons_interp;

    amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
				 cphysbc,0, fphysbc, 0, refRatio(lev-1),
				 mapper, bcs,0);
}

// utility to copy in data from phi_old and/or phi_new into another multifab
void AmrCoreAdv::GetData (int lev, Real time, Vector<MultiFab *> & data, Vector<Real> & datatime) {

    data.clear();
    datatime.clear();

    MultiFab* phi_new_stdptr=phi_new[lev].get();
    MultiFab* phi_old_stdptr=phi_old[lev].get();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps) {

        data.push_back(phi_new_stdptr);
        datatime.push_back(t_new[lev]);

    } else if (time > t_old[lev] - teps && time < t_old[lev] + teps) {

        data.push_back(phi_old_stdptr);
        datatime.push_back(t_old[lev]);

    } else {

        data.push_back(phi_old_stdptr);
        data.push_back(phi_new_stdptr);
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}



// advance a level by dt
// includes a recursive call for finer levels
void AmrCoreAdv::timeStep (int lev, Real time, int iteration)
{

    if (regrid_int > 0)  // We may need to regrid
    {

        // help keep track of whether a level was already regridded from a
        // coarser level call to regrid
        static Vector<int> last_regrid_step(max_level+1, 0);

        // regrid changes level "lev+1" so we don't regrid on max_level also
        // make sure we don't regrid fine levels again if it was taken care of
        // during a coarser regrid
        if (lev < max_level && istep[lev] > last_regrid_step[lev]) {
            if (istep[lev] % regrid_int == 0) {

                // regrid could add newly refine levels (if finest_level < max_level)
                // so we save the previous finest level index
                int old_finest = finest_level;
                regrid(lev, time);
                // mark that we have regridded this level already
                for (int k = lev; k <= finest_level; ++k) {
                    last_regrid_step[k] = istep[k];
                }

                // if there are newly created levels, set the time step
                for (int k = old_finest+1; k <= finest_level; ++k) {
                    dt[k] = dt[k-1] / MaxRefRatio(k-1);
                }
            }
        }
    }


    // advance a single level for a single time step, updates flux registers
    Advance(lev, time, dt[lev], iteration, nsubsteps[lev]);

    ++istep[lev];

    if (Verbose()) {
        amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
        amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
    }

    if (lev < finest_level)
    {
        // recursive call for next-finer level
        for (int i = 1; i <= nsubsteps[lev+1]; ++i)
        {
            timeStep(lev+1, time+(i-1)*dt[lev+1], i);

        }

        if (do_reflux)
        {
         // update lev based on coarse-fine flux mismatch
            flux_reg[lev+1]->Reflux(*phi_new[lev], 1.0, 0, 0, phi_new[lev]->nComp(), geom[lev]);
        }

        AverageDownTo(lev); // average lev+1 down to lev
    }

}


// advance a single level for a single time step, updates flux registers
void AmrCoreAdv::Advance (int lev, Real time, Real dt_lev, int iteration, int ncycle) {

    constexpr int num_grow = 3;

    std::swap(phi_old[lev], phi_new[lev]);
    t_old[lev]  = t_new[lev];
    t_new[lev] += dt_lev;

    const BoxArray & badp            = phi_new[lev]->boxArray();
    const DistributionMapping & dmdp = phi_new[lev]->DistributionMap();

    MultiFab ptSource(badp,dmdp,1,0);

    ptSource.setVal(0.);

    Vector<int> xloc;
    Vector<int> yloc;
    Vector<int> zloc;
    Real strength = 0.001;

    MultiFab &  S_new     = * phi_new[lev];
    MultiFab &  uface_lev = * uface[lev];
    MultiFab &  vface_lev = * vface[lev];
    MultiFab &  wface_lev = * wface[lev];
    iMultiFab & sloc_mf   = * source_loc;

    const Real * dx      = geom[lev].CellSize();
    const Real * prob_lo = geom[lev].ProbLo();

    MultiFab fluxes[BL_SPACEDIM];
    if (do_reflux) {
        for (int i = 0; i < BL_SPACEDIM; ++i) {
            BoxArray ba = grids[lev];
            ba.surroundingNodes(i);
            fluxes[i].define(ba, dmap[lev], S_new.nComp(), 0);
        }
    }

    // State with ghost cells
    MultiFab Sborder(grids[lev], dmap[lev], S_new.nComp(), num_grow);
    FillPatch(lev, time, Sborder, 0, Sborder.nComp());


    static Vector<Real> ib_init_pos;

    {
        ParmParse pp("mfix");
        int n = pp.countval("ib_init__pos");
        if (n > 0) pp.getarr("ib_init__pos",ib_init_pos, 0, n);
    }


#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox flux[BL_SPACEDIM];

        for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {
            const Box & bx = mfi.tilebox();

            const FArrayBox & statein =   Sborder[mfi];
            FArrayBox & stateout      =     S_new[mfi];
            FArrayBox & ptS           =  ptSource[mfi];
            IArrayBox & fabsl 	      =   sloc_mf[mfi];
            FArrayBox & uface_mf      = uface_lev[mfi];
            FArrayBox & vface_mf      = vface_lev[mfi];
            FArrayBox & wface_mf      = wface_lev[mfi];

            for (int i = 0; i < BL_SPACEDIM ; i++) {
                const Box& bxtmp = amrex::surroundingNodes(bx,i);
                flux[i].resize(bxtmp,S_new.nComp());
            }


            if (BL_SPACEDIM==2) {
                // fill the point source multifab from the tagged interface multifab
                get_ptsource_2d( bx.loVect(), bx.hiVect(),
                                 BL_TO_FORTRAN_3D(fabsl),
                                 BL_TO_FORTRAN_3D(ptS),
                                 & strength, dx, & ib_init_pos[0], & ib_init_pos[1],
                                 AMREX_ZFILL(prob_lo));

                // compute new state (stateout) and fluxes.
                advect_2d( & time, bx.loVect(), bx.hiVect(),
                           BL_TO_FORTRAN_3D(statein),
                           BL_TO_FORTRAN_3D(stateout),
                           BL_TO_FORTRAN_3D(ptS),
                           BL_TO_FORTRAN_3D(fabsl),
                           AMREX_D_DECL(BL_TO_FORTRAN_3D(uface_mf),
                                        BL_TO_FORTRAN_3D(uface_mf),
                                        BL_TO_FORTRAN_3D(uface_mf)),
                           AMREX_D_DECL(BL_TO_FORTRAN_3D(flux[0]),
                                        BL_TO_FORTRAN_3D(flux[1]),
                                        BL_TO_FORTRAN_3D(flux[2])),
                           dx, & dt_lev, &diffcoeff);
            } else {
                // fill the point source multifab from the tagged interface multifab
                get_ptsource_3d( bx.loVect(), bx.hiVect(),
                                 BL_TO_FORTRAN_3D(fabsl),
                                 BL_TO_FORTRAN_3D(ptS),
                                 & strength, dx, & ib_init_pos[0], & ib_init_pos[1],
                                 & ib_init_pos[2], AMREX_ZFILL(prob_lo));

                // compute new state (stateout) and fluxes.
                advect_3d(& time, bx.loVect(), bx.hiVect(),
                          BL_TO_FORTRAN_3D(statein),
                          BL_TO_FORTRAN_3D(stateout),
                          BL_TO_FORTRAN_3D(ptS),
                          BL_TO_FORTRAN_3D(fabsl),
                          AMREX_D_DECL(BL_TO_FORTRAN_3D(uface_mf),
                                       BL_TO_FORTRAN_3D(vface_mf),
                                       BL_TO_FORTRAN_3D(wface_mf)),
                          AMREX_D_DECL(BL_TO_FORTRAN_3D(flux[0]),
                                       BL_TO_FORTRAN_3D(flux[1]),
                                       BL_TO_FORTRAN_3D(flux[2])),
                          dx, & dt_lev, & diffcoeff);}


            if (do_reflux) {
                for (int i = 0; i < BL_SPACEDIM ; i++) {
                    fluxes[i][mfi].copy(flux[i],mfi.nodaltilebox(i));
                }
            }
        }
    }

    // After updating phi_new we compute the first derivatives
    MultiFab &  sx_mf       = * Dphi_x[lev];
    MultiFab &  sy_mf       = * Dphi_y[lev];
    MultiFab &  sz_mf       = * Dphi_z[lev];

    // phi_new including 1 ghost cell
    MultiFab S_new_fill(grids[lev], dmap[lev], S_new.nComp(), 1);
    S_new_fill.copy(S_new, 0, 0,1, 0, 0);
    S_new_fill.FillBoundary(geom[lev].periodicity());


#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter mfi(sx_mf, true); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();

            FArrayBox & stateout      =   S_new_fill[mfi];
            IArrayBox & fabsl         =      sloc_mf[mfi];
            FArrayBox & fabsx         =        sx_mf[mfi];
            FArrayBox & fabsy         =        sy_mf[mfi];
            FArrayBox & fabsz         =        sz_mf[mfi];

            // compute velocities on faces (prescribed function of space and time)
            if (BL_SPACEDIM==2) {
                get_phigrad_2d( bx.loVect(), bx.hiVect(),
                                BL_TO_FORTRAN_3D(stateout),
                                BL_TO_FORTRAN_3D(fabsx),
                                BL_TO_FORTRAN_3D(fabsy),
                                BL_TO_FORTRAN_3D(fabsl),
                                & ib_init_pos[0], & ib_init_pos[1],
                                dx, AMREX_ZFILL(prob_lo));
            } else {
                get_phigrad_3d( bx.loVect(), bx.hiVect(),
                                BL_TO_FORTRAN_3D(stateout),
                                BL_TO_FORTRAN_3D(fabsx),
                                BL_TO_FORTRAN_3D(fabsy),
                                BL_TO_FORTRAN_3D(fabsz),
                                BL_TO_FORTRAN_3D(fabsl),
                                & ib_init_pos[0], & ib_init_pos[1],
                                & ib_init_pos[2], dx, AMREX_ZFILL(prob_lo));
            }
        }
    }

    amrex::Print() << "simulated phi total"<< (phi_new[lev]->sum(0,false));
    amrex::Print() << "true phi total"<< ptSource.sum(0,false)*(time+dt[0])<< std::endl;

    // increment or decrement the flux registers by area and time-weighted
    // fluxes Note that the fluxes have already been scaled by dt and area In
    // this example we are solving phi_t = -div(+F) The fluxes contain, e.g.,
    // F_{i+1/2,j} = (phi*u)_{i+1/2,j} Keep this in mind when considering the
    // different sign convention for updating the flux registers from the coarse
    // or fine grid perspective NOTE: the flux register associated with
    // flux_reg[lev] is associated with the lev/lev-1 interface (and has grid
    // spacing associated with lev-1)

    if (do_reflux) {
        if (flux_reg[lev+1]) {
            for (int i = 0; i < BL_SPACEDIM; ++i) {
                // update the lev+1/lev flux register (index lev+1)
                flux_reg[lev+1]->CrseInit(fluxes[i],i,0,0,fluxes[i].nComp(), -1.0);
            }
        }
        if (flux_reg[lev]) {
            for (int i = 0; i < BL_SPACEDIM; ++i) {
                // update the lev/lev-1 flux register (index lev)
                flux_reg[lev]->FineAdd(fluxes[i],i,0,0,fluxes[i].nComp(), 1.0);
            }
        }
    }
}

// copy phi or one of it's first derivatives into a multifab mf, mf doesn't
// have to have the same dm and box array
void AmrCoreAdv::phi_new_copy(int lev, Vector<std::unique_ptr<MultiFab>> & MF,
                              int dcomp, int indicator) {
    // indicator gives which quantity is being copied into MF
    if (indicator==0)      MF[lev]->copy(* phi_new[lev], 0, dcomp,1, 0, 0);
    else if (indicator==1) MF[lev]->copy(* Dphi_x[lev], 0, dcomp,1, 0, 0);
    else if (indicator==2) MF[lev]->copy(* Dphi_y[lev], 0, dcomp,1, 0, 0);
    else if (indicator==3) MF[lev]->copy(* Dphi_z[lev], 0, dcomp,1, 0, 0);
    else amrex::Abort( "Incorrect indicator for copying information from AmrCoreAdv to Mfix" );
}


void AmrCoreAdv::ReadCheckpointFile ()
{
    amrex::Print() << "ReadCheckpointFile"<< std::endl;

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
        int nghost = 0;

        phi_old[lev]->define(grids[lev], dmap[lev], ncomp, nghost);

        phi_new[lev]->define(grids[lev], dmap[lev], ncomp, nghost);

        if (lev > 0 && do_reflux) {
            flux_reg[lev].reset(new FluxRegister(grids[lev], dmap[lev], refRatio(lev-1), lev, ncomp));
        }
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {

        VisMF::Read(*phi_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "phi"));
    }

}

// utility to skip to next line in Header
void
AmrCoreAdv::GotoNextLine (std::istream& is)
{

    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}
