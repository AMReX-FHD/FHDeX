
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_EB_levelset.H>

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif
#include <AmrCoreAdv_F.H>
#include <AmrCoreAdv.H>
#include <common_functions.H>
using namespace amrex;


/*******************************************************************************
 *                                                                             *
 * Initializize problem                                                        *
 *                                                                             *
 *******************************************************************************/

void AmrCoreAdv::Initialize( )
{

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    // initialize vectors found at the end of AmrCoreAdv.cpp
    con_pre.resize(nlevs_max);
    con_new.resize(nlevs_max);

    con_old.resize(nlevs_max);

    Dcon_x.resize(nlevs_max);
    Dcon_y.resize(nlevs_max);
    Dconc_x.resize(nlevs_max);
    Dconc_y.resize(nlevs_max);

    MagDcon.resize(nlevs_max);
    Dcon_z.resize(nlevs_max);
    Dconc_z.resize(nlevs_max);
    bcs.resize(1);

    // periodic boundaries
    int bc_con_lo[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    int bc_con_hi[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};

/*
    // walls (Neumann)
    int bc_con_lo[] = {amrex::BCType::foextrap, amrex::BCType::foextrap, amrex::BCType::foextrap};
    int bc_con_hi[] = {amrex::BCType::foextrap, amrex::BCType::foextrap, amrex::BCType::foextrap};
*/
//    bcs.resize(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        // lo-side BCs
        if (bc_con_lo[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_con_lo[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_con_lo[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs[0].setLo(idim, bc_con_lo[idim]);
        }
        else {
            amrex::Abort("Invalid bc_con_lo");
        }

        // hi-side BCSs
        if (bc_con_hi[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_con_hi[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_con_hi[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs[0].setHi(idim, bc_con_hi[idim]);
        }
        else {
            amrex::Abort("Invalid bc_con_hi");
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
 * a specified level called by advance in executable directory                                                          *
 *                                                                             *
 *******************************************************************************/

void AmrCoreAdv::EvolveChem(
        std::array<MultiFab, AMREX_SPACEDIM> & umac_pre,
        std::array<MultiFab, AMREX_SPACEDIM> & umac,
        const iMultiFab & iface_pre,
        const iMultiFab & iface,
        const iMultiFab & catalyst_pre,
        const iMultiFab & catalyst,
        const MultiFab & LevelSet_pre,
        const MultiFab & LevelSet,
        int lev, int nstep,
        Real dt_fluid, Real time, Real dc,
        const Vector<std::array<MultiFab, AMREX_SPACEDIM>> & face_coords, int corrector, Real source_strength)
{
    // if this is a predictor or corrector step
    Correct=corrector;
    // diffusion coefficent
    diffcoeff=dc;
    // " strength" of point sources
    strength= source_strength;
    // time step
    dt[lev] = dt_fluid;

    // initialize copies of velocities umac amd face center posiions face_coord, initalize surface gradient of concentration, note that _pre indicates that this quantity is from the predicter step. If this is the predicter step (Corrector=0)  then the predicter quanties are the same as the other quanity, ie uface =uface_pre
   uface.resize(max_level + 1);
   vface.resize(max_level + 1);

   uface_pre.resize(max_level + 1);
   vface_pre.resize(max_level + 1);
   xface.resize(max_level + 1);
   yface.resize(max_level + 1);
    DistributionMapping condm = con_new[lev]->DistributionMap();
    BoxArray conba            = con_new[lev]->boxArray();

   BoxArray x_face_ba = umac[0].boxArray();//conba;
   BoxArray y_face_ba = umac[1].boxArray(); //conba;

#if (AMREX_SPACEDIM>=3)

   wface.resize(max_level + 1);

   wface_pre.resize(max_level + 1);

   zface.resize(max_level + 1);
#endif


    BoxArray z_face_ba = umac[2].boxArray();//conba;

    for (lev = 0; lev <= finest_level; ++lev) {

        uface[lev].reset(new MultiFab(x_face_ba, condm, 1, 1));
        vface[lev].reset(new MultiFab(y_face_ba, condm, 1, 1));
        uface_pre[lev].reset(new MultiFab(x_face_ba, condm, 1, 1));
        vface_pre[lev].reset(new MultiFab(y_face_ba, condm, 1, 1));
        int mac_ncompx= (face_coords[lev])[0].nComp();
        int mac_ncompy= (face_coords[lev])[1].nComp();
        int mac_ngrowx= (face_coords[lev])[0].nGrow();
        int mac_ngrowy= (face_coords[lev])[1].nGrow();
        xface[lev].reset(new MultiFab(x_face_ba, condm, mac_ncompx, 0));
        yface[lev].reset(new MultiFab(y_face_ba, condm, mac_ncompy, 0));

         xface[lev]->setVal(0.);
        yface[lev]->setVal(0.);

        uface[lev]->setVal(0.);
        vface[lev]->setVal(0.);
        uface_pre[lev]->setVal(0.);
        vface_pre[lev]->setVal(0.);
        Dcon_x[lev].reset(new MultiFab(x_face_ba, condm, 1, 1));
        Dcon_y[lev].reset(new MultiFab(y_face_ba, condm, 1, 1));

        Dcon_x[lev]->setVal(0.);
        Dcon_y[lev]->setVal(0.);
        Dconc_x[lev].reset(new MultiFab(conba, condm, 1, 1));
        Dconc_y[lev].reset(new MultiFab(conba, condm, 1, 1));

        MagDcon[lev].reset(new MultiFab(conba, condm, 1, 0));

        Dconc_x[lev]->setVal(0.);
        Dconc_y[lev]->setVal(0.);
        uface[lev]->copy(umac[0], 0, 0, 1, 0, 1);
        vface[lev]->copy(umac[1], 0, 0, 1, 0, 1);
        uface_pre[lev]->copy(umac_pre[0], 0, 0, 1, 0, 1);
        vface_pre[lev]->copy(umac_pre[1], 0, 0, 1, 0, 1);
        xface[lev]->copy((face_coords[lev])[0], 0, 0, mac_ncompx , 0, 0);
        yface[lev]->copy((face_coords[lev])[1], 0, 0, mac_ncompy , 0, 0);
       uface_pre[lev]->FillBoundary(geom[lev].periodicity());
       vface_pre[lev]->FillBoundary(geom[lev].periodicity());
       uface[lev]->FillBoundary(geom[lev].periodicity());
       vface[lev]->FillBoundary(geom[lev].periodicity());

#if (AMREX_SPACEDIM>=3)
        wface[lev].reset(new MultiFab(z_face_ba, condm, 1, 1));

        wface_pre[lev].reset(new MultiFab(z_face_ba, condm, 1, 1));

        int mac_ncompz= (face_coords[lev])[2].nComp();
        int mac_ngrowz= (face_coords[lev])[2].nGrow();


        zface[lev].reset(new MultiFab(z_face_ba, condm, mac_ncompz, 0));

        wface[lev]->setVal(0.);

        wface_pre[lev]->setVal(0.);

        zface[lev]->setVal(0.);

        Dcon_z[lev].reset(new MultiFab(z_face_ba, condm, 1, 1));

        Dcon_z[lev]->setVal(0.);

        Dconc_z[lev].reset(new MultiFab(conba, condm, 1, 1));

        Dconc_z[lev]->setVal(0.);

        wface[lev]->copy(umac[2], 0, 0, 1, 0, 1);

        wface_pre[lev]->copy(umac_pre[2], 0, 0, 1, 0, 1);

        zface[lev]->copy((face_coords[lev])[2], 0, 0, mac_ncompz , 0, 0);

       wface_pre[lev]->FillBoundary(geom[lev].periodicity());

       wface[lev]->FillBoundary(geom[lev].periodicity());

#endif
    }
    // initialize copies of levelset, and the location of the interface and catalyst

    int ls_gst= LevelSet.nGrow();
    int ls_nc= LevelSet.nComp();

    interface_loc.reset(new iMultiFab(conba, condm, 1, 1));
    interface_loc->copy(iface, 0, 0, 1, 0, 1 );
    interface_loc->FillBoundary(geom[0].periodicity());

    interface_loc_pre.reset(new iMultiFab(conba, condm, 1, 1));
    interface_loc_pre->copy(iface_pre, 0, 0, 1, 0, 1 );
    interface_loc_pre->FillBoundary(geom[0].periodicity());

    source_loc.reset(new iMultiFab(conba, condm, 1, 1));
    source_loc->copy(catalyst, 0, 0, 1, 0, 1 );
    source_loc->FillBoundary(geom[0].periodicity());

    source_loc_pre.reset(new iMultiFab(conba, condm, 1, 1));
    source_loc_pre->copy(catalyst_pre, 0, 0, 1, 0, 1 );
    source_loc_pre->FillBoundary(geom[0].periodicity());

    DistributionMapping lsdm = LevelSet.DistributionMap();
    BoxArray lsba            = LevelSet.boxArray();
    levset.reset(new MultiFab(lsba, lsdm, 1, 1));
    levset->setVal(0.);
    levset->copy(LevelSet, 0, ls_nc-1, ls_nc, ls_gst, 1);
    levset->FillBoundary(geom[0].periodicity());

    levset_pre.reset(new MultiFab(lsba, lsdm, 1, 1));
    levset_pre->setVal(0.);
    levset_pre->copy(LevelSet, 0, ls_nc-1, ls_nc, ls_gst, 1);
    levset_pre->FillBoundary(geom[0].periodicity());

    /***************************************************************************
     * Evolve chemical field by integrating time step                          *
     ***************************************************************************/
    lev = 0;
    int iteration = 1;
    timeStep( lev, time, iteration);



#ifdef BL_MEM_PROFILING
    {
        std::ostringstream ss;
        ss << "[STEP " << step+1 << "]";
        MemProfiler::report(ss.str());
    }
#endif
}


/*******************************************************************************
 *                                                                             *
 * Initialize multilevel data                                                  *
 *                                                                             *
 *******************************************************************************/

void AmrCoreAdv::InitData ( BoxArray & ba, DistributionMapping & dm)
{
    Initialize();

    if (restart_chkfile == "") {

        // start simulation from the beginning
        const Real time = 0.0;

        // initialize Levels (done in AmrCore)
        finest_level = 0;

        const BoxArray& ba = MakeBaseGrids();
        DistributionMapping dm(ba);

        MakeNewLevelFromScratch(0, time, ba, dm);

        SetBoxArray(0, ba);
        SetDistributionMap(0, dm);

  //      InitFromScratch(time);
        //AverageDown();

        // DEBUG: write intial checkpoint
        // if (chk_int > 0) WriteCheckpointFile();

    }
    else {
        // restart from a checkpoint
        ReadCheckpointFile();
    }

   // initialize grad con
   for (int lev = 0; lev <= finest_level; ++lev) {
       SetBoxArray(lev, ba);
       SetDistributionMap(lev, dm);

       con_pre[lev]->setVal(0.);
       con_new[lev]->setVal(0.);
       con_old[lev]->setVal(0.);
       // fills in concentration
       MakeNewLevelFromScratch ( lev, 0., ba, dm);
       BoxArray x_face_ba=ba;
       BoxArray y_face_ba=ba;
       BoxArray z_face_ba=ba;
        x_face_ba.surroundingNodes(0);
        y_face_ba.surroundingNodes(1);
        z_face_ba.surroundingNodes(2);

       Dcon_x[lev].reset(new MultiFab(x_face_ba, dm, 1, 1));
       Dcon_y[lev].reset(new MultiFab(y_face_ba, dm, 1, 1));

       Dconc_x[lev].reset(new MultiFab(ba, dm, 1, 1));
       Dconc_y[lev].reset(new MultiFab(ba, dm, 1, 1));

       Dcon_x[lev]->setVal(0.);
       Dcon_y[lev]->setVal(0.);

       Dconc_x[lev]->setVal(0.);
       Dconc_y[lev]->setVal(0.);


       Dcon_z[lev].reset(new MultiFab(z_face_ba, dm, 1, 1));
       Dconc_z[lev].reset(new MultiFab(ba, dm, 1, 1));
       Dcon_z[lev]->setVal(0.);

       Dconc_z[lev]->setVal(0.);
       MagDcon[lev].reset(new MultiFab(ba, dm, 1, 0));


       MagDcon[lev]->setVal(0.);


   }
}


// Make a new level using provided BoxArray and DistributionMapping and fill
// with interpolated coarse level data. Note: overrides the pure virtual
// function in AmrCore
void AmrCoreAdv::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray & ba,
                                         const DistributionMapping & dm)
{

    const int ncomp  = con_new[lev-1]->nComp();
    const int nghost = con_new[lev-1]->nGrow();

    con_new[lev]->define(ba, dm, ncomp, nghost);
    con_old[lev]->define(ba, dm, ncomp, nghost);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }

    FillCoarsePatch(lev, time, *con_new[lev], 0, ncomp);
}


// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data. Note: overrides the pure virtual
// function in AmrCore
void AmrCoreAdv::RemakeLevel (int lev, Real time, const BoxArray & ba,
                              const DistributionMapping & dm)
{

    const int ncomp = con_new[lev]->nComp();
    const int nghost =con_new[lev]->nGrow();

    MultiFab new_state(ba, dm, ncomp, nghost);
    MultiFab old_state(ba, dm, ncomp, nghost);

    FillPatch(lev, time, new_state, 0, ncomp);

    std::swap(new_state, *con_new[lev]);
    std::swap(old_state, *con_old[lev]);

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

    con_new[lev]->clear();
    con_old[lev]->clear();
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

    con_new[lev].reset(new MultiFab(ba, dm, ncomp, nghost));
    con_pre[lev].reset(new MultiFab(ba, dm, ncomp, nghost));
    con_old[lev].reset(new MultiFab(ba, dm, ncomp, nghost));

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

//    if (lev > 0 && do_reflux) {
//        flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
//    }

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    Real cur_time = t_new[lev];

    MultiFab& state = *con_new[lev];

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
    static Vector<Real> conerr;

    // only do this during the first call to ErrorEst
    if (first)
    {
        first = false;
        // read in an array of "conerr", which is the tagging threshold
        // in this example, we tag values of "con" which are greater than conerr
        // for that particular level
        // in subroutine state_error, you could use more elaborate tagging, such
        // as more advanced logical expressions, or gradients, etc.
        ParmParse pp("adv");
        int n = pp.countval("conerr");
        if (n > 0) {
            pp.getarr("conerr", conerr, 0, n);
        }
    }

    if (lev >= conerr.size()) return;

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    const MultiFab& state = *con_new[lev];

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
                        AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo), &time, &conerr[lev]);
            //
            // Now update the tags in the TagBox in the tilebox region
            // to be equal to itags
            //
            tagfab.tags_and_untags(itags, tilebox);
        }
    }
}

// set covered coarse cells to be the average of overlying fine cells
void AmrCoreAdv::AverageDown ()
{

    for (int lev = finest_level-1; lev >= 0; --lev)
        amrex::average_down(
                * con_new[lev+1], * con_new[lev],
                     geom[lev+1],   geom[lev],
                0, con_new[lev]->nComp(), refRatio(lev)
            );
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void AmrCoreAdv::AverageDownTo (int crse_lev)
{

    amrex::average_down(
            * con_new[crse_lev+1], * con_new[crse_lev],
                 geom[crse_lev+1],   geom[crse_lev],
            0, con_new[crse_lev]->nComp(), refRatio(crse_lev)
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

// compute a new multifab by coping in con from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
AmrCoreAdv::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{

    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;

        GetData(0, time, smf, stime);


        BndryFuncArray bfunc(confill);

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

        BndryFuncArray bfunc(confill);
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
    BndryFuncArray bfunc(confill);
    PhysBCFunct<BndryFuncArray> cphysbc(geom[lev-1],bcs,bfunc);
    PhysBCFunct<BndryFuncArray> fphysbc(geom[lev  ],bcs,bfunc);

//    PhysBCFunct cphysbc(geom[lev-1],bcs,BndryFunctBase(confill));
//    PhysBCFunct fphysbc(geom[lev  ],bcs,BndryFunctBase(confill));

    Interpolater* mapper = &cell_cons_interp;

    amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
                                 cphysbc,0, fphysbc, 0, refRatio(lev-1),
                                 mapper, bcs,0);
}

// utility to copy in data from con_old and/or con_new into another multifab
void AmrCoreAdv::GetData (int lev, Real time, Vector<MultiFab *> & data, Vector<Real> & datatime) {

    data.clear();
    datatime.clear();

    MultiFab* con_new_stdptr=con_new[lev].get();
    MultiFab* con_old_stdptr=con_old[lev].get();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps) {

        data.push_back(con_new_stdptr);
        datatime.push_back(t_new[lev]);

    } else if (time > t_old[lev] - teps && time < t_old[lev] + teps) {

        data.push_back(con_old_stdptr);
        datatime.push_back(t_old[lev]);

    } else {

        data.push_back(con_old_stdptr);
        data.push_back(con_new_stdptr);
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}



// advance a level by dt
// includes a recursive call for finer levels
void AmrCoreAdv::timeStep (int lev, Real time, int iteration)
{
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

//        if (do_reflux)
//        {
//         // update lev based on coarse-fine flux mismatch
//            flux_reg[lev+1]->Reflux(*con_new[lev], 1.0, 0, 0, con_new[lev]->nComp(), geom[lev]);
//        }

//        AverageDownTo(lev); // average lev+1 down to lev
    }

}


// advance a single level for a single time step, updates flux registers
void AmrCoreAdv::Advance (int lev, Real time, Real dt_lev, int iteration, int ncycle) {

    // need 3 ghost cell to advance concentration
    constexpr int num_grow =3;
    int num_comp=1;
    t_old[lev]  = t_new[lev];
    t_new[lev] += dt_lev;

    const BoxArray & badp            = con_new[lev]->boxArray();
    const DistributionMapping & dmdp = con_new[lev]->DistributionMap();

    // initalize point source multifab
    MultiFab ptSource(badp,dmdp,1,0);
    ptSource.setVal(0.);
    MultiFab ptSource_pre(badp,dmdp,1,0);

    ptSource_pre.setVal(0.);

    // if we are in the correcting step the updated concetration will be con_new, otherwise it is con_pre, the predicted concentration
    MultiFab * S_new = NULL;
    if (Correct == 1){
        S_new     = con_new[lev].get();}
    else{
        S_new     =con_pre[lev].get();}

    // define fluid velocity multifab for this level
    MultiFab &  uface_lev_pre = * uface_pre[lev];
    MultiFab &  vface_lev_pre = * vface_pre[lev];

    MultiFab &  uface_lev = * uface[lev];
    MultiFab &  vface_lev = * vface[lev];

#if (AMREX_SPACEDIM>=3)

    MultiFab &  wface_lev_pre = * wface_pre[lev];
    MultiFab &  wface_lev = * wface[lev];
#endif
    // define location of interface for this level
    iMultiFab & iloc_mf   = * interface_loc;
    iMultiFab & iloc_mf_pre   = * interface_loc_pre;
    // define location of sources for this level
    iMultiFab & sloc_mf   = * source_loc;
    iMultiFab & sloc_mf_pre   = * source_loc_pre;
    int Num_loc=sloc_mf.sum(0,false);
    int Num_loc_Pre=sloc_mf_pre.sum(0,false);
    // problem set up
    const Real * dx      = geom[lev].CellSize();
    const Real * prob_lo = geom[lev].ProbLo();
    // initalize fluxes multifab, reflux is old code for AMR
    MultiFab fluxes[BL_SPACEDIM];
    if (do_reflux) {
        for (int i = 0; i < BL_SPACEDIM; ++i) {
            BoxArray ba = grids[lev];
            ba.surroundingNodes(i);
            fluxes[i].define(ba, dmap[lev], num_comp, 0);
        }
    }

    // Copy of con_old used to update con_pre or con_new
    MultiFab Sborder(badp,dmdp, num_comp, num_grow);
    Sborder.setVal(0.);
    Sborder.copy(*con_old[lev], 0,num_comp-1 , num_comp, 0, num_grow);
    Sborder.FillBoundary(geom[0].periodicity());

    MultiFab Sborder_pre(badp,dmdp, num_comp, num_grow);
    Sborder_pre.setVal(0.);

    // Copy of con_old/con_pre used to update con_pre or con_new

    if (Correct==0){
        Sborder_pre.copy(*con_old[lev], 0,(*con_old[lev]).nComp()-1 ,(*con_old[lev]).nComp(), 0, num_grow);
    }
    else{
        Sborder_pre.copy(*con_pre[lev], 0,(*con_pre[lev]).nComp()-1 ,(*con_pre[lev]).nComp(), 0, num_grow);
    }
    Sborder_pre.FillBoundary(geom[0].periodicity());
#ifdef _OPENMP
#pragma omp parallel
#endif
    {

        FArrayBox flux[BL_SPACEDIM];
        FArrayBox flux1[BL_SPACEDIM];
        FArrayBox flux2[BL_SPACEDIM];

        FArrayBox vel[BL_SPACEDIM];

        for (MFIter mfi(*S_new, true); mfi.isValid(); ++mfi) {
            const Box & bx = mfi.tilebox();

            const FArrayBox & statein_p =   Sborder_pre[mfi];
            const FArrayBox & statein =   Sborder[mfi];
            FArrayBox & stateout      =     (*S_new)[mfi];
            FArrayBox & ptS           =  ptSource[mfi];
            IArrayBox & fabsl         =   sloc_mf[mfi];
            IArrayBox & fabil         =   iloc_mf[mfi];
            FArrayBox & uface_mf      = uface_lev[mfi];
            FArrayBox & vface_mf      = vface_lev[mfi];

            FArrayBox & ptS_p           =  ptSource_pre[mfi];
            IArrayBox & fabsl_p         =   sloc_mf_pre[mfi];
            IArrayBox & fabil_p         =   iloc_mf_pre[mfi];
            FArrayBox & uface_mf_p      = uface_lev_pre[mfi];
            FArrayBox & vface_mf_p      = vface_lev_pre[mfi];

#if (AMREX_SPACEDIM>=3)
            FArrayBox & wface_mf      = wface_lev[mfi];
            FArrayBox & wface_mf_p      = wface_lev_pre[mfi];
#endif


            for (int i = 0; i < BL_SPACEDIM ; i++) {
                const Box& bxtmp = amrex::surroundingNodes(bx,i);
                flux[i].resize(bxtmp,S_new->nComp());
                flux1[i].resize(bxtmp,S_new->nComp());
                flux2[i].resize(bxtmp,S_new->nComp());
            }
#if (AMREX_SPACEDIM==2)
            // define previous and predicted pointsources
            get_ptsource_2d( bx.loVect(), bx.hiVect(),
                             BL_TO_FORTRAN_3D(fabsl),
                             BL_TO_FORTRAN_3D(ptS),
                             & strength, dx,
                             AMREX_ZFILL(prob_lo), & Num_loc);
            get_ptsource_2d( bx.loVect(), bx.hiVect(),
                             BL_TO_FORTRAN_3D(fabsl_p),
                             BL_TO_FORTRAN_3D(ptS_p),
                             & strength, dx,
                              AMREX_ZFILL(prob_lo),& Num_loc_Pre);


            // compute new state (stateout) and fluxes.
            advect_2d(& time, bx.loVect(), bx.hiVect(),
                      BL_TO_FORTRAN_3D(statein),
                      BL_TO_FORTRAN_3D(statein_p),
                      BL_TO_FORTRAN_3D(stateout),
                      BL_TO_FORTRAN_3D(ptS),
                      BL_TO_FORTRAN_3D(ptS_p),
                      BL_TO_FORTRAN_3D(fabil),
                      BL_TO_FORTRAN_3D(fabil_p),
                      AMREX_D_DECL(BL_TO_FORTRAN_3D(uface_mf),
                                   BL_TO_FORTRAN_3D(vface_mf),
                      AMREX_D_DECL(BL_TO_FORTRAN_3D(uface_mf_p),
                                   BL_TO_FORTRAN_3D(vface_mf_p),
                      AMREX_D_DECL(BL_TO_FORTRAN_3D(flux[0]),
                                   BL_TO_FORTRAN_3D(flux[1]),
                      AMREX_D_DECL(BL_TO_FORTRAN_3D(flux1[0]),
                                   BL_TO_FORTRAN_3D(flux1[1]),
                      AMREX_D_DECL(BL_TO_FORTRAN_3D(flux2[0]),
                                   BL_TO_FORTRAN_3D(flux2[1]),
                      dx, & dt_lev, & diffcoeff, & Correct);



#elif (AMREX_SPACEDIM>=3)

            // define previous and predicted pointsources

            get_ptsource_3d( bx.loVect(), bx.hiVect(),
                             BL_TO_FORTRAN_3D(fabsl),
                             BL_TO_FORTRAN_3D(ptS),
                             & strength, dx,
                             AMREX_ZFILL(prob_lo), & Num_loc);
            get_ptsource_3d( bx.loVect(), bx.hiVect(),
                             BL_TO_FORTRAN_3D(fabsl_p),
                             BL_TO_FORTRAN_3D(ptS_p),
                             & strength, dx,
                              AMREX_ZFILL(prob_lo), & Num_loc_Pre);


            // compute new state (stateout) and fluxes.
            advect_3d(& time, bx.loVect(), bx.hiVect(),
                      BL_TO_FORTRAN_3D(statein),
                      BL_TO_FORTRAN_3D(statein_p),
                      BL_TO_FORTRAN_3D(stateout),
                      BL_TO_FORTRAN_3D(ptS),
                      BL_TO_FORTRAN_3D(ptS_p),
                      BL_TO_FORTRAN_3D(fabil),
                      BL_TO_FORTRAN_3D(fabil_p),
                      AMREX_D_DECL(BL_TO_FORTRAN_3D(uface_mf),
                                   BL_TO_FORTRAN_3D(vface_mf),
                                   BL_TO_FORTRAN_3D(wface_mf)),
                      AMREX_D_DECL(BL_TO_FORTRAN_3D(uface_mf_p),
                                   BL_TO_FORTRAN_3D(vface_mf_p),
                                   BL_TO_FORTRAN_3D(wface_mf_p)),
                      AMREX_D_DECL(BL_TO_FORTRAN_3D(flux[0]),
                                   BL_TO_FORTRAN_3D(flux[1]),
                                   BL_TO_FORTRAN_3D(flux[2])),
                      AMREX_D_DECL(BL_TO_FORTRAN_3D(flux1[0]),
                                   BL_TO_FORTRAN_3D(flux1[1]),
                                   BL_TO_FORTRAN_3D(flux1[2])),
                      AMREX_D_DECL(BL_TO_FORTRAN_3D(flux2[0]),
                                   BL_TO_FORTRAN_3D(flux2[1]),
                                   BL_TO_FORTRAN_3D(flux2[2])),
                      dx, & dt_lev, & diffcoeff, & Correct);
#endif
            // old code to do AMR
            if (do_reflux) {
                for (int i = 0; i < BL_SPACEDIM ; i++) {
                    fluxes[i][mfi].copy(flux[i],mfi.nodaltilebox(i));
                }
            }
        }
    }
    S_new->FillBoundary(geom[lev].periodicity());
    con_pre[lev]->FillBoundary(geom[lev].periodicity());
    con_new[lev]->FillBoundary(geom[lev].periodicity());

    // if we are in a correcting step we update the concentration
    if (Correct==1)
        std::swap(con_old[lev], con_new[lev]);

    con_old[lev]->FillBoundary(geom[lev].periodicity());

    // After updating con_new we compute the first derivatives

    // tangential components of the derivatives, face centered
    MultiFab &  txf_mf       = * Dcon_x[lev];
    MultiFab &  tyf_mf       = * Dcon_y[lev];
    MultiFab &  txc_mf       = * Dconc_x[lev];
    MultiFab &  tyc_mf       = * Dconc_y[lev];
    const BoxArray & badpx           = Dcon_x[lev]->boxArray();
    const DistributionMapping & dmdpx =Dcon_x[lev]->DistributionMap();

    const BoxArray & badpy           = Dcon_y[lev]->boxArray();
    const DistributionMapping & dmdpy =Dcon_y[lev]->DistributionMap();

    txf_mf.setVal(0.);
    tyf_mf.setVal(0.);
    txc_mf.setVal(0.);
    tyc_mf.setVal(0.);

#if (AMREX_SPACEDIM>=3)
    MultiFab &  tzf_mf       = * Dcon_z[lev];

    MultiFab &  tzc_mf       = * Dconc_z[lev];
    const BoxArray & badpz           = Dcon_z[lev]->boxArray();
    const DistributionMapping & dmdpz =Dcon_z[lev]->DistributionMap();
    tzf_mf.setVal(0.);
    tzc_mf.setVal(0.);

#endif

    // level set
    MultiFab &  ls_mf       = * levset;
    // magnitude of the surface gradient
    MultiFab &  mdc_mf       = * MagDcon[lev];
    MultiFab S_new_fill(badp,dmdp,num_comp, 1);
    S_new_fill.setVal(0.);

    // if we are correcting we take the surface gradient of the correct concentration (recalling that above we swapped so it will be stored in con_old), otherwise we take the surface gradient of the predicted concentration
    if (Correct==1){
        S_new_fill.copy(*con_old[lev], 0,num_comp-1 ,num_comp,(*con_new[lev]).nGrow() , 1);}
    else{
        S_new_fill.copy(*con_pre[lev], 0,num_comp-1 ,num_comp,(*con_pre[lev]).nGrow() , 1);}
    S_new_fill.FillBoundary(geom[lev].periodicity());


#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter mfi(S_new_fill, true); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();

            FArrayBox & stateout      =   S_new_fill[mfi];
            IArrayBox & fabil         =      iloc_mf[mfi];
            FArrayBox & fabmdc        =      mdc_mf[mfi];
            FArrayBox & fabsls        =      ls_mf[mfi];

            FArrayBox & fabx         =       txc_mf[mfi];
            FArrayBox & faby         =       tyc_mf[mfi];
#if (AMREX_SPACEDIM>=3)
            FArrayBox & fabz         =       tzc_mf[mfi];
#endif
            // No support for 2D yet
#if (AMREX_SPACEDIM>=3)
            // compute surface gradient
            get_surfgrad_3d( bx.loVect(), bx.hiVect(),
                            BL_TO_FORTRAN_3D(stateout),
                            BL_TO_FORTRAN_3D(fabx),
                            BL_TO_FORTRAN_3D(faby),
                            BL_TO_FORTRAN_3D(fabz),
                            BL_TO_FORTRAN_3D(fabmdc),
                            BL_TO_FORTRAN_3D(fabsls),
                            BL_TO_FORTRAN_3D(fabil),
                            dx, AMREX_ZFILL(prob_lo));
#endif
        }
    }

    txc_mf.FillBoundary(geom[lev].periodicity());
    tyc_mf.FillBoundary(geom[lev].periodicity());

    // initialize array of face centered multifabs
    std::array< MultiFab, AMREX_SPACEDIM > Cxface_array;
    Cxface_array[0].define(badpx, dmdpx, 1, 0);
    Cxface_array[1].define(badpy, dmdpy, 1, 0);
    std::array< MultiFab, AMREX_SPACEDIM > Cyface_array;
    Cyface_array[0].define(badpx, dmdpx, 1, 0);
    Cyface_array[1].define(badpy, dmdpy, 1, 0);
    Cxface_array[0].setVal(0.);
    Cxface_array[1].setVal(0.);
    Cyface_array[0].setVal(0.);
    Cyface_array[1].setVal(0.);

#if (AMREX_SPACEDIM>=3)
    Cxface_array[2].define(badpz, dmdpz, 1, 0);
    tzc_mf.FillBoundary(geom[lev].periodicity());

    Cyface_array[2].define(badpz, dmdpz, 1, 0);

    std::array< MultiFab, AMREX_SPACEDIM > Czface_array;
    Czface_array[0].define(badpx, dmdpx, 1, 0);
    Czface_array[1].define(badpy, dmdpy, 1, 0);
    Czface_array[2].define(badpz, dmdpz, 1, 0);

    Cxface_array[2].setVal(0.);

    Cyface_array[2].setVal(0.);

    Czface_array[0].setVal(0.);
    Czface_array[1].setVal(0.);
    Czface_array[2].setVal(0.);
#endif

    // Average cell centered surface gradient to face centers
    AverageCCToFace(txc_mf,Cxface_array,0,1,SPEC_BC_COMP,geom);
    AverageCCToFace(tyc_mf,Cyface_array,0,1,SPEC_BC_COMP,geom);

    txf_mf.copy(Cxface_array[0], 0, 0,1, 0, 0);
    tyf_mf.copy(Cyface_array[1], 0, 0,1, 0, 0);
    txf_mf.FillBoundary(geom[lev].periodicity());
    tyf_mf.FillBoundary(geom[lev].periodicity());

#if (AMREX_SPACEDIM>=3)
    AverageCCToFace(tzc_mf,Czface_array,0,1,SPEC_BC_COMP,geom);
    tzf_mf.copy(Czface_array[2], 0, 0,1, 0, 0);
    tzf_mf.FillBoundary(geom[lev].periodicity());
#endif
    if( Correct==1){
        // Print out the total concentration in simulated domain vs the true total concentration
        amrex::Real SA=2*3.14*0.1*0.1;
        amrex::Print() << "time = "<< time+dt[0]<< " simulated con total grid "<< (con_old[lev]->sum(0,false)*(*dx)*(*dx)*(*dx));
        amrex::Print() << "simulated con surface "<< (ptSource.sum(0,false))*(time+dt[0])*(*dx)*(*dx)/4;
        amrex::Print() << "true con total"<< SA*strength/Num_loc*(time+dt[0])<< std::endl;
    }

}

void AmrCoreAdv::con_new_copy(int  lev, amrex::Vector<std::unique_ptr<MultiFab>> & MF, int indicator) {
    // indicator gives which quantity is being copied into MF
   // mf[lev].reset(new MultiFab(grids[lev], dmap[lev], ncomp, ngrow));

    DistributionMapping condm = con_new[lev]->DistributionMap();
    BoxArray conba            = con_new[lev]->boxArray();
   // BoxArray x_face_ba = conba;
   // BoxArray y_face_ba = conba;
   // x_face_ba.surroundingNodes(0);
   // y_face_ba.surroundingNodes(1);

#if (AMREX_SPACEDIM>=3)

   // BoxArray z_face_ba = conba;

   // z_face_ba.surroundingNodes(2);
#endif

    if (indicator==0){
        MF[lev].reset(new MultiFab(conba, condm, 1, 0));

        MF[lev]->setVal(0.);

        MF[lev]->copy(* con_old[lev], 0, 0,1, 0, 0);
    }
    // face centered surface gradients
    else if (indicator==1){
        int xng=Dcon_x[lev]->nGrow();
        BoxArray x_face_ba= Dcon_x[lev]->boxArray();
        MF[lev].reset(new MultiFab(x_face_ba, condm, 1, 0));

        MF[lev]->setVal(0.);

        MF[lev]->copy(* Dcon_x[lev], 0, 0,1, xng, 0);
    }
    else if (indicator==2){
        int yng=Dcon_y[lev]->nGrow();
        BoxArray y_face_ba= Dcon_y[lev]->boxArray();


        MF[lev].reset(new MultiFab(y_face_ba, condm, 1, 0));

        MF[lev]->setVal(0.);

        MF[lev]->copy(* Dcon_y[lev], 0, 0,1, yng, 0);
    }
    else if (indicator==3){
#if (AMREX_SPACEDIM>=3)

        int zng=Dcon_z[lev]->nGrow();
        BoxArray z_face_ba= Dcon_z[lev]->boxArray();

        MF[lev].reset(new MultiFab(z_face_ba, condm, 1, 0));

        MF[lev]->setVal(0.);
        MF[lev]->copy(* Dcon_z[lev], 0, 0,1, zng, 0);
#endif
    }
    else if (indicator==4){
        MF[lev].reset(new MultiFab(conba, condm, 1, 0));

        MF[lev]->setVal(0.);

        MF[lev]->copy(* MagDcon[lev], 0, 0,1, 0, 0);
    }
    else if (indicator==5){
        DistributionMapping condm = Dconc_x[lev]->DistributionMap();
        BoxArray conba            = Dconc_x[lev]->boxArray();
        int xng=Dconc_x[lev]->nGrow();

        MF[lev].reset(new MultiFab(conba, condm, 1, 0));

        MF[lev]->setVal(0.);

        MF[lev]->copy(* Dconc_x[lev], 0, 0,1, xng, 0);
    }
    // cell centered surface gradient
    else if (indicator==6){
        DistributionMapping condm = Dconc_y[lev]->DistributionMap();
        BoxArray conba            = Dconc_y[lev]->boxArray();
        int yng=Dconc_y[lev]->nGrow();


        MF[lev].reset(new MultiFab(conba, condm, 1, 0));

        MF[lev]->setVal(0.);

        MF[lev]->copy(* Dconc_y[lev], 0, 0,1, yng, 0);
    }
    else if (indicator==7){
#if (AMREX_SPACEDIM>=3)

        DistributionMapping condm = Dconc_z[lev]->DistributionMap();
        BoxArray conba            = Dconc_z[lev]->boxArray();
        int zng=Dconc_z[lev]->nGrow();

        MF[lev].reset(new MultiFab(conba, condm, 1, 0));

        MF[lev]->setVal(0.);
        MF[lev]->copy(* Dconc_z[lev], 0, 0,1, zng, 0);
#endif
    }

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

        con_old[lev]->define(grids[lev], dmap[lev], ncomp, nghost);

        con_new[lev]->define(grids[lev], dmap[lev], ncomp, nghost);

        if (lev > 0 && do_reflux) {
            flux_reg[lev].reset(new FluxRegister(grids[lev], dmap[lev], refRatio(lev-1), lev, ncomp));
        }
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {

        VisMF::Read(*con_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "con"));
    }

}

// utility to skip to next line in Header
void
AmrCoreAdv::GotoNextLine (std::istream& is)
{

    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}
