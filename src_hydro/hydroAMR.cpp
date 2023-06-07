#include "gmres_functions.H"
#include "common_functions.H"
#include "hydroAMR.H"
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

hydroAMR::hydroAMR(int ang, int * is_periodic, Real dt_in) : AmrCore() {

//    IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
//    IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
//    IntVect max_grid(AMREX_D_DECL(common::max_grid_size[0], common::max_grid_size[1], common::max_grid_size[2]));


    //domain.resize(nlevels);
    //ba.resize(nlevels);
    //dmap.resize(nlevels);
    pres.resize(nlevels);
    plot_mf.resize(nlevels);
    geom.resize(nlevels);
    alpha_fc.resize(nlevels);
    source.resize(nlevels);
    sourceTemp.resize(nlevels);
    sourceRFD.resize(nlevels);
    beta_ed.resize(nlevels);
    eta_ed.resize(nlevels);
    temp_ed.resize(nlevels);
    gamma.resize(nlevels);
    beta.resize(nlevels);
    eta_cc.resize(nlevels);
    temp_cc.resize(nlevels);
    
    cc_mask.resize(nlevels);
    fc_mask.resize(nlevels);

    umac.resize(nlevels);
    umacNew.resize(nlevels);    
    umacM.resize(nlevels);
    umacV.resize(nlevels);    
    alpha_fc.resize(nlevels);
    source.resize(nlevels);
    sourceTemp.resize(nlevels);
    sourceRFD.resize(nlevels);
    stochMfluxdiv.resize(nlevels);
    beta_ed.resize(nlevels);
    eta_ed.resize(nlevels);
    temp_ed.resize(nlevels);
    
    gmres_rhs_p.resize(nlevels);
    gmres_rhs_u.resize(nlevels);    

    ng = ang;
    dt = dt_in;

    int bc_lo[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    int bc_hi[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};

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


//    domain[0].setSmall(dom_lo);
//    domain[0].setBig(dom_hi);
    
//    ba.define(domain[0]);
//    ba.maxSize(max_grid);
//    dm.define(ba);
    

}

void hydroAMR::initData ()
{
    const Real time = 0.0;
    InitFromScratch(time);

}

void hydroAMR::updateAlpha(int lev, Real dt)
{

    Real dtinv = 1.0/dt;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        alpha_fc[lev][d].setVal(dtinv);
    }

}

void hydroAMR::updateGrid ()
{
    const Real time = 0.0;
    regrid(0,time);
    
    BoxArray bac = boxArray(0);
    BoxArray baf = boxArray(1);

    cc_mask[0].define(bac,DistributionMap(0), 1, 0);
    cc_mask[0].setVal(1);
    for(int d=0; d<AMREX_SPACEDIM; ++d) {

        fc_mask[0][d].define(convert(bac,nodal_flag_dir[d]), DistributionMap(0), 1, 0);
        fc_mask[0][d].setVal(1);       
    }
    
    if(nlevels>1)
    {
        cc_mask[1].define(baf.coarsen(2),DistributionMap(0), 1, 0);
        cc_mask[1].setVal(0);
        cc_mask[0].ParallelCopy(cc_mask[1], 0, 0, 1, 0, 0);
        cc_mask[1].setVal(1);
        for(int d=0; d<AMREX_SPACEDIM; ++d) {
            fc_mask[1][d].define(convert(baf.coarsen(2),nodal_flag_dir[d]), DistributionMap(0), 1, 0);
            fc_mask[1][d].setVal(0);
            fc_mask[0][d].ParallelCopy(fc_mask[1][d], 0, 0, 1, 0, 0);
            fc_mask[1][d].setVal(1);
        }

        
    }
    
   
    
}

void
hydroAMR::ErrorEst (int lev, TagBoxArray& tags, Real /*time*/, int /*ngrow*/)
{
    if(init)
    {
        init = false;
        mask.define(boxArray(0),DistributionMap(0),1,0);
        mask.setVal(0);
    }else
    {
        for (MFIter mfi(mask,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx  = mfi.tilebox();
            const auto tagfab  = tags.array(mfi);
            const auto maskfab = mask.array(mfi);

            const int   tagval = TagBox::SET;

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if(maskfab(i,j,k) > 0)
                {
                    tagfab(i,j,k) = tagval;
                }
            });
        }
    }
    
}

void hydroAMR::addPatch(IntVect patch_lo, IntVect patch_hi) {

        for (MFIter mfi(mask,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx  = mfi.tilebox();
            const auto maskfab = mask.array(mfi);

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if(i >= patch_lo[0] && i <= patch_hi[0] && j >= patch_lo[1] && j <= patch_hi[1] && k >= patch_lo[2] && k <= patch_hi[2])
                {
                    maskfab(i,j,k)++;
                }
            });
        }

}

void hydroAMR::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
                                          const DistributionMapping& dm)
{
    Print() << "Init level " << lev << std::endl;
    Print() << "dm: " << dm << std::endl;
    Print() << "ba: " << ba << std::endl;
    Print() << "size: " << umac.size() << std::endl;
    
     for (int d=0; d<AMREX_SPACEDIM; ++d) {
            umac[lev][d].define(convert(ba,nodal_flag_dir[d]), dm, 1, ng);
            umacNew[lev][d].define(convert(ba,nodal_flag_dir[d]), dm, 1, ng);
            source[lev][d].define(convert(ba,nodal_flag_dir[d]), dm, 1, ng);
            sourceRFD[lev][d].define(convert(ba,nodal_flag_dir[d]), dm, 1, ng);
            stochMfluxdiv[lev][d].define(convert(ba,nodal_flag_dir[d]), dm, 1, ng);
            sourceTemp[lev][d].define(convert(ba,nodal_flag_dir[d]), dm, 1, ng);
            umacM[lev][d].define(convert(ba,nodal_flag_dir[d]), dm, 1, 1);            
            umacV[lev][d].define(convert(ba,nodal_flag_dir[d]), dm, 1, 1);            
            gmres_rhs_u[lev][d].define(convert(ba,nodal_flag_dir[d]), dm, 1, 1);

            stochMfluxdiv[lev][d].setVal(0.);
            umac[lev][d].setVal(0.);
            source[lev][d].setVal(0.);
            sourceRFD[lev][d].setVal(0.);
            sourceTemp[lev][d].setVal(0.);                        
            umacM[lev][d].setVal(0.);
            umacV[lev][d].setVal(0.);
            gmres_rhs_u[lev][d].setVal(0.);
            //fc_mask[d].setVal(0);
    }
    
    pres[lev].define(ba,dm,1,1);
    pres[lev].setVal(0.);
    
    gmres_rhs_p[lev].define(ba,dm,1,1);
    gmres_rhs_p[lev].setVal(0.);
    
    plot_mf[lev].define(ba,dm,nplot,1);
    
//    RealBox realDomain({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
//                       {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});
//                       
//    geom[0].define(domain[0] ,&realDomain,CoordSys::cartesian,is_periodic);
    
    Real dtinv = 1.0/dt;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        alpha_fc[lev][d].define(convert(ba,nodal_flag_dir[d]), dm, 1, 1);
        alpha_fc[lev][d].setVal(dtinv);
    }
    
    // beta cell centred
    beta[lev].define(ba, dm, 1, 1);
    beta[lev].setVal(visc_coef);
    // cell-centered gamma
    gamma[lev].define(ba, dm, 1, 1);
    gamma[lev].setVal(0.);



    eta_cc[lev].define(ba, dm, 1, 1);
    temp_cc[lev].define(ba, dm, 1, 1);
    eta_cc[lev].setVal(visc_coef);
    temp_cc[lev].setVal(T_init[0]);


#if (AMREX_SPACEDIM == 2)
    beta_ed[lev][0].define(convert(ba,nodal_flag), dm, 1, 1);
    //ed_mask[lev].define(convert(ba,nodal_flag), dm, 1, 1);
#elif (AMREX_SPACEDIM == 3)
    beta_ed[lev][0].define(convert(ba,nodal_flag_xy), dm, 1, 1);
    beta_ed[lev][1].define(convert(ba,nodal_flag_xz), dm, 1, 1);
    beta_ed[lev][2].define(convert(ba,nodal_flag_yz), dm, 1, 1);
    //ed_mask[lev].define(convert(ba,nodal_flag_xy), dm, 1, 1);
    //ed_mask[lev].define(convert(ba,nodal_flag_xz), dm, 1, 1);
    //ed_mask[lev].define(convert(ba,nodal_flag_yz), dm, 1, 1);
#endif
    
#if (AMREX_SPACEDIM == 2)
    eta_ed[lev][0].define(convert(ba,nodal_flag),    dm, 1, 0);
    temp_ed[lev][0].define(convert(ba,nodal_flag),    dm, 1, 0);
#elif (AMREX_SPACEDIM == 3)
    eta_ed[lev][0].define(convert(ba,nodal_flag_xy), dm, 1, 0);
    eta_ed[lev][1].define(convert(ba,nodal_flag_xz), dm, 1, 0);
    eta_ed[lev][2].define(convert(ba,nodal_flag_yz), dm, 1, 0);
    temp_ed[lev][0].define(convert(ba,nodal_flag_xy), dm, 1, 0);
    temp_ed[lev][1].define(convert(ba,nodal_flag_xz), dm, 1, 0);
    temp_ed[lev][2].define(convert(ba,nodal_flag_yz), dm, 1, 0);
#endif                 

    for (int d=0; d<NUM_EDGE; ++d) {
        eta_ed[lev][d] .setVal(visc_coef);
        temp_ed[lev][d].setVal(T_init[0]);
        beta_ed[lev][d].setVal(visc_coef);
        //ed_mask[d].setVal(0);
    }
}


void
hydroAMR::FillFine()
{
    FillCoarsePatchCell(pres);
    FillCoarsePatchFace(umac);
}

void
hydroAMR::FillCoarsePatchCell(Vector<MultiFab>& mf)
{
    //BL_ASSERT(lev > 0);

    //Interpolater* mapper = &cell_cons_interp;
    CellConservativeLinear mapper;
    

    Real time=0;
    //const char br;

    if(Gpu::inLaunchRegion())
    {
        GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > cphysbc(geom[0],bcs,gpu_bndry_func);
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > fphysbc(geom[1],bcs,gpu_bndry_func);
        
        //PhysBCFunctNoOp cphysbc, fphysbc;

        amrex::InterpFromCoarseLevel(mf[1], time, mf[0], 0, 0, 1, geom[0], geom[1],
                                     cphysbc, 0, fphysbc, 0, refRatio(0),
                                     &mapper, bcs, 0);
    }
    else
    {
        CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[0],bcs,bndry_func);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[1],bcs,bndry_func);

        amrex::InterpFromCoarseLevel(mf[1], time, mf[0], 0, 0, 1, geom[0], geom[1],
                                     cphysbc, 0, fphysbc, 0, refRatio(0),
                                     &mapper, bcs, 0);
    }
}

void
hydroAMR::FillCoarsePatchFace(Vector<std::array< MultiFab, AMREX_SPACEDIM >>& mf)
{
    //BL_ASSERT(lev > 0);

    //Interpolater* mapper = &cell_cons_interp;
    Interpolater* mapper = &face_divfree_interp;
    //FaceDivFree mapper;
    
    std::array< MultiFab*, AMREX_SPACEDIM > mfCoarse;
    std::array< MultiFab*, AMREX_SPACEDIM > mfFine;
    
    mfCoarse = GetArrOfPtrs(mf[0]);
    mfFine = GetArrOfPtrs(mf[1]);
    
    Real time=0;

    if(Gpu::inLaunchRegion())
    {
        GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > cphysbc(geom[0],bcs,gpu_bndry_func);
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > fphysbc(geom[1],bcs,gpu_bndry_func);

        Array<PhysBCFunct<GpuBndryFuncFab<AmrCoreFill>>,AMREX_SPACEDIM> cp = {AMREX_D_DECL(cphysbc,cphysbc,cphysbc)};
        Array<PhysBCFunct<GpuBndryFuncFab<AmrCoreFill>>,AMREX_SPACEDIM> fp = {AMREX_D_DECL(fphysbc,fphysbc,fphysbc)};
        
        Array<Vector<BCRec>,AMREX_SPACEDIM> bcr = {AMREX_D_DECL(bcs,bcs,bcs)};
        
        amrex::InterpFromCoarseLevel(mfFine, time, mfCoarse, 0, 0, 1, geom[0], geom[1],
                                     cp, 0, fp, 0, refRatio(0),
                                     mapper, bcr, 0);

    }
    else
    {
        CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[0],bcs,bndry_func);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[1],bcs,bndry_func);
        
        Array<PhysBCFunct<CpuBndryFuncFab>,AMREX_SPACEDIM> cp = {AMREX_D_DECL(cphysbc,cphysbc,cphysbc)};
        Array<PhysBCFunct<CpuBndryFuncFab>,AMREX_SPACEDIM> fp = {AMREX_D_DECL(fphysbc,fphysbc,fphysbc)};
        
        Array<Vector<BCRec>,AMREX_SPACEDIM> bcr = {AMREX_D_DECL(bcs,bcs,bcs)};
        
        amrex::InterpFromCoarseLevel(mfFine, time, mfCoarse, 0, 0, 1, geom[0], geom[1],
                                     cp, 0, fp, 0, refRatio(0),
                                     mapper, bcr, 0);

    }


}

void
hydroAMR::advanceStokes()
{
    BL_PROFILE_VAR("hydroAMR::advanceStokes()",hydroadvance);
          
    Real theta_alpha = 0.;
    Real norm_pre_rhs;
    
    for(int lev=0;lev<nlevels;++lev)
    {
        gmres_rhs_p[lev].setVal(0.);
        
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
        gmres_rhs_u[lev][d].setVal(0.);

        }
    }


    //////////////////////////////////////////////////
    for(int lev=0;lev<nlevels;++lev)
    { 
    // add stochastic forcing to gmres_rhs_u
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Add(gmres_rhs_u[lev][d], stochMfluxdiv[lev][d], 0, 0, 1, 0);
            MultiFab::Add(gmres_rhs_u[lev][d], source[lev][d], 0, 0, 1, 0);
        }
    }

//    if (zero_net_force == 1)
//    {
//        Vector<Real> mean_stress_umac(AMREX_SPACEDIM);

//        SumStag(gmres_rhs_u,mean_stress_umac,true);

//        Print() << "correcting mean force: " << mean_stress_umac[0]*(prob_hi[0]-prob_lo[0])*(prob_hi[1]-prob_lo[1])*(prob_hi[2]-prob_lo[2]);

//        for (int d=1; d<AMREX_SPACEDIM; ++d) {
//            Print() << ", " << mean_stress_umac[d]*(prob_hi[0]-prob_lo[0])*(prob_hi[1]-prob_lo[1])*(prob_hi[2]-prob_lo[2]);
//        }

//        Print() << "\n";

//        Print() << "test force: " << sourceTerms[0].sum()*(prob_hi[0]-prob_lo[0])*(prob_hi[1]-prob_lo[1])*(prob_hi[2]-prob_lo[2]);

//        for (int d=1; d<AMREX_SPACEDIM; ++d) {
//            Print() << ", " << sourceTerms[d].sum()*(prob_hi[0]-prob_lo[0])*(prob_hi[1]-prob_lo[1])*(prob_hi[2]-prob_lo[2]);
//        }

//        Print() << "\n";

//        for (int d=0; d<AMREX_SPACEDIM; ++d) {
//            if (geom.isPeriodic(d)) {
//                gmres_rhs_u[d].plus(-mean_stress_umac[d],0,1,0);
//            }
//        }
//    }
//    
    //PrintMF(alpha_fc[0],0,0);
      
    // call GMRES
    GMRES gmres(grids,dmap,geom,nlevels);
    gmres.Solve(gmres_rhs_u,gmres_rhs_p,umac,pres,
                alpha_fc,beta,beta_ed,gamma,cc_mask,fc_mask,theta_alpha,geom,norm_pre_rhs);


//    GMRES gmres(grids[0],dmap[0],geom[0]);
//    gmres.Solve(gmres_rhs_u[0],gmres_rhs_p[0],umac[0],pres[0],
//                alpha_fc[0],beta[0],beta_ed[0],gamma[0],theta_alpha,geom[0],norm_pre_rhs);
                


    for (int i=0; i<AMREX_SPACEDIM; i++) {
        MultiFabPhysBCDomainVel(umac[0][i], geom[0], i);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umac[0][i], geom[0], i, is_inhomogeneous);
        umac[0][i].FillBoundary(geom[0].periodicity());
    }
}

void
hydroAMR::WritePlotFile(Real time, int step)
{
    int num_output_comp = 1 + AMREX_SPACEDIM;
    int num_levels = pres.size();
    //IntVect cc_flag = IntVect::TheZeroVector();
    Vector<std::unique_ptr<MultiFab> > output_cc(num_levels);
    
    for (int lev = 0; lev < nlevels; ++lev) {
        const BoxArray& cc_ba = pres[lev].boxArray();
        output_cc[lev].reset(new MultiFab(cc_ba,pres[lev].DistributionMap(), num_output_comp, 0));
        MultiFab::Copy(*output_cc[lev], pres[lev], 0, 0, 1, 0);
        AverageFaceToCC(umac[lev][0],*output_cc[lev],0,1);
        AverageFaceToCC(umac[lev][1],*output_cc[lev],1,2);
        AverageFaceToCC(umac[lev][2],*output_cc[lev],2,3);        
    }

    Vector<std::string> varnames;
    varnames.push_back("pressure");
    varnames.push_back("velx");
    varnames.push_back("vely");
    varnames.push_back("velz");

    Vector<int> level_steps;
    level_steps.push_back(step);
    level_steps.push_back(step);

    int output_levs = num_levels;

    Vector<IntVect> outputRR(output_levs);
    for (int lev = 0; lev < output_levs; ++lev) {
        outputRR[lev] = IntVect(AMREX_D_DECL(2, 2, 2));
    }

    const std::string& pltfile = amrex::Concatenate("pltML", step, 9);
    WriteMultiLevelPlotfile(pltfile, output_levs, GetVecOfConstPtrs(output_cc),
                            varnames, geom, 0.0, level_steps, outputRR);                                   
}


void
hydroAMR::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
                                    const DistributionMapping& dm)
{
    MakeNewLevelFromScratch (lev, time, ba, dm);
}


void
hydroAMR::RemakeLevel (int lev, Real time, const BoxArray& ba,
                         const DistributionMapping& dm)
{
    MakeNewLevelFromScratch (lev, time, ba, dm);
}

void
hydroAMR::ClearLevel (int lev)
{
       Print() << "CL\n";
}





