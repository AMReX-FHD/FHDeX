#include "common_functions.H"


// Divergence free coarsening of vel, assumes two levels with ratio 2

void FaceFillCoarse(Vector<std::array< MultiFab, AMREX_SPACEDIM >>& mf, int map)
{
    std::array< MultiFab*, AMREX_SPACEDIM > mfCoarse;
    std::array< MultiFab*, AMREX_SPACEDIM > mfFine;
    
    mfCoarse = GetArrOfPtrs(mf[0]);
    mfFine = GetArrOfPtrs(mf[1]);
    
    FaceFillCoarse(mfCoarse,mfFine,map);
}

void FaceFillCoarse(std::array<MultiFab, AMREX_SPACEDIM>* & mf, int map)
{

    std::array< MultiFab*, AMREX_SPACEDIM > mfCoarse;
    std::array< MultiFab*, AMREX_SPACEDIM > mfFine;
    
    mfCoarse = GetArrOfPtrs(mf[0]);
    mfFine = GetArrOfPtrs(mf[1]);
    
    FaceFillCoarse(mfCoarse,mfFine,map);
}

void FaceFillCoarse(std::array<MultiFab*, AMREX_SPACEDIM> & mfCoarse,std::array<MultiFab*, AMREX_SPACEDIM> & mfFine, int map)
{
    BL_PROFILE_VAR("FaceFillCoarse()",FaceFillCoarse);

    const IntVect ratio(AMREX_D_DECL(2,2,2));    
    
    if(map==0)
    {    
        const std::array<MultiFab*, AMREX_SPACEDIM> coarse = {AMREX_D_DECL(mfCoarse[0],mfCoarse[1],mfCoarse[2])};
        const std::array<const MultiFab*, AMREX_SPACEDIM> fine = {AMREX_D_DECL(mfFine[0],mfFine[1],mfFine[2])};
        amrex::average_down_faces(fine,coarse,ratio);
    }
    else if(map == 1)
    {
        
        std::array< MultiFab, AMREX_SPACEDIM > coarseTemp;
        const BoxArray& fine_BA = mfFine[0]->boxArray();
        BoxArray coarsened_BA = fine_BA;
        coarsened_BA.coarsen(2);
        
        for(int d=0;d<AMREX_SPACEDIM;d++)
        {
            
            coarseTemp[d].define(convert(coarsened_BA,nodal_flag_dir[d]), mfFine[d]->DistributionMap(), 1, mfFine[d]->nGrow());
            
            coarseTemp[d].setVal(0.0);
        }

        for ( MFIter mfi(coarseTemp[0],TilingIfNotGPU()); mfi.isValid(); ++mfi ) {


        AMREX_D_TERM(Box bx_x = mfi.tilebox(nodal_flag_x);,
                     Box bx_y = mfi.tilebox(nodal_flag_y);,
                     Box bx_z = mfi.tilebox(nodal_flag_z););

        AMREX_D_TERM(Array4<Real> const& phix_c_fab = coarseTemp[0].array(mfi);,
                     Array4<Real> const& phiy_c_fab = coarseTemp[1].array(mfi);,
                     Array4<Real> const& phiz_c_fab = coarseTemp[2].array(mfi););

        AMREX_D_TERM(Array4<Real const> const& phix_f_fab = mfFine[0]->array(mfi);,
                     Array4<Real const> const& phiy_f_fab = mfFine[1]->array(mfi);,
                     Array4<Real const> const& phiz_f_fab = mfFine[2]->array(mfi););


        amrex::ParallelFor(bx_x, bx_y, bx_z, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            phix_c_fab(i,j,k) = 0.125* ( phix_f_fab(2*i  ,2*j,2*k  ) + phix_f_fab(2*i  ,2*j+1,2*k  )
                                        +phix_f_fab(2*i  ,2*j,2*k+1) + phix_f_fab(2*i  ,2*j+1,2*k+1) )
                             + 0.0625* ( phix_f_fab(2*i+1,2*j,2*k  ) + phix_f_fab(2*i+1,2*j+1,2*k  )
                                        +phix_f_fab(2*i+1,2*j,2*k+1) + phix_f_fab(2*i+1,2*j+1,2*k+1) )
                             + 0.0625* ( phix_f_fab(2*i-1,2*j,2*k  ) + phix_f_fab(2*i-1,2*j+1,2*k  )
                                        +phix_f_fab(2*i-1,2*j,2*k+1) + phix_f_fab(2*i-1,2*j+1,2*k+1) );
        },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            phiy_c_fab(i,j,k) = 0.125* ( phiy_f_fab(2*i,2*j  ,2*k  ) + phiy_f_fab(2*i+1,2*j  ,2*k  )
                                        +phiy_f_fab(2*i,2*j  ,2*k+1) + phiy_f_fab(2*i+1,2*j  ,2*k+1) )
                             + 0.0625* ( phiy_f_fab(2*i,2*j+1,2*k  ) + phiy_f_fab(2*i+1,2*j+1,2*k  )
                                        +phiy_f_fab(2*i,2*j+1,2*k+1) + phiy_f_fab(2*i+1,2*j+1,2*k+1) )
                             + 0.0625* ( phiy_f_fab(2*i,2*j-1,2*k  ) + phiy_f_fab(2*i+1,2*j-1,2*k  )
                                        +phiy_f_fab(2*i,2*j-1,2*k+1) + phiy_f_fab(2*i+1,2*j-1,2*k+1) );
        },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            phiz_c_fab(i,j,k) = 0.125* ( phiz_f_fab(2*i,2*j  ,2*k  ) + phiz_f_fab(2*i+1,2*j  ,2*k  )
                                        +phiz_f_fab(2*i,2*j+1,2*k  ) + phiz_f_fab(2*i+1,2*j+1,2*k  ) )
                             + 0.0625* ( phiz_f_fab(2*i,2*j  ,2*k+1) + phiz_f_fab(2*i+1,2*j  ,2*k+1)
                                        +phiz_f_fab(2*i,2*j+1,2*k+1) + phiz_f_fab(2*i+1,2*j+1,2*k+1) )
                             + 0.0625* ( phiz_f_fab(2*i,2*j  ,2*k-1) + phiz_f_fab(2*i+1,2*j  ,2*k-1)
                                        +phiz_f_fab(2*i,2*j+1,2*k-1) + phiz_f_fab(2*i+1,2*j+1,2*k-1) );
        });

        }
        for(int d=0;d<AMREX_SPACEDIM;d++)
        {
            mfCoarse[d]->ParallelCopy(coarseTemp[d], 0, 0, 1, 0, 0);    
        }
    }   
}

void FaceFillFine(Vector<std::array< MultiFab, AMREX_SPACEDIM >>& mf, Vector<Geometry> geom,Vector<amrex::BCRec> bcs, int map)
{
    std::array< MultiFab*, AMREX_SPACEDIM > mfCoarse;
    std::array< MultiFab*, AMREX_SPACEDIM > mfFine;
    
    mfCoarse = GetArrOfPtrs(mf[0]);
    mfFine = GetArrOfPtrs(mf[1]);
    
    Geometry* geomp = &geom[0];
    
    FaceFillFine(mfCoarse,mfFine, geomp, bcs, map);
    
}

void FaceFillFine(std::array<MultiFab, AMREX_SPACEDIM>* & mf, Vector<Geometry> geom,Vector<amrex::BCRec> bcs, int map)
{

    std::array< MultiFab*, AMREX_SPACEDIM > mfCoarse;
    std::array< MultiFab*, AMREX_SPACEDIM > mfFine;
    
    mfCoarse = GetArrOfPtrs(mf[0]);
    mfFine = GetArrOfPtrs(mf[1]);
    Geometry* geomp = &geom[0];
    
    FaceFillFine(mfCoarse,mfFine, geomp, bcs, map);
}

void FaceFillFine(Vector<std::array< MultiFab, AMREX_SPACEDIM >>& mf, Vector<Geometry> geom, int map)
{
    std::array<MultiFab, AMREX_SPACEDIM>* mfp = &mf[0];
    Geometry* geomp = &geom[0];
    
    FaceFillFine(mfp,geomp,map);
}

void FaceFillFine(std::array<MultiFab, AMREX_SPACEDIM>* & mf, Geometry* geom, int map)
{

    std::array< MultiFab*, AMREX_SPACEDIM > mfCoarse;
    std::array< MultiFab*, AMREX_SPACEDIM > mfFine;
    
    mfCoarse = GetArrOfPtrs(mf[0]);
    mfFine = GetArrOfPtrs(mf[1]);
    
    amrex::Vector<amrex::BCRec> bcs(1);
    
    int bc_lo[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    int bc_hi[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};

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
    
    FaceFillFine(mfCoarse,mfFine, geom, bcs, map);
}

void FaceFillFine(std::array<MultiFab*, AMREX_SPACEDIM> & mfCoarse,std::array<MultiFab*, AMREX_SPACEDIM> & mfFine, Geometry* geom,Vector<amrex::BCRec>  bcs, int map)
{
    BL_PROFILE_VAR("VelFillFine()",VelFillFine);
    
    Interpolater* mapper;
    
    if(map==0)
    {
        mapper = &face_divfree_interp;
    }else
    {
        mapper = &face_linear_interp;
    }
    
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
                                     cp, 0, fp, 0, {2,2,2},
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
                                     cp, 0, fp, 0, {2,2,2},
                                     mapper, bcr, 0);

    }

}


void FaceFillGhost(Vector<std::array< MultiFab, AMREX_SPACEDIM >>& mf, Vector<Geometry> geom,Vector<amrex::BCRec> bcs, int map)
{
    std::array< MultiFab*, AMREX_SPACEDIM > mfCoarse;
    std::array< MultiFab*, AMREX_SPACEDIM > mfFine;
    
    mfCoarse = GetArrOfPtrs(mf[0]);
    mfFine = GetArrOfPtrs(mf[1]);

    Geometry* geomp = &geom[0];
    
    FaceFillGhost(mfCoarse,mfFine, geomp, bcs, map);
    
}

void FaceFillGhost(std::array<MultiFab, AMREX_SPACEDIM>* & mf, Vector<Geometry> geom,Vector<amrex::BCRec> bcs, int map)
{

    std::array< MultiFab*, AMREX_SPACEDIM > mfCoarse;
    std::array< MultiFab*, AMREX_SPACEDIM > mfFine;
    
    mfCoarse = GetArrOfPtrs(mf[0]);
    mfFine = GetArrOfPtrs(mf[1]);
    
    Geometry* geomp = &geom[0];
    
    FaceFillGhost(mfCoarse,mfFine, geomp, bcs, map);
}

void FaceFillGhost(Vector<std::array< MultiFab, AMREX_SPACEDIM >>& mf, Vector<Geometry> geom, int map)
{
    std::array<MultiFab, AMREX_SPACEDIM>* mfp = &mf[0];
    Geometry* geomp = &geom[0];
    
    FaceFillGhost(mfp,geomp,map);
}

void FaceFillGhost(std::array<MultiFab, AMREX_SPACEDIM>* & mf, Geometry* geom, int map)
{

    std::array< MultiFab*, AMREX_SPACEDIM > mfCoarse;
    std::array< MultiFab*, AMREX_SPACEDIM > mfFine;
    
    mfCoarse = GetArrOfPtrs(mf[0]);
    mfFine = GetArrOfPtrs(mf[1]);
    
    amrex::Vector<amrex::BCRec> bcs(1);
    
    int bc_lo[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    int bc_hi[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};

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
    
    FaceFillGhost(mfCoarse,mfFine, geom, bcs, map);
}

void FaceFillGhost(std::array<MultiFab*, AMREX_SPACEDIM> & mfCoarse,std::array<MultiFab*, AMREX_SPACEDIM> & mfFine, Geometry* geom,Vector<amrex::BCRec>  bcs, int map)
{
    BL_PROFILE_VAR("VelFillFine()",VelFillFine);

    std::array< MultiFab, AMREX_SPACEDIM > fineTemp;
    
    for(int d=0;d<AMREX_SPACEDIM;d++)
    {
        fineTemp[d].define(convert(mfFine[d]->boxArray(),nodal_flag_dir[d]), mfFine[d]->DistributionMap(), 1, mfFine[d]->nGrow());
        fineTemp[d].setVal(0.0);
    }
    
    std::array< MultiFab*, AMREX_SPACEDIM > tempPtr;
    tempPtr = GetArrOfPtrs(fineTemp);
    FaceFillFine(mfCoarse, tempPtr, geom, bcs, map);
    
    for(int d=0;d<AMREX_SPACEDIM;d++)
    {
        mfFine[d]->ParallelCopyToGhost(fineTemp[d],0,0,1,mfFine[d]->nGrowVect(),mfFine[d]->nGrowVect());
    }

}

void CellFillCoarse(Vector<MultiFab>& mf, Vector<Geometry> geom)
{

    MultiFab* mfp = &mf[0];
    Geometry* geomp = &geom[0];       
    CellFillCoarse(mfp,geomp);
}

void CellFillCoarse(Vector<MultiFab>& mf, Geometry* geom)
{
    MultiFab* mfp = &mf[0];    
       
    CellFillCoarse(mfp,geom);
}

void CellFillCoarse(MultiFab* & mf, Geometry* geom)
{
    BL_PROFILE_VAR("CellFillCoarse()",CellFillCoarse);

    const IntVect ratio(AMREX_D_DECL(2,2,2));    
    amrex::average_down(mf[1],mf[0], geom[1], geom[0],0,1,ratio);

}

void CellFillCoarse(MultiFab & mfFine, Geometry geomFine, MultiFab & mfCoarse, Geometry geomCoarse)
{
    BL_PROFILE_VAR("CellFillCoarse()",CellFillCoarse);

    const IntVect ratio(AMREX_D_DECL(2,2,2));    
    amrex::average_down(mfFine,mfCoarse, geomFine, geomCoarse,0,1,ratio);

}

void CellFillFine(Vector<MultiFab>& mf, Geometry* geom)
{   

    MultiFab* mfCoarse = &mf[0];    
    MultiFab* mfFine = &mf[1];
    
    CellFillFine(mfCoarse, mfFine, geom);
}    


void CellFillFine(MultiFab* & mfCoarse, MultiFab* & mfFine, Geometry*  geom)
{   
    amrex::Vector<amrex::BCRec> bcs(1);
    
    int bc_lo[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    int bc_hi[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};

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
    
    CellFillFine(mfCoarse, mfFine, geom, bcs);
}    

void CellFillFine(MultiFab* & mfCoarse, MultiFab* & mfFine, Geometry* geom, Vector<amrex::BCRec>  bcs)
{
    //BL_ASSERT(lev > 0);

    Interpolater* mapper = &cell_cons_interp;
    //CellConservativeLinear mapper;
    

    Real time=0;
    //const char br;

    if(Gpu::inLaunchRegion())
    {
        GpuBndryFuncFab<AmrCoreFill> gpu_bndry_func(AmrCoreFill{});
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > cphysbc(geom[0],bcs,gpu_bndry_func);
        PhysBCFunct<GpuBndryFuncFab<AmrCoreFill> > fphysbc(geom[1],bcs,gpu_bndry_func);
        
        //PhysBCFunctNoOp cphysbc, fphysbc;

        amrex::InterpFromCoarseLevel(mfFine[0], time, mfCoarse[0], 0, 0, 1, geom[0], geom[1],
                                     cphysbc, 0, fphysbc, 0, {2,2,2},
                                     mapper, bcs, 0);
    }
    else
    {
        CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[0],bcs,bndry_func);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[1],bcs,bndry_func);

        amrex::InterpFromCoarseLevel(mfFine[0], time, mfCoarse[0], 0, 0, 1, geom[0], geom[1],
                                     cphysbc, 0, fphysbc, 0, {2,2,2},
                                     mapper, bcs, 0);
        //Print() << "HERE!\n";
        //PrintMFLine(mfCoarse[0],0,0,5,5);
        //PrintMFLine(mfFine[0],0,0,12,12);
    }
}

void CellFillGhost(Vector<MultiFab> & mf, Vector<Geometry> & geom)
{

    //PrintMFLine(mf[0],0,0,5,5);

    MultiFab* mfCoarse = &mf[0];
    MultiFab* mfFine = &mf[1];
    
    Geometry* geomp = &geom[0];
    
    //PrintMFLine(mfCoarse[0],0,0,5,5);
        
    CellFillGhost(mfCoarse,mfFine, geomp);
}

void CellFillGhost(Vector<MultiFab> & mf, Geometry* geom)
{

    MultiFab* mfCoarse = &mf[0];
    MultiFab* mfFine = &mf[1];
        
    CellFillGhost(mfCoarse,mfFine, geom);
}

void CellFillGhost(MultiFab* & mf, Geometry* geom)
{

    MultiFab* mfCoarse = &mf[0];
    MultiFab* mfFine = &mf[1];
        
    CellFillGhost(mfCoarse,mfFine, geom);
}

void CellFillGhost(MultiFab* & mfCoarse, MultiFab* & mfFine, Geometry* geom)
{
    BL_PROFILE_VAR("CellFillGhost()",CellFillGhost);

    MultiFab fineTemp;


    fineTemp.define(mfFine->boxArray(), mfFine->DistributionMap(), 1, mfFine->nGrow());
    fineTemp.setVal(0.0);

    
    MultiFab* tempPtr = &fineTemp;
    
    CellFillFine(mfCoarse, tempPtr, geom);
   
    mfFine->ParallelCopyToGhost(fineTemp,0,0,1,mfFine->nGrowVect(),mfFine->nGrowVect());
}

