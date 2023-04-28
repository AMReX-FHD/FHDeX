#include "gmres_functions.H"
#include "common_functions.H"
#include "hydroAMR.H"

using namespace amrex;

hydroAMR::hydroAMR(int ang, int * is_periodic, Real dt) {

    IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
    IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));

    domain.setSmall(dom_lo);
    domain.setBig(dom_hi);
    
    ba.define(domain);
    ba.maxSize(IntVect(max_grid_size));
    dmap.define(ba);
    
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
            umac  [d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ang);
            umacNew  [d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ang);
            source[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ang);
            sourceRFD[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ang);
            sourceTemp[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ang);
            umacM [d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
            umac  [d].setVal(0.);
            source[d].setVal(0.);
            sourceRFD[d].setVal(0.);
            sourceTemp[d].setVal(0.);                        
            umacM [d].setVal(0.);
    }
    
    pres.define(ba,dmap,1,1);
    pres.setVal(0.);
    
    RealBox realDomain({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                       {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});
                       
    geom.define(domain ,&realDomain,CoordSys::cartesian,is_periodic);
    
    Real dtinv = 1.0/dt;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        alpha_fc[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
        mask_fc[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 1);
        alpha_fc[d].setVal(dtinv);
        mask_fc[d].setVal(0);
    }
    
    // beta cell centred
    beta.define(ba, dmap, 1, 1);
    beta.setVal(visc_coef);
    // cell-centered gamma
    gamma.define(ba, dmap, 1, 1);
    gamma.setVal(0.);

    eta_cc.define(ba, dmap, 1, 1);
    temp_cc.define(ba, dmap, 1, 1);
    eta_cc .setVal(visc_coef);
    temp_cc.setVal(T_init[0]);


#if (AMREX_SPACEDIM == 2)
    beta_ed[0].define(convert(ba,nodal_flag), dmap, 1, 1);
    mask_ed[0].define(convert(ba,nodal_flag), dmap, 1, 1);
#elif (AMREX_SPACEDIM == 3)
    beta_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
    beta_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
    beta_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
    mask_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
    mask_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
    mask_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
#endif
    
#if (AMREX_SPACEDIM == 2)
    eta_ed [0].define(convert(ba,nodal_flag),    dmap, 1, 0);
    temp_ed[0].define(convert(ba,nodal_flag),    dmap, 1, 0);
#elif (AMREX_SPACEDIM == 3)
    eta_ed [0].define(convert(ba,nodal_flag_xy), dmap, 1, 0);
    eta_ed [1].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
    eta_ed [2].define(convert(ba,nodal_flag_yz), dmap, 1, 0);
    temp_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 0);
    temp_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
    temp_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 0);
#endif                 

    for (int d=0; d<NUM_EDGE; ++d) {
        eta_ed[d] .setVal(visc_coef);
        temp_ed[d].setVal(T_init[0]);
        beta_ed[d].setVal(visc_coef);
        mask_ed[d].setVal(0);
    }
    
    npatches = 0;
    pg = 1;     

}

void hydroAMR::addPatch(IntVect patch_lo, IntVect patch_hi, Real dt) {

    npatches++;
        
    geom_patches.resize(npatches);
    domain_patches.resize(npatches);
    ba_patches.resize(npatches);
    dmap_patches.resize(npatches);
    
    domain_patches[npatches-1].setSmall(patch_lo);
    domain_patches[npatches-1].setBig(patch_hi);

    const Real* dx = geom.CellSize();

    RealBox patchRealDomain({AMREX_D_DECL(patch_lo[0]*dx[0],patch_lo[1]*dx[1],patch_lo[2]*dx[2])},
                       {AMREX_D_DECL(patch_hi[0]*dx[0],patch_hi[1]*dx[1],patch_hi[2]*dx[2])});

    Vector<int> patch_periodic  (AMREX_SPACEDIM,0);
                       
    geom_patches[npatches-1].define(domain_patches[npatches-1] ,&patchRealDomain,CoordSys::cartesian,patch_periodic.data());
       
    BoxList bl = ba.boxList();
    BoxList bl_patch_tmp = ba.boxList();
    BoxList bl_patch(domain_patches[npatches-1]);
    Box b_patch(domain_patches[npatches-1]);
    //Vector<Box> boxVec(bl.size());
    BoxList bl_patch_new;
    //Box b_temp;
    //Print() << bl_patch << std::endl;
    
//    for(int i=0; i<bl.size(); i++)
//    {
////        int k=0;
////        for(int j=0; j<ba.size(); j++)
////        {
//              Box b = bl.data()[i];
//              BoxList bl_temp(bl.data()[i]);
//              if(b.intersects(b_patch))
//              {
//                bl_temp.intersect(b_patch);
//                bl_patch_new.push_back(bl_temp.data()[0]);

//              }else
//              {
//                   // IntVect dom_zero(b.smallEnd());
//                   // Box null_box(dom_zero, dom_zero);
//                   // bl_patch_new.push_back(null_box);
//              }
//              
//              //if(
//              
//              //Print() << bl_temp << std::endl;
////            if(testBox.contains(testPatch))
////            {
////                dm[k]=dmap[j];
////                k++;
////                //dm.push_back(dmap[j]);
////            }
////        }
//    }
    
    bl_patch_tmp.intersect(bl_patch);
    //bl.join(bl);
    ba_patches[npatches-1].define(bl_patch);

//    Print() << ba_patches[npatches-1].size() << std::endl;

//    Vector<int> dm(ba.size(),0);
    Vector<int> dm;
//    
    for(int i=0; i<ba_patches[npatches-1].size(); i++)
    {
        //int k=0;
        for(int j=0; j<ba.size(); j++)
        {
            Box testBox = ba[j];
            Box testPatch = ba_patches[npatches-1][i];
            if(testBox.contains(testPatch))
            {
                //dm[k]=dmap[j];
                //k++;
                dm.push_back(dmap[j]);
            }
        }
    }
//    
    dmap_patches[npatches-1].define(dm);
//    dmap_patches[npatches-1] = dmap;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        for ( MFIter mfi(mask_fc[d]); mfi.isValid(); ++mfi ) {

            Box tile_box  = mfi.growntilebox(pg);

            Array4<int> mask_data = mask_fc[d].array(mfi);

            amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                mask_data(i,j,k,0) = npatches;        
            });

        }
    }
    
    ba_patches[npatches-1].refine(2);

    //mask_fc_patches.resize(npatches);
    alpha_fc_patches.resize(npatches);

    BoxArray ba_patch = ba_patches[npatches-1];
    DistributionMapping dmap_patch = dmap_patches[npatches-1];

    Real dtinv = 1.0/dt;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        alpha_fc_patches[npatches-1][d].define(convert(ba_patch,nodal_flag_dir[d]), dmap_patch, 1, 1);
       // mask_fc_patches[npatches-1][d].define(convert(ba_patch,nodal_flag_dir[d]), dmap_patch, 1, 1);
        alpha_fc_patches[npatches-1][d].setVal(dtinv);
        //mask_fc_patches[npatches-1][d].setVal(0);
    }



//    for (int d=0; d<AMREX_SPACEDIM; ++d) {
//        for ( MFIter mfi(mask_fc_patches[npatches-1][d]); mfi.isValid(); ++mfi ) {

//            Box tile_box  = mfi.growntilebox(pg);

//            Array4<int> mask_data = mask_fc_patches[npatches-1][d].array(mfi);

//            amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//            {
//                
//            });

//        }
//    }

    Print() << ba << std::endl;    
    Print() << ba_patches[npatches-1] << std::endl;
//    
//    Print() << convert(ba,nodal_flag_dir[0]) << std::endl;    
//    Print() << convert(ba_patches[npatches-1],nodal_flag_dir[0]) << std::endl;
//    
    Print() << dmap << std::endl;
    Print() << dmap_patches[npatches-1] << std::endl;

    umac_patches.resize(npatches);
    umacNew_patches.resize(npatches);
    source_patches.resize(npatches);
    sourceRFD_patches.resize(npatches);
    sourceTemp_patches.resize(npatches);
    umacM_patches.resize(npatches);

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umac_patches[npatches-1][d].define(convert(ba_patch,nodal_flag_dir[d]), dmap_patch, 1, 1);
        umacNew_patches[npatches-1][d].define(convert(ba_patch,nodal_flag_dir[d]), dmap_patch, 1, 1);
        source_patches[npatches-1][d].define(convert(ba_patch,nodal_flag_dir[d]), dmap_patch, 1, 1);
        sourceRFD_patches[npatches-1][d].define(convert(ba_patch,nodal_flag_dir[d]), dmap_patch, 1, 1);
        sourceTemp_patches[npatches-1][d].define(convert(ba_patch,nodal_flag_dir[d]), dmap_patch, 1, 1);
        umacM_patches[npatches-1][d].define(convert(ba_patch,nodal_flag_dir[d]), dmap_patch, 1, 1);
        umac_patches[npatches-1][d].setVal(0.);
        source_patches[npatches-1][d].setVal(0.);
        sourceRFD_patches[npatches-1][d].setVal(0.);
        sourceTemp_patches[npatches-1][d].setVal(0.);                        
        umacM_patches[npatches-1][d].setVal(0.);
    }

}




