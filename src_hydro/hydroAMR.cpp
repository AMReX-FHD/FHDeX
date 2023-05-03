#include "gmres_functions.H"
#include "common_functions.H"
#include "hydroAMR.H"

using namespace amrex;

hydroAMR::hydroAMR(int ang, int * is_periodic, Real dt_in) : AmrCore() {

//    IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
//    IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
//    IntVect max_grid(AMREX_D_DECL(common::max_grid_size[0], common::max_grid_size[1], common::max_grid_size[2]));


    //domain.resize(nlevels);
    //ba.resize(nlevels);
    //dmap.resize(nlevels);
    pres.resize(nlevels);
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

    umac.resize(nlevels);
    umacNew.resize(nlevels);    
    umacM.resize(nlevels);
    alpha_fc.resize(nlevels);
    source.resize(nlevels);
    sourceTemp.resize(nlevels);
    sourceRFD.resize(nlevels);
    beta_ed.resize(nlevels);
    eta_ed.resize(nlevels);
    temp_ed.resize(nlevels);

    ng = ang;
    dt = dt_in;


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

void hydroAMR::updateGrid ()
{
    const Real time = 0.0;
    regrid(0,time);
    
}

void
hydroAMR::ErrorEst (int lev, TagBoxArray& tags, Real /*time*/, int /*ngrow*/)
{
    if(init)
    {
        init = false;
        cc_mask.define(boxArray(0),DistributionMap(0),1,0);
        cc_mask.setVal(0);
    }else
    {
        for (MFIter mfi(cc_mask,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx  = mfi.tilebox();
            const auto tagfab  = tags.array(mfi);
            const auto maskfab = cc_mask.array(mfi);

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

        for (MFIter mfi(cc_mask,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx  = mfi.tilebox();
            const auto maskfab = cc_mask.array(mfi);

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
            sourceTemp[lev][d].define(convert(ba,nodal_flag_dir[d]), dm, 1, ng);
            umacM[lev][d].define(convert(ba,nodal_flag_dir[d]), dm, 1, 1);
            //fc_mask[d].define(convert(ba,nodal_flag_dir[d]), dm, 1, ng);

            umac[lev][d].setVal(0.);
            source[lev][d].setVal(0.);
            sourceRFD[lev][d].setVal(0.);
            sourceTemp[lev][d].setVal(0.);                        
            umacM[lev][d].setVal(0.);
            //fc_mask[d].setVal(0);
    }
    
    pres[lev].define(ba,dm,1,1);
    pres[lev].setVal(0.);
    
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

    //cc_mask.define(ba, dm, 1, 1);

    eta_cc[lev].define(ba, dm, 1, 1);
    temp_cc[lev].define(ba, dm, 1, 1);
    eta_cc[lev].setVal(visc_coef);
    temp_cc[lev].setVal(T_init[0]);


#if (AMREX_SPACEDIM == 2)
    beta_ed[lev][0].define(convert(ba,nodal_flag), dm, 1, 1);
    ed_mask[lev].define(convert(ba,nodal_flag), dm, 1, 1);
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



