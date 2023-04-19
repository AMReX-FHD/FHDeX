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
        alpha_fc[d].setVal(dtinv);
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
#elif (AMREX_SPACEDIM == 3)
    beta_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 1);
    beta_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 1);
    beta_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 1);
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
    }
    
    npatches = 0;     

}

void hydroAMR::addPatch(IntVect patch_lo, IntVect patch_hi) {

    npatches++;
        
    geom_patches.resize(npatches);
    domain_patches.resize(npatches);
    ba_patches.resize(npatches);
    dmap_patches.resize(npatches);
    
    domain_patches[npatches-1].setSmall(patch_lo);
    domain_patches[npatches-1].setBig(patch_hi);    
    
    umac_patches.resize(npatches); 

}



