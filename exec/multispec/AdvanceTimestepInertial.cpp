
#include "hydro_functions.H"

#include "common_functions.H"

#include "gmres_functions.H"

#include "multispec_functions.H"


#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>


// argv contains the name of the inputs file entered at the command line
void AdvanceTimestepInertial(std::array< MultiFab, AMREX_SPACEDIM >& umac,
                             MultiFab& rho_old, MultiFab& rho_new,
                             MultiFab& rhotot_old, MultiFab& rhotot_new,
                             MultiFab& pi, MultiFab& eta, 
                             std::array< MultiFab, NUM_EDGE >&  eta_ed,
                             MultiFab& kappa, MultiFab& Temp,
                             std::array< MultiFab, NUM_EDGE >& Temp_ed,
                             MultiFab& diff_mass_fluxdiv,
                             MultiFab& stoch_mass_fluxdiv,
                             std::array< MultiFab, AMREX_SPACEDIM >& stoch_mass_flux,
                             const Real& dt, const Real& time, const int& istep,
                             const Geometry& geom)
{
  
    BL_PROFILE_VAR("AdvanceTimestepInertial()",AdvanceTimestepInertial);

    BoxArray ba = rho_old.boxArray();
    DistributionMapping dmap = rho_old.DistributionMap();

    Vector<Real> weights;
    weights = {1.0};
    
    Real theta_alpha = 1./dt;

    MultiFab rho_update (ba,dmap,nspecies,0);
    MultiFab gmres_rhs_p(ba,dmap,       1,0);
    MultiFab dpi        (ba,dmap,       1,1);

    std::array< MultiFab, AMREX_SPACEDIM > mold;
    std::array< MultiFab, AMREX_SPACEDIM > mtemp;
    std::array< MultiFab, AMREX_SPACEDIM > adv_mom_fluxdiv;
    std::array< MultiFab, AMREX_SPACEDIM > diff_mom_fluxdiv;
    std::array< MultiFab, AMREX_SPACEDIM > gmres_rhs_v;
    std::array< MultiFab, AMREX_SPACEDIM > dumac;
    std::array< MultiFab, AMREX_SPACEDIM > gradpi;
    std::array< MultiFab, AMREX_SPACEDIM > rhotot_fc_old;
    std::array< MultiFab, AMREX_SPACEDIM > rhotot_fc_new;
    std::array< MultiFab, AMREX_SPACEDIM > rho_fc;
    std::array< MultiFab, AMREX_SPACEDIM > diff_mass_flux;
    std::array< MultiFab, AMREX_SPACEDIM > total_mass_flux;
    std::array< MultiFab, AMREX_SPACEDIM > stoch_mom_fluxdiv;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        mold[d]             .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 1);
        mtemp[d]            .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 1);
        adv_mom_fluxdiv[d]  .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        diff_mom_fluxdiv[d] .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        gmres_rhs_v[d]      .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        dumac[d]            .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 1);
        gradpi[d]           .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 0);
        rhotot_fc_old[d]    .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 1);
        rhotot_fc_new[d]    .define(convert(ba,nodal_flag_dir[d]), dmap,        1, 1);
        rho_fc[d]           .define(convert(ba,nodal_flag_dir[d]), dmap, nspecies, 0);
        diff_mass_flux[d]   .define(convert(ba,nodal_flag_dir[d]), dmap, nspecies, 0);
        total_mass_flux[d]  .define(convert(ba,nodal_flag_dir[d]), dmap, nspecies, 0);
        stoch_mom_fluxdiv[d].define(convert(ba,nodal_flag_dir[d]), dmap, nspecies, 0);
    }

    // make copies of old quantities
    // copy umac into umac_tmp if using bds

    //////////////////////////////////////////////
    /// Step 1 - Calculate Predictor Diffusive and Stochastic Fluxes
    //////////////////////////////////////////////

    // diff/stoch_mass_fluxdiv already contain F_i
    // this was already done in Step 0 (initialization) or Step 6 from the previous time step
    
    //////////////////////////////////////////////
    // Step 2 - Predictor Euler Step
    //////////////////////////////////////////////

    /*
    if (use_charged_fluid) {
        // compute old Lorentz force
    }
    */

    // average rho_old and rhotot_old to faces
    AverageCCToFace(rho_old,rho_fc,0,nspecies);
    AverageCCToFace(rhotot_old,rhotot_fc_old,0,1);

    // add D^n and St^n to rho_update
    MultiFab::Copy(rho_update,diff_mass_fluxdiv,0,0,nspecies,0);
    if (variance_coef_mass != 0.) {
        MultiFab::Add(rho_update,stoch_mass_fluxdiv,0,0,nspecies,0);
    }

    // add A^n to rho_update
    if (advection_type >= 1) {
        Abort("AdvanceTimestepInterial: bds not supported");
    }
    else {
        MkAdvSFluxdiv(umac,rho_fc,rho_update,geom,0,nspecies,true);
    }
   
    
    
    
}
