#include "common_functions.H"
#include "common_functions_F.H"
#include "common_params.H"

using namespace common;

void set_common_params() {

    prob_lo.resize(AMREX_SPACEDIM);
    prob_hi.resize(AMREX_SPACEDIM);
    n_cells.resize(AMREX_SPACEDIM);
    max_grid_size.resize(AMREX_SPACEDIM);
    plot_base_name.resize(128);
    chk_base_name.resize(128);
    grav.resize(AMREX_SPACEDIM);
    molmass.resize(MAX_SPECIES);
    rhobar.resize(MAX_SPECIES);
    u_init.resize(2);
    bc_lo.resize(AMREX_SPACEDIM);
    bc_hi.resize(AMREX_SPACEDIM);
    wallspeed_lo.resize((AMREX_SPACEDIM-1)*AMREX_SPACEDIM);
    wallspeed_hi.resize((AMREX_SPACEDIM-1)*AMREX_SPACEDIM);
    density_weights.resize(MAX_SPECIES);
    shift_cc_to_boundary.resize(AMREX_SPACEDIM*LOHI);
    
    copy_common_params_to_c(prob_lo.dataPtr(), prob_hi.dataPtr(), n_cells.dataPtr(),
                            max_grid_size.dataPtr(), &fixed_dt, &cfl, &max_step,
                            &plot_int, plot_base_name.c_str(), plot_base_name.size()+1,
                            &chk_int, chk_base_name.c_str(), chk_base_name.size()+1,
                            &prob_type, &restart, &print_int, &project_eos_int,
                            grav.dataPtr(), &nspecies, molmass.dataPtr(), rhobar.dataPtr(),
                            &rho0, &variance_coef_mom, &variance_coef_mass, &k_B, &Runiv,
                            &algorithm_type, &barodiffusion_type, &use_bl_rng, &seed,
                            &seed_momentum, &seed_diffusion, &seed_reaction, &seed_init_mass,
                            &seed_init_momentum, &visc_coef, &visc_type, &advection_type,
                            &filtering_width, &stoch_stress_form, u_init.dataPtr(),
                            &perturb_width, &smoothing_width, &initial_variance_mom,
                            &initial_variance_mass, bc_lo.dataPtr(), bc_hi.dataPtr(),
                            wallspeed_lo.dataPtr(), wallspeed_hi.dataPtr(), &histogram_unit,
                            density_weights.dataPtr(), shift_cc_to_boundary.dataPtr());

}
