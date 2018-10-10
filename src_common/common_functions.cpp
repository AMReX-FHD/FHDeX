#include "common_functions.H"
#include "common_functions_F.H"
#include "common_namespace.H"

using namespace common;

void InitializeCommonNamespace() {

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        nodal_flag[i] = 1;
        nodal_flag_x[i] = int(i==0);
        nodal_flag_y[i] = int(i==1);
        nodal_flag_z[i] = int(i==2);
        nodal_flag_xy[i] = int(i==0 || i==1);
        nodal_flag_xz[i] = int(i==0 || i==2);
        nodal_flag_yz[i] = int(i==1 || i==2);
    }

    prob_lo.resize(AMREX_SPACEDIM);
    prob_hi.resize(AMREX_SPACEDIM);
    n_cells.resize(AMREX_SPACEDIM);
    max_grid_size.resize(AMREX_SPACEDIM);
    plot_base_name.resize(128);
    chk_base_name.resize(128);
    grav.resize(AMREX_SPACEDIM);
    molmass.resize(MAX_SPECIES);
    diameter.resize(MAX_SPECIES);
    dof.resize(MAX_SPECIES);
    hcv.resize(MAX_SPECIES);
    hcp.resize(MAX_SPECIES);
    rhobar.resize(MAX_SPECIES);
    u_init.resize(2);
    T_init.resize(2);
    bc_lo.resize(AMREX_SPACEDIM);
    bc_hi.resize(AMREX_SPACEDIM);
    t_lo.resize(AMREX_SPACEDIM);
    t_hi.resize(AMREX_SPACEDIM);
    wallspeed_lo.resize((AMREX_SPACEDIM-1)*AMREX_SPACEDIM);
    wallspeed_hi.resize((AMREX_SPACEDIM-1)*AMREX_SPACEDIM);
    density_weights.resize(MAX_SPECIES);
    shift_cc_to_boundary.resize(AMREX_SPACEDIM*LOHI);
    
    initialize_common_namespace(prob_lo.dataPtr(), prob_hi.dataPtr(), n_cells.dataPtr(),
                                max_grid_size.dataPtr(), &cell_depth, &ngc, &nvars, &nprimvars, &fixed_dt, &cfl, &max_step,
                                &plot_int, plot_base_name.c_str(), plot_base_name.size()+1,
                                &chk_int, chk_base_name.c_str(), chk_base_name.size()+1,
                                &prob_type, &restart, &print_int, &project_eos_int,
                                grav.dataPtr(), &nspecies, molmass.dataPtr(), diameter.dataPtr(), dof.dataPtr(), hcv.dataPtr(), hcp.dataPtr(), 
                                rhobar.dataPtr(),
                                &rho0, &variance_coef_mom, &variance_coef_mass, &k_B, &Runiv,
                                T_init.dataPtr(),
                                &algorithm_type,  &advection_type,
                                &barodiffusion_type, &use_bl_rng, &seed,
                                &seed_momentum, &seed_diffusion, &seed_reaction, 
                                &seed_init_mass,
                                &seed_init_momentum, &visc_coef, &visc_type,
                                &filtering_width, &stoch_stress_form, u_init.dataPtr(),
                                &perturb_width, &smoothing_width, &initial_variance_mom,
                                &initial_variance_mass, bc_lo.dataPtr(), bc_hi.dataPtr(), t_lo.dataPtr(), t_hi.dataPtr(),
                                wallspeed_lo.dataPtr(), wallspeed_hi.dataPtr(),
                                &struct_fact_int, &n_steps_skip,
                                &histogram_unit,
                                density_weights.dataPtr(), shift_cc_to_boundary.dataPtr());

}

