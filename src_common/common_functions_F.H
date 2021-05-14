#ifndef _common_params_F_H_
#define _common_params_F_H_

#ifdef __cplusplus
extern "C"
{
#endif

    // these are interfaces to fortran subroutines

///////////////////////////
// in common_namelist.F90

    // read in fortran namelist into common_params_module
    void read_common_namelist (const char* inputs_file, const int length);

    // copy contents of common_params_module to C++ common namespace
    void initialize_common_namespace(amrex::Real* prob_lo, amrex::Real* prob_hi, int* n_cells,
                                     int* max_grid_size, int* max_grid_size_structfact,
                                     int* max_particle_tile_size,amrex::Real* cell_depth, int* ngc,
                                     int* nvars, int* nprimvars,
                                     int* membrane_cell, int* cross_cell, amrex::Real* transmission,
                                     amrex::Real* qval, int* pkernel_fluid, int* pkernel_es,
                                     int* eskernel_fluid, amrex::Real* eskernel_beta,
                                     amrex::Real* fixed_dt,
                                     amrex::Real* cfl, amrex::Real* rfd_delta, int* max_step, int* plot_int, int* plot_stag,
                                     const char* plot_base_name,
                                     const int plot_base_name_len,
                                     int* chk_int,
                                     const char* chk_base_name, const int chk_base_name_len,
                                     int* prob_type, int* restart, int* particle_restart, int* print_int,
                                     int* project_eos_int, amrex::Real* grav, int* nspecies,
                                     amrex::Real* molmass,amrex::Real* diameter,
                                     int* dof,amrex::Real* hcv,amrex::Real* hcp, amrex::Real* rhobar,
                                     amrex::Real* rho0, amrex::Real* variance_coef_mom,
                                     amrex::Real* variance_coef_mass, amrex::Real* k_B,
                                     amrex::Real* Runiv, amrex::Real* T_init,
                                     int* algorithm_type,
                                     int* advection_type,
                                     int* barodiffusion_type, int* use_bl_rng, int* seed,
                                     int* seed_momentum, int* seed_diffusion,
                                     int* seed_reaction,
                                     int* seed_init_mass, int* seed_init_momentum,
                                     amrex::Real* visc_coef, int* visc_type,
                                     int* filtering_width, int* stoch_stress_form,
                                     amrex::Real* u_init, amrex::Real* perturb_width,
                                     amrex::Real* smoothing_width,
                                     amrex::Real* initial_variance_mom,
                                     amrex::Real* initial_variance_mass,
				     amrex::Real* domega,
                                     int* bc_vel_lo, int* bc_vel_hi,
                                     int* bc_es_lo, int* bc_es_hi,
                                     int* bc_mass_lo, int* bc_mass_hi,
                                     int* bc_therm_lo, int* bc_therm_hi,
                                     amrex::Real * p_lo, amrex::Real * p_hi,
                                     amrex::Real* t_lo, amrex::Real* t_hi,
                                     amrex::Real* bc_Yk_x_lo, amrex::Real* bc_Yk_x_hi,
                                     amrex::Real* bc_Yk_y_lo, amrex::Real* bc_Yk_y_hi,
                                     amrex::Real* bc_Yk_z_lo, amrex::Real* bc_Yk_z_hi,
                                     amrex::Real* bc_Xk_x_lo, amrex::Real* bc_Xk_x_hi,
                                     amrex::Real* bc_Xk_y_lo, amrex::Real* bc_Xk_y_hi,
                                     amrex::Real* bc_Xk_z_lo, amrex::Real* bc_Xk_z_hi,
                                     amrex::Real* wallspeed_lo, amrex::Real* wallspeed_hi,
                                     amrex::Real* potential_lo, amrex::Real* potential_hi,
                                     int* struct_fact_int, int* radialdist_int,
                                     int* cartdist_int, int* n_steps_skip,
                                     amrex::Real* binSize, amrex::Real* searchDist,
				     int* project_dir, int* max_grid_projection,
                                     int* histogram_unit, amrex::Real* density_weights,
                                     int* shift_cc_to_boundary,
                                     int* particle_placement, int* particle_count, int* p_move_tog, int* p_force_tog, int* p_int_tog, int* p_int_tog_wall,
                                     amrex::Real* particle_neff, amrex::Real* particle_n0,
                                     amrex::Real* mass, amrex::Real* nfrac, amrex::Real* permittivity,
                                     int* wall_mob, amrex::Real* rmin, amrex::Real* rmax, amrex::Real* eepsilon, amrex::Real* sigma,
                                     amrex::Real* rmin_wall, amrex::Real* rmax_wall, amrex::Real* eepsilon_wall, amrex::Real* sigma_wall,
                                     int* poisson_verbose,
                                     int* poisson_bottom_verbose,
                                     int* poisson_max_iter, amrex::Real* poisson_rel_tol,
                                     amrex::Real* particle_grid_refine, amrex::Real* es_grid_refine,
                                     amrex::Real* diff, int* all_dry,
                                     int* fluid_tog, int* es_tog, int* drag_tog,
                                     int* move_tog, int* rfd_tog, int* dry_move_tog, int* sr_tog,
                                     int* graphene_tog, int* crange, int* thermostat_tog, int* zero_net_force, int* images,
                                     amrex::Real* eamp, amrex::Real* efreq, amrex::Real* ephase,
                                     int* plot_ascii, int* solve_chem, amrex::Real* diffcoeff,
                                     amrex::Real* scaling_factor, amrex::Real* source_strength,
                                     int* regrid_int, int* do_reflux, int* particle_motion,
                                     amrex::Real* turb_a, amrex::Real* turb_b, int* turbForcing);

    void set_max_step(int* max_step);
    void set_domega(amrex::Real* domega);
    
///////////////////////////

///////////////////////////
// in init_funcs.F90

    void init_vel(const int* lo, const int* hi,
                  amrex_real* data, const int* dlo, const int* dhi,
                  const amrex_real* dx, const amrex_real* prob_lo, const amrex_real* prob_hi, const int* di, const amrex::Real* real_lo, const amrex::Real* real_hi);

///////////////////////////

///////////////////////////
// in sqrt_mf.F90
    
  void sqrt_mf(const int* lo, const int* hi,
	       const amrex_real* mf, const int* mflo, const int* mfhi, const int* ncomp);

#ifdef __cplusplus
}
#endif

#endif
