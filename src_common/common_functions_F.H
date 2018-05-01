#ifndef _common_functions_F_H_
#define _common_functions_F_H_

#ifdef __cplusplus
extern "C" 
{
#endif

    // these are interfaces to fortran subroutines

    //////////////////////
    // in common_params.F90

    // read in fortran namelists into common_params_module
    void read_common_params (const char* inputs_file, const int length);

    // copy contents of common_params_module to C++ common namespace
    void copy_common_params_to_c (amrex::Real* prob_lo, amrex::Real* prob_hi, int* n_cells,
                                  int* max_grid_size, amrex::Real* fixed_dt, amrex::Real* cfl,
                                  int* max_step, int* plot_int, const char* plot_base_name,
                                  const int plot_base_name_len, int* chk_int,
                                  const char* chk_base_name, const int chk_base_name_len,
                                  int* prob_type, int* restart, int* print_int,
                                  int* project_eos_int, amrex::Real* grav, int* nspecies,
                                  amrex::Real* molmass, amrex::Real* rhobar,
                                  amrex::Real* rho0, amrex::Real* variance_coef_mom,
                                  amrex::Real* variance_coef_mass, amrex::Real* k_B,
                                  amrex::Real* Runiv, int* algorithm_type,
                                  int* barodiffusion_type, int* use_bl_rng, int* seed,
                                  int* seed_momentum, int* seed_diffusion, int* seed_reaction,
                                  int* seed_init_mass, int* seed_init_momentum,
                                  amrex::Real* visc_coef, int* visc_type, int* advection_type,
                                  int* filtering_width, int* stoch_stress_form,
                                  amrex::Real* u_init, amrex::Real* perturb_width,
                                  amrex::Real* smoothing_width,
                                  amrex::Real* initial_variance_mom,
                                  amrex::Real* initial_variance_mass,
                                  int* bc_lo, int* bc_hi,
                                  amrex::Real* wallspeed_lo, amrex::Real* wallspeed_hi,
                                  int* histogram_unit, amrex::Real* density_weights,
                                  int* shift_cc_to_boundary);
    //////////////////////

    //////////////////////
    // in common_functions.F90


    //////////////////////

#ifdef __cplusplus
}
#endif

#endif
