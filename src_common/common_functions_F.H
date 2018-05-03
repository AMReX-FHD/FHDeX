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
                                     int* max_grid_size, amrex::Real* fixed_dt, 
                                     amrex::Real* cfl, int* max_step, int* plot_int, 
                                     const char* plot_base_name, 
                                     const int plot_base_name_len, 
                                     int* chk_int,
                                     const char* chk_base_name, const int chk_base_name_len,
                                     int* prob_type, int* restart, int* print_int,
                                     int* project_eos_int, amrex::Real* grav, int* nspecies,
                                     amrex::Real* molmass, amrex::Real* rhobar,
                                     amrex::Real* rho0, amrex::Real* variance_coef_mom,
                                     amrex::Real* variance_coef_mass, amrex::Real* k_B,
                                     amrex::Real* Runiv, int* algorithm_type,
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
                                     int* bc_lo, int* bc_hi,
                                     amrex::Real* wallspeed_lo, amrex::Real* wallspeed_hi,
                                     int* histogram_unit, amrex::Real* density_weights,
                                     int* shift_cc_to_boundary);
///////////////////////////

///////////////////////////
// in convert_stag.F90

    void average_face_to_cc(const int* lo, const int* hi,
                            const amrex::Real* face, 
                            const int* f_lo, const int* f_hi, const int* nc_f,
                            amrex::Real* cc, 
                            const int* c_lo, const int* c_hi, const int* nc_c,
                            const int* face_comp, const int* cc_comp, const int* ncomp,
                            const int* av_dim);

    void average_cc_to_face(const int* lo, const int* hi,
                            const amrex::Real* cc, 
                            const int* c_lo, const int* c_hi, const int* nc_c,
                            amrex::Real* facex, 
                            const int* x_lo, const int* x_hi, const int* nc_x,
                            amrex::Real* facey, 
                            const int* y_lo, const int* y_hi, const int* nc_y,
#if (AMREX_SPACEDIM == 3)
                            amrex::Real* facez, 
                            const int* z_lo, const int* z_hi, const int* nc_z,
#endif
                            const int* cc_comp, const int* face_comp, const int* ncomp);

///////////////////////////

#ifdef __cplusplus
}
#endif

#endif
