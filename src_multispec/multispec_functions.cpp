#include "multispec_functions.H"

void InitializeMultispecNamespace() {

    Dbar.resize(MAX_ELEMENT);
    Dtherm.resize(MAX_SPECIES);
    H_offdiag.resize(MAX_ELEMENT);
    H_diag.resize(MAX_SPECIES);
    c_init.resize(2*MAX_SPECIES);
    c_bc.resize(AMREX_SPACEDIM*2*MAX_SPECIES);

    charge_per_mass.resize(MAX_SPECIES);
    Epot_wall_bc_type.resize(2*AMREX_SPACEDIM);
    Epot_wall.resize(2*AMREX_SPACEDIM);
    E_ext_value.resize(3);
    

    initialize_multispec_namespace( &inverse_type, &temp_type, 
				    &chi_iterations, &start_time, 
				    Dbar.dataPtr(), Dtherm.dataPtr(), 
				    H_offdiag.dataPtr(), H_diag.dataPtr(), 
				    &fraction_tolerance, &correct_flux, 
				    &print_error_norms,
				    &is_nonisothermal, &is_ideal_mixture,
				    &use_lapack, 
				    c_init.dataPtr(), 
				    c_bc.dataPtr(),
				    &midpoint_stoch_mass_flux_type, 
				    &avg_type, &mixture_type,
                                    &use_charged_fluid,
                                    &print_debye_len,
                                    &dielectric_const,
                                    &dielectric_type,
                                    charge_per_mass.dataPtr(),
                                    Epot_wall_bc_type.dataPtr(),
                                    Epot_wall.dataPtr(),
                                    &theta_pot,
                                    &num_pot_iters,
                                    &dpdt_factor,
                                    &relxn_param_charge,
                                    &E_ext_type,
                                    E_ext_value.dataPtr(),
                                    &electroneutral,
                                    &induced_charge_eo,
                                    &zero_eps_on_wall_type,
                                    &zero_charge_on_wall_type,
                                    &zero_eps_on_wall_left_end,
                                    &zero_eps_on_wall_right_start,
                                    &epot_mg_verbose,
                                    &epot_mg_abs_tol,
                                    &epot_mg_rel_tol,
                                    &bc_function_type,
                                    &L_pos,
                                    &L_trans,
                                    &L_zero);
    
}
