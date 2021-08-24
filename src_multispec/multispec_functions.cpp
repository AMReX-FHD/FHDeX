#include "multispec_functions.H"

int                                                         multispec::inverse_type;
int                                                         multispec::temp_type;
int                                                         multispec::chi_iterations;
amrex::Real                                                 multispec::start_time;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_ELEMENT> multispec::Dbar;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> multispec::Dtherm;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_ELEMENT> multispec::H_offdiag;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> multispec::H_diag;
amrex::Real                                                 multispec::fraction_tolerance;
int                                                         multispec::correct_flux;
int                                                         multispec::print_error_norms;
AMREX_GPU_MANAGED int                                       multispec::is_nonisothermal;
AMREX_GPU_MANAGED int                                       multispec::is_ideal_mixture;
int                                                         multispec::use_lapack;
AMREX_GPU_MANAGED int                                       multispec::use_multiphase;
AMREX_GPU_MANAGED amrex::Real                               multispec::kc_tension;
AMREX_GPU_MANAGED amrex::Real                               multispec::alpha_gex;
AMREX_GPU_MANAGED int                                       multispec::n_gex;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> multispec::c_init_1;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> multispec::c_init_2;
  
int                                                         multispec::midpoint_stoch_mass_flux_type;
AMREX_GPU_MANAGED int                                       multispec::avg_type;
int                                                         multispec::mixture_type;

// charged fluid
int                                                         multispec::use_charged_fluid;
int                                                         multispec::print_debye_len;
amrex::Real                                                 multispec::dielectric_const;
int                                                         multispec::dielectric_type;
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> multispec::charge_per_mass;
amrex::Real                                                 multispec::theta_pot;
int                                                         multispec::num_pot_iters;
amrex::Real                                                 multispec::dpdt_factor;
amrex::Real                                                 multispec::relxn_param_charge;
int                                                         multispec::E_ext_type;
amrex::Vector<amrex::Real>                                  multispec::E_ext_value;
int                                                         multispec::electroneutral;
int                                                         multispec::induced_charge_eo;
int                                                         multispec::zero_eps_on_wall_type;
int                                                         multispec::zero_charge_on_wall_type;
amrex::Real                                                 multispec::zero_eps_on_wall_left_end;
amrex::Real                                                 multispec::zero_eps_on_wall_right_start;
int                                                         multispec::bc_function_type;
amrex::Real                                                 multispec::L_pos;
amrex::Real                                                 multispec::L_trans;
amrex::Real                                                 multispec::L_zero;

void InitializeMultispecNamespace() {

    //Dtherm.resize(MAX_SPECIES);
    //H_offdiag.resize(MAX_ELEMENT);
    //H_diag.resize(MAX_SPECIES);

    E_ext_value.resize(3);
    
    initialize_multispec_namespace( &inverse_type, &temp_type, 
				    &chi_iterations, &start_time, 
				    Dbar.data(), Dtherm.data(), 
				    H_offdiag.data(), H_diag.data(), //HACK confirm changes
				    &fraction_tolerance, &correct_flux, 
				    &print_error_norms,
				    &is_nonisothermal, &is_ideal_mixture,
				    &use_lapack, &use_multiphase,
                                    &kc_tension, &alpha_gex, &n_gex,
				    c_init_1.data(),
                                    c_init_2.data(),
				    &midpoint_stoch_mass_flux_type, 
				    &avg_type, &mixture_type,
                                    &use_charged_fluid,
                                    &print_debye_len,
                                    &dielectric_const,
                                    &dielectric_type,
                                    charge_per_mass.data(),
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
                                    &bc_function_type,
                                    &L_pos,
                                    &L_trans,
                                    &L_zero);
    
}
