
#include "common_namespace.H" // included so we get MAX_SPECIES

namespace multispec {

    // Below are items that are read in from the multispec namelist

    extern int                        inverse_type;
    extern int                        temp_type;
    extern AMREX_GPU_MANAGED int      chi_iterations;
    extern amrex::Real                start_time;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_ELEMENT> Dbar;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> Dtherm;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_ELEMENT> H_offdiag;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> H_diag;
    extern AMREX_GPU_MANAGED amrex::Real fraction_tolerance;
    extern int                        correct_flux;
    extern int                        print_error_norms;
    extern AMREX_GPU_MANAGED int      is_nonisothermal;
    extern AMREX_GPU_MANAGED int      is_ideal_mixture;
    extern AMREX_GPU_MANAGED int      use_lapack;
    extern AMREX_GPU_MANAGED int      use_multiphase;
    extern AMREX_GPU_MANAGED int      use_flory_huggins;
    extern AMREX_GPU_MANAGED amrex::Array2D<amrex::Real,0,MAX_SPECIES-1,0,MAX_SPECIES-1>  fh_kappa;
    extern AMREX_GPU_MANAGED amrex::Array2D<amrex::Real,0,MAX_SPECIES-1,0,MAX_SPECIES-1>  fh_chi;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> fh_monomers;
    extern AMREX_GPU_MANAGED amrex::Real monomer_mass;
    extern AMREX_GPU_MANAGED amrex::Real kc_tension;
    extern AMREX_GPU_MANAGED amrex::Real alpha_gex;
    extern AMREX_GPU_MANAGED int      n_gex;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> c_init_1;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> c_init_2;

    extern int                        midpoint_stoch_mass_flux_type;
    extern AMREX_GPU_MANAGED int      avg_type;
    extern AMREX_GPU_MANAGED int      mixture_type;

    // charged fluid
    extern int                        use_charged_fluid;
    extern int                        print_debye_len;
    extern amrex::Real                dielectric_const;
    extern int                        dielectric_type;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> charge_per_mass;
    extern amrex::Real                theta_pot;
    extern int                        num_pot_iters;
    extern amrex::Real                dpdt_factor;
    extern amrex::Real                relxn_param_charge;
    extern int                        E_ext_type;
    extern amrex::Vector<amrex::Real> E_ext_value;
    extern int                        electroneutral;
    extern int                        induced_charge_eo;
    extern int                        zero_eps_on_wall_type;
    extern int                        zero_charge_on_wall_type;
    extern amrex::Real                zero_eps_on_wall_left_end;
    extern amrex::Real                zero_eps_on_wall_right_start;
    extern int                        bc_function_type;
    extern amrex::Real                L_pos;
    extern amrex::Real                L_trans;
    extern amrex::Real                L_zero;
}

