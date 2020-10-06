namespace multispec {

    // Below are items that are read in from the multispec namelist

    extern int                        inverse_type;
    extern int                        temp_type;
    extern int                        chi_iterations;
    extern amrex::Real                start_time;
    extern amrex::Vector<amrex::Real> Dbar;
    extern amrex::Vector<amrex::Real> Dtherm;
    extern amrex::Vector<amrex::Real> H_offdiag;
    extern amrex::Vector<amrex::Real> H_diag;
    extern amrex::Real                fraction_tolerance;
    extern int                        correct_flux;
    extern int                        print_error_norms;
    extern int                        is_nonisothermal;
    extern int                        is_ideal_mixture;
    extern int                        use_lapack;
    extern amrex::Vector<amrex::Real> c_init;
    extern amrex::Vector<amrex::Real> c_bc;
  
    extern int                        midpoint_stoch_mass_flux_type;
    extern int                        avg_type;
    extern int                        mixture_type;

    // charged fluid
    extern int                        use_charged_fluid;
    extern int                        print_debye_len;
    extern amrex::Real                dielectric_const;
    extern int                        dielectric_type;
    extern amrex::Vector<amrex::Real> charge_per_mass;
    extern amrex::Vector<amrex::Real> Epot_wall_bc_type;
    extern amrex::Vector<amrex::Real> Epot_wall;
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
    extern int                        epot_mg_verbose;
    extern amrex::Real                epot_mg_abs_tol;
    extern amrex::Real                epot_mg_rel_tol;
    extern int                        bc_function_type;
    extern amrex::Real                L_pos;
    extern amrex::Real                L_trans;
    extern amrex::Real                L_zero;
}

