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
    extern int                        plot_stag;
    extern int                        is_nonisothermal;
    extern int                        is_ideal_mixture;
    extern int                        use_lapack;
    extern amrex::Vector< amrex::Vector<amrex::Real> > c_init;
    extern amrex::Vector< amrex::Vector< amrex::Vector<amrex::Real> > > c_bc;
  
    extern amrex::Real                alpha1;
    extern amrex::Real                beta;
    extern amrex::Real                delta;
    extern amrex::Real                sigma;
  
    extern int                        midpoint_stoch_mass_flux_type;
    extern int                        avg_type;
    extern int                        mixture_type;

}

