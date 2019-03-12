int                        multispec::inverse_type;
int                        multispec::temp_type;
int                        multispec::chi_iterations;
amrex::Real                multispec::start_time;
amrex::Vector<amrex::Real> multispec::Dbar;
amrex::Vector<amrex::Real> multispec::Dtherm;
amrex::Vector<amrex::Real> multispec::H_offdiag;
amrex::Vector<amrex::Real> multispec::H_diag;
amrex::Real                multispec::fraction_tolerance;
int                        multispec::correct_flux;
int                        multispec::print_error_norms;
int                        multispec::plot_stag;
int                        multispec::is_nonisothermal;
int                        multispec::is_ideal_mixture;
int                        multispec::use_lapack;
amrex::Vector<amrex::Real> multispec::c_init;
amrex::Vector<amrex::Real> multispec::c_bc;
  
amrex::Real                multispec::alpha1;
amrex::Real                multispec::beta;
amrex::Real                multispec::delta;
amrex::Real                multispec::sigma;
  
int                        multispec::midpoint_stoch_mass_flux_type;
int                        multispec::avg_type;
int                        multispec::mixture_type;
