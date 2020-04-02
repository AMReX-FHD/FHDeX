#include "multispec_functions.H"


void InitializeMultispecNamespace() {

    int max_element = MAX_SPECIES*(MAX_SPECIES-1)/2;

    Dbar.resize(max_element);
    Dtherm.resize(MAX_SPECIES);
    H_offdiag.resize(max_element);
    H_diag.resize(MAX_SPECIES);

    c_init.resize(2*MAX_SPECIES);

    c_bc.resize(AMREX_SPACEDIM*2*MAX_SPECIES);
    
    // c_init.resize(2);
    // for (int i=0; i<2; ++i) {
    //   c_init[i].resize(MAX_SPECIES);
    // }

    // c_bc.resize(AMREX_SPACEDIM);
    // for (int i=0; i<AMREX_SPACEDIM; ++i) {
    //   c_bc[i].resize(2);
    //   for (int j=0; j<2; ++j) {
    // 	c_bc[i][j].resize(MAX_SPECIES);
    //   }
    // }

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
				    &avg_type, &mixture_type);
    
}
