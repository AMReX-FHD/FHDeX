#ifndef _multispec_functions_F_H_
#define _multispec_functions_F_H_

#ifdef __cplusplus
extern "C" 
{
#endif

    // these are interfaces to fortran subroutines

    ///////////////////////////
    // in multispec_namelist.F90

    // read in fortran namelist into multispec_params_module
    void read_multispec_namelist (const char* inputs_file, const int length);

    // copy contents of multispec_params_module to C++ multispec namespace
    void initialize_multispec_namespace (const int* inverse_type,
                                         const int* temp_type,
                                         const int* chi_iterations,
                                         const amrex_real* start_time, 
                                         const amrex_real* Dbar,
                                         const amrex_real* Dtherm,
                                         const amrex_real* H_offdiag,
                                         const amrex_real* H_diag,
                                         const amrex_real* fraction_tolerance,
                                         const int* correct_flux,
                                         const int* print_error_norms,
                                         const int* is_nonisothermal,
                                         const int* is_ideal_mixture,
                                         const int* use_lapack,
                                         const amrex_real* c_init,
                                         const amrex_real* c_bc,
                                         const int* midpoint_stoch_mass_flux_type,
                                         const int* avg_type,
                                         const int* mixture_type,
                                         const int* use_charged_fluid,
                                         const int* print_debye_len,
                                         const amrex_real* dielectric_const,
                                         const int* dielectric_type,
                                         const amrex_real* charge_per_mass,
                                         const amrex_real* Epot_wall_bc_type,
                                         const amrex_real* Epot_wall,
                                         const amrex_real* theta_pot,
                                         const int* num_pot_iters,
                                         const amrex_real* dpdt_factor,
                                         const amrex_real* relxn_param_charge,
                                         const int* E_ext_type,
                                         const amrex_real* E_ext_value,
                                         const int* electroneutral,
                                         const int* induced_charge_eo,
                                         const int* zero_eps_on_wall_type,
                                         const int* zero_charge_on_wall_type,
                                         const amrex_real* zero_eps_on_wall_left_end,
                                         const amrex_real* zero_eps_on_wall_right_start,
                                         const int* epot_mg_verbose,
                                         const amrex_real* epot_mg_abs_tol,
                                         const amrex_real* epot_mg_rel_tol,
                                         const int* bc_function_type,
                                         const amrex_real* L_pos,
                                         const amrex_real* L_trans,
                                         const amrex_real* L_zero);

    ///////////////////////////
    // in compute_mixture_properties.F90
  
    void mixture_properties_mass(const int* tlo, const int* thi,
                                 const amrex_real* rho, const int* rhlo, const int* rhhi,
                                 const amrex_real* rhotot, const int* rtlo, const int* rthi,
                                 const amrex_real* D_bar, const int* dblo, const int* dbhi,
                                 const amrex_real* D_therm, const int* dtlo, const int* dthi,
                                 const amrex_real* Hessian, const int* hslo, const int* hshi);
    ///////////////////////////
    // in matvec_mul.F90

    void matvec_mul(const int* tlo, const int* thi,
                    const amrex_real* xp, const int* xlo, const int* xhi, const int* nc,
                    const amrex_real* ap, const int* alo, const int* ahi, const int* nc2);

    ///////////////////////////
    // in correction_flux.F90
  
    void correction_flux(const int* lo, const int* hi,
                         const amrex_real* rho, const int* rhlo, const int* rhhi,
                         const amrex_real* rhotot, const int* rtlo, const int* rthi,
                         const amrex_real* flux_x, const int* fxlo, const int* fxhi, 
                         const amrex_real* flux_y, const int* fylo, const int* fyhi
#if (AMREX_SPACEDIM == 3)
                         ,const amrex_real* flux_z, const int* fzlo, const int* fzhi
#endif		       
                        );

    ///////////////////////////
    // in mass_flux_utilities.F90

    void compute_molconc_molmtot(const int* tlo, const int* thi,
                                 const amrex_real* rho, const int* rhlo, const int* rhhi,
                                 const amrex_real* rhotot, const int* rtlo, const int* rthi,
                                 const amrex_real* molarconc, const int* mclo, const int* mchi,
                                 const amrex_real* molmtot, const int* mtlo, const int* mthi);
  
    void compute_Gamma(const int* tlo, const int* thi,
                       const amrex_real* molarconc, const int* mlo, const int* mhi,
                       const amrex_real* Hessian, const int* hlo, const int* hhi, 
                       const amrex_real* Gamma, const int* glo, const int* ghi);
    
    void compute_rhoWchi(const int* tlo, const int* thi,
                         const amrex_real* rho, const int* rhlo, const int* rhhi,
                         const amrex_real* rhotot, const int* rtlo, const int* rthi,
                         const amrex_real* molarconc, const int* mclo, const int* mchi,
                         const amrex_real* rhoWchi, const int* rclo, const int* rchi,
                         const amrex_real* D_bar, const int* dblo, const int* dbhi);
    
    void compute_zeta_by_Temp(const int* tlo, const int* thi,
                              const amrex_real* molarconc, const int* mclo, const int* mchi,
                              const amrex_real* D_bar, const int* dblo, const int* dbhi,
                              const amrex_real* Temp, const int* tplo, const int* tphi,
                              const amrex_real* zeta_by_Temp, const int* ztlo, const int* zthi,
                              const amrex_real* D_therm, const int* dtlo, const int* dthi);

    void compute_sqrtLonsager_fc(const int* tlo, const int* thi,
                                 const amrex_real* rho,            const int* rhlo, const int* rhhi,
                                 const amrex_real* rhotot,         const int* rtlo, const int* rthi,
                                 const amrex_real* sqrtLonsager_x, const int* Lxlo, const int* Lxhi,
                                 const amrex_real* sqrtLonsager_y, const int* Lylo, const int* Lyhi,
#if (AMREX_SPACEDIM == 3)
                                 const amrex_real* sqrtLonsager_z, const int* Lzlo, const int* Lzhi,
#endif
                                 const amrex_real* dx);


#ifdef __cplusplus
}
#endif

#endif
