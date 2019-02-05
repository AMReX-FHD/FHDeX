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
    const int* plot_stag,
    const int* is_nonisothermal,
    const int* is_ideal_mixture,
    const int* use_lapack,
    const amrex_real*c_init,
    const amrex_real* c_bc,
    const amrex_real* alpha1,
    const amrex_real* beta,
    const amrex_real* delta,
    const amrex_real* sigma,
    const int* midpoint_stoch_mass_flux_type,
    const int* avg_type,
    const int* mixture_type);

///////////////////////////


///////////////////////////
// in mass_flux_utilities.F90

  void compute_molconc_molmtot(const int* tlo, const int* thi,
  			       const amrex_real* rho, const int* rlo, const int* rhi, const int* nspecies,
  			       const amrex_real* rhotot, const int* rtlo, const int* rthi,
  			       const amrex_real* molarconc, const int* mclo, const int* mchi,
  			       const amrex_real* molmtot, const int* mtlo, const int* mthi);
  
  void compute_Gamma(const int* tlo, const int* thi,
		     const amrex_real* molarconc, const int* mclo, const int* mchi, const int* nspecies,
		     const amrex_real* Hessian, const int* hlo, const int* hhi, 
		     const amrex_real* Gamma, const int* glo, const int* ghi);

  void compute_rhoWchi(const int* tlo, const int* thi,
		       const amrex_real* rho, const int* rlo, const int* rhi, const int* nspecies,
		       const amrex_real* rhotot, const int* rtlo, const int* rthi,
		       const amrex_real* molarconc, const int* mclo, const int* mchi,
		       const amrex_real* rhoWchi, const int* rwclo, const int* rwchi,
		       const amrex_real* D_bar, const int* dblo, const int* dbhi);

  void compute_zeta_by_Temp(const int* tlo, const int* thi,
			    const amrex_real* molarconc, const int* mclo, const int* mchi, const int* nspecies,
			    const amrex_real* D_bar, const int* dblo, const int* dbhi,
			    const amrex_real* Temp, const int* tplo, const int* tphi,
			    const amrex_real* zeta_by_Temp, const int* ztlo, const int* zthi,
			    const amrex_real* D_therm, const int* dtlo, const int* dthi);

///////////////////////////

///////////////////////////
// in compute_mixture_properties.F90
  
  void  mixture_properties_mass(const int* tlo, const int* thi,
				const amrex_real* rho, const int* rlo, const int* rhi, const int* nspecies,
				const amrex_real* rhotot, const int* rtlo, const int* rthi,
				const amrex_real* D_bar, const int* dblo, const int* dbhi,
				const amrex_real* D_therm, const int* dtlo, const int* dthi,
				const amrex_real* Hessian, const int* hslo, const int* hshi);
///////////////////////////
  
///////////////////////////
// in mass_flux_utilities.F90

  void matvec_mul(const int* tlo, const int* thi,
		  const amrex_real* xp, const int* xlo, const int* xhi,
                  const amrex_real* ap, const int* alo, const int* ahi, const int* nc);

///////////////////////////

///////////////////////////
// in correction_flux.F90
  
  void correction_flux(const int* lo, const int* hi,
		       const amrex_real* rho, const int* rlo, const int* rhi, const int* nspecies,
		       const amrex_real* rhotot, const int* rtlo, const int* rthi,
		       const amrex_real* flux_x, const int* fxlo, const int* fxhi, 
		       const amrex_real* flux_y, const int* fylo, const int* fyhi
#if (AMREX_SPACEDIM == 3)
		       ,const amrex_real* flux_z, const int* fzlo, const int* fzhi
#endif		       
		       );

///////////////////////////

#ifdef __cplusplus
}
#endif

#endif
