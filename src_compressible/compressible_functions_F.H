#ifndef _compressible_functions_F_H_
#define compressible_functions_F_H_

#include <AMReX_MultiFab.H>

#ifdef __cplusplus
extern "C" 
{
#endif

    // these are interfaces to fortran subroutines

///////////////////////////
// in compressible_namelist.F90

  // read in fortran namelist into compressible_params_module
  void read_compressible_namelist (const char* inputs_file, const int length);

  // copy contents of compressible_params_module to C++ compressible namespace
  void initialize_compressible_namespace(amrex::Real* Yk_lo, amrex::Real* Yk_hi, 
					 amrex::Real* Xk_lo, amrex::Real* Xk_hi);

///////////////////////////

///////////////////////////
// in diff_flux.F90

  void diff_flux(const int* lo, const int* hi,
		 const amrex::Real* cons,
		 const amrex::Real* prim,
		 const amrex::Real* eta,
		 const amrex::Real* zeta,
		 const amrex::Real* kappa,
		 const amrex::Real* chi,
		 const amrex::Real* Dij,
		 amrex::Real* xflux,
		 amrex::Real* yfluz,
#if (AMREX_SPACEDIM == 3)
		 amrex::Real* zflux,
#endif
		 amrex::Real* cornux,
		 amrex::Real* cornvx,
		 amrex::Real* cornwx,
		 amrex::Real* cornuy,
		 amrex::Real* cornvy,
		 amrex::Real* cornwy,
		 amrex::Real* cornuz,
		 amrex::Real* cornvz,
		 amrex::Real* cornwz,
		 amrex::Real* visccorn,
		 const amrex::Real* dx);

    void hyp_flux_prim(const int* lo, const int* hi,
		       const amrex::Real* cons,
		       const amrex::Real* prim,
		       amrex::Real* xflux,
		       amrex::Real* yfluz,
#if (AMREX_SPACEDIM == 3)
		       amrex::Real* zflux,
#endif
		       const amrex::Real* dx);

    void hyp_flux_cons(const int* lo, const int* hi,
		       const amrex::Real* cons,
		       const amrex::Real* prim,
		       amrex::Real* xflux,
		       amrex::Real* yfluz,
#if (AMREX_SPACEDIM == 3)
		       amrex::Real* zflux,
#endif
		       const amrex::Real* dx);

    void stoch_flux(const int* lo, const int* hi,
		    const amrex::Real* cons,
		    const amrex::Real* prim,
		    amrex::Real* xflux,
		    amrex::Real* yfluz,
#if (AMREX_SPACEDIM == 3)
		    amrex::Real* zflux,
#endif
		    amrex::Real* xsflux,
		    amrex::Real* ysfluz,
#if (AMREX_SPACEDIM == 3)
		    amrex::Real* zsflux,
#endif
		    const amrex::Real* rancorn,
		    const amrex::Real* eta,
		    const amrex::Real* zeta,
		    const amrex::Real* zappa,
		    const amrex::Real* chi,
		    const amrex::Real* Dij,
		    const amrex::Real* dx, const amrex::Real* dt);

  void trans_coeffs(const int* lo, const int* hi,
		    const amrex::Real* prim,
		    amrex::Real* eta,
		    amrex::Real* zeta,
		    amrex::Real* kappa);

  void makecoef(const int* lo, const int* hi,
		const amrex::Real* prim,
		amrex::Real* eta,
		amrex::Real* zeta,
		amrex::Real* kappa,
		amrex::Real* chi,
		amrex::Real* Dij);

  void rk3_stage1(const int* lo, const int* hi,
			      amrex::Real* cu,
			      amrex::Real* cup,
			      const amrex::Real* source,
			      const amrex::Real* xflux,
			      const amrex::Real* yfluz,
#if (AMREX_SPACEDIM == 3)
			      const amrex::Real* zflux,
#endif
                              const amrex::Real* dx,
                              const amrex::Real* dt);

  void rk3_stage2(const int* lo, const int* hi,
			      const amrex::Real* cu,
			      amrex::Real* cup,
			      const amrex::Real* cup2,
			      const amrex::Real* source,
			      const amrex::Real* xflux,
			      const amrex::Real* yfluz,
#if (AMREX_SPACEDIM == 3)
			      const amrex::Real* zflux,
#endif
                              const amrex::Real* dx,
                              const amrex::Real* dt);

  void rk3_stage3(const int* lo, const int* hi,
			      amrex::Real* cu,
			      const amrex::Real* cup,
			      const amrex::Real* cup2,
			      const amrex::Real* source,
			      const amrex::Real* xflux,
			      const amrex::Real* yfluz,
#if (AMREX_SPACEDIM == 3)
			      const amrex::Real* zflux,
#endif
                              const amrex::Real* dx,
                              const amrex::Real* dt);

  void cons_to_prim(const int* lo, const int* hi,
			      const amrex::Real* cons,
			      amrex::Real* prim);

  void prim_to_cons(const int* lo, const int* hi,
			      amrex::Real* cons,
			      const amrex::Real* prim);

  void set_bc(const int* lo, const int* hi,
	      amrex::Real* cons,
	      amrex::Real* prim,
	      amrex::Real* eta,
	      amrex::Real* zeta,
	      amrex::Real* kappa,
	      amrex::Real* chi,
	      amrex::Real* Dij);

  void get_energy(amrex::Real* energy, amrex::Real* massvec, const amrex::Real* temp);

  void get_pressure_gas(amrex::Real* pressure, amrex::Real* fracvec, amrex::Real* density, amrex::Real* temp);
  void get_density_gas(amrex::Real* pressure,amrex::Real* density, amrex::Real* temp);

  void get_hc_gas();

  void evaluate_means(const int* lo, const int* hi, const amrex::Real* cu, amrex::Real* cumeans, const amrex::Real* prim, amrex::Real* primmeans, const int* steps, amrex::Real* delHolder1, amrex::Real* delHolder2, amrex::Real* delHolder3, amrex::Real* delHolder4, amrex::Real* delHolder5, amrex::Real* delHolder6, amrex::Real* totalmass);
  void evaluate_corrs(const int* lo, const int* hi, const amrex::Real* cu, const amrex::Real* cumeans, amrex::Real* cuvars, const amrex::Real* prim, const amrex::Real* primmeans, amrex::Real* primvars, amrex::Real* spatialcross, const int* steps, amrex::Real* delHolder1, amrex::Real* delHolder2, amrex::Real* delHolder3, amrex::Real* delHolder4, amrex::Real* delHolder5, amrex::Real* delHolder6);


  void multifab_yzav(const int* lo, const int* hi, const amrex::Real* fabin, amrex::Real* fabout, const int* comps);

  void do_ssa(const int* lo, const int* hi,
			      const amrex::Real* cu,
			      const amrex::Real* prim,
			      amrex::Real* xflux,
                              const amrex::Real* dx,
                              const amrex::Real* dt);

  void do_langevin(const int* lo, const int* hi,
			      const amrex::Real* cu,
			      const amrex::Real* prim,
			      amrex::Real* xflux,
                              const amrex::Real* dx,
                              const amrex::Real* dt);

  void apply_effusion(const int* lo, const int* hi,
			      amrex::Real* cu,
			      const amrex::Real* xflux,
                              const amrex::Real* dx,
                              const amrex::Real* dt);

#ifdef __cplusplus
}
#endif

#endif
