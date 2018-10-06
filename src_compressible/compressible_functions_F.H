#ifndef _compressible_functions_F_H_
#define compressible_functions_F_H_

#include <AMReX_MultiFab.H>

#ifdef __cplusplus
extern "C" 
{
#endif

    // these are interfaces to fortran subroutines

///////////////////////////
// in diff_flux.F90

  void diff_flux(const int* lo, const int* hi,
			      const amrex::Real* cons,
			      const amrex::Real* prim,
			      const amrex::Real* eta,
			      const amrex::Real* zeta,
			      const amrex::Real* zappa,

			      amrex::Real* xflux,
			      amrex::Real* yfluz,
#if (AMREX_SPACEDIM == 3)
			      amrex::Real* zflux,
#endif
                              const amrex::Real* dx);

  void hyp_flux(const int* lo, const int* hi,
			      const amrex::Real* cons,
			      const amrex::Real* prim,
			      const amrex::Real* eta,
			      const amrex::Real* zeta,
			      const amrex::Real* zappa,

			      amrex::Real* xflux,
			      amrex::Real* yfluz,
#if (AMREX_SPACEDIM == 3)
			      amrex::Real* zflux,
#endif
                              const amrex::Real* dx);

  void trans_coeffs(const int* lo, const int* hi,
			      const amrex::Real* prim,
			      amrex::Real* eta,
			      amrex::Real* zeta,
			      amrex::Real* zappa);


  void rk3_stage1(const int* lo, const int* hi,
			      const amrex::Real* cu,
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

#ifdef __cplusplus
}
#endif

#endif
