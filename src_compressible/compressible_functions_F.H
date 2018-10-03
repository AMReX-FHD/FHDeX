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
                              const amrex::Real* dx, const int* nvars, const int* nprimvars);

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
                              const amrex::Real* dx, const int* nvars, const int* nprimvars);

#ifdef __cplusplus
}
#endif

#endif
