#ifndef _hydro_functions_F_H_
#define _hydro_functions_F_H_

#ifdef __cplusplus
extern "C" 
{
#endif

    // these are interfaces to fortran subroutines

///////////////////////////
// in mk_advective_m_fluxdiv.F90

  void mk_advective_m_fluxdiv(const int* lo, const int* hi,
			      const amrex_real* umac, const int* umaclo, const int* umachi,
			      const amrex_real* vmac, const int* vmaclo, const int* vmachi,
#if (AMREX_SPACEDIM == 3)
			      const amrex_real* wmac, const int* wmaclo, const int* wmachi,
#endif
			      const amrex_real* mx, const int* mxlo, const int* mxhi,
			      const amrex_real* my, const int* mylo, const int* myhi,
#if (AMREX_SPACEDIM == 3)
			      const amrex_real* mz, const int* mzlo, const int* mzhi,
#endif
			      const amrex_real* m_updatex, const int* mudxlo, const int* mudxhi,
			      const amrex_real* m_updatey, const int* mudylo, const int* mudyhi,
#if (AMREX_SPACEDIM == 3)
			      const amrex_real* m_updatez, const int* mudzlo, const int* mudzhi,
#endif
			      const amrex_real* dx, const bool* increment);

///////////////////////////



#ifdef __cplusplus
}
#endif

#endif
