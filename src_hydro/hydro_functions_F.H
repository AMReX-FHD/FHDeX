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
			      const amrex_real* dx, const int* increment=0);

///////////////////////////

///////////////////////////
// in mk_advective_s_fluxdiv.F90

  void mk_advective_s_fluxdiv(const int* lo, const int* hi,
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
			      const amrex_real* m, const int* mlo, const int* mhi,
			      const amrex_real* m_update, const int* mudlo, const int* mudhi,
			      const amrex_real* dx, const int* increment=0);

///////////////////////////


///////////////////////////
// in stoch_m_force.F90

  void stoch_m_force(const int* lo, const int* hi,
		     const amrex_real* flux_cc, const int* fluxcclo, const int* fluxcchi, const int* ncomp_cc,
#if (AMREX_SPACEDIM == 2)
		     const amrex_real* flux_nd, const int* fluxndlo, const int* fluxndhi, const int* ncomp_nd,		     	     
#endif
#if (AMREX_SPACEDIM == 3)
		     const amrex_real* flux_xy, const int* fluxxylo, const int* fluxxyhi, const int* ncomp_xy,
		     const amrex_real* flux_xz, const int* fluxxzlo, const int* fluxxzhi, const int* ncomp_xz,
		     const amrex_real* flux_yz, const int* fluxyzlo, const int* fluxyzhi, const int* ncomp_yz,		     	     
#endif
		     const amrex_real* divx, const int* divxlo, const int* divxhi,
		     const amrex_real* divy, const int* divylo, const int* divyhi,
#if (AMREX_SPACEDIM == 3)
		     const amrex_real* divz, const int* divzlo, const int* divzhi,
#endif
		     const amrex_real* dx, const int* increment=0);


  void mult_by_sqrt_eta_temp(const int* lo, const int* hi,
			     const amrex_real* flux_cc, const int* fluxcclo, const int* fluxcchi, const int* ncomp_cc,
#if (AMREX_SPACEDIM == 2)
			     const amrex_real* flux_nd, const int* fluxndlo, const int* fluxndhi, const int* ncomp_nd,		     	     
#endif
#if (AMREX_SPACEDIM == 3)
			     const amrex_real* flux_xy, const int* fluxxylo, const int* fluxxyhi, const int* ncomp_xy,
			     const amrex_real* flux_xz, const int* fluxxzlo, const int* fluxxzhi, const int* ncomp_xz,
			     const amrex_real* flux_yz, const int* fluxyzlo, const int* fluxyzhi, const int* ncomp_yz,		     	     
#endif
			     const amrex_real* eta, const int* etacclo, const int* etacchi,
#if (AMREX_SPACEDIM == 2)
			     const amrex_real* eta_nodal, const int* etandlo, const int* etandhi,
#endif
#if (AMREX_SPACEDIM == 3)
			     const amrex_real* eta_xy, const int* etaxylo, const int* etaxyhi,
			     const amrex_real* eta_xz, const int* etaxzlo, const int* etaxzhi,
			     const amrex_real* eta_yz, const int* etayzlo, const int* etayzhi,		     	     
#endif
			     const amrex_real* temperature, const int* tempcclo, const int* tempcchi,
#if (AMREX_SPACEDIM == 2)
			     const amrex_real* temperature_nodal, const int* tempndlo, const int* tempndhi);
#endif
#if (AMREX_SPACEDIM == 3)
  const amrex_real* temperature_xy, const int* tempxylo, const int* tempxyhi,
    const amrex_real* temperature_xz, const int* tempxzlo, const int* tempxzhi,
    const amrex_real* temperature_yz, const int* tempyzlo, const int* tempyzhi);
#endif

///////////////////////////

///////////////////////////
// in struct_fact.F90

  void fft_shift(const int* lo, const int* hi,
		 const amrex_real* dft, const int* dftlo, const int* dfthi, const int* ncomp);

  void sqrt_mf(const int* lo, const int* hi,
	       const amrex_real* mf, const int* mflo, const int* mfhi, const int* ncomp);

///////////////////////////

///////////////////////////
// in mac_solver.F90
void mac_solver_rhs(const int* lev, const int* lo, const int* hi,
		    amrex::Real* newrhs, const int* nrhs_lo, const int* nrhs_hi, 
		    const amrex::Real* oldrhs, const int* orhs_lo, const int* orhs_hi,
		    const amrex::Real* uedge, const int* u_lo, const int* u_hi, 
#if (AMREX_SPACEDIM >= 2)
		    const amrex::Real* vedge, const int* v_lo, const int* v_hi, 
#if (AMREX_SPACEDIM == 3)
		    const amrex::Real* wedge, const int* w_lo, const int* w_hi, 
#endif
#endif
		    const amrex::Real* dx);

void mac_bcoef_face(const int* lev, const int* lo, const int* hi,
		    amrex::Real* xface, const int* x_lo, const int* x_hi,
#if (AMREX_SPACEDIM >= 2)
		    amrex::Real* yface, const int* y_lo, const int* y_hi, 
#if (AMREX_SPACEDIM == 3)
		    amrex::Real* zface, const int* z_lo, const int* z_hi,
#endif
#endif
		    const amrex::Real* rhocc, const int* r_lo, const int* r_hi);
///////////////////////////

#ifdef __cplusplus
}
#endif

#endif
