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


    void stoch_flux(const int* lo, const int* hi,
                const amrex::Real* cons,
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
                const amrex::Real* dx, const amrex::Real* dt);

   void euler_step(const int* lo, const int* hi,
                  amrex::Real* cu,
                  const amrex::Real* xflux,
                  const amrex::Real* yfluz,
#if (AMREX_SPACEDIM == 3)
                  const amrex::Real* zflux,
#endif
                  const amrex::Real* dx,
                  const amrex::Real* dt);

  void rk3_stage1(const int* lo, const int* hi,
                  amrex::Real* cu,
                  amrex::Real* cup,
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
                  const amrex::Real* xflux,
                  const amrex::Real* yfluz,
#if (AMREX_SPACEDIM == 3)
                  const amrex::Real* zflux,
#endif
                  const amrex::Real* dx,
                  const amrex::Real* dt);

  void set_bc(const int* lo, const int* hi,
                  amrex::Real* cons,
                  amrex::Real* prim,
                  amrex::Real* eta,
                  amrex::Real* zeta,
                  amrex::Real* zappa);

  void evaluate_means(const int* lo, const int* hi, const amrex::Real* cu, amrex::Real* cumeans, const int* steps);
  void evaluate_corrs(const int* lo, const int* hi, const amrex::Real* cu, const amrex::Real* cumeans, amrex::Real* cuvars, const int* steps);



#ifdef __cplusplus
}
#endif

#endif
