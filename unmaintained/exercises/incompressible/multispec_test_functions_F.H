#ifndef _multispec_test_functions_F_H_
#define _multispec_test_functions_F_H_

#ifdef __cplusplus
extern "C"
{
#endif

    // these are interfaces to fortran subroutines

///////////////////////////
// in init_s_vel.F90

    void init_s_vel(const int* lo, const int* hi,
                    const amrex_real* vel, const int* vello, const int* velhi,
                    const amrex_real* dx,
                    const amrex_real* reallo, const amrex_real* realhi);

///////////////////////////

///////////////////////////
// in

    void init_rho_and_umac(const int* lo, const int* hi,
                           const amrex_real* rho, const int* rlo, const int* rhi,
                           const int* nspecies,
                           const amrex_real* u, const int* ulo, const int* uhi,
                           const amrex_real* v, const int* vlo, const int* vhi,
#if(AMREX_SPACEDIM == 3)
                           const amrex_real* w, const int* wlo, const int* whi,
#endif
                           const amrex_real* dx,
                           const amrex_real* prob_lo, const amrex_real* prob_hi,
                           const amrex_real* reallo, const amrex_real* realhi);

///////////////////////////


#ifdef __cplusplus
}
#endif

#endif