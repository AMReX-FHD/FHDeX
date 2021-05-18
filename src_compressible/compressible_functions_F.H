#ifndef _compressible_functions_F_H_
#define _compressible_functions_F_H_

#include <AMReX_MultiFab.H>

#ifdef __cplusplus
extern "C"
{
#endif

    // these are interfaces to fortran subroutines

///////////////////////////
// in membrane_funcs.F90

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
