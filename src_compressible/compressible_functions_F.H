#ifndef _compressible_functions_F_H_
#define _compressible_functions_F_H_

#include <AMReX_MultiFab.H>

#ifdef __cplusplus
extern "C"
{
#endif

    // these are interfaces to fortran subroutines

///////////////////////////
// in compressible_namelist.F90

    // read in fortran namelist into compressible_namelist_module
    void read_compressible_namelist (const char* inputs_file, const int length);

    // copy contents of compressible_params_module to C++ compressible namespace
    void initialize_compressible_namespace(int* plot_means, int* plot_vars, int* plot_covars, int* plot_cross);

///////////////////////////

///////////////////////////
// in diff_flux.F90

    void trans_coeffs(const int* lo, const int* hi,
                      const amrex::Real* prim,
                      amrex::Real* eta,
                      amrex::Real* zeta,
                      amrex::Real* kappa);

    void setup_bc();

    void setup_cwall();

    void get_hc_gas();

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
