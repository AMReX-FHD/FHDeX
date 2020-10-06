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
    void initialize_compressible_namespace(amrex::Real* bc_Yk, amrex::Real* bc_Xk,
                                           int* plot_means, int* plot_vars);

///////////////////////////

///////////////////////////
// in diff_flux.F90

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

    void set_bc(const int* lo, const int* hi,
                amrex::Real* cons,
                amrex::Real* prim);

    void setup_bc();

    void setup_cwall();

    void get_energy(amrex::Real* energy, amrex::Real* massvec, const amrex::Real* temp);

    void get_hc_gas();

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
