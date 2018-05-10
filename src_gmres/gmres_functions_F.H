#ifndef _gmres_functions_F_H_
#define _gmres_functions_F_H_

#ifdef __cplusplus
extern "C" 
{
#endif

    // these are interfaces to fortran subroutines

///////////////////////////
// in gmres_namelist.F90

    // read in fortran namelist into gmres_params_module
    void read_gmres_namelist (const char* inputs_file, const int length);

    // copy contents of gmres_params_module to C++ gmres namespace
    void initialize_gmres_namespace (int* precon_type, int* visc_schur_approx, 
                                     amrex::Real* p_norm_weight, amrex::Real* scale_factor,
                                     int* mg_verbose, int* cg_verbose, int* mg_max_vcycles,
                                     int* mg_minwidth, int* mg_bottom_solver,
                                     int* mg_nsmooths_down, int* mg_nsmooths_up,
                                     int* mg_nsmooths_bottom, int* mg_max_bottom_nlevels,
                                     amrex::Real* mg_rel_tol, amrex::Real* mg_abs_tol,
                                     int* stag_mg_verbosity, int* stag_mg_max_vcycles,
                                     int* stag_mg_minwidth, int* stag_mg_bottom_solver,
                                     int* stag_mg_nsmooths_down, int* stag_mg_nsmooths_up,
                                     int* stag_mg_nsmooths_bottom,
                                     int* stag_mg_max_bottom_nlevels, 
                                     amrex::Real* stag_mg_omega,
                                     int* stag_mg_smoother, amrex::Real* stag_mg_rel_tol,
                                     amrex::Real* gmres_rel_tol, amrex::Real* gmres_abs_tol,
                                     int* gmres_verbose, int* gmres_max_outer,
                                     int* gmres_max_inner, int* gmres_max_iter,
                                     int* gmres_min_iter, int* gmres_spatial_order);
///////////////////////////

///////////////////////////
// in stag_mg_solver.F90

    void cc_restriction (const int* lo_c, const int* hi_c,
                         const amrex::Real* phi_c,
                         const int* c_lo, const int* c_hi,
                         const amrex::Real* phi_f,
                         const int* f_lo, const int* f_hi);

    void stag_restriction (const int* lo_c, const int* hi_c,
                           const amrex::Real* phi_cx,
                           const int* cx_lo, const int* cx_hi,
                           const amrex::Real* phix_f,
                           const int* fx_lo, const int* fx_hi,
                           const amrex::Real* phi_cy,
                           const int* cy_lo, const int* cy_hi,
                           const amrex::Real* phiy_f,
                           const int* fy_lo, const int* fy_hi,
#if (AMREX_SPACEDIM == 3)
                           const amrex::Real* phi_cz,
                           const int* cz_lo, const int* cz_hi,
                           const amrex::Real* phiz_f,
                           const int* fz_lo, const int* fz_hi,
#endif
                           const int* simple_stencil);

    void nodal_restriction (const int* lo_c, const int* hi_c,
                            const amrex::Real* phi_c,
                            const int* c_lo, const int* c_hi,
                            const amrex::Real* phi_f,
                            const int* f_lo, const int* f_hi);

    void edge_restriction (const int* lo_c, const int* hi_c,
                           const amrex::Real* phi_cxy,
                           const int* cxy_lo, const int* cxy_hi,
                           const amrex::Real* phixy_f,
                           const int* fxy_lo, const int* fxy_hi,
                           const amrex::Real* phi_cxz,
                           const int* cxz_lo, const int* cxz_hi,
                           const amrex::Real* phixz_f,
                           const int* fxz_lo, const int* fxz_hi,
                           const amrex::Real* phi_cyz,
                           const int* cyz_lo, const int* cyz_hi,
                           const amrex::Real* phiyz_f,
                           const int* fyz_lo, const int* fyz_hi);

    void stag_prolongation (const int* lo_f, const int* hi_f,
                            const amrex::Real* phi_cx,
                            const int* cx_lo, const int* cx_hi,
                            const amrex::Real* phix_f,
                            const int* fx_lo, const int* fx_hi,
                            const amrex::Real* phi_cy,
                            const int* cy_lo, const int* cy_hi,
                            const amrex::Real* phiy_f,
                            const int* fy_lo, const int* fy_hi
#if (AMREX_SPACEDIM == 3)
                          , const amrex::Real* phi_cz,
                            const int* cz_lo, const int* cz_hi,
                            const amrex::Real* phiz_f,
                            const int* fz_lo, const int* fz_hi
#endif
                            );

///////////////////////////



#ifdef __cplusplus
}
#endif

#endif
