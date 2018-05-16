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

    void stag_mg_update (const int* lo, const int* hi,
                         const amrex::Real* phix,
                         const int* px_lo, const int* px_hi,
                         const amrex::Real* phiy,
                         const int* py_lo, const int* py_hi,
#if (AMREX_SPACEDIM == 3)
                         const amrex::Real* phiz,
                         const int* pz_lo, const int* pz_hi,
#endif
                         const amrex::Real* rhsx,
                         const int* rx_lo, const int* rx_hi,
                         const amrex::Real* rhsy,
                         const int* ry_lo, const int* ry_hi,
#if (AMREX_SPACEDIM == 3)
                         const amrex::Real* rhsz,
                         const int* rz_lo, const int* rz_hi,
#endif
                         const amrex::Real* Lpx,
                         const int* Lx_lo, const int* Lx_hi,
                         const amrex::Real* Lpy,
                         const int* Ly_lo, const int* Ly_hi,
#if (AMREX_SPACEDIM == 3)
                         const amrex::Real* Lpz,
                         const int* Lz_lo, const int* Lz_hi,
#endif
                         const amrex::Real* alphax,
                         const int* ax_lo, const int* ax_hi,
                         const amrex::Real* alphay,
                         const int* ay_lo, const int* ay_hi,
#if (AMREX_SPACEDIM == 3)
                         const amrex::Real* alphaz,
                         const int* az_lo, const int* az_hi,
#endif
                         const amrex::Real* beta,
                         const int* b_lo, const int* b_hi,
#if (AMREX_SPACEDIM == 2)
                         const amrex::Real* beta_ed,
                         const int* w_lo, const int* w_hi,
#elif (AMREX_SPACEDIM == 3)
                         const amrex::Real* beta_xy,
                         const int* w_lo, const int* w_hi,
                         const amrex::Real* beta_xz,
                         const int* x_lo, const int* x_hi,
                         const amrex::Real* beta_yz,
                         const int* y_lo, const int* y_hi,
#endif
                         const amrex::Real* gamma,
                         const int* g_lo, const int* g_hi,
                         const amrex::Real* dx,
                         const int* color);

    void stag_apply_op  (const int* lo, const int* hi,
                         const amrex_real* betacc, const int* betacclo, const int* betacchi,
                         const amrex_real* gammacc, const int* gammacclo, const int* gammacchi,
#if (AMREX_SPACEDIM == 2)
                         const amrex_real* betanodal, const int* betanodallo, const int* betanodalhi,
                         const amrex_real* gammanodal, const int* gammanodallo, const int* gammanodalhi,
#endif
#if (AMREX_SPACEDIM == 3)
                         const amrex_real* betaxy, const int* betaxylo, const int* betaxyhi,
                         const amrex_real* betaxz, const int* betaxzlo, const int* betaxzhi,
                         const amrex_real* betayz, const int* betayzlo, const int* betayzhi,
                         const amrex_real* gammaxy, const int* gammaxylo, const int* gammaxyhi,
                         const amrex_real* gammaxz, const int* gammaxzlo, const int* gammaxzhi,
                         const amrex_real* gammayz, const int* gammayzlo, const int* gammayzhi,

#endif
                         const amrex_real* velxin, const int* velxinlo, const int* velxinhi,
                         const amrex_real* velyin, const int* velyinlo, const int* velyinhi,
        #if (AMREX_SPACEDIM == 3)
                         const amrex_real* velzin, const int* velzinlo, const int* velzinhi,
        #endif
                         const amrex_real* velxout, const int* velxoutlo, const int* velxouthi,
                         const amrex_real* velyout, const int* velyoutlo, const int* velyouthi,
        #if (AMREX_SPACEDIM == 3)
                         const amrex_real* velzout, const int* velzoutlo, const int* velzouthi,
        #endif
                         const amrex_real* alphax, const int* alphaxlo, const int* alphaxhi,
                         const amrex_real* alphay, const int* alphaylo, const int* alphayhi,
        #if (AMREX_SPACEDIM == 3)
                         const amrex_real* alphaz, const int* alphazlo, const int* alphazhi,
        #endif
                         const amrex_real* dx,
                         const int* viscType);


///////////////////////////



#ifdef __cplusplus
}
#endif

#endif
