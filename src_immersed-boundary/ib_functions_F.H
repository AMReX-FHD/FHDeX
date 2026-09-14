#ifndef _ib_functions_F_H_
#define _ib_functions_F_H_

#include <AMReX_REAL.H>

#ifdef BL_SINGLE_PRECISION_PARTICLES
typedef float amrex_particle_real;
#else
typedef double amrex_particle_real;
#endif

#ifdef __cplusplus
extern "C"
{
#endif


    //__________________________________________________________________________
    // IN: src_F90/immbdy_namelist.F90

    // read fortran namelist into immbdy_namelist_module
    void read_immbdy_namelist(const char * inputs_file, const int length);

    void chlamy_flagellum_datafile_len(int * len_buffer);

    void flagellum_max_markers(int * max_n_makers);

    void initialize_immbdy_namespace(int * n_immbdy, int * contains_flagellum,
                                     int * contains_fourier, int * contains_colloid);

    void initialize_ib_flagellum_namespace(
            int          n_immbdy,   int *        n_markers, void *       offset_0,
            amrex_real * amplitude,  amrex_real * frequency, amrex_real * length,
            amrex_real * wavelength, amrex_real * k_spring,  amrex_real * k_driving,
            int * fourier_coef_len,  const char * chlamy_flagellum_datafile
        );

    void copy_ib_fourier_data(int i_marker, int i_immbdy,
                              int * N, void * a_coef, void * b_coef);

    void initialize_ib_colloid_namespace(
            int          n_immbdy, int *        n_markers, void *       center,
            amrex_real * radius,   amrex_real * rho,       amrex_real * k_spring
        );

    void test_interface(
            void * particle_info, int * np
        );

    void fill_levelset_ib(
            const int  * lo,            const int * hi,
            const void * particle_info, const int * np,
            amrex_real * phi,           const int * philo, const int * phihi,
            int        * tag,           const int * taglo, const int * taghi,
            amrex_real * vel,           const int * vello, const int * velhi,
            const amrex_real * dx
        );

    void fill_levelset_sphere(
            const int  * lo,            const int * hi,
            const void * particle_info,
            amrex_real * phi,           const int * philo, const int * phihi,
            const amrex_real * dx
        );
    void tag_catalyst_interface(
            const int  * lo,            const int * hi,
            const void * particle_info, int * np,
            const int  * iface,           const int * iflo, const int * ifhi,
            const int  * ctag,            const int * ctaglo, const int * ctaghi,
            const amrex_real * dx
        );

    void tag_interface_ib(
            int * iface,            const int * iflo,  const int * ifhi,
            const amrex_real * phi, const int * philo, const int * phihi,
            const int * tag,        const int * taglo, const int * taghi
        );

    void spread_markers(
            const int * lo,             const int * hi,
            const int * tile_lo,        const int * tile_hi,
                  amrex_real * mf_x,    const int * mfx_lo, const int * mfx_hi,
                  amrex_real * mf_y,    const int * mfy_lo, const int * mfy_hi,
                  amrex_real * mf_z,    const int * mfz_lo, const int * mfz_hi,
                  amrex_real * weightx, const int * wx_lo,  const int * wx_hi,
                  amrex_real * weighty, const int * wy_lo,  const int * wy_hi,
                  amrex_real * weightz, const int * wz_lo,  const int * wz_hi,
            const amrex_real * coord_x, const int * cx_lo,  const int * cx_hi,
            const amrex_real * coord_y, const int * cy_lo,  const int * cy_hi,
            const amrex_real * coord_z, const int * cz_lo,  const int * cz_hi,
            const void *       pos_marker,
            const void *       v_marker,
            const int *        n_marker,
            const amrex_real * dx, const int * ghost,
            const int * pkernel_fluid_in
        );

    void inv_interpolate_markers(
            const int * lo,             const int * hi,
            const int * tile_lo,        const int * tile_hi,
                  amrex_real * mf_x,    const int * mfx_lo, const int * mfx_hi,
                  amrex_real * mf_y,    const int * mfy_lo, const int * mfy_hi,
                  amrex_real * mf_z,    const int * mfz_lo, const int * mfz_hi,
            const amrex_real * coord_x, const int * cx_lo,  const int * cx_hi,
            const amrex_real * coord_y, const int * cy_lo,  const int * cy_hi,
            const amrex_real * coord_z, const int * cz_lo,  const int * cz_hi,
            const void *       pos_marker,
            const void *       v_marker,
            const int *        n_marker,
            const amrex_real * dx,
            const int * pkernel_fluid_in
        );

    void interpolate_markers(
            const int * lo,             const int * hi,
            const int * tile_lo,        const int * tile_hi,
            const amrex_real * mf_x,    const int * mfx_lo, const int * mfx_hi,
            const amrex_real * mf_y,    const int * mfy_lo, const int * mfy_hi,
            const amrex_real * mf_z,    const int * mfz_lo, const int * mfz_hi,
            const amrex_real * weightx, const int * wx_lo,  const int * wx_hi,
            const amrex_real * weighty, const int * wy_lo,  const int * wy_hi,
            const amrex_real * weightz, const int * wz_lo,  const int * wz_hi,
            const amrex_real * coord_x, const int * cx_lo,  const int * cx_hi,
            const amrex_real * coord_y, const int * cy_lo,  const int * cy_hi,
            const amrex_real * coord_z, const int * cz_lo,  const int * cz_hi,
            const void *       pos_marker,
                  void *       v_marker,
            const int *        n_marker,
            const amrex_real * dx,
            const int * pkernel_fluid_in
        );

    void inv_spread_markers(
            const int * lo,             const int * hi,
            const int * tile_lo,        const int * tile_hi,
            const amrex_real * mf_x,    const int * mfx_lo, const int * mfx_hi,
            const amrex_real * mf_y,    const int * mfy_lo, const int * mfy_hi,
            const amrex_real * mf_z,    const int * mfz_lo, const int * mfz_hi,
            const amrex_real * coord_x, const int * cx_lo,  const int * cx_hi,
            const amrex_real * coord_y, const int * cy_lo,  const int * cy_hi,
            const amrex_real * coord_z, const int * cz_lo,  const int * cz_hi,
            const void *       pos_marker,
                  void *       v_marker,
            const int *        n_marker,
            const amrex_real * dx,
            const int * pkernel_fluid_in
        );

    void interpolate_ib_staggered(
            const int * lo,         const int * hi,
                  amrex_real * u_d, const int * udlo,  const int * udhi,
                  amrex_real * v_d, const int * vdlo,  const int * vdhi,
                  amrex_real * w_d, const int * wdlo,  const int * wdhi,
            const amrex_real * u_s, const int * uslo,  const int * ushi,
            const amrex_real * v_s, const int * vslo,  const int * vshi,
            const amrex_real * w_s, const int * wslo,  const int * wshi,
                  int *        et,  const int * etlo,  const int * ethi,
            const amrex_real * phi, const int * philo, const int * phihi,
            const int *        tag, const int * taglo, const int * taghi,
            const amrex_real * vel, const int * vello, const int * velhi
        );

    void interpolate_ib_cc(
            const int * lo,         const int * hi,
                  amrex_real * u_d, const int * udlo,  const int * udhi,
            const amrex_real * u_s, const int * uslo,  const int * ushi,
                  int *        et,  const int * etlo,  const int * ethi,
            const amrex_real * phi, const int * philo, const int * phihi,
            const int *        tag, const int * taglo, const int * taghi,
            const amrex_real * vel, const int * vello, const int * velhi
        );

    void fill_fgds_ib(
            const int * lo,         const int * hi,
                  amrex_real * f_u, const int * fulo, const int * fuhi,
                  amrex_real * f_v, const int * fvlo, const int * fvhi,
                  amrex_real * f_w, const int * fwlo, const int * fwhi,
            const int *        et,  const int * etlo, const int * ethi
        );

    void fill_force_ib_staggered(
            const int * lo,         const int * hi,
                  amrex_real * f_u, const int * fulo, const int * fuhi,
                  amrex_real * f_v, const int * fvlo, const int * fvhi,
                  amrex_real * f_w, const int * fwlo, const int * fwhi,
            const amrex_real * u_g, const int * uglo, const int * ughi,
            const amrex_real * v_g, const int * vglo, const int * vghi,
            const amrex_real * w_g, const int * wglo, const int * wghi,
            const amrex_real * u_d, const int * udlo, const int * udhi,
            const amrex_real * v_d, const int * vdlo, const int * vdhi,
            const amrex_real * w_d, const int * wdlo, const int * wdhi,
            const int *        et,  const int * etlo, const int * ethi,
            const amrex_real * dt
        );

    void fill_force_ib_cc(
            const int * lo,         const int * hi,
                  amrex_real * f_u, const int * fulo, const int * fuhi,
            const amrex_real * u_g, const int * uglo, const int * ughi,
            const amrex_real * u_d, const int * udlo, const int * udhi,
            const int *        et,  const int * etlo, const int * ethi,
            const amrex_real * dt
        );


    void set_particle_properties(
            const int * pstate, const amrex_real * pradius, const amrex_real * pdensity,
            amrex_real * pvol,  amrex_real * pmass, amrex_real * omoi, amrex_real * omega
        );



#ifdef __cplusplus
}
#endif

#endif
