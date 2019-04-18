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

    void tag_interface_ib(
            int * iface,            const int * iflo,  const int * ifhi,
            const amrex_real * phi, const int * philo, const int * phihi,
            const int * tag,        const int * taglo, const int * taghi
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
