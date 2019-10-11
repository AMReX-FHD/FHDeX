#include <AMReX_REAL.H>

#ifdef __cplusplus

extern "C" 
{
#endif

///////////////////////////
// in rng_fucntions.F90

    void rng_initialize(int* fhd, int* particle, int* sel, int* theta, int* phi, int* general);
    void rng_init(int* fhd);
    void rng_checkpoint(int* step);
    void rng_restart(int* step, int* digits);
    void get_angles(double* costheta, double* sintheta, double* cosphi, double* sinphi);
    double get_uniform_func();
    double get_particle_normal_func();
    double get_fhd_normal_func();

///////////////////////////
// in multifab_fill_random.F90

    void multifab_fill_random(const int* lo, const int* hi,
			      const amrex_real* mf,
			      const int* mflo, const int* mfhi,
			      const int* ncomp, const int* comp);

    void multifab_fill_random_hack(const int* lo, const int* hi,
			      const amrex_real* mf,
			      const int* mflo, const int* mfhi,
			      const int* ncomp, const int* comp, const double* mytop, const double* mybottom);

#ifdef __cplusplus
}
#endif

