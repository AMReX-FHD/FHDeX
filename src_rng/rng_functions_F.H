

#ifdef __cplusplus

extern "C" 
{
#endif

// in rng_fucntions.F90

    void rng_initialize(int* fhd, int* particle, int* sel, int* theta, int* phi, int* general);
    void get_angles(double* costheta, double* sintheta, double* cosphi, double* sinphi);
    double get_uniform_func();
    double get_particle_normal_func();

#ifdef __cplusplus
}
#endif

