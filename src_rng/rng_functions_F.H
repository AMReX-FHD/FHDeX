

#ifdef __cplusplus

extern "C" 
{
#endif

// in rng_fucntions.F90

    void rng_initialize(int* fhd, int* particle, int* sel, int* theta, int* phi, int* general);
    void get_angles(double* costheta, double* sintheta, double* cosphi, double* sinphi);

#ifdef __cplusplus
}
#endif

