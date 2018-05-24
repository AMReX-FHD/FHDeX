#ifndef _particle_functions_F_H_
#define _particle_functions_F_H_

#ifdef __cplusplus
extern "C" 
{
#endif

// in particle_fucntions.F90

    void update_particles(void* particles, const int* np, const amrex::Real* dt, const amrex::Real* dx,
                         const int* index_lo, const int* index_hi, const amrex::Real* real_lo, const amrex::Real* real_hi,
                         const amrex_real* velx, const int* velxlo, const int* velxhi,
                         const amrex_real* vely, const int* velylo, const int* velyhi,
#if (AMREX_SPACEDIM == 3)
                         const amrex_real* velz, const int* velzlo, const int* velzhi,
#endif
#if (AMREX_SPACEDIM == 2)
                         const amrex_real* betanodal, const int* betanodallo, const int* betanodalhi,
#endif
#if (AMREX_SPACEDIM == 3)
                         const amrex_real* betaxy, const int* betaxylo, const int* betaxyhi,
                         const amrex_real* betaxz, const int* betaxzlo, const int* betaxzhi,
                         const amrex_real* betayz, const int* betayzlo, const int* betayzhi,
#endif
                         const amrex_real* rho, const int* rholo, const int* rhohi,

                         const amrex_real* sourcex, const int* sourcexlo, const int* sourcexhi,
                         const amrex_real* sourcey, const int* sourceylo, const int* sourceyhi
#if (AMREX_SPACEDIM == 3)
                         , const amrex_real* sourcez, const int* sourcezlo, const int* sourcezhi
#endif
                        );



#ifdef __cplusplus
}
#endif

#endif


