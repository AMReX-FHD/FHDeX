#ifndef _particle_functions_F_H_
#define _particle_functions_F_H_

#ifdef __cplusplus
extern "C" 
{
#endif

// in particle_fucntions.F90

    void move_particles(void* particles, const int* np,
                              const amrex::Real* dt,
                              const amrex::Real* problo,
                              const amrex::Real* problhi);



#ifdef __cplusplus
}
#endif

#endif
