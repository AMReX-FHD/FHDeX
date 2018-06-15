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
                         const amrex_real* coordsx, const int* coordsxlo, const int* coordsxhi,
                         const amrex_real* coordsy, const int* coordsylo, const int* coordsyhi,
#if (AMREX_SPACEDIM == 3)
                         const amrex_real* coordsz, const int* coordszlo, const int* coordszhi,
#endif
                         const amrex_real* beta, const int* betalo, const int* betahi,
                         const amrex_real* rho, const int* rholo, const int* rhohi,

                         const amrex_real* sourcex, const int* sourcexlo, const int* sourcexhi,
                         const amrex_real* sourcey, const int* sourceylo, const int* sourceyhi,
#if (AMREX_SPACEDIM == 3)
                         const amrex_real* sourcez, const int* sourcezlo, const int* sourcezhi,
#endif
                         const int* cellmembers, const int* cmlo, const int* cmhi,
                         const int* celllists, const int* cllo, const int* clhi, const amrex::Real* cdx, const int* hivect, const int* ppc
                        );


    void init_particles(void* particles, const int* np, const amrex::Real* dx,
                        const amrex::Real* real_lo, const amrex::Real* real_hi,
                        const int* cellmembers, const int* cmlo, const int* cmhi,
                        const int* celllists, const int* cllo, const int* clhi, const int* ppc);

    void insert_particles(void* particles, const int* np,
                        const int* cellmembers, const int* cmlo, const int* cmhi,
                        const int* celllists, const int* cllo, const int* clhi, const int* ppc);



#ifdef __cplusplus
}
#endif

#endif


