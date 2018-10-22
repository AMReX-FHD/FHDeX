#ifndef _particle_functions_F_H_
#define _particle_functions_F_H_

#ifdef __cplusplus

extern "C" 
{
#endif

// in particle_fucntions.F90


    void move_particles_dsmc(void* particles,
                        const int *np,
                        const int *lo,
                        const int *hi,
                        int** c_vectors,
                        int* sizes,
                        const int* clo,
                        const int* chi,
                        const amrex_real* plo,
                        const amrex_real* phi,
                        const amrex_real* dx,
                        const amrex_real* dt,
                        const void* surfaces,
                        const int* ns);

    void move_particles_fhd(void* particles,
                        const int *np,
                        const int *lo,
                        const int *hi,
                        int** c_vectors,
                        int* sizes,
                        const int* clo,
                        const int* chi,
                        const amrex_real* plo,
                        const amrex_real* phi,
                        const amrex_real* dx,
                        const amrex_real* dt,
                        const amrex_real* plof,
                        const amrex_real* dxf,
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
                        const void* surfaces,
                        const int* ns);


    void init_particles(void* particles, const int* np, const amrex::Real* dx,
                        const amrex::Real* real_lo, const amrex::Real* real_hi,
                        const int* cellmembers, const int* cmlo, const int* cmhi,
                        const int* celllists, const int* cllo, const int* clhi, const int* ppc);


// in collision_functions.F90

    void init_cells(void* particles,
                        const int* lo, const int* hi,
                        int** c_vectors, int* sizes, const int* clo, const int* chi,
                        amrex_real* cellpairs, const int* cplo, const int* cphi, 
                        amrex_real* cellfactor, const int* cflo, const int* cfhi,
                        amrex_real* cellVols, const int* cvlo, const int* cvhi,
                        const int* np, const amrex_real* neff, const amrex_real* cp, const amrex_real* d, const amrex_real* delt);

    void collide_cells(void* particles,
                        const int* lo, const int* hi,
                        int** c_vectors, int* sizes, const int* clo, const int* chi,
                        amrex_real* cellpairs, const int* cplo, const int* cphi, 
                        amrex_real* cellfactor, const int* cflo, const int* cfhi,
                        amrex_real* cellVols, const int* cvlo, const int* cvhi,
                        const int* np, const amrex_real* neff, const amrex_real* cp, const amrex_real* d, const amrex_real* delt);

   void evaluate_fields(void* particles,
                        const int* lo, const int* hi,
                        int** c_vectors, int* sizes, const int* clo, const int* chi,
                        amrex_real* instant, const int* ilo, const int* ihi,
                        amrex_real* cellVols, const int* cvlo, const int* cvhi,
                        const amrex_real* neff, const int* np, amrex_real* del1, amrex_real* del2, amrex_real* tm, amrex_real* te);

   void initialize_fields(void* particles,
                        const int* lo, const int* hi,
                        int** c_vectors, int* sizes, const int* clo, const int* chi,
                        amrex_real* instant, const int* plo, const int* phi,
                        amrex_real* cellVols, const int* cvlo, const int* cvhi,
                        const amrex_real* neff, const int* np, const amrex_real* r, const amrex_real* t);



   void evaluate_means(void* particles,
                        const int* lo, const int* hi,
                        int** c_vectors, int* sizes, const int* clo, const int* chi,
                        amrex_real* instant, const int* ilo, const int* ihi,
                        amrex_real* means, const int* mlo, const int* mhi,
                        amrex_real* vars, const int* vlo, const int* vhi,

                        amrex_real* membraneflux, const int* mflo, const int* mfhi,

                        amrex_real* cellVols, const int* cvlo, const int* cvhi,
                        const int* np, const amrex_real* neff, const amrex_real* n0, const amrex_real* T0, const amrex_real* delt, int* steps, amrex_real* del1, amrex_real* del2);

   void evaluate_corrs(void* particles,
                        const int* lo, const int* hi,
                        int** c_vectors, int* sizes, const int* clo, const int* chi,
                        amrex_real* instant, const int* ilo, const int* ihi,
                        amrex_real* means, const int* mlo, const int* mhi,
                        amrex_real* vars, const int* vlo, const int* vhi,

                        amrex_real* cellVols, const int* cvlo, const int* cvhi,
                        const int* np, const amrex_real* neff, const amrex_real* n0, const amrex_real* T0, const amrex_real* delt, int* steps, amrex_real* del1, amrex_real* del2, amrex_real* del3);


#ifdef __cplusplus
}
#endif

#endif


