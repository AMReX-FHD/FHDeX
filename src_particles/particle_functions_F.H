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
                        const amrex_real* dx,
                        const amrex_real* dt);

    void move_particles_fhd(void* particles,
                        const int *np,
                        const int *lo,
                        const int *hi,
                        int** c_vectors,
                        int* sizes,
                        const int* clo,
                        const int* chi,
                        const amrex_real* plo,
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
                         const amrex_real* sourcey, const int* sourceylo, const int* sourceyhi
#if (AMREX_SPACEDIM == 3)
                         , const amrex_real* sourcez, const int* sourcezlo, const int* sourcezhi
#endif
                        );


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
                        amrex_real* members, const int* mlo, const int* mhi,
                        amrex_real* density, const int* dlo, const int* dhi,
                        amrex_real* velx, const int* velxlo, const int* velxhi,
                        amrex_real* vely, const int* velylo, const int* velyhi,
                        amrex_real* velz, const int* velzlo, const int* velzhi,
                        amrex_real* temp, const int* templo, const int* temphi,
                        amrex_real* px, const int* pxlo, const int* pxhi, 
                        amrex_real* py, const int* pylo, const int* pyhi, 
                        amrex_real* pz, const int* pzlo, const int* pzhi, 
                        amrex_real* energy, const int* energylo, const int* energyhi,    
                        amrex_real* cellVols, const int* cvlo, const int* cvhi,
                        const amrex_real* neff, const int* np);

   void initialize_fields(void* particles,
                        const int* lo, const int* hi,
                        int** c_vectors, int* sizes, const int* clo, const int* chi,
                        amrex_real* members, const int* mlo, const int* mhi,
                        amrex_real* density, const int* dlo, const int* dhi,
                        amrex_real* velx, const int* velxlo, const int* velxhi,
                        amrex_real* vely, const int* velylo, const int* velyhi,
                        amrex_real* velz, const int* velzlo, const int* velzhi,
                        amrex_real* temp, const int* templo, const int* temphi,
                        amrex_real* cellVols, const int* cvlo, const int* cvhi,
                        const amrex_real* neff, const int* np, const amrex_real* r, const amrex_real* t);



   void evaluate_stats(void* particles,
                        const int* lo, const int* hi,
                        int** c_vectors, int* sizes, const int* clo, const int* chi,
                        amrex_real* members, const int* mlo, const int* mhi,
                        amrex_real* density, const int* dlo, const int* dhi,
                        amrex_real* velx, const int* velxlo, const int* velxhi,
                        amrex_real* vely, const int* velylo, const int* velyhi,
                        amrex_real* velz, const int* velzlo, const int* velzhi,
                        amrex_real* temp, const int* templo, const int* temphi,
                        amrex_real* px, const int* pxlo, const int* pxhi, 
                        amrex_real* py, const int* pylo, const int* pyhi, 
                        amrex_real* pz, const int* pzlo, const int* pzhi, 
                        amrex_real* energy, const int* energylo, const int* energyhi,    


                        amrex_real* membersmean, const int* mmlo, const int* mmhi,
                        amrex_real* densitymean, const int* dmlo, const int* dmhi,
                        amrex_real* velxmean, const int* velxmlo, const int* velxmhi,
                        amrex_real* velymean, const int* velymlo, const int* velymhi,
                        amrex_real* velzmean, const int* velzmlo, const int* velzmhi,
                        amrex_real* tempmean, const int* tempmlo, const int* tempmhi,
                        amrex_real* pxmean, const int* pxmlo, const int* pxmhi, 
                        amrex_real* pymean, const int* pymlo, const int* pymhi, 
                        amrex_real* pzmean, const int* pzmlo, const int* pzmhi, 
                        amrex_real* energymean, const int* energymlo, const int* energymhi,    


                        amrex_real* membersvar, const int* mvlo, const int* mvhi,
                        amrex_real* densityvar, const int* dvlo, const int* dvhi,
                        amrex_real* velxvar, const int* velxvlo, const int* velxvhi,
                        amrex_real* velyvar, const int* velyvlo, const int* velyvhi,
                        amrex_real* velzvar, const int* velzvlo, const int* velzvhi,
                        amrex_real* tempvar, const int* tempvlo, const int* tempvhi,
                        amrex_real* pxvar, const int* pxvlo, const int* pxvhi, 
                        amrex_real* pyvar, const int* pyvlo, const int* pyvhi, 
                        amrex_real* pzvar, const int* pzvlo, const int* pzvhi, 
                        amrex_real* energyvar, const int* energyvlo, const int* energyvhi, 

                        amrex_real* gvar, const int* gvlo, const int* gvhi, 
                        amrex_real* kgcross, const int* kgclo, const int* kgchi, 
                        amrex_real* krcross, const int* krclo, const int* krchi, 
                        amrex_real* rgcross, const int* rgclo, const int* rgchi, 

                        amrex_real* cellVols, const int* cvlo, const int* cvhi,
                        const int* np, const amrex_real* neff, const amrex_real* n0, const amrex_real* T0, const amrex_real* delt, int* steps);


#ifdef __cplusplus
}
#endif

#endif


