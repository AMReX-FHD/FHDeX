#ifndef _particle_functions_F_H_
#define _particle_functions_F_H_

#ifdef __cplusplus

extern "C" 
{
#endif

// in particle_fucntions.F90

		void repulsive_force(void* part1,
			void* part2,
			const amrex_real* dx,
			const amrex_real* dr);

		void force_function2(void* part1,
			void* part2,
                        const amrex_real* domsize);

                void near_wall_check(void* part, const int* dir, int* near_wall); 

                void calc_im_charge_loc(void* part, const int* dir, amrex_real* im_charge_pos); 

                void compute_p3m_force_mag(const amrex_real* r, amrex_real* mag, const amrex_real* dx); 

		void calculate_force(void* particles,
                        const int *np,
                        const int *lo,
                        const int *hi,
                        int** c_vectors,
                        int* sizes,
                        const int* clo,
                        const int* chi,
                        const amrex_real* plo,
                        const amrex_real* phi,
                        const amrex_real* dx);

                void amrex_compute_p3m_sr_correction_nl(void* particles, const int* np,
                        const void* ghosts,    const int* nn,
                        const int* neighbor_list, const int* size, amrex_real* rcount, const amrex_real* charge, const int* chargelo,
                        const int* chargehi, const amrex_real* coords, const int* coordslo, const int* coordshi, const int* lo, const int* hi, const amrex_real* dx);

                void amrex_compute_forces_nl(void* particles, const int* np,
                        const void* ghosts,    const int* nn,
                        const int* neighbor_list, const int* size, amrex_real* rcount);


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
                        const void* paramplanes,
			      const int* ns,
			      amrex_real* time,
                        int* flux);

    void move_ions_fhd(void* particles,
                       const int *np_tile,
                       const amrex_real *rejected_tile,
                       const amrex_real *moves_tile,
                       const amrex_real *maxspeed_tile,
                       const amrex_real *maxdist_tile,
                       const amrex_real *diffinst_tile,                   
                       const int *lo, const int *hi,
                       int** c_vectors,
                       int* sizes,
                       const int* clo,
                       const int* chi,
                       const amrex_real* plo, const amrex_real* phi,
                       const amrex_real* dx,
                       const amrex_real* dt,
                       const amrex_real* plof,
                       const amrex_real* dxf,
                       const amrex_real* dxe,
                       const amrex_real* velx, const int* velxlo, const int* velxhi,
                       const amrex_real* vely, const int* velylo, const int* velyhi,
#if (AMREX_SPACEDIM == 3)
                       const amrex_real* velz, const int* velzlo, const int* velzhi,
#endif
                       const amrex_real* efx, const int* efxlo, const int* efxhi,
                       const amrex_real* efy, const int* efylo, const int* efyhi,
#if (AMREX_SPACEDIM == 3)
                       const amrex_real* efz, const int* efzlo, const int* efzhi,
#endif
                       const amrex_real* coordsx, const int* coordsxlo, const int* coordsxhi,
                       const amrex_real* coordsy, const int* coordsylo, const int* coordsyhi,
#if (AMREX_SPACEDIM == 3)
                       const amrex_real* coordsz, const int* coordszlo, const int* coordszhi,
#endif
                       const amrex_real* sourcex, const int* sourcexlo, const int* sourcexhi,
                       const amrex_real* sourcey, const int* sourceylo, const int* sourceyhi,
#if (AMREX_SPACEDIM == 3)
                       const amrex_real* sourcez, const int* sourcezlo, const int* sourcezhi,
#endif
                       const amrex_real* mobility, const int* mlo, const int* mhi,
                       const void* paramplanes,
                       const int* ns, amrex_real * kinetic, int* sw);

    void spread_ions_fhd(void* particles,
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
                        const amrex_real* dxe,
                         const amrex_real* velx, const int* velxlo, const int* velxhi,
                         const amrex_real* vely, const int* velylo, const int* velyhi,
#if (AMREX_SPACEDIM == 3)
                         const amrex_real* velz, const int* velzlo, const int* velzhi,
#endif
                         const amrex_real* efx, const int* efxlo, const int* efxhi,
                         const amrex_real* efy, const int* efylo, const int* efyhi,
#if (AMREX_SPACEDIM == 3)
                         const amrex_real* efz, const int* efzlo, const int* efzhi,
#endif
                         const amrex_real* charge, const int* chargelo, const int* chargehi,
                         const amrex_real* coordsx, const int* coordsxlo, const int* coordsxhi,
                         const amrex_real* coordsy, const int* coordsylo, const int* coordsyhi,
#if (AMREX_SPACEDIM == 3)
                         const amrex_real* coordsz, const int* coordszlo, const int* coordszhi,
#endif
                         const amrex_real* cellcenters, const int* cellcenterslo, const int* cellcentershi,
                         const amrex_real* sourcex, const int* sourcexlo, const int* sourcexhi,
                         const amrex_real* sourcey, const int* sourceylo, const int* sourceyhi,
#if (AMREX_SPACEDIM == 3)
                        const amrex_real* sourcez, const int* sourcezlo, const int* sourcezhi,
#endif
                        const void* paramplanes,
                        const int* ns, amrex_real* potential, int* sw);

    void do_rfd(void* particles,
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
                        const amrex_real* dxe,
                         const amrex_real* velx, const int* velxlo, const int* velxhi,
                         const amrex_real* vely, const int* velylo, const int* velyhi,
#if (AMREX_SPACEDIM == 3)
                         const amrex_real* velz, const int* velzlo, const int* velzhi,
#endif
                         const amrex_real* efx, const int* efxlo, const int* efxhi,
                         const amrex_real* efy, const int* efylo, const int* efyhi,
#if (AMREX_SPACEDIM == 3)
                         const amrex_real* efz, const int* efzlo, const int* efzhi,
#endif
                         const amrex_real* coordsx, const int* coordsxlo, const int* coordsxhi,
                         const amrex_real* coordsy, const int* coordsylo, const int* coordsyhi,
#if (AMREX_SPACEDIM == 3)
                         const amrex_real* coordsz, const int* coordszlo, const int* coordszhi,
#endif
                         const amrex_real* cellcenters, const int* cellcenterslo, const int* cellcentershi,
                         const amrex_real* sourcex, const int* sourcexlo, const int* sourcexhi,
                         const amrex_real* sourcey, const int* sourceylo, const int* sourceyhi,
#if (AMREX_SPACEDIM == 3)
                        const amrex_real* sourcez, const int* sourcezlo, const int* sourcezhi,
#endif
                        const void* paramplanes,
                        const int* ns, int* sw);

   void collect_charge(void* particles,
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
                        const amrex_real* ploes,
                        const amrex_real* dxes,
                        const amrex_real* cellcenters, const int* cellcenterslo, const int* cellcentershi,

                         const amrex_real* charge, const int* chargelo, const int* chargehi);

   void compute_dry_mobility( const int *lo,
                        const int *hi,
                        amrex_real* mobility,
                        const int *mlo,
                        const int *mhi,
                        const amrex_real* dx,
                        const amrex_real* plo,
                        const amrex_real* phi,
                        const int* ngc, const species* speciesInfo);


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

   void evaluate_fields_mem(void* particles,
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



   void evaluate_means_mem(void* particles,
                        const int* lo, const int* hi,
                        int** c_vectors, int* sizes, const int* clo, const int* chi,
                        amrex_real* instant, const int* ilo, const int* ihi,
                        amrex_real* means, const int* mlo, const int* mhi,
                        amrex_real* vars, const int* vlo, const int* vhi,

                        amrex_real* membraneflux, const int* mflo, const int* mfhi,

                        amrex_real* cellVols, const int* cvlo, const int* cvhi,
                        const int* np, const amrex_real* neff, const amrex_real* n0, const amrex_real* T0, const amrex_real* delt, int* steps, amrex::Real* delHolder1, amrex::Real* delHolder2, amrex::Real* delHolder3, amrex::Real* delHolder4, amrex::Real* delHolder5, amrex::Real* delHolder6, amrex::Real* totalmass);

   void evaluate_corrs_mem(void* particles,
                        const int* lo, const int* hi,
                        int** c_vectors, int* sizes, const int* clo, const int* chi,
                        amrex_real* instant, const int* ilo, const int* ihi,
                        amrex_real* means, const int* mlo, const int* mhi,
                        amrex_real* vars, const int* vlo, const int* vhi,

                        amrex_real* cellVols, const int* cvlo, const int* cvhi,
                        const int* np, const amrex_real* neff, const amrex_real* n0, const amrex_real* T0, const amrex_real* delt, int* steps, amrex::Real* delHolder1, amrex::Real* delHolder2, amrex::Real* delHolder3, amrex::Real* delHolder4, amrex::Real* delHolder5, amrex::Real* delHolder6);


   void evaluate_fields(void* particles,
                        const int* lo, const int* hi,
                        int** c_vectors, int* sizes, const int* clo, const int* chi,
                        amrex_real* instant, const int* ilo, const int* ihi,
                        amrex_real* cellVols, const int* cvlo, const int* cvhi,
                        const amrex_real* neff, const int* np);



   void evaluate_means(const int* lo, const int* hi,
                       amrex_real* instant, const int* ilo, const int* ihi,
                       amrex_real* means, const int* mlo, const int* mhi,
                       amrex_real* vars, const int* vlo, const int* vhi,
                       amrex_real* cellVols, const int* cvlo, const int* cvhi,
                       const amrex_real* n0, const amrex_real* T0,
                       const amrex_real* delt, int* steps,
                       amrex_real* avcurrent);

   void evaluate_corrs(void* particles,
                       const int* lo, const int* hi,
                       int** c_vectors, int* sizes, const int* clo, const int* chi,
                       amrex_real* instant, const int* ilo, const int* ihi,
                       amrex_real* means, const int* mlo, const int* mhi,
                       amrex_real* vars, const int* vlo, const int* vhi,
                       amrex_real* cellVols, const int* cvlo, const int* cvhi,
                       const int* np, const amrex_real* neff, const amrex_real* n0, const amrex_real* T0, const amrex_real* delt, int* steps,
                       amrex_real* varcurrent);

   void sync_particles(amrex_real* spec3xPos, amrex_real* spec3yPos, amrex_real* spec3zPos, amrex_real* spec3xForce, amrex_real* spec3yForce,
                       amrex_real* spec3zForce, void* particles, const int* np, const int* length);

   void force_particles(amrex_real* spec3xPos, amrex_real* spec3yPos, amrex_real* spec3zPos, amrex_real* spec3xForce, amrex_real* spec3yForce,
                        amrex_real* spec3zForce, void* particles, const int* np, const int* length);


   void user_force_calc_init(const char* inputs_file, const int length);

   void user_force_calc(const amrex_real* spec3xPos, const amrex_real* spec3yPos, const amrex_real* spec3zPos, amrex_real* spec3xForce, amrex_real* spec3yForce, amrex_real* spec3zForce, const int* length, const int* step, const species* particleInfo);

   void user_force_calc_destroy();

   void get_explicit_mobility(amrex_real* mob, void* part, const amrex_real* plo, const amrex_real* phi);

   void dry(const amrex_real*dt,void* part,amrex_real* dry_terms, amrex_real* mb);

   void remove_particle_from_cell(int* cell_parts, int * cell_np, int* new_np, int* i);


#ifdef __cplusplus
}
#endif

#endif


