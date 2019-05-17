#ifndef _GL_functions_F_H_
#define GL_functions_F_H_

#include <AMReX_MultiFab.H>

#ifdef __cplusplus
extern "C" 
{
#endif

  void rk2_stage1(const int* lo, const int* hi,
			      amrex::Real* phi,
			      amrex::Real* phin,
			      const amrex::Real* rannums,
                              amrex::Real* integral,
                              amrex::Real* energy,
                              amrex::Real* teng,
                              const amrex::Real* dx,
                              const amrex::Real* dt);

  void rk2_stage2(const int* lo, const int* hi,
			      amrex::Real* phi,
			      amrex::Real* phin,
			      const amrex::Real* rannums,
                              const amrex::Real* integral,
                              const amrex::Real* dx,
                              const amrex::Real* dt);

  void read_GL_namelist( const char* inputs_file, const int length);

  void integrate(const int* lo, const int* hi,
                              amrex::Real* phi,
                              const amrex::Real* dx,
                              amrex::Real* integral);

   void  setdt ( const amrex::Real* dx, amrex::Real* dt );

   void  inc_phi0 ( const int* step );

   void  initphi (const int* lo, const int* hi,
                              amrex::Real* phi,
                              const amrex::Real* dx);


#ifdef __cplusplus
}
#endif

#endif
