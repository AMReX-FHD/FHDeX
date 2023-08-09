#ifndef _main_driver_F_H_
#define _main_driver_F_H_

#ifdef __cplusplus
extern "C"
{
#endif

// these are interfaces to fortran subroutines

///////////////////////////
// in init_s_vel.F90

void init_s_vel(const int * lo, const int * hi,
                const amrex_real * vel, const int * vello, const int * velhi,
                const amrex_real * dx,
                const amrex_real * reallo, const amrex_real * realhi);

///////////////////////////


///////////////////////////
// in init_funcs.F90

void init_vel(const int* lo, const int* hi,
              amrex_real* data, const int* dlo, const int* dhi,
              const amrex_real* dx, const amrex_real* prob_lo,
              const amrex_real* prob_hi, const int* di,
              const amrex::Real* real_lo, const amrex::Real* real_hi,
              const int* prob_type);

#ifdef __cplusplus
}
#endif

#endif
