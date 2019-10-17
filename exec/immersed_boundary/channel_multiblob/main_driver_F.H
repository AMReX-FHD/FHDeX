#ifndef _main_driver_F_H_
#define _main_driver_F_H_

#ifdef __cplusplus
extern "C"
{
#endif

    // these are interfaces to fortran subroutines

///////////////////////////
// in init_s_vel.F90

void init_s_vel(const int * lo,         const int * hi,
                const amrex_real * vel, const int * vello, const int * velhi,
                const amrex_real * dx,
                const amrex_real * reallo, const amrex_real * realhi);

///////////////////////////



#ifdef __cplusplus
}
#endif

#endif
