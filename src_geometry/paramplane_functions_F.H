#ifndef _geometry_F_H_
#define _geometry_F_H_


#ifdef __cplusplus
extern "C"
{
#endif


    void find_inter(void* part, const amrex_real* delt, const void* paramplanes, const int* ns, int* intsurf,
                        amrex_real* inttime, int* intside, const amrex_real* phi, const amrex_real* plo);

    void app_bc(const void* surf, const void* part, int* intside, amrex_real* domsize, int* push, const amrex_real* time,
                    const amrex_real* inttime);

#ifdef __cplusplus
}
#endif

#endif


