#ifndef _INS_functions_H_
#define _INS_functions_H_

#include "common_functions.H"
#include "gmres_functions.H"
#include "FhdParticleContainer.H"

///////////////////////////
// in WritePlotFile.cpp

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   const MultiFab& particleInstant,
                   const MultiFab& particleMeans,
                   const MultiFab& particleVars,
                   const MultiFab& cellVols,
                   FhdParticleContainer& particles);

///////////////////////////
// in particleGen.cpp


    void getCellVols(MultiFab & vols, const Geometry & Geom, int samples);


///////////////////////////
// in cellvols.F90
#ifdef __cplusplus
extern "C"
{
#endif
    void get_cell_vols(
                    amrex_real* vols,
                    const int *vlo,
                    const int *vhi,
                    const amrex_real* dx,
                    const int* samples,const amrex_real* plo);

#ifdef __cplusplus
}
#endif


///////////////////////////
// in thermostat.F90
#ifdef __cplusplus
extern "C"
{
#endif
    void thermostat(
                    void* particles,
                    const int* lo, const int* hi,
                    int** c_vectors, int* sizes, const int* clo, const int* chi,
                    amrex_real* cellVols, const int* cvlo, const int* cvhi,
                    const amrex_real* neff, const int* np,
                    const surface* surfaceList, const int* surfaceCount,
                    const Real* tL, const Real* tR);

#ifdef __cplusplus
}
#endif

#endif