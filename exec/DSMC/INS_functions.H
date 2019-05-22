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

#endif
