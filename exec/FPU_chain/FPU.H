#ifndef _FPU_H_
#define _FPU_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

void init_p(MultiFab& state,
            const Real& beta);

#endif
