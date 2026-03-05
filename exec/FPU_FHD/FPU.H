#ifndef _FPU_H_
#define _FPU_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

void ComputePhiFromState(MultiFab& phi,
                         const Real& r_eq,
                         const Real& p_eq,
                         const Real& e_eq,
                         const Real& R_00,
                         const Real& R_01,
                         const Real& R_02,
                         const Real& R_10,
                         const Real& R_11,
                         const Real& R_12,
                         const Real& R_20,
                         const Real& R_21,
                         const Real& R_22);

void ComputeCalphaalpha(MultiFab& C_alphaalpha,
                        const MultiFab& phi,
                        const MultiFab& phi0);

#endif
