#ifndef _FPU_H_
#define _FPU_H_

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include "rng_functions.H"

#include "FPU_namespace.H"

#include <chrono>
#include <cmath>

using namespace amrex;
using namespace FPU;

void ComputePhiFromState(MultiFab& phi);

void ComputeCalphaalpha(MultiFab& C_alphaalpha,
                        const MultiFab& phi,
                        const MultiFab& phi0);

#endif
