#ifndef _electrostatic_H_
#define _electrostatic_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_MLPoisson.H>

using namespace amrex;

void poissonSolve(MultiFab& solution, const MultiFab& rhs);

#endif
