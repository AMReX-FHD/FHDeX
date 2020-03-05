#ifndef _rng_functions_H_
#define _rng_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

#include "rng_functions_F.H"

using namespace amrex;

///////////////////////////
// in MultiFABFillRandom.cpp

void MultiFABFillRandom(MultiFab& mf, const int& comp, const Real& variance, const Geometry& geom);

void MultiFABFillRandomHack(MultiFab& mf, const int& comp, const Real& variance, const Geometry& geom);

#endif
