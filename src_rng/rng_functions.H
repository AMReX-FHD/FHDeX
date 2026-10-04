#ifndef _rng_functions_H_
#define _rng_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

using namespace amrex;

///////////////////////////
// in MultiFabFillRandom.cpp

// Fill component comp of mf with N(0,variance) samples in the valid region plus
// ng ghost cells (ng <= mf.nGrow()). Interior and periodic ghost cells are then
// overwritten by FillBoundary. Ghost cells outside a non-periodic boundary keep
// their own samples, so two boxes adjacent along a wall hold independent values
// for the same out-of-domain ghost cell.
void MultiFabFillRandom(MultiFab& mf, const int& comp, const Real& variance, const Geometry& geom, const int& ng=0);

void MultiFabFillRandomNormal(MultiFab& mf, const int& scomp, const int& ncomp,
                              const amrex::Real& mean, const amrex::Real& variance,
                              const Geometry& geom,
                              bool overridesync = true, bool fillboundary = true);

void MultiFabFillRandomUniform(MultiFab& mf, const int& scomp, const int& ncomp,
                               const Geometry& geom,
                               bool overridesync = true, bool fillboundary = true);

#endif
