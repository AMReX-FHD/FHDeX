#ifndef _electrostatic_H_
#define _electrostatic_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_MLPoisson.H>


using namespace amrex;

void esSolve(MultiFab& potential, MultiFab& charge, std::array< MultiFab, AMREX_SPACEDIM >& efield, const std::array< MultiFab, AMREX_SPACEDIM >& external, const Geometry geom);

void calculateField(MultiFab& potential, const Geometry geom);

#endif
