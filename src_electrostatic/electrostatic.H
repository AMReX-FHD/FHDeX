#ifndef _electrostatic_H_
#define _electrostatic_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_MLPoisson.H>


using namespace amrex;

void esSolve(MultiFab& potential, MultiFab& charge, std::array< MultiFab, AMREX_SPACEDIM >& efieldCC,
	     const std::array< MultiFab, AMREX_SPACEDIM >& permittivity_fc,
             const std::array< MultiFab, AMREX_SPACEDIM >& beta_es,
             const std::array< MultiFab, AMREX_SPACEDIM >& externalCC,
             const std::array< MultiFab, AMREX_SPACEDIM >& externalFC, 
	     const Geometry geom);

void calculateField(MultiFab& potential, const Geometry geom);

#endif
