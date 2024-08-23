#ifndef _reactdiff_functions_H_
#define _reactdiff_functions_H_

#include <AMReX.H>

#include "common_functions.H"
#include "common_namespace.H"
#include "reactDiff_namespace.H"

using namespace reactDiff;
using namespace amrex;
using namespace common;

////////////////////////
// In reactDiff_functions.cpp
////////////////////////
void InitializeReactDiffNamespace();

////////////////////////
// In AdvanceDiffusion.cpp
////////////////////////
void AdvanceDiffusion(const MultiFab& n_old,
                      MultiFab& n_new,
                      const MultiFab& ext_src,
                      const Real& dt,
                      const Real& time,
                      const Geometry& geom);

////////////////////////
// In AdvanceTimestep.cpp
////////////////////////
void AdvanceTimestep(const MultiFab& n_old,
                     MultiFab& n_new,
                     const Real& dt,
                     const Real& time,
                     const Geometry& geom);

////////////////////////
// In WritePlotFile.cpp
////////////////////////

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry& geom,
		   const MultiFab& n_in);

#endif
