#ifndef _reactdiff_functions_H_
#define _reactdiff_functions_H_

#include <AMReX.H>

#include "common_functions.H"
#include "common_namespace.H"
#include "reactDiff_namespace.H"

using namespace reactDiff;
using namespace amrex;
using namespace common;

void InitializeReactDiffNamespace();

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry& geom,
		   const MultiFab& n_in);
#endif
