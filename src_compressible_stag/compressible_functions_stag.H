#ifndef _compressible_functions_stag_H_
#define _compressible_functions_stag_H_

#define LOHI 2

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

#include "compressible_functions_F.H"
#include "compressible_namespace.H"

#include "common_functions.H"

using namespace compressible;
using namespace amrex;

void WritePlotFileStag(int step,
                       const Real time,
                       const Geometry geom,
                       const MultiFab& cu,
                       const MultiFab& cuMeans,
                       const MultiFab& cuVars,
                       const MultiFab& prim,
                       const MultiFab& primMeans,
                       const MultiFab& primVars,
                       const MultiFab& spatialCross,
                       const MultiFab& eta, 
                       const MultiFab& kappa);

#endif
