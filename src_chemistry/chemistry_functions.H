#ifndef _chemistry_functions_H_
#define _chemistry_functions_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <math.h>
#include "common_functions.H"
#include "common_namespace.H"
#include "chemistry_namespace.H"

using namespace chemistry;
using namespace amrex;
using namespace common;

void InitializeChemistryNamespace();

// used in compressible code only
void compute_compressible_chemistry_source_CLE(amrex::Real dt, amrex::Real dV,
                                               MultiFab& prim, MultiFab& source, MultiFab& ranchem);
#endif
