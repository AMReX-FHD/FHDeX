#ifndef _common_functions_H_
#define _common_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

#define MAX_SPECIES 10
#define LOHI 2

void main_driver ();

// copy contents of common_params_module to C++ common namespace
void set_common_params();

#endif
