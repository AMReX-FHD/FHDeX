#ifndef _gmres_functions_H_
#define _gmres_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

using namespace amrex;

// copy contents of gmres_params_module to C++ gmres namespace
void InitializeGmresNamespace();

#endif
