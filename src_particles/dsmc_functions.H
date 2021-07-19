#ifndef _dsmc_functions_H_
#define _dsmc_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include "dsmc_namespace.H"


using namespace dsmc;
using namespace amrex;

// copy contents of gmres_params_module to C++ gmres namespace
void InitializeDsmcNamespace();

#endif
