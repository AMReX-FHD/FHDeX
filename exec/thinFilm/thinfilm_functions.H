#ifndef _thinfilm_functions_H_
#define _thinfilm_functions_H_

#include <AMReX_MultiFab.H>

#include "thinfilm_namespace.H"

#ifdef AMREX_USE_CUDA
#include <cufft.h>
#endif

using namespace thinfilm;
using namespace amrex;

void InitializeThinfilmNamespace();

#ifdef AMREX_USE_CUDA
std::string cufftErrorToString (const cufftResult& err);
#endif

#endif
