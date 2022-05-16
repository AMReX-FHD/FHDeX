#ifndef _MFsurfchem_functions_H_
#define _MFsurfchem_functions_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <math.h>
#include "common_functions.H"
#include "common_namespace.H"
#include "MFsurfchem_namespace.H"

#define AVONUM 6.02214076e23

using namespace MFsurfchem;
using namespace amrex;
using namespace common;

void InitializeMFSurfchemNamespace();

void init_surfcov(MultiFab& surfcov);

void sample_MFsurfchem(MultiFab& cu, MultiFab& prim, MultiFab& surfcov, MultiFab& dNadsdes, const amrex::Real* dx, const amrex::Real dt);

void update_MFsurfchem(MultiFab& cu, MultiFab& surfcov, MultiFab& dNadsdes, const amrex::Real* dx, const amrex::Real dt);

#endif
