#ifndef _surfchem_mui_functions_H_
#define _surfchem_mui_functions_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <math.h>
#include "common_functions.H"
#include "common_namespace.H"

#include <mui.h>
#include "surfchem_mui_namespace.H"

#define AVONUM 6.02214076e23

using namespace surfchem_mui;
using namespace mui;
using namespace amrex;
using namespace common;

void InitializeSurfChemMUINamespace();

void mui_push(MultiFab& cu, MultiFab& prim, const amrex::Real* dx, mui::uniface2d &uniface, const int step);

void mui_fetch(MultiFab& cu, MultiFab& prim, const amrex::Real* dx, mui::uniface2d &uniface, const int step);

void mui_announce_send_recv_span(mui::uniface2d &uniface, MultiFab& mf, const Real* dx);

#endif
