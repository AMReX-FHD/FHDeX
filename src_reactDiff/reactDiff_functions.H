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
void AdvanceDiffusion(MultiFab& n_old,
                      MultiFab& n_new,
                      const MultiFab& ext_src,
                      const Real& dt,
                      const Real& time,
                      const Geometry& geom);

void DiffusiveNFluxdiv(MultiFab& n_in,
                       MultiFab& diff_fluxdiv,
                       const Geometry& geom,
                       const Real& time);

////////////////////////
// In AdvanceReaction.cpp
////////////////////////
void AdvanceReaction(MultiFab& n_old,
                     MultiFab& n_new,
                     const MultiFab& ext_src,
                     const Real& dt,
                     const Real& time,
                     const Geometry& geom);

////////////////////////
// In AdvanceTimestep.cpp
////////////////////////
void AdvanceTimestep(MultiFab& n_old,
                     MultiFab& n_new,
                     const Real& dt,
                     const Real& time,
                     const Geometry& geom);

////////////////////////
// In InitN.cpp
////////////////////////
void InitN(MultiFab& n_in,
           const Geometry& geom,
           const Real& time);

////////////////////////
// In WritePlotFile.cpp
////////////////////////

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry& geom,
		   const MultiFab& n_in);

#endif
