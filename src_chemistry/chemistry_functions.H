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

void ChemicalRates(const MultiFab& n_cc, MultiFab& chem_rate, const amrex::Geometry& geom, const amrex::Real& dt,
                   const MultiFab& n_interm, Vector<Real>& lin_comb_coef_in, Real volume_factor_in=1.);

AMREX_GPU_HOST_DEVICE void compute_reaction_rates(GpuArray<Real,MAX_SPECIES>& n_in,
                                                  GpuArray<Real,MAX_REACTION>& reaction_rates,
                                                  const amrex::Real& dv);

AMREX_GPU_HOST_DEVICE void sample_num_reactions(GpuArray<Real,MAX_SPECIES>& n_in,
                                                GpuArray<Real,MAX_REACTION>& num_reactions,
                                                GpuArray<Real,MAX_REACTION>& avg_num_reactions,
                                                const amrex::RandomEngine& engine);
#endif
