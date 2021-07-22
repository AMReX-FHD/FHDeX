#ifndef _chemistry_functions_H_
#define _chemistry_functions_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_MultiFab.H>
#include <math.h>
#include "common_functions.H"
#include "common_functions_F.H"
#include "common_namespace.H"
#include "chemistry_namespace.H"

using namespace chemistry;
using namespace amrex;
using namespace common;

void InitializeChemistryNamespace();

void compute_reaction_rates(amrex::Real n_dens[MAX_SPECIES], amrex::Real a_r[MAX_REACTION]);
void compute_Omega(MultiFab& rho, MultiFab& Omega);
void advance_reaction_SSA_cell(amrex::Real n_old[MAX_SPECIES],amrex::Real n_new[MAX_SPECIES],amrex::Real dt,amrex::Real dV,RandomEngine const& engine);
void advance_reaction_det_cell(amrex::Real n_old[MAX_SPECIES],amrex::Real n_new[MAX_SPECIES],amrex::Real dt);
void advance_reaction_CLE_cell(amrex::Real n_old[MAX_SPECIES],amrex::Real n_new[MAX_SPECIES],amrex::Real dt,amrex::Real dV,RandomEngine const& engine);

#endif