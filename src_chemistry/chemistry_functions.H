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

void compute_chemistry_source_CLE(amrex::Real dt, amrex::Real dV,
                                  MultiFab& prim, MultiFab& source, MultiFab& ranchem);

AMREX_GPU_HOST_DEVICE void compute_reaction_rates(GpuArray<amrex::Real,MAX_SPECIES>& n_dens,
                                                  GpuArray<amrex::Real,MAX_REACTION>& a_r);

void compute_chemistry_source_CLE_1(amrex::Real dt, amrex::Real dV, MultiFab& mf_in, int startComp_in,
                                    MultiFab& source, int startComp_out);

void compute_chemistry_source_CLE_2(amrex::Real dt, amrex::Real dV, MultiFab& mf_in, int startComp_in,
                                    MultiFab& source, int startComp_out, MultiFab& ranchem);

void EMstep_chem_only(MultiFab& rho_old, MultiFab& rho_new,
                      const amrex::Geometry geom, const amrex::Real dt);

void RK3step_chem_only(MultiFab& rho_old, MultiFab& rho_new,
                       const amrex::Geometry geom, const amrex::Real dt);

AMREX_GPU_HOST_DEVICE void advance_reaction_det_cell(GpuArray<amrex::Real,MAX_SPECIES>& n_old,
                                                     GpuArray<amrex::Real,MAX_SPECIES>& n_new,amrex::Real dt);

AMREX_GPU_HOST_DEVICE void advance_reaction_CLE_cell(GpuArray<amrex::Real,MAX_SPECIES>& n_old,
                                                     GpuArray<amrex::Real,MAX_SPECIES>& n_new,
                                                     amrex::Real dt,amrex::Real dV,RandomEngine const& engine);

AMREX_GPU_HOST_DEVICE void advance_reaction_SSA_cell(GpuArray<amrex::Real,MAX_SPECIES>& n_old,
                                                     GpuArray<amrex::Real,MAX_SPECIES>& n_new,
                                                     amrex::Real dt,amrex::Real dV,RandomEngine const& engine);
#endif
