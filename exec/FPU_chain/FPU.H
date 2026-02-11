#ifndef _FPU_H_
#define _FPU_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

void init(MultiFab& state,
          const Real& beta,
          const Real& pressure,
          const Real& a,
          const Real& b,
          const Real& c,
          const Real& r0,
          const int& burn,
          const Real& step_size,
          const int& n_particles,
          const int& n_ensembles,
          const Geometry& geom);

void FPU_RK4(MultiFab& state,
             const Real& a,
             const Real& b,
             const Real& c,
             const Real& dt,
             const int& n_particles,
             const int& n_ensembles,
             const Geometry& geom);

void rhs(MultiFab& rhs,
         MultiFab& state,
         const Real& a,
         const Real& b,
         const Real& c,
         const Geometry& geom);

void compute_means(MultiFab& state,
                   const int& n_particles,
                   const int& n_ensembles,
                   const int& step);

void compute_energy(MultiFab& state,
                    const Real& a,
                    const Real& b,
                    const Real& c);

void compute_S_alphaalpha(const MultiFab& state,
                          const MultiFab& g_alpha_zero,
                          MultiFab& Salphaalpha);

void ComputePhiFromState(MultiFab& phi,
                         const Real& r_eq,
                         const Real& p_eq,
                         const Real& e_eq,
                         const Real& R_00,
                         const Real& R_01,
                         const Real& R_02,
                         const Real& R_10,
                         const Real& R_11,
                         const Real& R_12,
                         const Real& R_20,
                         const Real& R_21,
                         const Real& R_22);

void ComputeCalphaalpha(MultiFab& C_alphaalpha,
                        const MultiFab& phi,
                        const MultiFab& phi0);

#endif
