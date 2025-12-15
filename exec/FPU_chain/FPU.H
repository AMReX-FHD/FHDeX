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

void compute_mean_stretch_momentum(MultiFab& state,
                                   const int& n_particles,
                                   const int& n_ensembles);

void compute_energy(MultiFab& state,
                    const Real& a,
                    const Real& b,
                    const Real& c);

void compute_mean_energy(MultiFab& energy,
                         const int& n_particles,
                         const int& n_ensembles);
            

#endif
