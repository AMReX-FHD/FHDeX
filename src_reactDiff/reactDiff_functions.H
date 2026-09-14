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

void GenerateStochasticFluxdivCorrector(MultiFab& n_old,
                                        MultiFab& n_new,
                                        MultiFab& stoch_fluxdiv,
                                        const std::array< MultiFab, AMREX_SPACEDIM >& diff_coef_face,
                                        const Real& dt,
                                        const Real& time,
                                        const Geometry& geom);

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
// In AdvanceReactionDiffusion.cpp
////////////////////////
void AdvanceReactionDiffusion(MultiFab& n_old,
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
// In ImplicitDiffusion.cpp
////////////////////////
void ImplicitDiffusion(MultiFab& n_old,
                       MultiFab& n_new,
                       const MultiFab& rhs,
                       const std::array< MultiFab, AMREX_SPACEDIM >& diff_coef_face,
                       const Geometry& geom,
                       const Real& dt_fac,
                       const Real& time);

////////////////////////
// In DiffusiveNFluxdiv.cpp
////////////////////////
void DiffusiveNFluxdiv(MultiFab& n_in,
                       MultiFab& diff_fluxdiv,
                       const std::array< MultiFab, AMREX_SPACEDIM >& diff_coef_face,
                       const Geometry& geom,
                       const Real& time);

////////////////////////
// In InitN.cpp
////////////////////////
void InitN(MultiFab& n_in,
           const Geometry& geom,
           const Real& time);

////////////////////////
// In MultinomialDiffusion.cpp
////////////////////////
void MultinomialDiffusion(MultiFab& n_old,
                          MultiFab& n_new,
                          const std::array< MultiFab, AMREX_SPACEDIM >& diff_coef_face,
                          const Geometry& geom,
                          const Real& dt,
                          const Real& time);

AMREX_GPU_HOST_DEVICE void multinomial_rng(GpuArray<Real,2*AMREX_SPACEDIM>& samples,
                                           const int& N,
                                           GpuArray<Real,2*AMREX_SPACEDIM>& p,
                                           const amrex::RandomEngine& engine);

////////////////////////
// In StochasticNFluxdiv.cpp
////////////////////////
void StochasticNFluxdiv(MultiFab& n_in,
                        MultiFab& stoch_fluxdiv,
                        const std::array< MultiFab, AMREX_SPACEDIM >& diff_coef_face,
                        const Geometry& geom,
                        const Real& dt,
                        const Real& time,
                        int increment_div=0);

AMREX_GPU_HOST_DEVICE Real average_to_faces(const Real& value1,
                                            const Real& value2,
                                            const Real& dv);

////////////////////////
// In WritePlotFile.cpp
////////////////////////

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry& geom,
                   const MultiFab& n_in);

#endif