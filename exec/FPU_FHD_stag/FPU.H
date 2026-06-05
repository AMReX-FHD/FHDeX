#ifndef _FPU_H_
#define _FPU_H_

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include "common_functions.H"
#include "FPU_namespace.H"

using namespace amrex;
using namespace FPU;

void InitializeNamespace();

//void ComputePhiFromState(MultiFab& phi);
//
//void ComputeCalphaalpha(MultiFab& C_alphaalpha,
//                        const MultiFab& phi,
//                        const MultiFab& phi0);

void ReadCheckPoint(int& step,
                    Real& time,
                    int& statsCount,
                    Geometry& geom,
                    const amrex::Box& domain,
                    MultiFab& cu,
                    MultiFab& cuMeans,
                    MultiFab& cuVars,
                    MultiFab& cumom,
                    MultiFab& cumomMeans,
                    MultiFab& cumomVars,
                    MultiFab& coVars,
                    BoxArray& ba, DistributionMapping& dmap);

void WriteCheckPoint(int step,
                     const Real time,
                     int statsCount,
                     const Geometry& geom,
                     const MultiFab& cu,
                     const MultiFab& cuMeans,
                     const MultiFab& cuVars,
                     const MultiFab& cumom,
                     const MultiFab& cumomMeans,
                     const MultiFab& cumomVars,
                     const MultiFab& coVars);

void InitConsVarStag(MultiFab& cu,
                     MultiFab& cumom,
                     Geometry& geom);

void WritePlotFile(int step,
                       const Real time,
                       const Geometry& geom,
                       const MultiFab& cu,
                       const MultiFab& cuMeans,
                       const MultiFab& cuVars,
                       const MultiFab& cumom,
                       const MultiFab& cumomMeans,
                       const MultiFab& cumomVars,
                       const MultiFab& coVars);

void RK3step(MultiFab& cu,
             MultiFab& cumom,
             MultiFab& faceflux,
             MultiFab& cenflux,
             const Geometry& geom, 
             const Real dt, 
             const int step);

void calculateFlux(const MultiFab& cu,
                   const MultiFab& cumom,
                   MultiFab& faceflux,
                   MultiFab& cenflux,
                   const MultiFab& stochface,
                   const MultiFab& stochcen,
                   const Geometry& geom);
                     
void evaluateStats(const MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
                         const MultiFab& cumom,
                         MultiFab& cumomMean,
                         MultiFab& cumomVar,
                         MultiFab& coVar,
                         const int steps,
                         const Geometry& geom);

void GotoNextLine(std::istream& is);

void WritePlotFilesSF_1D(const amrex::MultiFab& mag, const amrex::MultiFab& realimag,
                         const int step, const Real time, const amrex::Vector< std::string >& names, std::string plotfile_base);

//void refreshCellCenteredMomentum(MultiFab& cu,
//                                  const MultiFab& cumom,
//                                  const Geometry& geom);
//
//void initializeState(MultiFab& cu,
//                     MultiFab& cumom,
//                     const Geometry& geom);

#endif
