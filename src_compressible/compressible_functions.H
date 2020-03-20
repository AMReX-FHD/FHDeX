#ifndef _compressible_functions_H_
#define _compressible_functions_H_

#define MAX_SPECIES 10
#define LOHI 2

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

#include "compressible_functions_F.H"
#include "compressible_namespace.H"

using namespace compressible;
using namespace amrex;

// copy contents of compressible_params_module to C++ compressible namespace
void InitializeCompressibleNamespace();

void calculateFlux(const MultiFab& cons, const MultiFab& prim,
                   const MultiFab& eta, const MultiFab& zeta, const MultiFab& kappa,
		   const MultiFab& chi, const MultiFab& D,
                   std::array<MultiFab, AMREX_SPACEDIM>& flux,
                   std::array<MultiFab, AMREX_SPACEDIM>& stochFlux,
                   std::array<MultiFab, AMREX_SPACEDIM>& cornx,
                   std::array<MultiFab, AMREX_SPACEDIM>& corny,
                   std::array<MultiFab, AMREX_SPACEDIM>& cornz,
                   MultiFab& visccorn,
                   MultiFab& rancorn,
                   const Geometry geom,
		   const Vector< Real >& stoch_weights,
	           const Real* dx, const Real dt);

void calculateTransportCoeffs(const MultiFab& prim, 
			      MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
			      MultiFab& chi, MultiFab& D);

void RK3step(MultiFab& cu, MultiFab& cup, MultiFab& cup2, MultiFab& cup3,
             MultiFab& prim, MultiFab& source,
             MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
             MultiFab& chi, MultiFab& D,
             std::array<MultiFab, AMREX_SPACEDIM>& flux,
             std::array<MultiFab, AMREX_SPACEDIM>& stochFlux, 
             std::array<MultiFab, AMREX_SPACEDIM>& cornx,
             std::array<MultiFab, AMREX_SPACEDIM>& corny,
             std::array<MultiFab, AMREX_SPACEDIM>& cornz,
             MultiFab& visccorn, MultiFab& rancorn,
             const Geometry geom, const Real* dx, const Real dt);

void conservedToPrimitive(MultiFab& prim, const MultiFab& cons);

void primitiveToConserved(const MultiFab& prim, MultiFab& cons);

void setBC(MultiFab& prim, MultiFab& cons);

void evaluateStats(const MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
                   const MultiFab& prim, MultiFab& primMean, MultiFab& primVar,
                   MultiFab& spatialCross,
                   Real* delHolder1, Real* delHolder2, Real* delHolder3,
                   Real* delHolder4, Real* delHolder5, Real* delHolder6,
                   const int steps, const Real* dx);

void doMembrane(MultiFab& cons, MultiFab& prim, std::array<MultiFab, AMREX_SPACEDIM>& flux,
                const Geometry geom, const Real* dx, const Real dt);

void WritePlotFile(int step,
                   const Real time,
                   const Geometry geom,
	           const MultiFab& cu,
	           const MultiFab& cuMeans,
	           const MultiFab& cuVars,
	           const MultiFab& prim,
	           const MultiFab& primMeans,
	           const MultiFab& primVars,
                   const MultiFab& spatialCross,
                   const MultiFab& eta, 
		   const MultiFab& kappa);

void yzAverage(const MultiFab& consMean, 
               const MultiFab& consVar, 
               const MultiFab& primMean, 
               const MultiFab& primVar, 
               const MultiFab& spatialCross, 
               MultiFab& consMeanAv, 
               MultiFab& consVarAv, 
               MultiFab& primMeanAv, 
               MultiFab& primVarAv, 
               MultiFab& spatialCrossAv);

void WriteCheckPoint(int step,
                     const Real time,
                     int statsCount,
                     const Geometry geom,
                     const MultiFab& cu,
                     const MultiFab& cuMeans,
                     const MultiFab& cuVars,
                     const MultiFab& prim,
                     const MultiFab& primMeans,
                     const MultiFab& primVars,
                     const MultiFab& spatialCross, 
                     const MultiFab& eta, 
                     const MultiFab& kappa);

void ReadCheckPoint(int& step,
                     Real& time,
                     int& statsCount,
                     Geometry geom,
                     MultiFab& cu,
                     MultiFab& cuMeans,
                     MultiFab& cuVars,
                     MultiFab& prim,
                     MultiFab& primMeans,
                     MultiFab& primVars,
                     MultiFab& spatialCross,
                     MultiFab& eta,
                     MultiFab& kappa);

#endif
