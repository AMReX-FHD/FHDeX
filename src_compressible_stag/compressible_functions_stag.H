#ifndef _compressible_functions_stag_H_
#define _compressible_functions_stag_H_

#define LOHI 2

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

#include "compressible_functions_F.H"
#include "compressible_namespace.H"

#include "common_functions.H"

using namespace compressible;
using namespace amrex;

void WritePlotFileStag(int step,
                       const Real time,
                       const Geometry geom,
                       const MultiFab& cu,
                       const MultiFab& cuMeans,
                       const MultiFab& cuVars,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumomMeans,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumomVars,
                       const MultiFab& prim,
                       const MultiFab& primMeans,
                       const MultiFab& primVars,
                       const std::array<MultiFab, AMREX_SPACEDIM>& vel,
                       const std::array<MultiFab, AMREX_SPACEDIM>& velMeans,
                       const std::array<MultiFab, AMREX_SPACEDIM>& velVars,
//                       const MultiFab& spatialCross,
                       const MultiFab& eta, 
                       const MultiFab& kappa);

void conservedToPrimitiveStag(MultiFab& prim_in, std::array<MultiFab, AMREX_SPACEDIM>& velStag_in,
                              MultiFab& cons_in, const std::array<MultiFab, AMREX_SPACEDIM>& momStag_in);

void InitConsVarStag(MultiFab& cons, std::array< MultiFab, AMREX_SPACEDIM >& momStag, 
                     const MultiFab& prim, const amrex::Geometry geom);

void setBCStag(MultiFab& prim_in, MultiFab& cons_in,
                 std::array< MultiFab, AMREX_SPACEDIM >& cumom_in,
                 std::array< MultiFab, AMREX_SPACEDIM >& vel_in,
                 const amrex::Geometry geom);

void BCMassTempPress(MultiFab& prim_in,const amrex::Geometry geom,int dim);

void BCMomNormal(MultiFab& mom_in, MultiFab& vel_in,
                 const amrex::Geometry geom, int dim);

void BCMomTrans(MultiFab& mom_in, MultiFab& vel_in,
                 const amrex::Geometry geom, int dim);


void BCRhoRhoE(MultiFab& cons_in, MultiFab& prim_in, 
               std::array< MultiFab, AMREX_SPACEDIM >& cumom_in, 
               const amrex::Geometry geom);

void StochFluxStag(std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in, std::array< MultiFab, 2 >& edgeflux_x_in,
                   std::array< MultiFab, 2 >& edgeflux_y_in, std::array< MultiFab, 2 >& edgeflux_z_in,
                   const amrex::Geometry geom);

void RK3stepStag(MultiFab& cu, 
                 std::array< MultiFab, AMREX_SPACEDIM >& cumom,
                 MultiFab& prim, std::array< MultiFab, AMREX_SPACEDIM >& facevel,
                 MultiFab& source,
                 MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
                 MultiFab& chi, MultiFab& D,
                 std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                 std::array< MultiFab, 2 >& edgeflux_x,
                 std::array< MultiFab, 2 >& edgeflux_y,
                 std::array< MultiFab, 2 >& edgeflux_z,
                 std::array< MultiFab, AMREX_SPACEDIM>& cenflux,
                 const amrex::Geometry geom, const amrex::Real dt);

void calculateFluxStag(const MultiFab& cons_in, const std::array< MultiFab, AMREX_SPACEDIM >& momStag_in, 
                       const MultiFab& prim_in, const std::array< MultiFab, AMREX_SPACEDIM >& velStag_in,
                       const MultiFab& eta_in, const MultiFab& zeta_in, const MultiFab& kappa_in,
                       const MultiFab& chi_in, const MultiFab& D_in,
                       std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in,
                       std::array< MultiFab, 2 >& edgeflux_x_in,
                       std::array< MultiFab, 2 >& edgeflux_y_in,
                       std::array< MultiFab, 2 >& edgeflux_z_in,
                       std::array< MultiFab, AMREX_SPACEDIM>& cenflux_in,
                       std::array< MultiFab, AMREX_SPACEDIM>& stochface_in,
                       std::array< MultiFab, 2 >& stochedge_x_in,
                       std::array< MultiFab, 2 >& stochedge_y_in,
                       std::array< MultiFab, 2 >& stochedge_z_in,
                       std::array< MultiFab, AMREX_SPACEDIM>& stochcen_in,
                       const amrex::Geometry geom,
		                   const amrex::Vector< amrex::Real >& stoch_weights,
		                   const amrex::Real dt);

void evaluateStatsStag(const MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
                       const MultiFab& prim_in, MultiFab& primMean, MultiFab& primVar,
                       const std::array<MultiFab, AMREX_SPACEDIM>& vel, 
                       std::array<MultiFab, AMREX_SPACEDIM>& velMean, 
                       std::array<MultiFab, AMREX_SPACEDIM>& velVar, 
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                       std::array<MultiFab, AMREX_SPACEDIM>& cumomMean,
                       std::array<MultiFab, AMREX_SPACEDIM>& cumomVar,
                       const int steps, const amrex::Real* dx);

//void yzAverageStag(const MultiFab& consMean, 
//                   const MultiFab& consVar, 
//                   const MultiFab& primMean, 
//                   const MultiFab& primVar, 
//                   const MultiFab& spatialCross, 
//                   MultiFab& consMeanAv, 
//                   MultiFab& consVarAv, 
//                   MultiFab& primMeanAv, 
//                   MultiFab& primVarAv, 
//                   MultiFab& spatialCrossAv);
#endif
