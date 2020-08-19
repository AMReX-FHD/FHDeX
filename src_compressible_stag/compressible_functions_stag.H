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
                       const MultiFab& prim,
                       const MultiFab& primMeans,
                       const MultiFab& primVars,
                       const MultiFab& spatialCross,
                       const MultiFab& eta, 
                       const MultiFab& kappa);

void conservedToPrimitiveStag(MultiFab& prim_in, std::array<MultiFab, AMREX_SPACEDIM>& velStag_in,
                              MultiFab& cons_in, const std::array<MultiFab, AMREX_SPACEDIM>& momStag_in);

void InitConsVarStag(MultiFab& cons, std::array< MultiFab, AMREX_SPACEDIM >& momStag, 
                     const MultiFab& prim, const amrex::Geometry geom);

void RK3stepStag(MultiFab& cu, MultiFab& cup, MultiFab& cup2,
                 std::array< MultiFab, AMREX_SPACEDIM >& momStag,
                 std::array< MultiFab, AMREX_SPACEDIM >& mompStag,
                 std::array< MultiFab, AMREX_SPACEDIM >& momp2Stag,
                 MultiFab& prim, std::array< MultiFab, AMREX_SPACEDIM >& velStag,
                 MultiFab& source,
                 MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
                 MultiFab& chi, MultiFab& D,
                 std::array<MultiFab, AMREX_SPACEDIM>& flux,
                 std::array<MultiFab, AMREX_SPACEDIM>& stochFlux,
                 std::array<MultiFab, AMREX_SPACEDIM>& cornx,
                 std::array<MultiFab, AMREX_SPACEDIM>& corny,
                 std::array<MultiFab, AMREX_SPACEDIM>& cornz,
                 MultiFab& visccorn, MultiFab& rancorn,
                 const Geometry geom, const Real* dxp, const amrex::Real dt);

void calculateFluxStag(const MultiFab& cons_in, const std::array< MultiFab, AMREX_SPACEDIM >& momStag_in, 
                       const MultiFab& prim_in, const std::array< MultiFab, AMREX_SPACEDIM >& velStag_in,
                       const MultiFab& eta_in, const MultiFab& zeta_in, const MultiFab& kappa_in,
                       const MultiFab& chi_in, const MultiFab& D_in,
                       std::array<MultiFab, AMREX_SPACEDIM>& flux_in,
                       std::array<MultiFab, AMREX_SPACEDIM>& stochFlux_in,
                       std::array<MultiFab, AMREX_SPACEDIM>& cornx_in,
                       std::array<MultiFab, AMREX_SPACEDIM>& corny_in,
                       std::array<MultiFab, AMREX_SPACEDIM>& cornz_in,
                       MultiFab& visccorn_in,
                       MultiFab& rancorn_in,
                       const amrex::Geometry geom,
		                   const amrex::Vector< amrex::Real >& stoch_weights,
		                   const amrex::Real* dx, const amrex::Real dt);

#endif
