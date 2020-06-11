#ifndef _compressible_functions_H_
#define _compressible_functions_H_

#define LOHI 2

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

#include "compressible_functions_F.H"
#include "compressible_namespace.H"

#include "common_functions.H"

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

void calculateTransportCoeffs(const MultiFab& prim_in,
			      MultiFab& eta_in, MultiFab& zeta_in, MultiFab& kappa_in,
			      MultiFab& chi_in, MultiFab& Dij_in);

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

void conservedToPrimitive(MultiFab& prim_in, const MultiFab& cons_in);

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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GetMolfrac ( int i, int j, int k,
                  const Array4<Real>& Yk,
                  GpuArray<Real,MAX_SPECIES>& Xk,
                  int const nspecies,
                  const GpuArray<Real,MAX_SPECIES>& molmass)
{
    Real molmix = 0.;

    for (int n=0; n<nspecies; ++n) {
        molmix += Yk(i,j,k,n)/molmass[n];
    }
    molmix = 1./molmix;
    for (int n=0; n<nspecies; ++n) {
        Xk[n] = Yk(i,j,k,n)*molmix/molmass[n];
    }

};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GetTemperature ( int i, int j, int k,
                      Real const energy,
                      const Array4<Real>& Yk,
                      Real& temp,
                      int const nspecies,
                      const GpuArray<Real,MAX_SPECIES>& hcv)
{
    Real cvmix = 0.;
    Real e0 = 0.;

    for (int n=0; n<nspecies; ++n) {
        cvmix += Yk(i,j,k,n)*hcv[n];
    }

    temp = (energy-e0)/cvmix;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GetPressureGas ( int i, int j, int k,
                      Real& pressure,
                      const Array4<Real>& Yk,
                      Real const density,
                      Real const temp,
                      int const nspecies,
                      Real const Runiv,
                      const GpuArray<Real,MAX_SPECIES>& molmass)
{
    Real molmix = 0.;

    for (int n=0; n<nspecies; ++n) {
        molmix += Yk(i,j,k,n)/molmass[n];
    }
    molmix = 1./molmix;

    pressure = density*(Runiv/molmix)*temp;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void IdealMixtureTransport ( int i, int j, int k,
                             Real const density,
                             Real const temperature,
                             Real const pressure,
                             const Array4<Real>& Yk,
                             const Array4<Real>& Xk,
                             Real& eta,
                             Real& kappa,
                             Real& zeta,
                             const Array4<Real>& diff_ij,
                             const Array4<Real>& chitil,
                             const Array4<Real>& xxtr,
                             int const nspecies,
                             Real const * const AMREX_RESTRICT molmass)
{

    // mole fracions correction - EGLIB
    for (int ii=0; ii<nspecies; ++ii) {
        Real sumXk = 0.;
        for (int jj=0; jj<nspecies; ++jj) {
            sumXk += Xk(i,j,k,jj);
        }
        // GM: Why this factor of 1E-15???
        xxtr(i,j,k,ii) = Xk(i,j,k,ii) + 1.e-15*(sumXk/nspecies - Xk(i,j,k,ii));
    }
    
    // molecular weight of mixture - EGLIB
    Real MWmix = 0.;
    for (int ii=0; ii<nspecies; ++ii) {
        MWmix += xxtr(i,j,k,ii)*molmass[ii];
    }
};
                      
#endif
