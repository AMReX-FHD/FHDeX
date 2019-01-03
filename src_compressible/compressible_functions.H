#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

using namespace amrex;

void calculateFlux(const MultiFab& cons, const MultiFab& prim,
                   const MultiFab& eta, const MultiFab& zeta, const MultiFab& kappa,
                   std::array<MultiFab, AMREX_SPACEDIM>& flux,
                   std::array<MultiFab, AMREX_SPACEDIM>& stochFlux,
                   std::array<MultiFab, AMREX_SPACEDIM>& cornx,
                   std::array<MultiFab, AMREX_SPACEDIM>& corny,
                   std::array<MultiFab, AMREX_SPACEDIM>& cornz,
                   MultiFab& visccorn,
                   MultiFab& rancorn,
                   const amrex::Geometry geom,
	           const amrex::Real* dx, const amrex::Real dt);

void calculateTransportCoeffs(const MultiFab& prim, MultiFab& eta, MultiFab& zeta, MultiFab& kappa);

void RK3step(MultiFab& cu, MultiFab& cup, MultiFab& cup2, MultiFab& cup3, MultiFab& prim, MultiFab& source, MultiFab& eta, MultiFab& zeta, MultiFab& kappa, std::array<MultiFab, AMREX_SPACEDIM>& flux, std::array<MultiFab, AMREX_SPACEDIM>& stochFlux, 
                                                 std::array<MultiFab, AMREX_SPACEDIM>& cornx, std::array<MultiFab, AMREX_SPACEDIM>& corny, std::array<MultiFab, AMREX_SPACEDIM>& cornz, MultiFab& visccorn, MultiFab& rancorn, const amrex::Geometry geom, const amrex::Real* dx, const amrex::Real dt);

void conservedToPrimitive(MultiFab& prim, const MultiFab& cons);

void primitiveToConserved(const MultiFab& prim, MultiFab& cons);

void setBC(MultiFab& prim, MultiFab& cons, MultiFab& eta, MultiFab& zeta, MultiFab& kappa);

void evaluateStats(const MultiFab& cons, MultiFab& consMean, MultiFab& consVar, const MultiFab& prim, MultiFab& primMean, MultiFab& primVar, MultiFab& spatialCross, MultiFab& eta, MultiFab& etaMean, MultiFab& kappa, MultiFab& kappaMean, Real* delHolder1, Real* delHolder2, Real* delHolder3, Real* delHolder4, Real* delHolder5, Real* delHolder6, const int steps, const amrex::Real* dx);

void doSSA(MultiFab& cons, MultiFab& prim, std::array<MultiFab, AMREX_SPACEDIM>& flux, const amrex::Geometry geom, const amrex::Real* dx, const amrex::Real dt);

void yzAverage(const MultiFab& consMean, const MultiFab& consVar, const MultiFab& primMean, const MultiFab& primVar, const MultiFab& spatialCross, const MultiFab& etaMean, const MultiFab& kappaMean, MultiFab& consMeanAv, MultiFab& consVarAv, MultiFab& primMeanAv, MultiFab& primVarAv, MultiFab& spatialCrossAv, MultiFab& etaMeanAv, MultiFab& kappaMeanAv);

