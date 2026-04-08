#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

using namespace amrex;

void calculateFlux(const MultiFab& cons,
                   std::array<MultiFab, AMREX_SPACEDIM>& flux,
                   std::array<MultiFab, AMREX_SPACEDIM>& stochFlux,
       const amrex::Real* dx, const amrex::Real dt);

void eulerStep(MultiFab& cu, std::array<MultiFab,AMREX_SPACEDIM>& flux,
               std::array<MultiFab, AMREX_SPACEDIM>& stochFlux,
               const amrex::Geometry geom, const amrex::Real* dx,
               const amrex::Real dt);

void RK3step(MultiFab& cu, MultiFab& cup, MultiFab& cup2, MultiFab& cup3,
    std::array<MultiFab,AMREX_SPACEDIM>& flux, std::array<MultiFab,
    AMREX_SPACEDIM>& stochFlux, const amrex::Geometry geom,
    const amrex::Real* dx, const amrex::Real dt);

void setBC(MultiFab& prim, MultiFab& cons, MultiFab& eta, MultiFab& zeta, MultiFab& kappa);

void evaluateStats(const MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
    const MultiFab& prim, MultiFab& primMean, MultiFab& primVar,
    MultiFab& spatialCross, MultiFab& eta, MultiFab& etaMean,
    MultiFab& kappa, MultiFab& kappaMean, Real* delHolder1,
    Real* delHolder2, Real* delHolder3, Real* delHolder4,
    Real* delHolder5, Real* delHolder6, const int steps,
    const amrex::Real* dx);