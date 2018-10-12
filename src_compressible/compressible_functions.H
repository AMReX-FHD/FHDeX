#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

using namespace amrex;

void calculateFlux(const MultiFab& cons, const MultiFab& prim,
                   const MultiFab& eta, const MultiFab& zeta, const MultiFab& kappa,
                   std::array<MultiFab, AMREX_SPACEDIM>& flux, std::array<MultiFab, AMREX_SPACEDIM>& stochFlux,const amrex::Geometry geom,
	        	   const amrex::Real* dx, const amrex::Real dt);

void calculateTransportCoeffs(const MultiFab& prim, MultiFab& eta, MultiFab& zeta, MultiFab& kappa);

void RK3step(MultiFab& cu, MultiFab& cup, MultiFab& cup2, MultiFab& cup3, MultiFab& prim, MultiFab& source, MultiFab& eta, MultiFab& zeta, MultiFab& kappa, std::array<MultiFab, AMREX_SPACEDIM>& flux, std::array<MultiFab, AMREX_SPACEDIM>& stochFlux,
	                   const amrex::Geometry geom, const amrex::Real* dx, const amrex::Real dt);

void conservedToPrimitive(MultiFab& prim, const MultiFab& cons);

void primitiveToConserved(const MultiFab& prim, MultiFab& cons);

void setBC(MultiFab& prim, MultiFab& cons, MultiFab& eta, MultiFab& zeta, MultiFab& kappa);

