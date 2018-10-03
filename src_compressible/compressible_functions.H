#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

using namespace amrex;

void calculateFlux(const MultiFab& cons, const MultiFab& prim,
                   const MultiFab& eta, const MultiFab& zeta, const MultiFab& kappa,
                   std::array<MultiFab, AMREX_SPACEDIM>& flux,
	        	   const amrex::Real* dx, const int nvars, const int nprimvars);

