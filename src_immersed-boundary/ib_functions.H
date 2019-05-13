#ifndef _ib_functions_H_
#define _ib_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>


#include <IBParticleContainer.H>

// In IBGMRES.cpp
void IBGMRES(std::array<MultiFab, AMREX_SPACEDIM> & b_u,
             const MultiFab & b_p,
             std::array<MultiFab, AMREX_SPACEDIM> & x_u,
             MultiFab & x_p,
             const std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
             MultiFab & beta,
             std::array<MultiFab, NUM_EDGE> & beta_ed,
             MultiFab & gamma,
             Real theta_alpha,
             const IBParticleContainer & ib_pc,
             const Geometry & geom,
             Real & norm_pre_rhs);

#endif
