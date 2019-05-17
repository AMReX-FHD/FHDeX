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


void IBMPrecon(const std::array<MultiFab, AMREX_SPACEDIM> & b_u, const MultiFab & b_p,
               std::array<MultiFab, AMREX_SPACEDIM> & x_u, MultiFab & x_p,
               const std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
               const MultiFab & beta, const std::array<MultiFab, NUM_EDGE> & beta_ed,
               const MultiFab & gamma, const Real & theta_alpha,
               const IBParticleContainer & ib_pc,
               const Vector<std::pair<int, int>> & pindex_list,
               std::map<std::pair<int, int>, Vector<RealVect>> & marker_forces,
               const std::map<std::pair<int, int>, Vector<RealVect>> & marker_W,
               const Geometry & geom);


void ApplyIBM(      std::array<MultiFab, AMREX_SPACEDIM>            & b_u,
                    std::map<std::pair<int, int>, Vector<RealVect>> & b_lambda,
              const std::array<MultiFab, AMREX_SPACEDIM>            & x_u,
              const IBParticleContainer                             & ib_pc,
              const Vector<std::pair<int, int>>                     & part_indices,
              const std::map<std::pair<int, int>, Vector<RealVect>> & x_lambda,
              int ib_grow, int ibpc_lev, const Geometry & geom );


void MarkerInvSub(Vector<RealVect> & a, const Vector<RealVect> & b);
void MarkerInvSub(const Vector<std::pair<int, int>> & part_indices,
                        std::map<std::pair<int, int>, Vector<RealVect>> & a,
                  const std::map<std::pair<int, int>, Vector<RealVect>> & b);


void MarkerInnerProd(const Vector<RealVect> & a, const Vector<RealVect> & b, Real & v);


void MarkerL2Norm(const Vector<RealVect> & markers, Real & norm_l2);
void MarkerL2Norm(const Vector<std::pair<int, int>> & part_indices,
                  const std::map<std::pair<int, int>, Vector<RealVect>> & b, Real & v);


void MarkerCopy(Vector<RealVect> & a, const Vector<RealVect> & b);


#endif
