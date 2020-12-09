#ifndef _ib_functions_H_
#define _ib_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>


#include <IBParticleContainer.H>

void InitializeImmbdyNamespace();
void InitializeIBFlagellumNamespace();
void InitializeIBColloidNamespace();


// In IBGMRES.cpp
void IBGMRES(std::array<MultiFab, AMREX_SPACEDIM> & b_u,
             const MultiFab & b_p,
             std::array<MultiFab, AMREX_SPACEDIM> & x_u,
             MultiFab & x_p,
             std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
             MultiFab & beta,
             std::array<MultiFab, NUM_EDGE> & beta_ed,
             MultiFab & gamma,
             Real theta_alpha,
             const IBParticleContainer & ib_pc,
             const Geometry & geom,
             Real & norm_pre_rhs);


void IBMPrecon(const std::array<MultiFab, AMREX_SPACEDIM> & b_u, const MultiFab & b_p,
               std::array<MultiFab, AMREX_SPACEDIM> & x_u, MultiFab & x_p,
               std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
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


void MarkerAdd(Vector<RealVect> & a, const Vector<RealVect> & b);
void MarkerAdd(const Vector<std::pair<int, int>> & part_indices,
                     std::map<std::pair<int, int>, Vector<RealVect>> & a,
               const std::map<std::pair<int, int>, Vector<RealVect>> & b);
void MarkerAdd(const Vector<std::pair<int, int>> & part_indices, int comp,
                     std::map<std::pair<int, int>,        Vector<RealVect>>  & a,
               const std::map<std::pair<int, int>, Vector<Vector<RealVect>>> & b);


void MarkerSub(Vector<RealVect> & a, const Vector<RealVect> & b);
void MarkerSub(const Vector<std::pair<int, int>> & part_indices,
                     std::map<std::pair<int, int>, Vector<RealVect>> & a,
               const std::map<std::pair<int, int>, Vector<RealVect>> & b);


void MarkerInvSub(Vector<RealVect> & a, const Vector<RealVect> & b);
void MarkerInvSub(const Vector<std::pair<int, int>> & part_indices,
                        std::map<std::pair<int, int>, Vector<RealVect>> & a,
                  const std::map<std::pair<int, int>, Vector<RealVect>> & b);


void MarkerMult(Real factor, Vector<RealVect> & a);
void MarkerMult(const Vector<std::pair<int, int>> & part_indices, Real factor,
                std::map<std::pair<int, int>, Vector<RealVect>> & a);
void MarkerMult(const Vector<std::pair<int, int>> & part_indices, int comp, Real factor,
                std::map<std::pair<int, int>, Vector<Vector<RealVect>>> & a);


void MarkerInnerProd(const Box & bx, const Geometry & geom, const Vector<RealVect> & marker_pos,
                     const Vector<RealVect> & a, const Vector<RealVect> & b, Real & v);
void MarkerInnerProd(const Vector<std::pair<int, int>> & part_indices,
                     const MultiFab & cc_iter, const Geometry & geom,
                     const std::map<std::pair<int, int>, Vector<RealVect>> & marker_pos,
                     const std::map<std::pair<int, int>, Vector<RealVect>> & a,
                     const std::map<std::pair<int, int>, Vector<RealVect>> & b,
                     Real & v);
void MarkerInnerProd(const Vector<std::pair<int, int>> & part_indices, int comp,
                     const MultiFab & cc_iter, const Geometry & geom,
                     const std::map<std::pair<int, int>, Vector<RealVect>> & marker_pos,
                     const std::map<std::pair<int, int>,        Vector<RealVect>>  & a,
                     const std::map<std::pair<int, int>, Vector<Vector<RealVect>>> & b,
                     Real & v);


void MarkerL2Norm(const Box & bx, const Geometry & geom, const Vector<RealVect> & marker_pos,
                  const Vector<RealVect> & markers, Real & norm_l2);
void MarkerL2Norm(const Vector<std::pair<int, int>> & part_indices,
                  const MultiFab & cc_iter, const Geometry & geom,
                  const std::map<std::pair<int, int>, Vector<RealVect>> & marker_pos,
                  const std::map<std::pair<int, int>, Vector<RealVect>> & b, Real & v);


void MarkerCopy(Vector<RealVect> & a, const Vector<RealVect> & b);
void MarkerCopy(const Vector<std::pair<int, int>> & part_indices,
                      std::map<std::pair<int, int>, Vector<RealVect>> & a,
                const std::map<std::pair<int, int>, Vector<RealVect>> & b);
void MarkerCopy(const Vector<std::pair<int, int>> & part_indices, int comp,
                      std::map<std::pair<int, int>, Vector<Vector<RealVect>>> & a,
                const std::map<std::pair<int, int>,        Vector<RealVect>>  & b);
void MarkerCopy(const Vector<std::pair<int, int>> & part_indices, int comp,
                      std::map<std::pair<int, int>,        Vector<RealVect>>  & a,
                const std::map<std::pair<int, int>, Vector<Vector<RealVect>>> & b);


void UpdateSolIBM(const Vector<std::pair<int, int>>                       & part_indices,
                  std::array<MultiFab, AMREX_SPACEDIM>                    & x_u,
                  MultiFab                                                & x_p,
                  std::map<std::pair<int, int>,        Vector<RealVect>>  & x_lambda,
                  std::array<MultiFab, AMREX_SPACEDIM>                    & V_u,
                  MultiFab                                                & V_p,
                  std::map<std::pair<int, int>, Vector<Vector<RealVect>>> & V_lambda,
                  const Vector<Real>                                      & y,
                  int i );

#endif
