#ifndef _gmres_functions_H_
#define _gmres_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include "gmres_namespace.H"
#include "GMRES.H"
#include "MacProj.H"
#include "StagMGSolver.H"
#include "Precon.H"

using namespace gmres;
using namespace amrex;

// copy contents of gmres_params_module to C++ gmres namespace
void InitializeGmresNamespace();

// In ApplyMatrix.cpp
void ApplyMatrix(std::array<MultiFab, AMREX_SPACEDIM> & b_u,
                 MultiFab & b_p,
                 std::array<MultiFab, AMREX_SPACEDIM> & x_u,
                 MultiFab & x_p,
                 std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
                 MultiFab & beta,
                 std::array<MultiFab, NUM_EDGE> & beta_ed,
                 MultiFab & gamma,
                 Real & theta_alpha,
                 Geometry & geom,
                 int is_inhomogeneous = 0);
                 
void ApplyMatrix(Vector<std::array< MultiFab*, AMREX_SPACEDIM >> & b_u,
                 Vector<MultiFab*>                               & b_p,
                 Vector<std::array< MultiFab*, AMREX_SPACEDIM >> & x_u,
                 Vector<MultiFab*>                               & x_p,
                 Vector<std::array< MultiFab*, AMREX_SPACEDIM >> & alpha_fc,
                 const Vector<MultiFab*>                         & beta,
                 const Vector<std::array< MultiFab*, NUM_EDGE >> & beta_ed,
                 const Vector<MultiFab*>                         & gamma,
                 const Real                                      & theta_alpha,
                 const Vector<Geometry>                          & geom,
                 int                                             is_inhomogeneous = 0,
                 int nlev = 0);

void UpdateSol(std::array<MultiFab, AMREX_SPACEDIM> & x_u,
               MultiFab & x_p,
               std::array<MultiFab, AMREX_SPACEDIM> & V_u,
               MultiFab & V_p,
               Vector<Real> & y,
               int i);

void LeastSquares(int i,
                  Vector<Vector<Real>> & H,
                  Vector<Real> & cs,
                  Vector<Real> & sn,
                  Vector<Real> & s);

void RotMat(Real a,    Real b,
            Real & cs, Real & sn);

void SolveUTriangular(int k,
                      Vector<Vector<Real>> & H,
                      Vector<Real> & s,
                      Vector<Real> & y);

// In Utility.cpp
void SubtractWeightedGradP(std::array<MultiFab, AMREX_SPACEDIM> & x_u,
                           const std::array<MultiFab, AMREX_SPACEDIM> & alphainv_fc,
                           MultiFab & phi,
                           std::array<MultiFab, AMREX_SPACEDIM> & gradp,
                           const Geometry & geom);

void CCApplyNegLap(MultiFab & phi,
                   MultiFab & Lphi,
                   const std::array<MultiFab, AMREX_SPACEDIM> & beta_fc,
                   const Geometry & geom);

// In StagApplyOp.cpp
void StagApplyOp(const Geometry & geom,
                 const MultiFab & beta_cc,
                 const MultiFab & gamma_cc,
                 const std::array<MultiFab, NUM_EDGE> & beta_ed,
                 const std::array<MultiFab, AMREX_SPACEDIM> & umacIn,
                 std::array<MultiFab, AMREX_SPACEDIM> & umacOut,
                 const std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
                 const Real * dx,
                 const Real & theta_alpha,
                 const int & color=0);

#endif
