#ifndef _gmres_functions_H_
#define _gmres_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include "gmres_functions_F.H"
#include "gmres_namespace.H"
#include "GMRES.H"
#include "StagMGSolver.H"
#include "Precon.H"

using namespace gmres;
using namespace amrex;

#ifndef _num_edge_
#define _num_edge_
#if (AMREX_SPACEDIM == 2)
    const int NUM_EDGE = 1;
#elif (AMREX_SPACEDIM == 3)
    const int NUM_EDGE = 3;
#endif
#endif

// copy contents of gmres_params_module to C++ gmres namespace
void InitializeGmresNamespace();

// In ApplyMatrix.cpp
void ApplyMatrix(std::array<MultiFab, AMREX_SPACEDIM> & b_u,
                 MultiFab & b_p,
                 std::array<MultiFab, AMREX_SPACEDIM> & x_u,
                 MultiFab & x_p,
                 std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
                 const MultiFab & beta,
                 const std::array<MultiFab, NUM_EDGE> & beta_ed,
                 const MultiFab & gamma,
                 const Real & theta_alpha,
                 const Geometry & geom);

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

// In MacProj.cpp
void MacProj(const std::array<MultiFab, AMREX_SPACEDIM> & alphainv_fc,
             MultiFab & mac_rhs,
             MultiFab & phi,
             const Geometry & geom,
             bool full_solve=false);

void SubtractWeightedGradP(std::array<MultiFab, AMREX_SPACEDIM> & x_u,
                           const std::array<MultiFab, AMREX_SPACEDIM> & alphainv_fc,
                           MultiFab & phi,
                           std::array<MultiFab, AMREX_SPACEDIM> & gradp,
                           const Geometry & geom);

void CCApplyNegLap(MultiFab & phi,
                   MultiFab & Lphi,
                   const std::array<MultiFab, AMREX_SPACEDIM> & beta_fc,
                   const Geometry & geom);

// In NormInnerProduct.cpp
void SumStag(const Geometry& geom,
             const std::array<MultiFab, AMREX_SPACEDIM> & m1,
	     const int & comp,
	     Vector<Real> & sum,
	     const bool & divide_by_ncells=false);

void SumCC(const MultiFab & m1,
	   const int & comp,
	   Real & sum,
	   const bool & divide_by_ncells=false);

void StagInnerProd(const Geometry& geom,
                   const std::array<MultiFab, AMREX_SPACEDIM> & m1,
                   const int & comp1,
                   const std::array<MultiFab, AMREX_SPACEDIM> & m2,
                   const int & comp2,
                   std::array<MultiFab, AMREX_SPACEDIM>& mscr,
                   Vector<Real> & prod_val);

void CCInnerProd(const MultiFab & m1,
		 const int & comp1,
		 const MultiFab & m2,
		 const int & comp2,
                 amrex::MultiFab& mscr,
		 Real & prod_val);

void StagL2Norm(const Geometry& geom,
                const std::array<MultiFab, AMREX_SPACEDIM> & m1,
		const int & comp,
                std::array<MultiFab, AMREX_SPACEDIM>& mscr,
		Real & norm_l2);

void CCL2Norm(const MultiFab & m1,
	      const int & comp,
              amrex::MultiFab& mscr,
	      Real & norm_l2);

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
