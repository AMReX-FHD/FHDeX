#ifndef _gmres_functions_H_
#define _gmres_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>


#ifndef _num_edge_
#define _num_edge_
#if (AMREX_SPACEDIM == 2)
    const int NUM_EDGE = 1;
#elif (AMREX_SPACEDIM == 3)
    const int NUM_EDGE = 3;
#endif
#endif


using namespace amrex;

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

// In ApplyPrecon.cpp
void ApplyPrecon(const std::array<MultiFab, AMREX_SPACEDIM> & b_u,
                 const MultiFab & b_p,
                 std::array<MultiFab, AMREX_SPACEDIM> & x_u,
                 MultiFab & x_p,
                 const std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
                 const MultiFab & beta,
                 const std::array<MultiFab, NUM_EDGE> & beta_ed,
                 const MultiFab & gamma,
                 const Real & theta_alpha,
                 const Geometry & geom);

// In GMRES.cpp
void GMRES(std::array<MultiFab, AMREX_SPACEDIM> & b_u,
           const MultiFab & b_p,
           std::array<MultiFab, AMREX_SPACEDIM> & x_u,
           MultiFab & x_p,
           std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
           MultiFab & beta,
           std::array<MultiFab, NUM_EDGE> & beta_ed,
           MultiFab & gamma,
           Real theta_alpha,
           const Geometry & geom,
           Real & norm_pre_rhs);

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
                           const Geometry & geom);

void CCApplyOp(MultiFab & phi,
               MultiFab & Lphi,
               const MultiFab & alpha,
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
                   Vector<Real> & prod_val);

void CCInnerProd(const MultiFab & m1,
		 const int & comp1,
		 const MultiFab & m2,
		 const int & comp2,
		 Real & prod_val);

void StagL2Norm(const Geometry& geom,
                const std::array<MultiFab, AMREX_SPACEDIM> & m1,
		const int & comp,
		Real & norm_l2);

void CCL2Norm(const MultiFab & m1,
	      const int & comp,
	      Real & norm_l2);

// solve "(theta*alpha*I - L) phi = rhs" using multigrid with Jacobi relaxation
// if abs(visc_type) = 1, L = div beta grad
// if abs(visc_type) = 2, L = div [ beta (grad + grad^T) ]
// if abs(visc_type) = 3, L = div [ beta (grad + grad^T) + I (gamma - (2/3)*beta) div ]
// if visc_type > 1 we assume constant coefficients
// if visc_type < 1 we assume variable coefficients
// beta_cc, and gamma_cc are cell-centered
// alpha_fc, phi_fc, and rhs_fc are face-centered
// beta_ed is nodal (2d) or edge-centered (3d)
// phi_fc must come in initialized to some value, preferably a reasonable guess
void StagMGSolver(const std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
                  const MultiFab & beta_cc,
                  const std::array<MultiFab, NUM_EDGE> & beta_ed,
                  const MultiFab & gamma_cc,
                  std::array<MultiFab, AMREX_SPACEDIM> & phi_fc,
                  const std::array<MultiFab, AMREX_SPACEDIM> & phiorig_fc,
                  const Real & theta,
                  const Geometry & geom);

// compute the number of multigrid levels assuming minwidth is the length of the
// smallest dimension of the smallest grid at the coarsest multigrid level
int ComputeNlevsMG(const BoxArray & ba);

void CCRestriction(MultiFab & phi_c, const MultiFab & phi_f,
                   const Geometry & geom_c);

void StagRestriction(std::array<MultiFab, AMREX_SPACEDIM> & phi_c,
                     const std::array<MultiFab, AMREX_SPACEDIM > & phi_f,
                     int simple_stencil=0);

void NodalRestriction(MultiFab & phi_c, const MultiFab & phi_f);

void EdgeRestriction(std::array<MultiFab, NUM_EDGE> & phi_c,
                      const std::array<MultiFab, NUM_EDGE> & phi_f);

void StagProlongation(const std::array<MultiFab, AMREX_SPACEDIM> & phi_c,
                      std::array<MultiFab, AMREX_SPACEDIM> & phi_f);

void StagMGUpdate(std::array<MultiFab, AMREX_SPACEDIM> & phi_fc,
                  const std::array<MultiFab, AMREX_SPACEDIM> & rhs_fc,
                  const std::array<MultiFab, AMREX_SPACEDIM> & Lphi_fc,
                  const std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
                  const MultiFab & beta_cc,
                  const std::array<MultiFab, NUM_EDGE> & beta_ed,
                  const MultiFab & gamma_cc,
                  const Real * dx,
                  const int & color=0);

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
