#ifndef _gmres_functions_H_
#define _gmres_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#if (AMREX_SPACEDIM == 2)
const int NUM_EDGE = 1;
#elif (AMREX_SPACEDIM == 3)
const int NUM_EDGE = 3;
#endif

// copy contents of gmres_params_module to C++ gmres namespace
void InitializeGmresNamespace();

// In ApplyPrecon.cpp
void ApplyPrecon(const std::array<amrex::MultiFab, AMREX_SPACEDIM>& b_u,
                 const amrex::MultiFab& b_p,
                 std::array<amrex::MultiFab, AMREX_SPACEDIM>& x_u,
                 amrex::MultiFab& x_p,
                 const std::array<amrex::MultiFab, AMREX_SPACEDIM>& alpha_fc,
                 const amrex::MultiFab& beta,
                 const std::array<amrex::MultiFab, NUM_EDGE>& beta_ed,
                 const amrex::MultiFab& gamma,
                 const amrex::Real& theta_alpha,
                 const amrex::Real* dx,
                 const amrex::Geometry geom);

// In NormInnerProduct.cpp
void StagInnerProd(const std::array<amrex::MultiFab, AMREX_SPACEDIM>& m1,
                   const int& comp1, 
                   const std::array<amrex::MultiFab, AMREX_SPACEDIM>& m2,
                   const int& comp2,
                   amrex::Vector<amrex::Real> prod_val);

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
void StagMGSolver(const std::array< amrex::MultiFab, AMREX_SPACEDIM >& alpha_fc,
                  const amrex::MultiFab& beta_cc,
                  const std::array< amrex::MultiFab, NUM_EDGE >& beta_ed,
                  const amrex::MultiFab& gamma_cc,
                  std::array< amrex::MultiFab, AMREX_SPACEDIM >& phi_fc,
                  const std::array< amrex::MultiFab, AMREX_SPACEDIM >& phiorig_fc,
                  const amrex::Real& theta,
                  const amrex::Geometry& geom);

// compute the number of multigrid levels assuming minwidth is the length of the
// smallest dimension of the smallest grid at the coarsest multigrid level
int ComputeNlevsMG(const amrex::BoxArray& ba);

void CCRestriction(amrex::MultiFab& phi_c, const amrex::MultiFab& phi_f,
                   const amrex::Geometry& geom_c);

void StagRestriction(std::array< amrex::MultiFab, AMREX_SPACEDIM >& phi_c, 
                     const std::array< amrex::MultiFab, AMREX_SPACEDIM >& phi_f,
                     int simple_stencil=0);

void NodalRestriction(amrex::MultiFab& phi_c, const amrex::MultiFab& phi_f); 

void EdgeRestriction(std::array< amrex::MultiFab, NUM_EDGE >& phi_c, 
                      const std::array< amrex::MultiFab, NUM_EDGE >& phi_f);

void StagProlongation(const std::array< amrex::MultiFab, AMREX_SPACEDIM >& phi_c, 
                      std::array< amrex::MultiFab, AMREX_SPACEDIM >& phi_f);

void StagMGUpdate(std::array< amrex::MultiFab, AMREX_SPACEDIM >& phi_fc,
                  const std::array< amrex::MultiFab, AMREX_SPACEDIM >& rhs_fc,
                  const std::array< amrex::MultiFab, AMREX_SPACEDIM >& Lphi_fc,
                  const std::array< amrex::MultiFab, AMREX_SPACEDIM >& alpha_fc,
                  const amrex::MultiFab& beta_cc,
                  const std::array< amrex::MultiFab, NUM_EDGE >& beta_ed,
                  const amrex::MultiFab& gamma_cc,
                  const amrex::Real* dx,
                  const int& color = 0);

// In StagApplyOp.cpp
void StagApplyOp(const amrex::MultiFab& betaCC, 
                 const amrex::MultiFab& gammaCC, 
                 const std::array<amrex::MultiFab, NUM_EDGE>& betaEdge, 
                 const std::array<amrex::MultiFab, AMREX_SPACEDIM>& umacIn, 
                 std::array<amrex::MultiFab, AMREX_SPACEDIM>& umacOut,
                 const std::array<amrex::MultiFab, AMREX_SPACEDIM>& alpha,
                 const amrex::Real* dx,
                 const int& color=0);

#endif
