#ifndef _gmres_functions_H_
#define _gmres_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>

// copy contents of gmres_params_module to C++ gmres namespace
void InitializeGmresNamespace();

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
#if (AMREX_SPACEDIM == 2)
                  const std::array< amrex::MultiFab, 1 >& beta_ed,
#elif (AMREX_SPACEDIM == 3)
                  const std::array< amrex::MultiFab, 3 >& beta_ed,
#endif
                  const amrex::MultiFab& gamma_cc,
                  const amrex::Real& theta,
                  const amrex::Geometry& geom);

// compute the number of multigrid levels assuming minwidth is the length of the
// smallest dimension of the smallest grid at the coarsest multigrid level
int ComputeNlevsMG(const amrex::BoxArray& ba);

void CCRestriction(amrex::MultiFab& phi_c, const amrex::MultiFab& phi_f);

void StagRestriction(amrex::MultiFab& phi_c, const amrex::MultiFab& phi_f);

void EdgeRestriction(amrex::MultiFab& phi_c, const amrex::MultiFab& phi_f);

void NodalRestriction(amrex::MultiFab& phi_c, const amrex::MultiFab& phi_f);

#endif
