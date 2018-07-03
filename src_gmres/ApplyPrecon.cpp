#include "common_functions.H"
#include "common_functions_F.H"
#include "common_namespace.H"

#include "gmres_functions.H"
#include "gmres_functions_F.H"
#include "gmres_namespace.H"

using namespace common;
using namespace gmres;

void ApplyPrecon(const std::array<MultiFab, AMREX_SPACEDIM>& b_u,
                 const MultiFab& b_p,
                 std::array<MultiFab, AMREX_SPACEDIM>& x_u,
                 MultiFab& x_p,
                 const std::array<MultiFab, AMREX_SPACEDIM>& alpha_fc,
                 const MultiFab& beta,
                 const std::array<MultiFab, NUM_EDGE>& beta_ed,
                 const MultiFab& gamma,
                 const Real& theta_alpha,
                 const Real* dx)
{

    // projection preconditioner
    if (abs(precon_type) == 1) {

        /*
          STEP 1: Solve for an intermediate state, x_u^star, using an implicit viscous solve
                  x_u^star = A^{-1} b_u
        */

        // x_u^star = A^{-1} b_u
        //
        //

        /*
          STEP 2: Construct RHS for pressure Poisson problem
        */

        // set mac_rhs = D(x_u^star)
        //
        //

        // add b_p to mac_rhs
        //
        //

        /*
          STEP 3: Compute x_u
        */

        // use multigrid to solve for Phi
        // x_u^star is only passed in to get a norm for absolute residual criteria
        //
        //

        // x_u = x_u^star - (alpha I)^-1 grad Phi
        //
        //

        /*
          STEP 4: Compute x_p by applying the Schur complement approximation
        */

        if (visc_schur_approx == 0) {
            // if precon_type = +1, or theta_alpha=0 then x_p = theta_alpha*Phi - c*beta*(mac_rhs)
            // if precon_type = -1                   then x_p = theta_alpha*Phi - c*beta*L_alpha Phi
           
            if (precon_type == 1 || theta_alpha == 0) {
                // first set x_p = -mac_rhs 
                //
                //
            }
            else {
                // first set x_p = -L_alpha Phi
                //
                //
            }

            if ( abs(visc_type) == 1 || abs(visc_type) == 2) {
                // multiply x_p by beta; x_p = -beta L_alpha Phi
                //
                //

                if (abs(visc_type) == 2) {
                    // multiply by c=2; x_p = -2*beta L_alpha Phi
                    //
                    //
                }
            }
            else if (abs(visc_type) == 3) {

                // multiply x_p by gamma, use mac_rhs a temparary to save x_p 
                //
                //
                // multiply x_p by beta; x_p = -beta L_alpha Phi
                //
                //
                // multiply by c=4/3; x_p = -(4/3) beta L_alpha Phi
                //
                //
                // x_p = -(4/3) beta L_alpha Phi - gamma L_alpha Phi
                //
                //
            }

            // multiply Phi by theta_alpha
            //
            //

            // add theta_alpha*Phi to x_p
            //
            //
        }
        else {
            Abort("StagApplyOp: visc_schur_approx != 0 not supported");
        }
    }
    else {
        Abort("StagApplyOp: unsupposed precon_type");
    }

}
