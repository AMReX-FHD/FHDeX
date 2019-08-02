
#include "common_functions.H"
#include "common_functions_F.H"

#include "common_namespace.H"

#include "gmres_functions.H"
#include "gmres_functions_F.H"

#include "gmres_namespace.H"

using namespace common;

//Takes cell centred and nodal viscosity multifabs, and face centred velocity
//multifab, and outputs to face-centered velocity multifab.

void StagApplyOp(const MultiFab& beta_cc,
                 const MultiFab& gamma_cc,
                 const std::array<MultiFab, NUM_EDGE>& beta_ed,
                 const std::array<MultiFab, AMREX_SPACEDIM>& phi,
                 std::array<MultiFab, AMREX_SPACEDIM>& Lphi,
                 std::array<MultiFab, AMREX_SPACEDIM>& alpha_fc,
                 const Real* dx,
                 const amrex::Real& theta_alpha,
                 const int& color)
{

    BL_PROFILE_VAR("StagApplyOp()",StagApplyOp);

    AMREX_D_DECL(bool do_x=true,
                      do_y=true,
                      do_z=true);

    int offset = 1;

    if (color == 1 || color == 2) {
        AMREX_D_TERM(do_x = true;,
                     do_y = false;,
                     do_z = false;);
        offset = 2;
    }
    else if (color == 3 || color == 4) {
        AMREX_D_TERM(do_x = false;,
                     do_y = true;,
                     do_z = false;);
        offset = 2;
    }
#if (AMREX_SPACEDIM == 3)
    else if (color == 4 || color == 5) {
        AMREX_D_TERM(do_x = false;,
                     do_y = false;,
                     do_z = true;);
        offset = 2;
    }
#endif    

    Real dxsqinv = 1./(dx[0]*dx[0]);
    Real dysqinv = 1./(dx[1]*dx[1]);
    Real dxdyinv = 1./(dx[0]*dx[1]);
#if (AMREX_SPACEDIM == 3)
    Real dzsqinv = 1./(dx[2]*dx[2]);
    Real dxdzinv = 1./(dx[0]*dx[2]);
    Real dydzinv = 1./(dx[1]*dx[2]);
#endif
    
    // multiply alpha by theta_alpha
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        alpha_fc[d].mult(theta_alpha);
    }

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(beta_cc); mfi.isValid(); ++mfi) {

        const Box & bx = mfi.validbox();
        
        Array4<Real const> const& beta_cc_fab = beta_cc.array(mfi);
        Array4<Real const> const& gamma_cc_fab = gamma_cc.array(mfi);

        Array4<Real const> const& beta_xy_fab = beta_ed[0].array(mfi);
#if (AMREX_SPACEDIM == 3)
        Array4<Real const> const& beta_xz_fab = beta_ed[1].array(mfi);
        Array4<Real const> const& beta_yz_fab = beta_ed[2].array(mfi);
#endif
        
        AMREX_D_TERM(Array4<Real const> const& phix_fab = phi[0].array(mfi);,
                     Array4<Real const> const& phiy_fab = phi[1].array(mfi);,
                     Array4<Real const> const& phiz_fab = phi[2].array(mfi););
        
        AMREX_D_TERM(Array4<Real> const& Lphix_fab = Lphi[0].array(mfi);,
                     Array4<Real> const& Lphiy_fab = Lphi[1].array(mfi);,
                     Array4<Real> const& Lphiz_fab = Lphi[2].array(mfi););
        
        AMREX_D_TERM(Array4<Real> const& alphax_fab = alpha_fc[0].array(mfi);,
                     Array4<Real> const& alphay_fab = alpha_fc[1].array(mfi);,
                     Array4<Real> const& alphaz_fab = alpha_fc[2].array(mfi););

        if (visc_type == -1) {

        }
        else if (visc_type == 1) {

        }
        else if (visc_type == -2) {

        }
        else if (visc_type == 2) {

        }
        else {
            Abort("StagApplyOp.cpp: unsupported visc_type");
        }

        stag_apply_op(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                      BL_TO_FORTRAN_ANYD(beta_cc[mfi]),
                      BL_TO_FORTRAN_ANYD(gamma_cc[mfi]),
                      BL_TO_FORTRAN_ANYD(beta_ed[0][mfi]),
#if (AMREX_SPACEDIM == 3)
                      BL_TO_FORTRAN_ANYD(beta_ed[1][mfi]),
                      BL_TO_FORTRAN_ANYD(beta_ed[2][mfi]),
#endif
                      BL_TO_FORTRAN_ANYD(phi[0][mfi]),
                      BL_TO_FORTRAN_ANYD(phi[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                      BL_TO_FORTRAN_ANYD(phi[2][mfi]),
#endif
                      BL_TO_FORTRAN_ANYD(Lphi[0][mfi]),
                      BL_TO_FORTRAN_ANYD(Lphi[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                      BL_TO_FORTRAN_ANYD(Lphi[2][mfi]),
#endif
                      BL_TO_FORTRAN_ANYD(alpha_fc[0][mfi]),
                      BL_TO_FORTRAN_ANYD(alpha_fc[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                      BL_TO_FORTRAN_ANYD(alpha_fc[2][mfi]),
#endif
                      dx, &color);

    }

    // divide alpha by theta_alpha
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        alpha_fc[d].mult(1./theta_alpha);
    }
}

