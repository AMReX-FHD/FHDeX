
#include "common_functions.H"
#include "common_functions_F.H"

#include "common_namespace.H"

#include "gmres_functions.H"
#include "gmres_functions_F.H"

#include "gmres_namespace.H"

using namespace common;

//Takes cell centred and nodal viscosity multifabs, and face centred velocity
//multifab, and outputs to face-centered velocity multifab.

void StagApplyOp(const MultiFab& beta_cc, const MultiFab& gamma_cc,
                 const std::array<MultiFab, NUM_EDGE>& beta_ed,
                 const std::array<MultiFab, AMREX_SPACEDIM>& umacIn,
                 std::array<MultiFab, AMREX_SPACEDIM>& umacOut,
                 const std::array<MultiFab, AMREX_SPACEDIM>& alpha_fc,
                 const Real* dx,
                 const amrex::Real& theta_alpha,
                 const int& color)
{

    BL_PROFILE_VAR("StagApplyOp()",StagApplyOp);

    const BoxArray& ba = beta_cc.boxArray();
    const DistributionMapping& dmap = beta_cc.DistributionMap();

    // alpha_fc_temp arrays
    std::array< MultiFab, AMREX_SPACEDIM > alpha_fc_temp;
    AMREX_D_TERM(alpha_fc_temp[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 alpha_fc_temp[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 alpha_fc_temp[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););
    AMREX_D_TERM(MultiFab::Copy(alpha_fc_temp[0],alpha_fc[0],0,0,1,0);,
		 MultiFab::Copy(alpha_fc_temp[1],alpha_fc[1],0,0,1,0);,
		 MultiFab::Copy(alpha_fc_temp[2],alpha_fc[2],0,0,1,0););
    AMREX_D_TERM(alpha_fc_temp[0].mult(theta_alpha,1);,
		 alpha_fc_temp[1].mult(theta_alpha,1);,
		 alpha_fc_temp[2].mult(theta_alpha,1););


    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(beta_cc); mfi.isValid(); ++mfi) {

        const Box& validBox = mfi.validbox();

        stag_apply_op(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                      BL_TO_FORTRAN_ANYD(beta_cc[mfi]),
                      BL_TO_FORTRAN_ANYD(gamma_cc[mfi]),
                      BL_TO_FORTRAN_ANYD(beta_ed[0][mfi]),
#if (AMREX_SPACEDIM == 3)
                      BL_TO_FORTRAN_ANYD(beta_ed[1][mfi]),
                      BL_TO_FORTRAN_ANYD(beta_ed[2][mfi]),
#endif
                      BL_TO_FORTRAN_ANYD(umacIn[0][mfi]),
                      BL_TO_FORTRAN_ANYD(umacIn[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                      BL_TO_FORTRAN_ANYD(umacIn[2][mfi]),
#endif
                      BL_TO_FORTRAN_ANYD(umacOut[0][mfi]),
                      BL_TO_FORTRAN_ANYD(umacOut[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                      BL_TO_FORTRAN_ANYD(umacOut[2][mfi]),
#endif
                      BL_TO_FORTRAN_ANYD(alpha_fc_temp[0][mfi]),
                      BL_TO_FORTRAN_ANYD(alpha_fc_temp[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                      BL_TO_FORTRAN_ANYD(alpha_fc_temp[2][mfi]),
#endif
                      dx, &color);

    }

}
