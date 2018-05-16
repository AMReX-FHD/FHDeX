#include "INS_functions.H"
#include "common_functions.H"
#include "common_functions_F.H"

#include "gmres_functions.H"


//Explicit Euler
void eulerStep(const MultiFab& betaCC, const MultiFab& gammaCC,
#if (AMREX_SPACEDIM == 2) 
                 const MultiFab& betaNodal, const MultiFab& gammaNodal,
#endif
#if(AMREX_SPACEDIM == 3)
                 std::array<MultiFab, AMREX_SPACEDIM>& betaEdge,
                 std::array<MultiFab, AMREX_SPACEDIM>& gammaEdge,
#endif
                 std::array<MultiFab, AMREX_SPACEDIM>& umacIn, 
                 std::array<MultiFab, AMREX_SPACEDIM>& umacOut,
                 std::array<MultiFab, AMREX_SPACEDIM>& umacNew,
                 std::array<MultiFab, AMREX_SPACEDIM>& alpha,
                 const Geometry geom,
                 int viscType, Real* dt)
{

    StagApplyOp(betaCC, gammaCC,
#if (AMREX_SPACEDIM == 2)
                betaNodal, gammaNodal,
#endif
#if (AMREX_SPACEDIM == 3)
                betaEdge, gammaEdge,
#endif
                umacIn, umacOut, alpha, geom, viscType);

    const int xOff[3] = {1,0,0};
    const int yOff[3] = {0,1,0};
    const int zOff[3] = {0,0,1};

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(betaCC); mfi.isValid(); ++mfi) 
    {
        const Box& validBox = mfi.validbox();

        AMREX_D_TERM(euler_step_stag(BL_TO_FORTRAN_BOX(validBox),
                           BL_TO_FORTRAN_ANYD(umacIn[0][mfi]),
                           BL_TO_FORTRAN_ANYD(umacOut[0][mfi]),
                           BL_TO_FORTRAN_ANYD(umacNew[0][mfi]),
                           xOff, dt);,
        euler_step_stag(BL_TO_FORTRAN_BOX(validBox),
                           BL_TO_FORTRAN_ANYD(umacIn[1][mfi]),
                           BL_TO_FORTRAN_ANYD(umacOut[1][mfi]),
                           BL_TO_FORTRAN_ANYD(umacNew[1][mfi]),
                           yOff, dt);,
        euler_step_stag(BL_TO_FORTRAN_BOX(validBox),
                           BL_TO_FORTRAN_ANYD(umacIn[2][mfi]),
                           BL_TO_FORTRAN_ANYD(umacOut[2][mfi]),
                           BL_TO_FORTRAN_ANYD(umacNew[2][mfi]),
                           zOff, dt);
                );
    }
}

