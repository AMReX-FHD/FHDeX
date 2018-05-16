#include "common_functions.H"
#include "common_functions_F.H"

#include "gmres_functions.H"
#include "gmres_functions_F.H"
#include "gmres_namespace.H"

//Takes cell centred and nodal viscosity multifabs, and face centred velcoity multifab, and outputs to 

void StagApplyOp(const MultiFab& betaCC, const MultiFab& gammaCC,
#if (AMREX_SPACEDIM == 2) 
                 const MultiFab& betaNodal, const MultiFab& gammaNodal, 
#endif
#if (AMREX_SPACEDIM == 3) 
                 std::array<MultiFab, AMREX_SPACEDIM>& betaEdge, 
                 std::array<MultiFab, AMREX_SPACEDIM>& gammaEdge, 
#endif
                 std::array<MultiFab, AMREX_SPACEDIM>& umacIn, 
                 std::array<MultiFab, AMREX_SPACEDIM>& umacOut,
                 std::array<MultiFab, AMREX_SPACEDIM>& alpha,
                 const Geometry geom,
                 int viscType)
{
    const int v = viscType;
    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(betaCC); mfi.isValid(); ++mfi) {

        const Box& validBox = mfi.validbox();

        stag_apply_op(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                           BL_TO_FORTRAN_ANYD(betaCC[mfi]),
                           BL_TO_FORTRAN_ANYD(gammaCC[mfi]),
#if (AMREX_SPACEDIM == 2)
                           BL_TO_FORTRAN_ANYD(betaNodal[mfi]),
                           BL_TO_FORTRAN_ANYD(gammaNodal[mfi]),
#endif
#if (AMREX_SPACEDIM == 3)
                           BL_TO_FORTRAN_ANYD(betaEdge[0][mfi]),
                           BL_TO_FORTRAN_ANYD(betaEdge[1][mfi]),
                           BL_TO_FORTRAN_ANYD(betaEdge[1][mfi]),
                           BL_TO_FORTRAN_ANYD(gammaEdge[0][mfi]),
                           BL_TO_FORTRAN_ANYD(gammaEdge[1][mfi]),
                           BL_TO_FORTRAN_ANYD(gammaEdge[1][mfi]),

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

                           BL_TO_FORTRAN_ANYD(alpha[0][mfi]),
                           BL_TO_FORTRAN_ANYD(alpha[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                           BL_TO_FORTRAN_ANYD(alpha[2][mfi]),
#endif

                           geom.CellSize(),&v);

    }

}







