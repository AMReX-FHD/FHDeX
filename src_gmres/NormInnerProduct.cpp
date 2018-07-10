#include "gmres_functions.H"
#include "gmres_functions_F.H"
#include "gmres_namespace.H"

#include "AMReX_ArrayLim.H"
#include "AMReX_ParallelDescriptor.H"

using namespace amrex;

void StagInnerProd(const std::array<MultiFab, AMREX_SPACEDIM>& m1,
                   const int& comp1, 
                   const std::array<MultiFab, AMREX_SPACEDIM>& m2,
                   const int& comp2,
                   amrex::Vector<amrex::Real>& prod_val)
{
  std::fill(prod_val.begin(), prod_val.end(), 0.);

  // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
  for (MFIter mfi(m1[0]); mfi.isValid(); ++mfi) {

    const Box& validBox_cc = enclosedCells(mfi.validbox());

    stag_inner_prod(ARLIM_3D(validBox_cc.loVect()), ARLIM_3D(validBox_cc.hiVect()),
		    BL_TO_FORTRAN_N_ANYD(m1[0][mfi],comp1),
		    BL_TO_FORTRAN_N_ANYD(m1[1][mfi],comp1),
#if (AMREX_SPACEDIM == 3)
		    BL_TO_FORTRAN_N_ANYD(m1[2][mfi],comp1),
#endif
		    BL_TO_FORTRAN_N_ANYD(m2[0][mfi],comp2),
		    BL_TO_FORTRAN_N_ANYD(m2[1][mfi],comp2),
#if (AMREX_SPACEDIM == 3)
		    BL_TO_FORTRAN_N_ANYD(m2[2][mfi],comp2),
#endif
		    prod_val.dataPtr());

  }

  ParallelDescriptor::ReduceRealSum(prod_val.dataPtr(),AMREX_SPACEDIM);

}
