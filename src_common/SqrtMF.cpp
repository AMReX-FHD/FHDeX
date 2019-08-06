#include "common_functions.H"
#include "common_functions_F.H"

void SqrtMF(MultiFab& mf) {
  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {

    // Note: Make sure that multifab is cell-centered
    const Box& validBox_cc = enclosedCells(mfi.validbox());

    sqrt_mf(ARLIM_3D(validBox_cc.loVect()), ARLIM_3D(validBox_cc.hiVect()),
	    BL_TO_FORTRAN_FAB(mf[mfi]));
  }
}
