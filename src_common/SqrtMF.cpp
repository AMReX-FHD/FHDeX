#include "common_functions.H"
#include "common_functions_F.H"

void SqrtMF(MultiFab& mf) {
    
    BL_PROFILE_VAR("SqrtMF()",SqrtMF);
    
    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();
      
      sqrt_mf(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
              BL_TO_FORTRAN_FAB(mf[mfi]));
    }
    
}
