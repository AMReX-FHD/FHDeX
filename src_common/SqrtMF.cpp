#include "common_functions.H"
#include "common_functions_F.H"

void SqrtMF(MultiFab& mf) {
    
    BL_PROFILE_VAR("SqrtMF()",SqrtMF);

    int ncomp = mf.nComp();
    
    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();
      
      const Array4<Real> & mf_fab = mf.array(mfi);
        
      AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, n,
      {
          mf_fab(i,j,k,n) = sqrt(mf_fab(i,j,k,n));
      });
    }
}
