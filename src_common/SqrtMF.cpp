#include "common_functions.H"

void SqrtMF(MultiFab& mf) {

    BL_PROFILE_VAR("SqrtMF()",SqrtMF);

    int ncomp = mf.nComp();

    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();

      const Array4<Real> & mf_fab = mf.array(mfi);

      amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
          mf_fab(i,j,k,n) = sqrt(mf_fab(i,j,k,n));
      });
    }
}
