#include "multispec_functions.H"

void ProjectOntoEOS(MultiFab& rho_in)
{
    if (algorithm_type == 4 || algorithm_type == 6) {

        if (nspecies == 1) {
            rho_in.setVal(rho0);
            return;
        }

        for ( MFIter mfi(rho_in,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            const Box& bx = mfi.tilebox();
            const Array4<Real>& rho = rho_in.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real sum = 0.;
                for (int n=0; n<nspecies-1; ++n) {
                    sum += rho(i,j,k,n);
                }
                rho(i,j,k,nspecies-1) = rho0 - sum;
            });

        }

    } else {

        Abort("ProjectOntoEOS algorithm_type = 0, 2, 3, 5 not written");

    }

}
