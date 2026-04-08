#include "multispec_functions.H"

/**
 * Performs matrix vector multiplication for vectors of length nspecies,
 * and matrix size nspecies x nspecies.
 *
 * \param[in,out] x_in vector at each i,j,k location in multifab x_in
 * \param[in] A_in Matrix at each i,j,k location in multifab A_in
 *
 *
 */

void MatvecMul(MultiFab& x_in,
    const MultiFab& A_in)
{

    BL_PROFILE_VAR("MatvecMul()",MatvecMul);

    // Loop over boxes
    for (MFIter mfi(x_in); mfi.isValid(); ++mfi) {

        // Create a box that matches the NODALITY of MultiFab
        const Box& validBox = mfi.validbox();

        const Array4<      Real>& x = x_in.array(mfi);
        const Array4<const Real>& A = A_in.array(mfi);


        amrex::ParallelFor(validBox, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            GpuArray<Real,MAX_SPECIES> x_temp;
            Real sum;

            for (int m=0; m<nspecies; ++m){
                x_temp[m] = x(i,j,k,m);
            }

            for (int m=0; m<nspecies; ++m){
                sum = 0.0;
                for (int n=0; n<nspecies; ++n){
                    sum += A(i,j,k,n*nspecies+m) * x_temp[n];
                }

                x(i,j,k,m) = sum;
            }


        });
    }

}