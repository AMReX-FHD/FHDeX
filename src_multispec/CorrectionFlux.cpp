#include "multispec_functions.H"

void CorrectionFlux(const MultiFab& rho, const MultiFab& rhotot,
    std::array< MultiFab, AMREX_SPACEDIM >& flux)
{

    BL_PROFILE_VAR("CorrectionFlux()",CorrectionFlux);

    // Loop over boxes
    for (MFIter mfi(rho); mfi.isValid(); ++mfi) {


        AMREX_D_TERM(const Array4<Real> & flux_x = flux[0].array(mfi);,
                     const Array4<Real> & flux_y = flux[1].array(mfi);,
                     const Array4<Real> & flux_z = flux[2].array(mfi););

        AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                     const Box & bx_y = mfi.nodaltilebox(1);,
                     const Box & bx_z = mfi.nodaltilebox(2););

        amrex::ParallelFor(bx_x, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real sum = 0.;
            Real corr = 0.;

            // sum the fluxes upto nspecies-1
            for (int n=0; n<nspecies-1; ++n) {
                sum += flux_x(i,j,k,n);
            }

            // caculate corr and print error if not zero
            corr = flux_x(i,j,k,nspecies-1) + sum;

            // correct flux for last species
            flux_x(i,j,k,nspecies-1) = -sum;
        });


        amrex::ParallelFor(bx_y, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real sum = 0.;
            Real corr = 0.;

            // sum the fluxes upto nspecies-1
            for (int n=0; n<nspecies-1; ++n) {
                sum += flux_y(i,j,k,n);
            }

            // caculate corr and print error if not zero
            corr = flux_y(i,j,k,nspecies-1) + sum;

            // correct flux for last species
            flux_y(i,j,k,nspecies-1) = -sum;

        });

#if (AMREX_SPACEDIM == 3)
        amrex::ParallelFor(bx_z, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real sum = 0.;
            Real corr = 0.;

            // sum the fluxes upto nspecies-1
            for (int n=0; n<nspecies-1; ++n) {
                sum += flux_z(i,j,k,n);
            }

            // caculate corr and print error if not zero
            corr = flux_z(i,j,k,nspecies-1) + sum;

            // correct flux for last species
            flux_z(i,j,k,nspecies-1) = -sum;
        });
#endif
    } // end MFIter

}