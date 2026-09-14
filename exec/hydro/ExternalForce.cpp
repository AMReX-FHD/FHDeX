#include "common_functions.H"

void ExternalForce(std::array< MultiFab, AMREX_SPACEDIM >& gmres_rhs_u,
                   MultiFab& gmres_rhs_p) {

    // Loop over boxes (note that mfi takes a cell-centered multifab as an argument)
    for ( MFIter mfi(gmres_rhs_p,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        AMREX_D_TERM(const Array4<Real> & gru = gmres_rhs_u[0].array(mfi);,
                     const Array4<Real> & grv = gmres_rhs_u[1].array(mfi);,
                     const Array4<Real> & grw = gmres_rhs_u[2].array(mfi););

        AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                     const Box & bx_y = mfi.nodaltilebox(1);,
                     const Box & bx_z = mfi.nodaltilebox(2););

        amrex::ParallelFor(bx_x, bx_y,
#if (AMREX_SPACEDIM == 3)
                           bx_z,
#endif
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // modify gmres_rhs_u[0] here
                // note for 2D runs, k=0
                /*
                if (i == 15 && j == 15) {
                    gru(i,j,k) += 1000.;
                }
                */

            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // modify gmres_rhs_u[1] here

            }
#if (AMREX_SPACEDIM == 3)
          , [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // modify gmres_rhs_u[2] here

            }
#endif
            );
    }
}
