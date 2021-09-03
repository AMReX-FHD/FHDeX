#include "multispec_functions.H"

void ComputeMixtureProperties(const MultiFab& rho_in,
			      const MultiFab& rhotot_in,
			      MultiFab& D_bar_in,
			      MultiFab& D_therm_in,
			      MultiFab& Hessian_in)
{

    BL_PROFILE_VAR("ComputeMixtureProperties()",ComputeMixtureProperties);

    int ng = D_bar_in.nGrow();
    
    // Loop over boxes
    for (MFIter mfi(D_bar_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& bx = mfi.growntilebox(ng);

        const Array4<const Real>& rho_n = rho_in.array(mfi);
        const Array4<const Real>& rhotot = rhotot_in.array(mfi);
        const Array4<      Real>& D_bar_nn = D_bar_in.array(mfi);
        const Array4<      Real>& D_therm_n = D_therm_in.array(mfi);
        const Array4<      Real>& Hessian_nn = Hessian_in.array(mfi); 
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            GpuArray<Real, MAX_SPECIES> Rho;
            Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> DBar;
            GpuArray<Real, MAX_SPECIES> DTherm;
            Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> Hessian;

            //Read in multifab data
            for (int n=0; n<nspecies; ++n){
                Rho[n] = rho_n(i,j,k,n);
                DTherm[n] = D_therm_n(i,j,k,n);
                for (int m=0; m<nspecies; ++m){
                    DBar(m+1, n+1) = D_bar_nn(i,j,k,n*nspecies+m);
                    Hessian(m+1, n+1) = Hessian_nn(i,j,k,n*nspecies+m);
                }
            }

            MixturePropsMassLocal(Rho, rhotot(i,j,k), DBar, DTherm, Hessian);

            //write back to multifab 
            for (int n=0; n<nspecies; ++n ){
                D_therm_n(i,j,k,n) = DTherm[n];
                for (int m=0; m<nspecies; ++m){
                    D_bar_nn(i,j,k,n*nspecies+m) = DBar(m+1, n+1);
                    Hessian_nn(i,j,k,n*nspecies+m) = Hessian(m+1, n+1);
                }
            }

        });
    }

}
