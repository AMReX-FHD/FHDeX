#include "compressible_functions.H"

#include "common_functions.H"

using namespace common;

void calculateTransportCoeffs(const MultiFab& prim_in, 
			      MultiFab& eta_in, MultiFab& zeta_in, MultiFab& kappa_in,
			      MultiFab& chi_in, MultiFab& Dij_in)
{
    BL_PROFILE_VAR("calculateTransportCoeffs()",calculateTransportCoeffs);

    // see comments in conservedPrimitiveConversions.cpp regarding alternate ways of declaring
    // thread shared and thread private arrays on GPUs
    // if the size is not known at compile time, alternate approaches are required
    // here we know the size at compile time
    
    // Loop over boxes
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        // grow the box by ngc
        const Box& bx = amrex::grow(mfi.tilebox(), ngc);

        const Array4<const Real>& prim = prim_in.array(mfi);

        const Array4<Real>& eta   =   eta_in.array(mfi);
        const Array4<Real>& zeta  =  zeta_in.array(mfi);
        const Array4<Real>& kappa = kappa_in.array(mfi);
        const Array4<Real>& chi   =   chi_in.array(mfi);
        const Array4<Real>& Dij   =   Dij_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
        
            GpuArray<Real,MAX_SPECIES> Yk_fixed;
            GpuArray<Real,MAX_SPECIES> Xk_fixed;
            
            Real sumYk = 0.;
            for (int n=0; n<nspecies; ++n) {
                Yk_fixed[n] = amrex::max(0.,amrex::min(1.,prim(i,j,k,6+n)));
                sumYk += Yk_fixed[n];
            }

            for (int n=0; n<nspecies; ++n) {
                Yk_fixed[n] /= sumYk;
            }

            // compute mole fractions from mass fractions
            GetMolfrac(Yk_fixed, Xk_fixed);

            IdealMixtureTransport(i,j,k, prim(i,j,k,0), prim(i,j,k,4), prim(i,j,k,5),
                                  Yk_fixed, Xk_fixed, eta(i,j,k), kappa(i,j,k), zeta(i,j,k),
                                  Dij, chi);

            // want this multiplied by rho for all times
            for (int kk=0; kk<nspecies; ++kk) {
                for (int ll=0; ll<nspecies; ++ll) {
                    int n = kk*nspecies + ll;
                    Dij(i,j,k,n) *= prim(i,j,k,0);
                }
            }

        });
    }
}
