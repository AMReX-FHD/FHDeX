#include "compressible_functions.H"
#include "common_functions.H"

using namespace common;
using namespace compressible;

void calculateTransportCoeffs(const MultiFab& prim_in,
                              MultiFab& eta_in, MultiFab& zeta_in, MultiFab& kappa_in,
                              MultiFab& chi_in, MultiFab& Dij_in)
{
    BL_PROFILE_VAR("calculateTransportCoeffs()",calculateTransportCoeffs);

    // see comments in conservedPrimitiveConversions.cpp regarding alternate ways of declaring
    // thread shared and thread private arrays on GPUs
    // if the size is not known at compile time, alternate approaches are required
    // here we know the size at compile time

    amrex::IntVect ng_temp;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        ng_temp[d] = (do_1D == 1) ? 1 : ngc[d];
    }

    // Loop over boxes
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        // grow the box by ngc
        const Box& bx = amrex::grow(mfi.tilebox(), ng_temp);

        const Array4<const Real>& prim = prim_in.array(mfi);

        const Array4<Real>& eta   =   eta_in.array(mfi);
        const Array4<Real>& zeta  =  zeta_in.array(mfi);
        const Array4<Real>& kappa = kappa_in.array(mfi);
        const Array4<Real>& chi   =   chi_in.array(mfi);
        const Array4<Real>& Dij   =   Dij_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            GpuArray<Real,MAX_SPECIES> Yk;
            for (int n=0; n<nspecies; ++n) {
                Yk[n] = prim(i,j,k,6+n);
            }

            amrex::GpuArray<amrex::Real,MAX_SPECIES*MAX_SPECIES> Dloc;
            amrex::GpuArray<amrex::Real,MAX_SPECIES> chiloc;

            TransportCoeffs(prim(i,j,k,0), prim(i,j,k,4), prim(i,j,k,5),
                            Yk, eta(i,j,k), kappa(i,j,k), zeta(i,j,k),
                            Dloc, chiloc);

            for (int kk=0; kk<nspecies; ++kk) {
                chi(i,j,k,kk) = chiloc[kk];
            }
            // want this multiplied by rho for all times (rho*D_tilde = rho*Y*D)
            for (int kk=0; kk<nspecies; ++kk) {
                for (int ll=0; ll<nspecies; ++ll) {
                    int n = kk*nspecies + ll;
                    Dij(i,j,k,n) = Dloc[n];
                }
            }

        });
    }
}

