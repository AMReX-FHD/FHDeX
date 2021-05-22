#include "multispec_test_functions.H"

using namespace amrex;

void InitRhoUmac(std::array< MultiFab, AMREX_SPACEDIM >& umac,
                 MultiFab& rho_in,
                 const Geometry& geom)
{
    
    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    GpuArray<Real,AMREX_SPACEDIM> center;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        center[d] = 0.5*(prob_hi[d]-prob_lo[d]);
    }

    BoxArray ba = rho_in.boxArray();
    DistributionMapping dmap = rho_in.DistributionMap();

    MultiFab conc(ba,dmap,nspecies,0);

    // set velocity to zero; overwrite below if needed
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umac[d].setVal(0.);
    }
    
    for (MFIter mfi(rho_in); mfi.isValid(); ++mfi ) {

        Box bx = mfi.tilebox();

        const Array4<Real>& c = conc.array(mfi);

        if (prob_type == 1) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
            
            });
        } else if (prob_type == 2) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
            
            });
        } else if (prob_type == 4) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
            
            });
        } else if (prob_type == 12) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
            
            });
        } else if (prob_type == 15) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
            
            });
        } else {
            Abort("Init.cpp: Invalid prob_type");
        }
    }

    for (MFIter mfi(rho_in); mfi.isValid(); ++mfi ) {

        Box bx = mfi.tilebox();

        const Array4<Real>& c = conc.array(mfi);
        const Array4<Real>& rho = rho_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // set final c_i such that sumtot(c_i) = 1 to within roundoff
            Real sumtot = 0.;
            for (int n=0; n<nspecies-1; ++n) {
                sumtot = sumtot + c(i,j,k,n);
            }

            c(i,j,k,nspecies-1) = 1. - sumtot;

            Real rho_total;
            
            // calculate rho_total from eos
            if (algorithm_type == 6) {
                rho_total = rho0;
            } else {
                sumtot = 0.;
                for (int n=0; n<nspecies; ++n) {
                    // sumtot represents rhoinv
                    sumtot = sumtot + c(i,j,k,n)/rhobar[n];
                }
                rho_total = 1./sumtot;
            }

            // calculate rho_i
            for (int n=0; n<nspecies; ++n) {
                rho(i,j,k,n) = rho_total*c(i,j,k,n);
            }
                    
        });
     }
}
