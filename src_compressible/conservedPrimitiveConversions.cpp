#include "common_functions.H"
#include "compressible_functions.H"

void conservedToPrimitive(MultiFab& prim_in, const MultiFab& cons_in)
{
    BL_PROFILE_VAR("conservedToPrimitive()",conservedToPrimitive);

    // Loop over boxes
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.tilebox();

#if 0        
        cons_to_prim(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
                       cons_in[mfi].dataPtr(),  
                       prim_in[mfi].dataPtr());
#endif

        const Array4<const Real>& cons = cons_in.array(mfi);
        const Array4<      Real>& prim = prim_in.array(mfi);

        amrex::Gpu::ManagedVector<amrex::Real> Yk_vect(nspecies);
        amrex::Gpu::ManagedVector<amrex::Real> Yk_fixed_vect(nspecies);
        amrex::Gpu::ManagedVector<amrex::Real> Xk_vect(nspecies);

        Real * const AMREX_RESTRICT Yk = Yk_vect.dataPtr();
        Real * const AMREX_RESTRICT Yk_fixed = Yk_fixed_vect.dataPtr();
        Real * const AMREX_RESTRICT Xk = Xk_vect.dataPtr();

        int nspecies_gpu = nspecies;
        Real Runiv_gpu = Runiv;

        amrex::Gpu::ManagedVector<amrex::Real> molmass_vect(nspecies);
        for (int n=0; n<nspecies; ++n) {
            molmass_vect[n] = molmass[n];
        }
        Real const * const AMREX_RESTRICT molmass_gpu = molmass_vect.dataPtr();

        amrex::Gpu::ManagedVector<amrex::Real> hcv_vect(nspecies);
        for (int n=0; n<nspecies; ++n) {
            hcv_vect[n] = hcv[n];
        }
        Real const * const AMREX_RESTRICT hcv_gpu = hcv_vect.dataPtr();
                
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            prim(i,j,k,0) = cons(i,j,k,0);
            prim(i,j,k,1) = cons(i,j,k,1)/cons(i,j,k,0);
            prim(i,j,k,2) = cons(i,j,k,2)/cons(i,j,k,0);
            prim(i,j,k,3) = cons(i,j,k,3)/cons(i,j,k,0);

            Real vsqr = prim(i,j,k,1)*prim(i,j,k,1) + prim(i,j,k,2)*prim(i,j,k,2) + prim(i,j,k,3)*prim(i,j,k,3);
            Real intenergy = cons(i,j,k,4)/cons(i,j,k,0) - 0.5*vsqr;

            Real sumYk = 0.;
            for (int n=0; n<nspecies_gpu; ++n) {
                Yk[n] = cons(i,j,k,5+n)/cons(i,j,k,0);
                Yk_fixed[n] = std::max(0.,std::min(1.,Yk[n]));
                sumYk += Yk_fixed[n];
            }
            
            for (int n=0; n<nspecies_gpu; ++n) {
                Yk_fixed[n] = Yk_fixed[n] / sumYk;
            }

            // update temperature in-place using internal energy
            GetTemperature(intenergy,Yk_fixed,prim(i,j,k,4),nspecies_gpu,hcv_gpu);

            // compute mole fractions from mass fractions
            GetMolfrac(Yk,Xk,nspecies_gpu,molmass_gpu);

            // mass fractions
            for (int n=0; n<nspecies_gpu; ++n) {
                prim(i,j,k,6+n) = Yk[n];
                prim(i,j,k,6+nspecies_gpu+n) = Xk[n];
            }

            GetPressureGas(prim(i,j,k,5),Yk,prim(i,j,k,0),prim(i,j,k,4),nspecies_gpu,Runiv_gpu,molmass_gpu);
            
        });
    }
}
