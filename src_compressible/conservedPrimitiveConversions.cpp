#include "common_functions.H"
#include "compressible_functions.H"

void conservedToPrimitive(MultiFab& prim_in, const MultiFab& cons_in)
{
    BL_PROFILE_VAR("conservedToPrimitive()",conservedToPrimitive);

    // from namelist
    int nspecies_gpu = nspecies;

    // from namelist
    Real Runiv_gpu = Runiv;

    // from namelist
    amrex::Vector<amrex::Real> hcv_vect_host(nspecies);
    for (int n=0; n<nspecies; ++n) {
        hcv_vect_host[n] = hcv[n];
    }
    
    amrex::Gpu::DeviceVector<amrex::Real> hcv_vect(nspecies);
    Real const * const AMREX_RESTRICT hcv_gpu = hcv_vect.dataPtr();
    amrex::Gpu::copy(amrex::Gpu::HostToDevice,hcv_vect_host.begin(),hcv_vect_host.end(),
                     hcv_vect.begin());


    // FIX
    /*
    // from namelist
    amrex::Gpu::DeviceVector<amrex::Real> molmass_vect(nspecies);
    for (int n=0; n<nspecies; ++n) {
        molmass_vect[n] = molmass[n];
    }
    Real const * const AMREX_RESTRICT molmass_gpu = molmass_vect.dataPtr();
    */
    
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

        // FIX
        /*
        amrex::Gpu::ManagedVector<amrex::Real> Yk_vect(nspecies);
        amrex::Gpu::ManagedVector<amrex::Real> Yk_fixed_vect(nspecies);
        amrex::Gpu::ManagedVector<amrex::Real> Xk_vect(nspecies);

        Real * const AMREX_RESTRICT Yk = Yk_vect.dataPtr();
        Real * const AMREX_RESTRICT Yk_fixed = Yk_fixed_vect.dataPtr();
        Real * const AMREX_RESTRICT Xk = Xk_vect.dataPtr();
        */

        // this is allocated on the DEVICE (no page faults)
        FArrayBox Yk_fab(bx,nspecies);
        const Array4<Real>& Yk = Yk_fab.array(mfi);
        // make sure Yk_fab doesn't go out of scope once the CPU finishes and GPU isn't done
        auto Yk_eli = Yk_fab.elixir();
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // option if MAX_SPECIES is a compile-time constant
//            Real Yk[MAX_SPECIES];
            
            prim(i,j,k,0) = cons(i,j,k,0);
            prim(i,j,k,1) = cons(i,j,k,1)/cons(i,j,k,0);
            prim(i,j,k,2) = cons(i,j,k,2)/cons(i,j,k,0);
            prim(i,j,k,3) = cons(i,j,k,3)/cons(i,j,k,0);

            Real vsqr = prim(i,j,k,1)*prim(i,j,k,1) + prim(i,j,k,2)*prim(i,j,k,2) + prim(i,j,k,3)*prim(i,j,k,3);
            Real intenergy = cons(i,j,k,4)/cons(i,j,k,0) - 0.5*vsqr;

            Real sumYk = 0.;
            for (int n=0; n<nspecies_gpu; ++n) {
                Yk(i,j,k,n) = cons(i,j,k,5+n)/cons(i,j,k,0);
                Yk_fixed[n] = std::max(0.,std::min(1.,Yk(i,j,k,n)));
                sumYk += Yk_fixed[n];
            }
            
            for (int n=0; n<nspecies_gpu; ++n) {
                Yk_fixed[n] = Yk_fixed[n] / sumYk;
            }

            // update temperature in-place using internal energy
            GetTemperature(i,j,k,intenergy,Yk_fixed,prim(i,j,k,4),nspecies_gpu,hcv_gpu);

            // compute mole fractions from mass fractions
            GetMolfrac(Yk,Xk,nspecies_gpu,molmass_gpu);

            // mass fractions
            for (int n=0; n<nspecies_gpu; ++n) {
                prim(i,j,k,6+n) = Yk(i,j,k,n);
                prim(i,j,k,6+nspecies_gpu+n) = Xk[n];
            }

            GetPressureGas(prim(i,j,k,5),Yk,prim(i,j,k,0),prim(i,j,k,4),nspecies_gpu,Runiv_gpu,molmass_gpu);
            
        });
    } // end MFIter
}
