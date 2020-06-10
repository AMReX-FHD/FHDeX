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
    Vector<Real> hcv_vect_host(nspecies); // create a vector on the host and copy the values in
    for (int n=0; n<nspecies; ++n) {
        hcv_vect_host[n] = hcv[n];
    }
    Gpu::DeviceVector<Real> hcv_vect(nspecies); // create vector on GPU and copy values over
    Gpu::copy(Gpu::hostToDevice,
              hcv_vect_host.begin(),hcv_vect_host.end(),
              hcv_vect.begin());
    
    Real const * const AMREX_RESTRICT hcv_gpu = hcv_vect.dataPtr(); // pointer to data

    // from namelist

    /* OLD WAY
    Vector<Real> molmass_vect_host(nspecies); // create a vector on the host and copy the values in
    for (int n=0; n<nspecies; ++n) {
        molmass_vect_host[n] = molmass[n];
    }
    Gpu::DeviceVector<Real> molmass_vect(nspecies); // create vector on GPU and copy values over
    Gpu::copy(Gpu::hostToDevice,
              molmass_vect_host.begin(),molmass_vect_host.end(),
              molmass_vect.begin());
    Real const * const AMREX_RESTRICT molmass_gpu = molmass_vect.dataPtr();  // pointer to data
    */
    GpuArray<Real,MAX_SPECIES> molmass_gpu;
    for (int n=0; n<nspecies; ++n) {
        molmass_gpu[n] = molmass[n];
    }











    
    
    
    // Loop over boxes
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.tilebox();

        const Array4<const Real>& cons = cons_in.array(mfi);
        const Array4<      Real>& prim = prim_in.array(mfi);

        // this is allocated on the DEVICE (no page faults)
        FArrayBox Yk_fab(bx,nspecies);
        const Array4<Real>& Yk = Yk_fab.array();
        // make sure Yk_fab doesn't go out of scope once the CPU finishes and GPU isn't done
        auto Yk_eli = Yk_fab.elixir();

        // this is allocated on the DEVICE (no page faults)
        FArrayBox Yk_fixed_fab(bx,nspecies);
        const Array4<Real>& Yk_fixed = Yk_fixed_fab.array();
        // make sure Yk_fixed_fab doesn't go out of scope once the CPU finishes and GPU isn't done
        auto Yk_fixed_eli = Yk_fixed_fab.elixir();
        
        // this is allocated on the DEVICE (no page faults)
        FArrayBox Xk_fab(bx,nspecies);
        const Array4<Real>& Xk = Xk_fab.array();
        // make sure Xk_fab doesn't go out of scope once the CPU finishes and GPU isn't done
        auto Xk_eli = Xk_fab.elixir();
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // option if MAX_SPECIES is a compile-time constant
            // Real Yk[MAX_SPECIES];

            prim(i,j,k,0) = cons(i,j,k,0);
            prim(i,j,k,1) = cons(i,j,k,1)/cons(i,j,k,0);
            prim(i,j,k,2) = cons(i,j,k,2)/cons(i,j,k,0);
            prim(i,j,k,3) = cons(i,j,k,3)/cons(i,j,k,0);

            Real vsqr = prim(i,j,k,1)*prim(i,j,k,1) + prim(i,j,k,2)*prim(i,j,k,2) + prim(i,j,k,3)*prim(i,j,k,3);
            Real intenergy = cons(i,j,k,4)/cons(i,j,k,0) - 0.5*vsqr;

            Real sumYk = 0.;
            for (int n=0; n<nspecies_gpu; ++n) {
                Yk(i,j,k,n) = cons(i,j,k,5+n)/cons(i,j,k,0);
                Yk_fixed(i,j,k,n) = std::max(0.,std::min(1.,Yk(i,j,k,n)));
                sumYk += Yk_fixed(i,j,k,n);
            }
            
            for (int n=0; n<nspecies_gpu; ++n) {
                Yk_fixed(i,j,k,n) /= sumYk;
            }

            // update temperature in-place using internal energy
            GetTemperature(i,j,k, intenergy, Yk_fixed, prim(i,j,k,4), nspecies_gpu, hcv_gpu);

            // compute mole fractions from mass fractions
            GetMolfrac(i,j,k, Yk, Xk, nspecies_gpu, molmass_gpu);

            // mass fractions
            for (int n=0; n<nspecies_gpu; ++n) {
                prim(i,j,k,6+n) = Yk(i,j,k,n);
                prim(i,j,k,6+nspecies_gpu+n) = Xk(i,j,k,n);
            }

            GetPressureGas(i,j,k, prim(i,j,k,5), Yk, prim(i,j,k,0), prim(i,j,k,4), nspecies_gpu, Runiv_gpu, molmass_gpu);
        });
        
    } // end MFIter
}
