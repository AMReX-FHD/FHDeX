#include "common_functions.H"
#include "compressible_functions_stag.H"
#include "compressible_functions.H"

void conservedToPrimitiveStag(MultiFab& prim_in, std::array<MultiFab, AMREX_SPACEDIM>& velStag_in,
                              MultiFab& cons_in, const std::array<MultiFab, AMREX_SPACEDIM>& momStag_in)
{
    BL_PROFILE_VAR("conservedToPrimitiveStag()",conservedToPrimitiveStag);

    // from namelist
    int nspecies_gpu = nspecies;

    // from namelist
    Real Runiv_gpu = Runiv;

    // from namelist
    /* 
    // method 1 to create a thread shared array
    // must use if the size of the array is not known at compile time
    // note when passing this into a function, you need to use the type,
    // "Real const * const AMREX_RESTRICT"
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

    // method 2 to create a thread shared array
    // can use when the size of the array is known at compile-time
    GpuArray<Real,MAX_SPECIES> molmass_gpu;
    for (int n=0; n<nspecies; ++n) {
        molmass_gpu[n] = molmass[n];
    }
    
    // from namelist
    GpuArray<Real,MAX_SPECIES> hcv_gpu;
    for (int n=0; n<nspecies; ++n) {
        hcv_gpu[n] = hcv[n];
    }
    
    // Loop over boxes
    for ( MFIter mfi(prim_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.tilebox();

        //const Array4<const Real>& cons = cons_in.array(mfi);
        const Array4<      Real>& cons = cons_in.array(mfi);
        const Array4<      Real>& prim = prim_in.array(mfi);

        AMREX_D_TERM(Array4<Real const> const& momx = momStag_in[0].array(mfi);,
                     Array4<Real const> const& momy = momStag_in[1].array(mfi);,
                     Array4<Real const> const& momz = momStag_in[2].array(mfi););

        AMREX_D_TERM(const Array4<Real>& velx = velStag_in[0].array(mfi);,
                     const Array4<Real>& vely = velStag_in[1].array(mfi);,
                     const Array4<Real>& velz = velStag_in[2].array(mfi););

        /*
        // method 1 to create a thread private array
        // must use if the size of the array is not known at compile time
        // note when passing this into a function, you need to use the type,
        // const Array4<Real>&
        // this is allocated on the DEVICE (no page faults)
        FArrayBox Xk_fab(bx,nspecies);
        const Array4<Real>& Xk = Xk_fab.array();
        // make sure Xk_fab doesn't go out of scope once the CPU finishes and GPU isn't done
        auto Xk_eli = Xk_fab.elixir();
        */
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // method 2 to create a thread private array
            // can use when the size of the array is known at compile-time
            GpuArray<Real,MAX_SPECIES> Xk;
            GpuArray<Real,MAX_SPECIES> Yk;
            GpuArray<Real,MAX_SPECIES> Yk_fixed;
            

            prim(i,j,k,0) = cons(i,j,k,0);

            cons(i,j,k,1) = 0.5*(momx(i+1,j,k) + momx(i,j,k));
            cons(i,j,k,2) = 0.5*(momy(i,j+1,k) + momy(i,j,k));
            cons(i,j,k,3) = 0.5*(momz(i,j,k+1) + momz(i,j,k));

            prim(i,j,k,1) = cons(i,j,k,1)/cons(i,j,k,0);
            prim(i,j,k,2) = cons(i,j,k,2)/cons(i,j,k,0);
            prim(i,j,k,3) = cons(i,j,k,3)/cons(i,j,k,0);

            velx(i,j,k) = 2*momx(i,j,k)/(prim(i,j,k,0) + prim(i+1,j,k,0));
            vely(i,j,k) = 2*momy(i,j,k)/(prim(i,j,k,0) + prim(i,j+1,k,0));
            velz(i,j,k) = 2*momz(i,j,k)/(prim(i,j,k,0) + prim(i,j,k+1,0));

            // Do we need to calculate staggered velocities here as well? (from rho averaged to all faces) -- Ishan

            Real vsqr = prim(i,j,k,1)*prim(i,j,k,1) + prim(i,j,k,2)*prim(i,j,k,2) + prim(i,j,k,3)*prim(i,j,k,3);
            Real intenergy = cons(i,j,k,4)/cons(i,j,k,0) - 0.5*vsqr;

            Real sumYk = 0.;
            for (int n=0; n<nspecies_gpu; ++n) {
                Yk[n] = cons(i,j,k,5+n)/cons(i,j,k,0);
                Yk_fixed[n] = amrex::max(0.,amrex::min(1.,Yk[n]));
                sumYk += Yk_fixed[n];
            }
            
            for (int n=0; n<nspecies_gpu; ++n) {
                Yk_fixed[n] /= sumYk;
            }

            // update temperature in-place using internal energy
            GetTemperature(intenergy, Yk_fixed, prim(i,j,k,4), nspecies_gpu, hcv_gpu);

            // compute mole fractions from mass fractions
            GetMolfrac(i,j,k, Yk, Xk, nspecies_gpu, molmass_gpu);

            // mass fractions
            for (int n=0; n<nspecies_gpu; ++n) {
                prim(i,j,k,6+n) = Yk[n];
                prim(i,j,k,6+nspecies_gpu+n) = Xk[n];
            }

            GetPressureGas(prim(i,j,k,5), Yk, prim(i,j,k,0), prim(i,j,k,4),
                           nspecies_gpu, Runiv_gpu, molmass_gpu);
        });
        
    } // end MFIter
}
