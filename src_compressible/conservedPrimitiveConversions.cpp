#include "common_functions.H"
#include "compressible_functions.H"

void conservedToPrimitive(MultiFab& prim_in, const MultiFab& cons_in)
{
    BL_PROFILE_VAR("conservedToPrimitive()",conservedToPrimitive);

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

    // Loop over boxes
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<const Real>& cons = cons_in.array(mfi);
        const Array4<      Real>& prim = prim_in.array(mfi);

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
            prim(i,j,k,1) = cons(i,j,k,1)/cons(i,j,k,0);
            prim(i,j,k,2) = cons(i,j,k,2)/cons(i,j,k,0);
            prim(i,j,k,3) = cons(i,j,k,3)/cons(i,j,k,0);

            Real vsqr = prim(i,j,k,1)*prim(i,j,k,1) + prim(i,j,k,2)*prim(i,j,k,2) + prim(i,j,k,3)*prim(i,j,k,3);
            Real intenergy = cons(i,j,k,4)/cons(i,j,k,0) - 0.5*vsqr;

            Real sumYk = 0.;
            for (int n=0; n<nspecies; ++n) {
                Yk[n] = cons(i,j,k,5+n)/cons(i,j,k,0);
                Yk_fixed[n] = amrex::max(0.,amrex::min(1.,Yk[n]));
                sumYk += Yk_fixed[n];
            }

            for (int n=0; n<nspecies; ++n) {
                Yk_fixed[n] /= sumYk;
            }

            // update temperature in-place using internal energy
            GetTemperature(intenergy, Yk_fixed, prim(i,j,k,4));

            // compute mole fractions from mass fractions
            GetMolfrac(Yk, Xk);

            // mass fractions
            for (int n=0; n<nspecies; ++n) {
                prim(i,j,k,6+n) = Yk[n];
                prim(i,j,k,6+nspecies+n) = Xk[n];
            }

            GetPressureGas(prim(i,j,k,5), Yk, prim(i,j,k,0), prim(i,j,k,4));
        });

    } // end MFIter
}
