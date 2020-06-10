#include "compressible_functions.H"

#include "common_functions.H"

using namespace common;

void calculateTransportCoeffs(const MultiFab& prim_in, 
			      MultiFab& eta_in, MultiFab& zeta_in, MultiFab& kappa_in,
			      MultiFab& chi_in, MultiFab& Dij_in)
{
    BL_PROFILE_VAR("calculateTransportCoeffs()",calculateTransportCoeffs);

    /*
    // nspecies from namelist
    int nspecies_gpu = nspecies;
    
    // molmass from namelist
    Vector<Real> molmass_vect_host(nspecies); // create a vector on the host and copy the values in
    for (int n=0; n<nspecies; ++n) {
        molmass_vect_host[n] = molmass[n];
    }
    Gpu::DeviceVector<Real> molmass_vect(nspecies); // create vector on GPU and copy values over
    Gpu::copy(Gpu::hostToDevice,
              molmass_vect_host.begin(),molmass_vect_host.end(),
              molmass_vect.begin());
    Real const * const AMREX_RESTRICT molmass_gpu = molmass_vect.dataPtr();  // pointer to data

    // compute molecular_mass by dividing molmass by Avogadro's
    Vector<Real> molecular_mass_vect_host(nspecies); // create a vector on the host and compute values
    for (int n=0; n<nspecies; ++n) {
        molecular_mass_vect_host[n] = molmass[n]*(k_B/Runiv);;
    }
    Gpu::DeviceVector<Real> molecular_mass_vect(nspecies); // create vector on GPU and copy values over
    Gpu::copy(Gpu::hostToDevice,
              molecular_mass_vect_host.begin(),molecular_mass_vect_host.end(),
              molecular_mass_vect.begin());
    Real const * const AMREX_RESTRICT molecular_mass_gpu = molecular_mass_vect.dataPtr();  // pointer to data
    */
    
    // Loop over boxes
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        /*
        // grow the box by ngc
        const Box& bx = amrex::grow(mfi.tilebox(), ngc);

        const Array4<const Real>& prim = prim_in.array(mfi);

        const Array4<Real>& eta   =   eta_in.array(mfi);
        const Array4<Real>& zeta  =  zeta_in.array(mfi);
        const Array4<Real>& kappa = kappa_in.array(mfi);
        const Array4<Real>& chi   =   chi_in.array(mfi);
        const Array4<Real>& Dij   =   Dij_in.array(mfi);

        // this is allocated on the DEVICE (no page faults)
        FArrayBox Yk_fixed_fab(bx,nspecies);
        const Array4<Real>& Yk_fixed = Yk_fixed_fab.array();
        // make sure Yk_fixed_fab doesn't go out of scope once the CPU finishes and GPU isn't done
        auto Yk_fixed_eli = Yk_fixed_fab.elixir();

        // this is allocated on the DEVICE (no page faults)
        FArrayBox Xk_fixed_fab(bx,nspecies);
        const Array4<Real>& Xk_fixed = Xk_fixed_fab.array();
        // make sure Xk_fixed_fab doesn't go out of scope once the CPU finishes and GPU isn't done
        auto Xk_fixed_eli = Xk_fixed_fab.elixir();

        // this is allocated on the DEVICE (no page faults)
        FArrayBox xxtr_fab(bx,nspecies);
        const Array4<Real>& xxtr = xxtr_fab.array();
        // make sure xxtr_fab doesn't go out of scope once the CPU finishes and GPU isn't done
        auto xxtr_eli = xxtr_fab.elixir();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real sumYk = 0.;
            for (int n=0; n<nspecies_gpu; ++n) {
                Yk_fixed(i,j,k,n) = std::max(0.,std::min(1.,prim(i,j,k,6+n)));
                sumYk += Yk_fixed(i,j,k,n);
            }

            for (int n=0; n<nspecies_gpu; ++n) {
                Yk_fixed(i,j,k,n) /= sumYk;
            }

            // compute mole fractions from mass fractions
            GetMolfrac(i,j,k, Yk_fixed, Xk_fixed, nspecies_gpu, molmass_gpu);

            IdealMixtureTransport(i,j,k, prim(i,j,k,0), prim(i,j,k,4), prim(i,j,k,5),
                                  Yk_fixed, Xk_fixed, eta(i,j,k), kappa(i,j,k), zeta(i,j,k),
                                  Dij, chi, xxtr, nspecies_gpu, molmass_gpu);

            // want this multiplied by rho for all times
            for (int kk=0; kk<nspecies_gpu; ++kk) {
                for (int ll=0; ll<nspecies_gpu; ++ll) {
                    int n = kk*nspecies_gpu + ll;
                    Dij(i,j,k,n) *= prim(i,j,k,0);
                }
            }

        });
        */        

        const Box& bx = mfi.validbox();

        makecoef(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
		 prim_in [mfi].dataPtr(),  
		 eta_in  [mfi].dataPtr(),  
		 zeta_in [mfi].dataPtr(),  
		 kappa_in[mfi].dataPtr(),
		 chi_in  [mfi].dataPtr(),
		 Dij_in  [mfi].dataPtr());
    }

    

    /*
    // Loop over boxes
    for ( MFIter mfi(prim); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.validbox();

        trans_coeffs(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
	 	     prim[mfi].dataPtr(),  
                     eta[mfi].dataPtr(),  
                     zeta[mfi].dataPtr(),  
                     kappa[mfi].dataPtr());
    }
    */
}
