#include "multispec_functions.H"

void ComputeMolconcMolmtot(const MultiFab& rho_in,
			   const MultiFab& rhotot_in,
			   MultiFab& molarconc_in,
			   MultiFab& molmtot_in)
{

    BL_PROFILE_VAR("ComputeMolconcMolmtot()",ComputeMolconcMolmtot);

    int ng = molarconc_in.nGrow();
        
    // Loop over boxes
    for (MFIter mfi(rho_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& bx = mfi.growntilebox(ng);

        const Array4<const Real>& rho = rho_in.array(mfi);
        const Array4<const Real>& rhotot = rhotot_in.array(mfi);
        const Array4<      Real>& molarconc = molarconc_in.array(mfi);
        const Array4<      Real>& molmtot = molmtot_in.array(mfi);
        

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            GpuArray<Real, MAX_SPECIES> RhoN;
            GpuArray<Real, MAX_SPECIES> MolarConcN;

            for (int n=0; n<nspecies; ++n ){
                RhoN[n] = rho(i,j,k,n);
            }


            ComputeMolconcMolmtotLocal(nspecies, molmass, 
                            RhoN, rhotot(i,j,k),          
                            MolarConcN, molmtot(i,j,k));

            for (int n=0; n<nspecies; ++n ){
                molarconc(i,j,k,n) = MolarConcN[n] ;
            }

        });


    }

}

void ComputeGamma(const MultiFab& molarconc_in,
		      const MultiFab& Hessian_in,
          MultiFab& Gamma_in)
{
  
    BL_PROFILE_VAR("ComputeGamma()",ComputeGamma);

    int ng = Gamma_in.nGrow();
    
    // Loop over boxes
    for (MFIter mfi(Gamma_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& bx = mfi.growntilebox(ng);

        const Array4<const Real>& molarconc = molarconc_in.array(mfi);
        const Array4<const Real>& Hessian = Hessian_in.array(mfi); 
        const Array4<      Real>& Gamma = Gamma_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            GpuArray<Real, MAX_SPECIES> MolarConcN;
            Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> GammaN; 
            Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> HessianN;


            // Read MultiFab data into arrays
            for (int n=0; n<nspecies; ++n ){
                MolarConcN[n] = molarconc(i,j,k,n);

                for (int m=0; m<nspecies; ++m){ 
                    GammaN(m+1,n+1) = Gamma(i,j,k,n*nspecies+m);  
                    HessianN(m+1,n+1) = Hessian(i,j,k,n*nspecies+m); 
                } 
            }
        
            ComputeGammaLocal(MolarConcN, HessianN, GammaN, nspecies);

            // Write back to MultiFab
            for (int n=0; n<nspecies; ++n ){
                for (int m=0; m<nspecies; ++m){ 
                    Gamma(i,j,k,n*nspecies+m) = GammaN(m+1,n+1);  
                } 
            }
        });
    }
}

void ComputeRhoWChi(const MultiFab& rho,
		    const MultiFab& rhotot,
		    const MultiFab& molarconc,
		    MultiFab& rhoWchi,
		    const MultiFab& D_bar)
{
    BL_PROFILE_VAR("ComputeRhoWChi()",ComputeRhoWChi);

    int ng = rhoWchi.nGrow();
    
    // Loop over boxes
    for (MFIter mfi(rhoWchi,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& bx = mfi.growntilebox(ng);

        compute_rhoWchi(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			BL_TO_FORTRAN_ANYD(rho[mfi]),
			BL_TO_FORTRAN_ANYD(rhotot[mfi]),
			BL_TO_FORTRAN_ANYD(molarconc[mfi]),
			BL_TO_FORTRAN_ANYD(rhoWchi[mfi]),
			BL_TO_FORTRAN_ANYD(D_bar[mfi]));
    }

}

void ComputeZetaByTemp(const MultiFab& molarconc,
 		       const MultiFab& D_bar,
 		       const MultiFab& Temp,
 		       MultiFab& zeta_by_Temp,
 		       const MultiFab& D_therm)
{
    BL_PROFILE_VAR("ComputeZetaByTemp()",ComputeZetaByTemp);

    int ng = zeta_by_Temp.nGrow();

    // Loop over boxes
    for (MFIter mfi(zeta_by_Temp,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& bx = mfi.growntilebox(ng);

        compute_zeta_by_Temp(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                             BL_TO_FORTRAN_ANYD(molarconc[mfi]),
                             BL_TO_FORTRAN_ANYD(D_bar[mfi]),
                             BL_TO_FORTRAN_ANYD(Temp[mfi]),
                             BL_TO_FORTRAN_ANYD(zeta_by_Temp[mfi]),
                             BL_TO_FORTRAN_ANYD(D_therm[mfi]));
    }
}

void ComputeSqrtLonsagerFC(const MultiFab& rho, const MultiFab& rhotot,
                           std::array< MultiFab, AMREX_SPACEDIM >& sqrtLonsager_fc,
                           const Geometry& geom)
{
    BL_PROFILE_VAR("ComputeSqrtLonsagerFC()",ComputeSqrtLonsagerFC);

    const Real* dx = geom.CellSize();

    // for GPU later
    // const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    
    // Loop over boxes
    for (MFIter mfi(rho); mfi.isValid(); ++mfi) {

        // note: tiling or GPU-ing requires nodal tileboxes and changes to
        // loop indices in fortran
        
        // Create cell-centered box
        const Box& validBox = mfi.validbox();

        compute_sqrtLonsager_fc(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                                BL_TO_FORTRAN_ANYD(rho[mfi]),
                                BL_TO_FORTRAN_ANYD(rhotot[mfi]),
                                BL_TO_FORTRAN_ANYD(sqrtLonsager_fc[0][mfi]),
                                BL_TO_FORTRAN_ANYD(sqrtLonsager_fc[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                                BL_TO_FORTRAN_ANYD(sqrtLonsager_fc[2][mfi]),
#endif
                                dx);
    }

}
