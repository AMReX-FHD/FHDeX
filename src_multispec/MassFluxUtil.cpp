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

// EP
void ComputeGamma(const MultiFab& molarconc_in,
		      const MultiFab& Hessian_in,
          MultiFab& Gamma_in)
//void ComputeGamma(const MultiFab& molarconc,
//		  const MultiFab& Hessian,
//		  MultiFab& Gamma)
{
  
    BL_PROFILE_VAR("ComputeGamma()",ComputeGamma);

    //int ng = Gamma.nGrow();
    int ng = Gamma_in.nGrow();
    
    // Loop over boxes
    //for (MFIter mfi(Gamma,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    for (MFIter mfi(Gamma_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& bx = mfi.growntilebox(ng);

        //compute_Gamma(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		    //  BL_TO_FORTRAN_ANYD(molarconc[mfi]),
		    //  BL_TO_FORTRAN_ANYD(Hessian[mfi]),
		    //  BL_TO_FORTRAN_ANYD(Gamma[mfi]));

        //EP-Starts here
        const Array4<const Real>& molarconc = molarconc_in.array(mfi);
        const Array4<const Real>& Hessian = Hessian_in.array(mfi); //What if I want a 6-dim array?
        const Array4<      Real>& Gamma = Gamma_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            GpuArray<Real, MAX_SPECIES> MolarConcN;
            //GpuArray<Real, MAX_SPECIES> HessianN;
            //GpuArray<Real, MAX_SPECIES> GammaN;
            Array2D<Real, MAX_SPECIES, MAX_SPECIES*MAX_SPECIES> GammaN;
            Array2D<Real, MAX_SPECIES, MAX_SPECIES*MAX_SPECIES> HessianN;


            for (int n=0; n<nspecies; ++n ){
                MolarConcN[n] = molarconc(i,j,k,n);
                
                for (int m=0; m<nspecies*nspecies; ++m){ 
                    GammaN(n+1,m+1) = Gamma(i,j,k,n*nspecies+m+1);  //Gamma's index starts at 1 right?
                    HessianN(n+1,m+1) = Hessian(i,j,k,n*nspecies+m+1); 
                } 
                //HessianN[n] = Hessian(i,j,k,n);   
                //GammaN[n] = Gamma(i,j,k,n);
            }
        
            


            //ComputeGammaLocal(MolarConcN, HessianN, GammaN, nspecies);
            Array2D<Reak, MAX_SPECIES, MAX_SPECIES> I;
            Array2D<Reak, MAX_SPECIES, MAX_SPECIES> X_xxt;

            if ((use_multiphase == 1) && (nspecies == 2)){ 
                            
                Real w1 = MolarConcN[0];
                Real w2 = MolarConcN[1];

                if (abs(w1+w2-1.0) > 1e-14){ 
                    //verify this is the correct print statment
                    Print() << " mole fractions do not add up in gamma computation"; 
                }
                if (w1 < 0){ 
                    w1 = 0.0;
                    w2 = 1.0;
                }
                if (w2 < 0){ 
                    w1 = 1.0;
                    w2 = 0.0;
                }

                GammaN(1,2) = w1 * n_gex * n_gex * alpha_gex * pow(w1, n_gex-1) * pow(w2, n_gex-1);
                GammaN(2,1) = w2 * n_gex * n_gex * alpha_gex * pow(w1, n_gex-1) * pow(w2, n_gex-1);
                GammaN(1,1) = 1.0 + w1 + n_gex * (n_gex-1) * alpha_gex * pow(w1, n_gex-2) * pow(w2, n_gex); 
                GammaN(2,2) = 1.0 + w1 + n_gex * (n_gex-1) * alpha_gex * pow(w2, n_gex-2) * pow(w1, n_gex); 


            } else {
            ////construct identity matrix
            //    for (int n=1; n<=nspecies; ++n ){
            //        I(n,n) = 1.0;
            //    }

                //populate X_xxt
                if (is_ideal_mixture == 1){
                    for (int n=1; n<=nspecies; ++n ){
                       for (int m=1; m<=nspecies; ++m ){
                           X_xxt(n,m) = 0.0;
                       }
                    }
                } else {
                    for (int row=1; row<=nspecies;++row ){
                        // diagonal entries
                        X_xxt(row,row) = MolarConcN[row-1] - pow(MolarConcN[row-1],2);
                        for (int column=1; column<=nspecies; ++column ){
                            // form x*traspose(x) off diagonals -- is x the MolarConc vectorT?
                            X_xxt(row,column) = -MolarConcN[row-1]*MolarConcN[column-1];
                            //symmetric
                            X_xxt(column,row) = X_xxt(column,row);
                        }
                    }
                }
            }
            
              
            for (int row=1; row<=nspecies; ++row){
                for (int column=1; column<=nspecies; ++column){
                    if (row == column) {
                        GammaN(row,column) = 1.0;   // add the identity matrix
                    } else {
                        GammaN(row,column) = 0.0;   // intialize off-diagonals to 0
                    }
                    for (int n=1; n<=nspecies; ++n){
                        GammaN(row, column) += X_xxt(row,n) * HessianN(n,column);
                    }
                }
            }

            

            //This isn't going to work anymore
            //for (int n=0; n<nspecies; ++n ){
            //    Gamma(i,j,k,n) = GammaN[n];
            //}
            for (int n=0; n<nspecies; ++n ){
                for (int m=0; m<nspecies*nspecies; ++m){ 
                    Gamma(i,j,k,n*nspecies+m+1) = GammaN(n+1,m+1);  
                } 
            }

        });
        //EP-Ends here
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
