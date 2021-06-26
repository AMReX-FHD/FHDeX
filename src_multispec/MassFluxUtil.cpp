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

            //Print() << ">>>> *MolarConcN[0] " << MolarConcN[0] << std::endl;
            //Print() << ">>>> *MolarConcN[1] " << MolarConcN[1] << std::endl;
            //Print() << ">>>> *MolarConcN[2] " << MolarConcN[2] << std::endl;
            //Print() << ">>>> *molmtot(i,j,k) " << molmtot(i,j,k) << std::endl;
            //Abort();
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
        const Array4<const Real>& Hessian = Hessian_in.array(mfi); 
        const Array4<      Real>& Gamma = Gamma_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            GpuArray<Real, MAX_SPECIES> MolarConcN;
            Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> GammaN; // verified that Array2D are ok to use
            Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> HessianN;

           //Print() << "Sigfault sandwhich 0 " << std::endl; 

            for (int n=0; n<nspecies; ++n ){
                MolarConcN[n] = molarconc(i,j,k,n);

                for (int m=0; m<nspecies; ++m){ 

                    //Print() << nspecies << " n = " << n << " m = " << m << " iter = " << n*nspecies+m+1 << std::endl;


                    GammaN(n+1,m+1) = Gamma(i,j,k,n*nspecies+m);  //Gamma's index starts at 1 right?
                    HessianN(n+1,m+1) = Hessian(i,j,k,n*nspecies+m); 
                } 
                //HessianN[n] = Hessian(i,j,k,n);   
                //GammaN[n] = Gamma(i,j,k,n);
            }
        
           //Print() << "Sigfault sandwhich 1 " << std::endl; 

            //ComputeGammaLocal(MolarConcN, HessianN, GammaN, nspecies);
            // Fill this in later


            Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> I;
            Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> X_xxt;

            //if ((use_multiphase == 1) && (nspecies == 2)){ 
            //if (true){ 
            GammaN(1,1) = 3.0;
            GammaN(1,2) = 5.0;
            GammaN(2,1) = 7.0;
            GammaN(2,2) = 11.0;
            n_gex = 3.3;
            alpha_gex = 0.72;

            Print() << GammaN(1,1) << " " << GammaN(1,2) << std::endl;
            Print() << GammaN(2,1) << " " << GammaN(2,2) << std::endl;


            //Print() << "w1 " << w1 << std::endl;
            //Print() << "w2 " << w2 << std::endl;
            Print() << "n_gex " << n_gex << " alpha_gex " << alpha_gex << std::endl;

            if ((use_multiphase == 1) && (nspecies == 2)){ 
                            
                Real w1 = MolarConcN[0];
                Real w2 = MolarConcN[1];



                if (std::abs(w1+w2-1.0) > 1e-14){  //Tested and working
                    Print() << "Mole fractions do not add up in gamma computation" << std::endl; 
                }
                if (w1 < 0){ 
                    w1 = 0.0;
                    w2 = 1.0;
                }
                if (w2 < 0){ 
                    w1 = 1.0;
                    w2 = 0.0;
                }

                //These calculations were tested -- working 
                GammaN(1,2) = w1 * n_gex * n_gex * alpha_gex * std::pow(w1,n_gex-1) * std::pow(w2,n_gex-1);
                GammaN(2,1) = w2 * n_gex * n_gex * alpha_gex * std::pow(w2,n_gex-1) * std::pow(w1,n_gex-1);
                GammaN(1,1) = 1.0 + w1 * n_gex * (n_gex-1) * alpha_gex * std::pow(w1,n_gex-2) * std::pow(w2,n_gex); 
                GammaN(2,2) = 1.0 + w2 * n_gex * (n_gex-1) * alpha_gex * std::pow(w2,n_gex-2) * std::pow(w1,n_gex); 


            } else {
            ////construct identity matrix
            //    for (int n=1; n<=nspecies; ++n ){
            //        I(n,n) = 1.0;
            //    }

                //populate X_xxt
                //if (is_ideal_mixture == 1){
                if (is_ideal_mixture == -1){
                    for (int n=1; n<=nspecies; ++n ){
                       for (int m=1; m<=nspecies; ++m ){
                           X_xxt(n,m) = 0.0;
                       }
                    }
                } else {
                    for (int row=1; row<=nspecies;++row ){
                            Print() << "row " << row << std::endl;
                        // diagonal entries
                        X_xxt(row,row) = MolarConcN[row-1] - std::pow(MolarConcN[row-1],2);
                        for (int column=1; column<=row-1; ++column ){
                            Print() << "column " << column << std::endl;
                            // form x*traspose(x) off diagonals -- is x the MolarConc vectorT?
                            X_xxt(row,column) = -MolarConcN[row-1]*MolarConcN[column-1];
                            //symmetric
                            X_xxt(column,row) = X_xxt(row,column);
                        }
                    }
                }
            }
            //
                //HACK HACK HACK
            Print() << ">>>> MolarConcN[0] " << MolarConcN[0] << std::endl;
            Print() << ">>>> MolarConcN[1] " << MolarConcN[1] << std::endl;
            Print() << ">>>> MolarConcN[2] " << MolarConcN[2] << std::endl;
            Print() << std::endl;
            Print() << ">>>> X_xxt(1,1) " << X_xxt(1,1) << std::endl;
            Print() << ">>>> X_xxt(1,2) " << X_xxt(1,2) << std::endl;
            Print() << ">>>> X_xxt(1,3) " << X_xxt(1,3) << std::endl;
            Print() << ">>>> X_xxt(2,1) " << X_xxt(2,1) << std::endl;
            Print() << ">>>> X_xxt(2,2) " << X_xxt(2,2) << std::endl;
            Print() << ">>>> X_xxt(2,3) " << X_xxt(2,3) << std::endl;
            Print() << ">>>> X_xxt(3,1) " << X_xxt(3,1) << std::endl;
            Print() << ">>>> X_xxt(3,2) " << X_xxt(3,2) << std::endl;
            Print() << ">>>> X_xxt(3,3) " << X_xxt(3,3) << std::endl;
            Abort();
            
           //Print() << "Sigfault sandwhich 2 " << std::endl; 
              
            for (int row=1; row<=nspecies; ++row){
                for (int column=1; column<=nspecies; ++column){

                    //Print() << "Row = " << row << " Column = " << column << std::endl;

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
             
           //Print() << "Sigfault sandwhich 3 " << std::endl; 
            

            for (int n=0; n<nspecies; ++n ){
                for (int m=0; m<nspecies; ++m){ 
                    Gamma(i,j,k,n*nspecies+m) = GammaN(n+1,m+1);  
                } 
            }


            //Print() << GammaN(1,1) << " " << GammaN(1,2) << " " << GammaN(1,3) << std::endl;
            //Print() << GammaN(2,1) << " " << GammaN(2,2) << " " << GammaN(2,3) << std::endl;
            //Print() << GammaN(3,1) << " " << GammaN(3,2) << " " << GammaN(3,3) << std::endl;

        });

        //Print() << "end parfor" << std::endl;

        //EP-Ends here
    }
    //Abort();

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
