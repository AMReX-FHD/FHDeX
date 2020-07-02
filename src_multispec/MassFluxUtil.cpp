#include "multispec_functions.H"

void ComputeMolconcMolmtot(const MultiFab& rho,
			   const MultiFab& rhotot,
			   MultiFab& molarconc,
			   MultiFab& molmtot)
{

    BL_PROFILE_VAR("ComputeMolconcMolmtot()",ComputeMolconcMolmtot);

    int ng = molarconc.nGrow();
        
    // Loop over boxes
    for (MFIter mfi(rho,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& bx = mfi.growntilebox(ng);

        compute_molconc_molmtot(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				BL_TO_FORTRAN_ANYD(rho[mfi]),
				BL_TO_FORTRAN_ANYD(rhotot[mfi]),
				BL_TO_FORTRAN_ANYD(molarconc[mfi]),
				BL_TO_FORTRAN_ANYD(molmtot[mfi]));
    }

}

void ComputeGamma(const MultiFab& molarconc,
		  const MultiFab& Hessian,
		  MultiFab& Gamma)
{
  
    BL_PROFILE_VAR("ComputeGamma()",ComputeGamma);

    int ng = Gamma.nGrow();
    
    // Loop over boxes
    for (MFIter mfi(Gamma,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& bx = mfi.growntilebox(ng);

        compute_Gamma(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		      BL_TO_FORTRAN_ANYD(molarconc[mfi]),
		      BL_TO_FORTRAN_ANYD(Hessian[mfi]),
		      BL_TO_FORTRAN_ANYD(Gamma[mfi]));
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
