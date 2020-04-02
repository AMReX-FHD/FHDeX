#include "multispec_functions.H"

void ComputeMolconcMolmtot(const MultiFab& rho,
			   const MultiFab& rhotot,
			   MultiFab& molarconc,
			   MultiFab& molmtot)
{

    BL_PROFILE_VAR("ComputeMolconcMolmtot()",ComputeMolconcMolmtot);

    // Loop over boxes
    for (MFIter mfi(rho); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& validBox = mfi.validbox();

        compute_molconc_molmtot(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
				BL_TO_FORTRAN_FAB(rho[mfi]),
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

    // Loop over boxes
    for (MFIter mfi(molarconc); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& validBox = mfi.validbox();

        compute_Gamma(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
		      BL_TO_FORTRAN_FAB(molarconc[mfi]),
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

    // Loop over boxes
    for (MFIter mfi(rho); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& validBox = mfi.validbox();

        compute_rhoWchi(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
			BL_TO_FORTRAN_FAB(rho[mfi]),
			BL_TO_FORTRAN_ANYD(rhotot[mfi]),
			BL_TO_FORTRAN_ANYD(molarconc[mfi]),
			BL_TO_FORTRAN_ANYD(rhoWchi[mfi]),
			BL_TO_FORTRAN_ANYD(D_bar[mfi]));
    }

}

// void ComputeZetaByTemp(const MultiFab& molarconc,
// 		       const MultiFab& D_bar,
// 		       const MultiFab& Temp,
// 		       MultiFab& zeta_by_Temp,
// 		       const MultiFab& D_therm)
// {
  
//     BL_PROFILE_VAR("ComputeZetaByTemp()",ComputeZetaByTemp);

//     // Loop over boxes
//     for (MFIter mfi(molarconc); mfi.isValid(); ++mfi) {

//         // Create cell-centered box
//         const Box& validBox = mfi.validbox();

//         compute_zeta_by_Temp(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
// 			     BL_TO_FORTRAN_FAB(molarconc[mfi]),
// 			     BL_TO_FORTRAN_ANYD(D_bar[mfi]),
// 			     BL_TO_FORTRAN_ANYD(Temp[mfi]),
// 			     BL_TO_FORTRAN_ANYD(zeta_by_Temp[mfi]),
// 			     BL_TO_FORTRAN_ANYD(D_therm[mfi]));
//     }

// }
