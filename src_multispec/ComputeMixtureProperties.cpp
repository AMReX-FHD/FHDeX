#include "multispec_functions.H"

void ComputeMixtureProperties(const MultiFab& rho,
			      const MultiFab& rhotot,
			      MultiFab& D_bar,
			      MultiFab& D_therm,
			      MultiFab& Hessian)
{

    BL_PROFILE_VAR("ComputeMixtureProperties()",ComputeMixtureProperties);

    int ng = D_bar.nGrow();
    
    // Loop over boxes
    for (MFIter mfi(D_bar,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& bx = mfi.growntilebox(ng);

        mixture_properties_mass(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				BL_TO_FORTRAN_ANYD(rho[mfi]),
				BL_TO_FORTRAN_ANYD(rhotot[mfi]),
				BL_TO_FORTRAN_ANYD(D_bar[mfi]),
				BL_TO_FORTRAN_ANYD(D_therm[mfi]),
				BL_TO_FORTRAN_ANYD(Hessian[mfi]));
    }

}
