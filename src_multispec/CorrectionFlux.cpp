#include "multispec_functions.H"
#include "multispec_functions_F.H"

#include "common_functions.H"

#include "multispec_namespace.H"

using namespace multispec;
using namespace amrex;

void CorrectionFlux(const MultiFab& rho, const MultiFab& rhotot,
		    std::array< MultiFab, AMREX_SPACEDIM >& flux)
{

  BL_PROFILE_VAR("CorrectionFlux()",CorrectionFlux);

      // Loop over boxes
    for (MFIter mfi(rho); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& validBox = mfi.validbox();

        correction_flux(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
			BL_TO_FORTRAN_FAB(rho[mfi]),
			BL_TO_FORTRAN_ANYD(rhotot[mfi]),
			BL_TO_FORTRAN_ANYD(flux[0][mfi]),
			BL_TO_FORTRAN_ANYD(flux[1][mfi])
#if (AMREX_SPACEDIM == 3)
			,BL_TO_FORTRAN_ANYD(flux[2][mfi])
#endif
			);
    }

}
