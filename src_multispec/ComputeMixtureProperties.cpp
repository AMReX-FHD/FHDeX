#include "multispec_functions.H"
#include "multispec_functions_F.H"

#include "common_functions.H"

#include "multispec_namespace.H"

using namespace multispec;
using namespace amrex;

void ComputeMixtureProperties(const MultiFab& rho,
			      const MultiFab& rhotot,
			      MultiFab& D_bar,
			      MultiFab& D_therm,
			      MultiFab& Hessian)
{

    BL_PROFILE_VAR("ComputeMixtureProperties()",ComputeMixtureProperties);

    // Loop over boxes
    for (MFIter mfi(rho); mfi.isValid(); ++mfi) {

        // Create cell-centered box
        const Box& validBox = mfi.validbox();

        mixture_properties_mass(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
				BL_TO_FORTRAN_FAB(rho[mfi]),
				BL_TO_FORTRAN_ANYD(rhotot[mfi]),
				BL_TO_FORTRAN_ANYD(D_bar[mfi]),
				BL_TO_FORTRAN_ANYD(D_therm[mfi]),
				BL_TO_FORTRAN_ANYD(Hessian[mfi]));
    }

}
