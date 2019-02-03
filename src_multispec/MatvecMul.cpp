#include "multispec_functions.H"
#include "multispec_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "multispec_namespace.H"
#include "common_namespace.H"

using namespace multispec;
using namespace common;
using namespace amrex;

void MatvecMul(MultiFab& x,
	       const MultiFab& A,
	       const int& nc)
{

    BL_PROFILE_VAR("MatvecMul()",MatvecMul);

    // Loop over boxes
    for (MFIter mfi(x); mfi.isValid(); ++mfi) {

        // Create a box that matches the NODALITY of MultiFab
        const Box& validBox = mfi.validbox();

        matvec_mul(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
		   BL_TO_FORTRAN_ANYD(x[mfi]),
		   BL_TO_FORTRAN_FAB(A[mfi]));
    }

}
