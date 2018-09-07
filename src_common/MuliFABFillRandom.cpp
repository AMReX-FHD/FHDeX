#include "common_functions.H"
#include "common_functions_F.H"
#include "common_namespace.H"

using namespace common;
using namespace amrex;

void MultiFABFillRandom(MultiFab& mf, const int& comp, const amrex::Real& variance, const Geometry& geom)
{

    BL_PROFILE_VAR("MultiFABFillRandom()",MultiFABFillRandom);

    // Print() << "C++ hack:" << comp << "/" << mf.nComp() << "\n";

    // Loop over boxes
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {

        const Box& validBox = mfi.validbox();

	multifab_fill_random(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
			     BL_TO_FORTRAN_FAB(mf[mfi]), &comp);
    }
    
    // Scale standard gaussian samples by standard deviation
    mf.mult(sqrt(variance), comp, 1, 0);

    // Enforce boundary conditions on nodal boundaries & ghost cells
    mf.OverrideSync(geom.periodicity());
    mf.FillBoundary(geom.periodicity());
}
