#include "common_functions.H"
#include "common_functions_F.H"
#include "common_namespace.H"

using namespace common;
using namespace amrex;

void MultiFABFillRandom(MultiFab& mf, const int& comp, const amrex::Real& variance, const Geometry& geom)
{

    // Loop over boxes
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {

        const Box& validBox = mfi.validbox();

	multifab_fill_random(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
			     BL_TO_FORTRAN_FAB(mf[mfi]), &comp, &comp);
    }

    mf.mult(sqrt(variance), comp, mf.nComp(), 0);

    mf.OverrideSync(geom.periodicity());
    mf.FillBoundary(geom.periodicity());
}
