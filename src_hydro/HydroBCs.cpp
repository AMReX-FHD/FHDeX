#include "hydro_functions.H"
#include "hydro_functions_F.H"

#include "common_functions.H"


void SetPressureBC(MultiFab & p0, const Geometry & geom) {

    Box dom(geom.Domain());

    for (MFIter mfi(p0); mfi.isValid(); ++mfi) {

        const Box & bx = mfi.validbox();
        set_pressure_bc(BL_TO_FORTRAN_BOX(bx),
                        BL_TO_FORTRAN_BOX(dom),
                        BL_TO_FORTRAN_FAB(p0[mfi]), p0.nGrow());

    }

    // Don't do this because `set_pressure_bc` fills only the ghost cells
    // p0.FillBoundary(geom.periodicity());
}
