#include "compressible_functions.H"
#include "compressible_functions_F.H"

void setBC(MultiFab& prim, MultiFab& cons)
{
    // Loop over boxes
    for ( MFIter mfi(prim); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        set_bc(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                       cons[mfi].dataPtr(),
                       prim[mfi].dataPtr());

    }
}

