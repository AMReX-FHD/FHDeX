#include "compressible_functions.H"

void setBC(MultiFab& prim, MultiFab& cons, MultiFab& eta, MultiFab& zeta, MultiFab& kappa)
{

    // Loop over boxes
    for ( MFIter mfi(prim); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        set_bc(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                       cons[mfi].dataPtr(),
                       prim[mfi].dataPtr(),
                       eta[mfi].dataPtr(),
                       zeta[mfi].dataPtr(),
                       kappa[mfi].dataPtr());

    }

}

