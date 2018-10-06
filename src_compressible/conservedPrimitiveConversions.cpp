#include "compressible_functions.H"
#include "compressible_functions_F.H"

void conservedToPrimative(MultiFab& prim, const MultiFab& cons)
{

    // Loop over boxes
    for ( MFIter mfi(prim); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

//        trans_coeffs(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
//                       prim[mfi].dataPtr(),  
//                       eta[mfi].dataPtr(),  
//                       zeta[mfi].dataPtr(),  
//                       kappa[mfi].dataPtr());
    }

}

void primativeToConserved(const MultiFab& prim, MultiFab& cons)
{

    // Loop over boxes
    for ( MFIter mfi(prim); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

//        trans_coeffs(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
//                       prim[mfi].dataPtr(),  
//                       eta[mfi].dataPtr(),  
//                       zeta[mfi].dataPtr(),  
//                       kappa[mfi].dataPtr());
    }

}
