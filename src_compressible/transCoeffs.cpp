#include "compressible_functions.H"
#include "compressible_functions_F.H"

void calculateTransportCoeffs(const MultiFab& prim, MultiFab& eta, MultiFab& zeta, MultiFab& kappa)
{

    // Loop over boxes
    for ( MFIter mfi(prim); mfi.isValid(); ++mfi)
    {
        
        // amrex::Print() << "Hack: Got here \n";
        
        const Box& bx = mfi.validbox();

        trans_coeffs(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
                       prim[mfi].dataPtr(),  
                       eta[mfi].dataPtr(),  
                       zeta[mfi].dataPtr(),  
                       kappa[mfi].dataPtr());
    }

}
