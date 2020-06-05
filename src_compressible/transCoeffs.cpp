#include "compressible_functions.H"

void calculateTransportCoeffs(const MultiFab& prim, 
			      MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
			      MultiFab& chi, MultiFab& D)
{
    BL_PROFILE_VAR("calculateTransportCoeffs()",calculateTransportCoeffs);
    
    // Loop over boxes
    for ( MFIter mfi(prim); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.validbox();

        makecoef(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
		 prim[mfi].dataPtr(),  
		 eta[mfi].dataPtr(),  
		 zeta[mfi].dataPtr(),  
		 kappa[mfi].dataPtr(),
		 chi[mfi].dataPtr(),
		 D[mfi].dataPtr());

        // trans_coeffs(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
	// 	     prim[mfi].dataPtr(),  
	// 	     eta[mfi].dataPtr(),  
	// 	     zeta[mfi].dataPtr(),  
	// 	     kappa[mfi].dataPtr());
    }
}
