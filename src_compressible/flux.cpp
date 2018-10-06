#include "compressible_functions.H"
#include "compressible_functions_F.H"

void calculateFlux(const MultiFab& cons, const MultiFab& prim,
                   const MultiFab& eta, const MultiFab& zeta, const MultiFab& kappa,
                   std::array<MultiFab, AMREX_SPACEDIM>& flux,
	           const amrex::Real* dx)
{

    AMREX_D_TERM(flux[0].setVal(0.0);,
                 flux[1].setVal(0.0);,
                 flux[2].setVal(0.0););

    // Loop over boxes
    for ( MFIter mfi(cons); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        diff_flux(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                       cons[mfi].dataPtr(),  
                       prim[mfi].dataPtr(),  
                       eta[mfi].dataPtr(),  
                       zeta[mfi].dataPtr(),  
                       kappa[mfi].dataPtr(),  
        		       flux[0][mfi].dataPtr(),
        		       flux[1][mfi].dataPtr(),
#if (AMREX_SPACEDIM == 3)
        		       flux[2][mfi].dataPtr(),
#endif
    			       ZFILL(dx));

//        hyp_flux(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
//                       cons[mfi].dataPtr(),  
//                       prim[mfi].dataPtr(),  
//                       eta[mfi].dataPtr(),  
//                       zeta[mfi].dataPtr(),  
//                       kappa[mfi].dataPtr(),  
//        		       flux[0][mfi].dataPtr(),
//        		       flux[1][mfi].dataPtr(),
//#if (AMREX_SPACEDIM == 3)
//        		       flux[2][mfi].dataPtr(),
//#endif
//    			       ZFILL(dx), &nvars, &nprimvars);
    }

}
