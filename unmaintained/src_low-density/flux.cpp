#include "compressible_functions.H"
#include "common_functions.H"



void calculateFlux(const MultiFab& cons,
                   std::array<MultiFab, AMREX_SPACEDIM>& flux,
                   std::array<MultiFab, AMREX_SPACEDIM>& stochFlux,
    const amrex::Real* dx, const amrex::Real dt)
{

    AMREX_D_TERM(flux[0].setVal(0);,
                 flux[1].setVal(0);,
                 flux[2].setVal(0););

    // Loop over boxes
    for ( MFIter mfi(cons); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        //NOTE: Must do stoch. fluxes first,
        //      because fluxes at boundaries are weighted according to BCs

        stoch_flux(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                   cons[mfi].dataPtr(),
                   flux[0][mfi].dataPtr(),
                   flux[1][mfi].dataPtr(),
#if (AMREX_SPACEDIM == 3)
                   flux[2][mfi].dataPtr(),
#endif
                   stochFlux[0][mfi].dataPtr(),
                   stochFlux[1][mfi].dataPtr(),
#if (AMREX_SPACEDIM == 3)
                   stochFlux[2][mfi].dataPtr(),
#endif
                   ZFILL(dx), &dt);


//    diff_flux(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
//        cons[mfi].dataPtr(),
//        prim[mfi].dataPtr(),
//        eta[mfi].dataPtr(),
//        zeta[mfi].dataPtr(),
//        kappa[mfi].dataPtr(),
//        chi[mfi].dataPtr(),
//        D[mfi].dataPtr(),
//        flux[0][mfi].dataPtr(),
//        flux[1][mfi].dataPtr(),
//#if (AMREX_SPACEDIM == 3)
//        flux[2][mfi].dataPtr(),
//#endif
//        cornx[0][mfi].dataPtr(),
//        cornx[1][mfi].dataPtr(),
//        cornx[2][mfi].dataPtr(),
//        corny[0][mfi].dataPtr(),
//        corny[1][mfi].dataPtr(),
//        corny[2][mfi].dataPtr(),
//        cornz[0][mfi].dataPtr(),
//        cornz[1][mfi].dataPtr(),
//        cornz[2][mfi].dataPtr(),
//        visccorn[mfi].dataPtr(),
//        ZFILL(dx));

    }

}