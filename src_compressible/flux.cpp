#include "compressible_functions.H"
#include "compressible_functions_F.H"
#include "common_functions.H"


void calculateFlux(const MultiFab& cons, const MultiFab& prim,
                   const MultiFab& eta, const MultiFab& zeta, const MultiFab& kappa,
                   std::array<MultiFab, AMREX_SPACEDIM>& flux,
                   std::array<MultiFab, AMREX_SPACEDIM>& stochFlux,const amrex::Geometry geom,
	           const amrex::Real* dx, const amrex::Real dt)
{

    AMREX_D_TERM(flux[0].setVal(0.0);,
                 flux[1].setVal(0.0);,
                 flux[2].setVal(0.0););

    for(int i=0;i<6;i++)
    {
        AMREX_D_TERM(MultiFABFillRandom(stochFlux[0], i, 1, geom);,
                     MultiFABFillRandom(stochFlux[1], i, 1, geom);,
                     MultiFABFillRandom(stochFlux[2], i, 1, geom););
    }

    // Loop over boxes
    for ( MFIter mfi(cons); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        //Must do stoch first

        stoch_flux(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                       cons[mfi].dataPtr(),  
                       prim[mfi].dataPtr(),    
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
                       eta[mfi].dataPtr(),  
                       zeta[mfi].dataPtr(),  
                       kappa[mfi].dataPtr(),
    			       ZFILL(dx), &dt);

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

        hyp_flux(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                       cons[mfi].dataPtr(),  
                       prim[mfi].dataPtr(),    
        		       flux[0][mfi].dataPtr(),
        		       flux[1][mfi].dataPtr(),
#if (AMREX_SPACEDIM == 3)
        		       flux[2][mfi].dataPtr(),
#endif
    			       ZFILL(dx));
   
    }

}
