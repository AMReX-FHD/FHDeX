#include "compressible_functions.H"
#include "compressible_functions_F.H"
#include "common_functions.H"


void calculateFlux(const MultiFab& cons, const MultiFab& prim,
                   const MultiFab& eta, const MultiFab& zeta, const MultiFab& kappa,
                   std::array<MultiFab, AMREX_SPACEDIM>& flux,
                   std::array<MultiFab, AMREX_SPACEDIM>& stochFlux,
                   std::array<MultiFab, AMREX_SPACEDIM>& cornx,
                   std::array<MultiFab, AMREX_SPACEDIM>& corny,
                   std::array<MultiFab, AMREX_SPACEDIM>& cornz,
                   MultiFab& visccorn,
                   MultiFab& rancorn,
                   const amrex::Geometry geom,
	               const amrex::Real* dx, const amrex::Real dt)
{

    AMREX_D_TERM(flux[0].setVal(0);,
                 flux[1].setVal(0);,
                 flux[2].setVal(0););

    for(int d=0;d<AMREX_SPACEDIM;d++)
      {
    	for(int i=1;i<5;i++)
    	  {
    	    MultiFABFillRandom(stochFlux[d], i, 1, geom);
    	  }
      }

//    for(int i=0;i<2;i++)
//    {
//          MultiFABFillRandom(rancorn, i, 1, geom);
//    }

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
                       rancorn[mfi].dataPtr(),
                       eta[mfi].dataPtr(),  
                       zeta[mfi].dataPtr(),  
                       kappa[mfi].dataPtr(),
    			       ZFILL(dx), &dt);


//         diff_flux(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
//                        cons[mfi].dataPtr(),  
//                        prim[mfi].dataPtr(),  
//                        eta[mfi].dataPtr(),  
//                        zeta[mfi].dataPtr(),  
//                        kappa[mfi].dataPtr(),  
//         		       flux[0][mfi].dataPtr(),
//         		       flux[1][mfi].dataPtr(),
// #if (AMREX_SPACEDIM == 3)
//         		       flux[2][mfi].dataPtr(),
// #endif
//     			       ZFILL(dx));

	diff_flux_sym(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
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
                              cornx[0][mfi].dataPtr(),
                              cornx[1][mfi].dataPtr(),
                              cornx[2][mfi].dataPtr(),
                              corny[0][mfi].dataPtr(),
                              corny[1][mfi].dataPtr(),
                              corny[2][mfi].dataPtr(),
                              cornz[0][mfi].dataPtr(),
                              cornz[1][mfi].dataPtr(),
                              cornz[2][mfi].dataPtr(),
                              visccorn[mfi].dataPtr(),
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
