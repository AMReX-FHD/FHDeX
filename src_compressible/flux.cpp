#include "compressible_functions.H"
#include "compressible_functions_F.H"
#include "common_functions.H"

#include "common_namespace.H"

using namespace common;

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
		   const amrex::Vector< amrex::Real >& stoch_weights,
		   const amrex::Real* dx, const amrex::Real dt)
{

    AMREX_D_TERM(flux[0].setVal(0);,
                 flux[1].setVal(0);,
                 flux[2].setVal(0););

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

	if (abs(visc_type) > 1) {

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

	} else {

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

	}

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
