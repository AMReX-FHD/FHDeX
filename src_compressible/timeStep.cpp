#include "compressible_functions.H"
#include "compressible_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

void RK3step(MultiFab& cu, MultiFab& cup, MultiFab& cup2, MultiFab& cup3, MultiFab& prim, MultiFab& source, MultiFab& eta, MultiFab& zeta, MultiFab& kappa, std::array<MultiFab, AMREX_SPACEDIM>& flux, std::array<MultiFab, AMREX_SPACEDIM>& stochFlux,
	                                         const amrex::Geometry geom, const amrex::Real* dx, const amrex::Real dt)
{

    for ( MFIter mfi(cu); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        rk3_stage1(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                   cu[mfi].dataPtr(),  
                   cup[mfi].dataPtr(),  
                   source[mfi].dataPtr(),
      		       flux[0][mfi].dataPtr(),
       		       flux[1][mfi].dataPtr(),
#if (AMREX_SPACEDIM == 3)
       		       flux[2][mfi].dataPtr(),
#endif
      	           ZFILL(dx), &dt);   
    }

    conservedToPrimitive(prim, cup);
    cup.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());

    calculateTransportCoeffs(prim, eta, zeta, kappa);

    eta.FillBoundary(geom.periodicity());
    zeta.FillBoundary(geom.periodicity());
    kappa.FillBoundary(geom.periodicity());

    setBC(prim, cup, eta, zeta, kappa);

    calculateFlux(cup, prim, eta, zeta, kappa, flux, stochFlux, geom, ZFILL(dx), dt);


    for ( MFIter mfi(cu); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        rk3_stage2(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                   cu[mfi].dataPtr(),  
                   cup[mfi].dataPtr(),  
                   cup2[mfi].dataPtr(), 
                   source[mfi].dataPtr(),
      		       flux[0][mfi].dataPtr(),
       		       flux[1][mfi].dataPtr(),
#if (AMREX_SPACEDIM == 3)
       		       flux[2][mfi].dataPtr(),
#endif
      	           ZFILL(dx), &dt);
    }

    conservedToPrimitive(prim, cup2);
    cup2.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());

    calculateTransportCoeffs(prim, eta, zeta, kappa);

    eta.FillBoundary(geom.periodicity());
    zeta.FillBoundary(geom.periodicity());
    kappa.FillBoundary(geom.periodicity());

    setBC(prim, cup2, eta, zeta, kappa);

    calculateFlux(cup2, prim, eta, zeta, kappa, flux, stochFlux, geom, dx, dt);


    for ( MFIter mfi(cu); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        rk3_stage3(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                   cu[mfi].dataPtr(),  
                   cup[mfi].dataPtr(),
                   cup2[mfi].dataPtr(), 
                   source[mfi].dataPtr(),
      		       flux[0][mfi].dataPtr(),
       		       flux[1][mfi].dataPtr(),
#if (AMREX_SPACEDIM == 3)
       		       flux[2][mfi].dataPtr(),
#endif
      	           ZFILL(dx), &dt);
    
    }

    conservedToPrimitive(prim, cu);
    cu.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());

    doSSA(cu, prim, flux, geom, ZFILL(dx), dt);

    calculateTransportCoeffs(prim, eta, zeta, kappa);

    eta.FillBoundary(geom.periodicity());
    zeta.FillBoundary(geom.periodicity());
    kappa.FillBoundary(geom.periodicity());

    setBC(prim, cu, eta, zeta, kappa);

    calculateFlux(cu, prim, eta, zeta, kappa, flux, stochFlux, geom, dx, dt);

}

