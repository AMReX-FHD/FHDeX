#include "compressible_functions.H"
#include "compressible_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "common_namespace.H"

using namespace common;

void RK3step(MultiFab& cu, MultiFab& cup, MultiFab& cup2, MultiFab& cup3, MultiFab& prim, MultiFab& source, MultiFab& eta, MultiFab& zeta, MultiFab& kappa, std::array<MultiFab, AMREX_SPACEDIM>& flux, std::array<MultiFab, AMREX_SPACEDIM>& stochFlux, 
                                                 std::array<MultiFab, AMREX_SPACEDIM>& cornx, std::array<MultiFab, AMREX_SPACEDIM>& corny, std::array<MultiFab, AMREX_SPACEDIM>& cornz, MultiFab& visccorn, MultiFab& rancorn, const amrex::Geometry geom, const amrex::Real* dx, const amrex::Real dt)
{

    amrex::Vector< amrex::Real > stoch_weights;
    amrex::Real swgt1, swgt2;
    swgt1 = 1.0;

    calculateTransportCoeffs(prim, eta, zeta, kappa);

    conservedToPrimitive(prim, cu);
    cu.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());

    calculateTransportCoeffs(prim, eta, zeta, kappa);

    eta.FillBoundary(geom.periodicity());
    zeta.FillBoundary(geom.periodicity());
    kappa.FillBoundary(geom.periodicity());
    
    // Impose membrane BCs
    setBC(prim, cu, eta, zeta, kappa);

    // Set stochastic weights
    swgt2 = ( 2.0*std::sqrt(2.0) + 1.0*std::sqrt(3.0) ) / 5.0;
    stoch_weights = {swgt1, swgt2};

    calculateFlux(cu, prim, eta, zeta, kappa, flux, stochFlux, cornx, corny, cornz, visccorn, rancorn, geom, stoch_weights, dx, dt);

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
    
    // Set stochastic weights
    swgt2 = ( -4.0*std::sqrt(2.0) + 3.0*std::sqrt(3.0) ) / 5.0;
    stoch_weights = {swgt1, swgt2};

    calculateFlux(cup, prim, eta, zeta, kappa, flux, stochFlux, cornx, corny, cornz, visccorn, rancorn, geom, stoch_weights, dx, dt);

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
    
    // Set stochastic weights
    swgt2 = ( 1.0*std::sqrt(2.0) - 2.0*std::sqrt(3.0) ) / 10.0;
    stoch_weights = {swgt1, swgt2};

    calculateFlux(cup2, prim, eta, zeta, kappa, flux, stochFlux, cornx, corny, cornz, visccorn, rancorn, geom, stoch_weights, dx, dt);

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

}

