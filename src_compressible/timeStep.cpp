#include "compressible_functions.H"
#include "compressible_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "rng_functions.H"
#include "rng_functions_F.H"

#include "common_namespace.H"

using namespace common;

void RK3step(MultiFab& cu, MultiFab& cup, MultiFab& cup2, MultiFab& cup3, MultiFab& prim, MultiFab& source, MultiFab& eta, MultiFab& zeta, MultiFab& kappa, MultiFab& chi, MultiFab& D, std::array<MultiFab, AMREX_SPACEDIM>& flux, std::array<MultiFab, AMREX_SPACEDIM>& stochFlux, 
                                                 std::array<MultiFab, AMREX_SPACEDIM>& cornx, std::array<MultiFab, AMREX_SPACEDIM>& corny, std::array<MultiFab, AMREX_SPACEDIM>& cornz, MultiFab& visccorn, MultiFab& rancorn, const amrex::Geometry geom, const amrex::Real* dx, const amrex::Real dt)
{
    /////////////////////////////////////////////////////
    // Initialize white noise fields

    amrex::Vector< amrex::Real > stoch_weights;
    amrex::Real swgt1, swgt2;
    swgt1 = 1.0;

    // Temp. stoch. fluxes

    // field "A"
    std::array< MultiFab, AMREX_SPACEDIM > stochFlux_A;
    AMREX_D_TERM(stochFlux_A[0].define(stochFlux[0].boxArray(), stochFlux[0].DistributionMap(), nvars, 0);,
                 stochFlux_A[1].define(stochFlux[1].boxArray(), stochFlux[1].DistributionMap(), nvars, 0);,
                 stochFlux_A[2].define(stochFlux[2].boxArray(), stochFlux[2].DistributionMap(), nvars, 0););

    MultiFab rancorn_A;
    rancorn_A.define(rancorn.boxArray(), rancorn.DistributionMap(), 1, 0);
    
    // field "B"
    std::array< MultiFab, AMREX_SPACEDIM > stochFlux_B;
    AMREX_D_TERM(stochFlux_B[0].define(stochFlux[0].boxArray(), stochFlux[0].DistributionMap(), nvars, 0);,
                 stochFlux_B[1].define(stochFlux[1].boxArray(), stochFlux[1].DistributionMap(), nvars, 0);,
                 stochFlux_B[2].define(stochFlux[2].boxArray(), stochFlux[2].DistributionMap(), nvars, 0););

    MultiFab rancorn_B;
    rancorn_B.define(rancorn.boxArray(), rancorn.DistributionMap(), 1, 0);

    AMREX_D_TERM(stochFlux_A[0].setVal(0.0);,
                 stochFlux_A[1].setVal(0.0);,
                 stochFlux_A[2].setVal(0.0););
    rancorn_A.setVal(0.0);

    AMREX_D_TERM(stochFlux_B[0].setVal(0.0);,
                 stochFlux_B[1].setVal(0.0);,
                 stochFlux_B[2].setVal(0.0););
    rancorn_B.setVal(0.0);

    // Fill random (only momentum and energy)
    for(int d=0;d<AMREX_SPACEDIM;d++)
      {
    	for(int i=1;i<nvars;i++)
    	  {
    	    MultiFABFillRandom(stochFlux_A[d], i, 1.0, geom);
	    MultiFABFillRandom(stochFlux_B[d], i, 1.0, geom);
	  }
      }

    MultiFABFillRandom(rancorn_A, 0, 1.0, geom);
    MultiFABFillRandom(rancorn_B, 0, 1.0, geom);

    /////////////////////////////////////////////////////

    conservedToPrimitive(prim, cu);
    cu.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());

    calculateTransportCoeffs(prim, eta, zeta, kappa, chi, D);

    eta.FillBoundary(geom.periodicity());
    zeta.FillBoundary(geom.periodicity());
    kappa.FillBoundary(geom.periodicity());
    
    // Impose membrane BCs
    setBC(prim, cu, eta, zeta, kappa);

    ///////////////////////////////////////////////////////////
    // Perform weighting of white noise fields

    // Set stochastic weights
    swgt2 = ( 2.0*std::sqrt(2.0) + 1.0*std::sqrt(3.0) ) / 5.0;
    stoch_weights = {swgt1, swgt2};

    AMREX_D_TERM(stochFlux[0].setVal(0.0);,
                 stochFlux[1].setVal(0.0);,
                 stochFlux[2].setVal(0.0););
    rancorn.setVal(0.0);

    // apply weights (only momentum and energy)
    for(int d=0;d<AMREX_SPACEDIM;d++)
      {

	MultiFab::LinComb(stochFlux[d], 
			  stoch_weights[0], stochFlux_A[d], 1, 
			  stoch_weights[1], stochFlux_B[d], 1,
			  1, nvars-1, 0);

      }

    MultiFab::LinComb(rancorn, 
		      stoch_weights[0], rancorn_A, 0, 
		      stoch_weights[1], rancorn_B, 0,
		      0, 1, 0);

    ///////////////////////////////////////////////////////////

    calculateFlux(cu, prim, eta, zeta, kappa, chi, D, flux, stochFlux, cornx, corny, cornz, visccorn, rancorn, geom, stoch_weights, dx, dt);

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

    calculateTransportCoeffs(prim, eta, zeta, kappa, chi, D);

    eta.FillBoundary(geom.periodicity());
    zeta.FillBoundary(geom.periodicity());
    kappa.FillBoundary(geom.periodicity());

    setBC(prim, cup, eta, zeta, kappa);

    ///////////////////////////////////////////////////////////
    // Perform weighting of white noise fields

    // Set stochastic weights
    swgt2 = ( -4.0*std::sqrt(2.0) + 3.0*std::sqrt(3.0) ) / 5.0;
    stoch_weights = {swgt1, swgt2};

    AMREX_D_TERM(stochFlux[0].setVal(0.0);,
                 stochFlux[1].setVal(0.0);,
                 stochFlux[2].setVal(0.0););
    rancorn.setVal(0.0);

    // apply weights (only momentum and energy)
    for(int d=0;d<AMREX_SPACEDIM;d++)
      {

	MultiFab::LinComb(stochFlux[d], 
			  stoch_weights[0], stochFlux_A[d], 1, 
			  stoch_weights[1], stochFlux_B[d], 1,
			  1, nvars-1, 0);

      }

    MultiFab::LinComb(rancorn, 
		      stoch_weights[0], rancorn_A, 0, 
		      stoch_weights[1], rancorn_B, 0,
		      0, 1, 0);

    ///////////////////////////////////////////////////////////

    calculateFlux(cup, prim, eta, zeta, kappa, chi, D, flux, stochFlux, cornx, corny, cornz, visccorn, rancorn, geom, stoch_weights, dx, dt);

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

    calculateTransportCoeffs(prim, eta, zeta, kappa, chi, D);

    eta.FillBoundary(geom.periodicity());
    zeta.FillBoundary(geom.periodicity());
    kappa.FillBoundary(geom.periodicity());

    setBC(prim, cup2, eta, zeta, kappa);

    ///////////////////////////////////////////////////////////
    // Perform weighting of white noise fields

    // Set stochastic weights
    swgt2 = ( 1.0*std::sqrt(2.0) - 2.0*std::sqrt(3.0) ) / 10.0;
    stoch_weights = {swgt1, swgt2};
    
    AMREX_D_TERM(stochFlux[0].setVal(0.0);,
                 stochFlux[1].setVal(0.0);,
                 stochFlux[2].setVal(0.0););
    rancorn.setVal(0.0);

    // apply weights (only momentum and energy)
    for(int d=0;d<AMREX_SPACEDIM;d++)
      {

	MultiFab::LinComb(stochFlux[d], 
			  stoch_weights[0], stochFlux_A[d], 1, 
			  stoch_weights[1], stochFlux_B[d], 1,
			  1, nvars-1, 0);

      }

    MultiFab::LinComb(rancorn, 
		      stoch_weights[0], rancorn_A, 0, 
		      stoch_weights[1], rancorn_B, 0,
		      0, 1, 0);

    ///////////////////////////////////////////////////////////

    calculateFlux(cup2, prim, eta, zeta, kappa, chi, D, flux, stochFlux, cornx, corny, cornz, visccorn, rancorn, geom, stoch_weights, dx, dt);

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

