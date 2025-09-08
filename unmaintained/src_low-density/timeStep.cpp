#include "compressible_functions.H"

#include "common_functions.H"

#include "rng_functions.H"



void eulerStep(MultiFab& cu, std::array<MultiFab,AMREX_SPACEDIM>& flux,
               std::array<MultiFab, AMREX_SPACEDIM>& stochFlux,
               const amrex::Geometry geom, const amrex::Real* dx,
               const amrex::Real dt)
{
    /////////////////////////////////////////////////////
    //fill ghost cells

    cu.FillBoundary(geom.periodicity());

    // Impose membrane BCs
    // setBC(prim, cu, eta, zeta, kappa);

    //zero out stoch fluxes
    AMREX_D_TERM(stochFlux[0].setVal(0.0);,
                 stochFlux[1].setVal(0.0);,
                 stochFlux[2].setVal(0.0););


    //get the flux
    calculateFlux(cu, flux, stochFlux, dx, dt);

    //euler step
    for ( MFIter mfi(cu); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        euler_step(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                   cu[mfi].dataPtr(),
                   flux[0][mfi].dataPtr(),
                   flux[1][mfi].dataPtr(),
#if (AMREX_SPACEDIM == 3)
                   flux[2][mfi].dataPtr(),
#endif
                   ZFILL(dx), &dt);
    }


}

void RK3step(MultiFab& cu, MultiFab& cup, MultiFab& cup2, MultiFab& cup3,
             std::array<MultiFab,AMREX_SPACEDIM>& flux, std::array<MultiFab,
             AMREX_SPACEDIM>& stochFlux, const amrex::Geometry geom,
             const amrex::Real* dx, const amrex::Real dt)
{
    /////////////////////////////////////////////////////
    // Initialize white noise fields

    amrex::Vector< amrex::Real > stoch_weights;
    amrex::Real swgt1, swgt2;
    swgt1 = 1.0;

    cu.FillBoundary(geom.periodicity());

    // Impose membrane BCs
    // setBC(prim, cu, eta, zeta, kappa);

    ///////////////////////////////////////////////////////////
    // Perform weighting of white noise fields

    // Set stochastic weights
    swgt2 = ( 2.0*std::sqrt(2.0) + 1.0*std::sqrt(3.0) ) / 5.0;
    stoch_weights = {swgt1, swgt2};

    AMREX_D_TERM(stochFlux[0].setVal(0.0);,
                 stochFlux[1].setVal(0.0);,
                 stochFlux[2].setVal(0.0););


    calculateFlux(cu, flux, stochFlux, dx, dt);

    for ( MFIter mfi(cu); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        rk3_stage1(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                   cu[mfi].dataPtr(),
                   cup[mfi].dataPtr(),
                   flux[0][mfi].dataPtr(),
                   flux[1][mfi].dataPtr(),
#if (AMREX_SPACEDIM == 3)
                   flux[2][mfi].dataPtr(),
#endif
                   ZFILL(dx), &dt);
    }

    cup.FillBoundary(geom.periodicity());

    calculateFlux(cup, flux, stochFlux, dx, dt);

    for ( MFIter mfi(cu); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        rk3_stage2(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                   cu[mfi].dataPtr(),
                   cup[mfi].dataPtr(),
                   cup2[mfi].dataPtr(),
                   flux[0][mfi].dataPtr(),
                   flux[1][mfi].dataPtr(),
#if (AMREX_SPACEDIM == 3)
                   flux[2][mfi].dataPtr(),
#endif
                   ZFILL(dx), &dt);
    }

    cup2.FillBoundary(geom.periodicity());

    calculateFlux(cup2, flux, stochFlux, dx, dt);

    for ( MFIter mfi(cu); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        rk3_stage3(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                   cu[mfi].dataPtr(),
                   cup[mfi].dataPtr(),
                   cup2[mfi].dataPtr(),
                   flux[0][mfi].dataPtr(),
                   flux[1][mfi].dataPtr(),
#if (AMREX_SPACEDIM == 3)
                   flux[2][mfi].dataPtr(),
#endif
                   ZFILL(dx), &dt);

    }

}