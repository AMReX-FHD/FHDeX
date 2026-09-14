#include "compressible_functions.H"

#include "common_functions.H"

#include "chemistry_functions.H"

#include "rng_functions.H"



void RK3step(MultiFab& cu, MultiFab& cup, MultiFab& cup2, MultiFab& /*cup3*/,
             MultiFab& prim, MultiFab& source,
             MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
             MultiFab& chi, MultiFab& D,
             std::array<MultiFab, AMREX_SPACEDIM>& flux,
             std::array<MultiFab, AMREX_SPACEDIM>& stochFlux,
             std::array<MultiFab, AMREX_SPACEDIM>& cornx,
             std::array<MultiFab, AMREX_SPACEDIM>& corny,
             std::array<MultiFab, AMREX_SPACEDIM>& cornz,
             MultiFab& visccorn, MultiFab& rancorn, MultiFab& ranchem,
             const amrex::Geometry& geom, const amrex::Real dt)
{
    BL_PROFILE_VAR("RK3step()",RK3step);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    /////////////////////////////////////////////////////
    // Initialize white noise fields

    // weights for stochastic fluxes; swgt2 changes each stage
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

    // chemistry
    MultiFab ranchem_A;
    MultiFab ranchem_B;
    if (nreaction>0)
    {
        ranchem_A.define(ranchem.boxArray(), ranchem.DistributionMap(), nreaction, 0);
        ranchem_B.define(ranchem.boxArray(), ranchem.DistributionMap(), nreaction, 0);
    }

    // fill random numbers (can skip density component 0)
    for(int d=0;d<AMREX_SPACEDIM;d++) {
        for(int i=1;i<nvars;i++) {
            Real variance;
            if (i>=1 && i <= 3) {
                variance = variance_coef_mom;
            } else if (i == 4) {
                variance = variance_coef_ener;
            } else {
                variance = variance_coef_mass;
            }
            MultiFabFillRandom(stochFlux_A[d], i, variance*variance, geom);
            MultiFabFillRandom(stochFlux_B[d], i, variance*variance, geom);
        }
    }

    MultiFabFillRandom(rancorn_A, 0, variance_coef_mom*variance_coef_mom, geom);
    MultiFabFillRandom(rancorn_B, 0, variance_coef_mom*variance_coef_mom, geom);

    if (nreaction>0) {
        for (int m=0;m<nreaction;m++) {
            MultiFabFillRandom(ranchem_A, m, 1.0, geom);
            MultiFabFillRandom(ranchem_B, m, 1.0, geom);
        }
    }

    /////////////////////////////////////////////////////

    // Compute transport coefs after setting BCs
    calculateTransportCoeffs(prim, eta, zeta, kappa, chi, D);

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
    for(int d=0;d<AMREX_SPACEDIM;d++) {
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

    calculateFlux(cu, prim, eta, zeta, kappa, chi, D, flux, stochFlux, cornx, corny, cornz,
                  visccorn, rancorn, geom, stoch_weights, dt);

    if (nreaction>0)
    {
        MultiFab::LinComb(ranchem,
                    stoch_weights[0], ranchem_A, 0,
                    stoch_weights[1], ranchem_B, 0,
                    0, nreaction, 0);

        compute_compressible_chemistry_source_CLE(dt,dx[0]*dx[1]*dx[2],prim,source,ranchem);
    }

    for ( MFIter mfi(cu,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real> & cu_fab = cu.array(mfi);
        const Array4<Real> & cup_fab = cup.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);
        AMREX_D_TERM(Array4<Real const> const& xflux_fab = flux[0].array(mfi);,
                     Array4<Real const> const& yflux_fab = flux[1].array(mfi);,
                     Array4<Real const> const& zflux_fab = flux[2].array(mfi););

        // for the box,
        // for loop for GPUS
        amrex::ParallelFor(bx, nvars, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept // <- just leave it
        {
                // nth component of cell i,j,k
                // nvars -> nspecies*nspecies
            cup_fab(i,j,k,n) = cu_fab(i,j,k,n) - dt *
                ( AMREX_D_TERM(  (xflux_fab(i+1,j,k,n) - xflux_fab(i,j,k,n)) / dx[0],
                               + (yflux_fab(i,j+1,k,n) - yflux_fab(i,j,k,n)) / dx[1],
                               + (zflux_fab(i,j,k+1,n) - zflux_fab(i,j,k,n)) / dx[2])
                                                                                       )
                + dt*source_fab(i,j,k,n);
        });

        //        amrex::ParallelFor(bx, nvars, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        //        {
        //            cup_fab(i,j,k,n) = cu_fab(i,j,k,n) - dt *
        //                (  (xflux_fab(i+1,j,k,n) - xflux_fab(i,j,k,n)) / dx[0] )
        //                + dt*source_fab(i,j,k,n);
        //        });

        //  what to do about tests

        //              if(cup(i,j,k,1) .lt. 0) then
        //                print *, "Aborting. Negative density at", i,j,k
        //                call exit(0)
        //              endif
        //              if(cup(i,j,k,5) .lt. 0) then
        //                print *, "Aborting. Negative energy at", i,j,k
        //                call exit(0)
        //              endif
        //


        amrex::ParallelFor(bx, 3, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            cup_fab(i,j,k,n+1) += dt * cu_fab(i,j,k,0)*grav[n];
        });

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cup_fab(i,j,k,4) += dt * (  grav[0]*cu_fab(i,j,k,1)
                                      + grav[1]*cu_fab(i,j,k,2)
                                     + grav[2]*cu_fab(i,j,k,3)
                                                                        );
        });
    }

    conservedToPrimitive(prim, cup);

    // Set BC: 1) fill boundary 2) physical
    cup.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());
    setBC(prim, cup);

    // Compute transport coefs after setting BCs
    calculateTransportCoeffs(prim, eta, zeta, kappa, chi, D);

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
    for(int d=0;d<AMREX_SPACEDIM;d++) {
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

    calculateFlux(cup, prim, eta, zeta, kappa, chi, D, flux, stochFlux, cornx, corny, cornz,
                  visccorn, rancorn, geom, stoch_weights, dt);

    if (nreaction>0)
    {
        MultiFab::LinComb(ranchem,
                    stoch_weights[0], ranchem_A, 0,
                    stoch_weights[1], ranchem_B, 0,
                    0, nreaction, 0);

        compute_compressible_chemistry_source_CLE(dt,dx[0]*dx[1]*dx[2],prim,source,ranchem);
    }

    for ( MFIter mfi(cu,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real> & cu_fab = cu.array(mfi);
        const Array4<Real> & cup_fab = cup.array(mfi);
        const Array4<Real> & cup2_fab = cup2.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);
        AMREX_D_TERM(Array4<Real const> const& xflux_fab = flux[0].array(mfi);,
                     Array4<Real const> const& yflux_fab = flux[1].array(mfi);,
                     Array4<Real const> const& zflux_fab = flux[2].array(mfi););

        amrex::ParallelFor(bx, nvars, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            cup2_fab(i,j,k,n) = 0.25*( 3.0* cu_fab(i,j,k,n) + cup_fab(i,j,k,n) - dt *
                                       ( AMREX_D_TERM(  (xflux_fab(i+1,j,k,n) - xflux_fab(i,j,k,n)) / dx[0],
                                                      + (yflux_fab(i,j+1,k,n) - yflux_fab(i,j,k,n)) / dx[1],
                                                      + (zflux_fab(i,j,k+1,n) - zflux_fab(i,j,k,n)) / dx[2])
                                                                                                                )
                                       +dt*source_fab(i,j,k,n)  );
        });

        //        amrex::ParallelFor(bx, nvars, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        //        {
        //            cup2_fab(i,j,k,n) = 0.25*( 3.0* cu_fab(i,j,k,n) + cup_fab(i,j,k,n) - dt *
        //                                       (  (xflux_fab(i+1,j,k,n) - xflux_fab(i,j,k,n)) / dx[0])
        //                                       +dt*source_fab(i,j,k,n)  );
        //        });

        //  what to do about tests

        //              if(cup(i,j,k,1) .lt. 0) then
        //                print *, "Aborting. Negative density at", i,j,k
        //                call exit(0)
        //              endif
        //              if(cup(i,j,k,5) .lt. 0) then
        //                print *, "Aborting. Negative energy at", i,j,k
        //                call exit(0)
        //              endif
        //


        amrex::ParallelFor(bx, 3, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            cup2_fab(i,j,k,n+1) += 0.25* dt * cup_fab(i,j,k,0)*grav[n];
        });

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cup2_fab(i,j,k,4) += 0.25 * dt * (  grav[0]*cup_fab(i,j,k,1)
                                              + grav[1]*cup_fab(i,j,k,2)
                                              + grav[2]*cup_fab(i,j,k,3));
        });
    }

    conservedToPrimitive(prim, cup2);

    // Set BC: 1) fill boundary 2) physical
    cup2.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());
    setBC(prim, cup2);

    // Compute transport coefs after setting BCs
    calculateTransportCoeffs(prim, eta, zeta, kappa, chi, D);

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
    for(int d=0;d<AMREX_SPACEDIM;d++) {
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

    calculateFlux(cup2, prim, eta, zeta, kappa, chi, D, flux, stochFlux, cornx, corny, cornz,
                  visccorn, rancorn, geom, stoch_weights, dt);

    if (nreaction>0)
    {
        MultiFab::LinComb(ranchem,
                    stoch_weights[0], ranchem_A, 0,
                    stoch_weights[1], ranchem_B, 0,
                    0, nreaction, 0);

        compute_compressible_chemistry_source_CLE(dt,dx[0]*dx[1]*dx[2],prim,source,ranchem);
    }

    for ( MFIter mfi(cu,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real> & cu_fab = cu.array(mfi);
        const Array4<Real> & cup2_fab = cup2.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);
        AMREX_D_TERM(Array4<Real const> const& xflux_fab = flux[0].array(mfi);,
                     Array4<Real const> const& yflux_fab = flux[1].array(mfi);,
                     Array4<Real const> const& zflux_fab = flux[2].array(mfi););

        amrex::ParallelFor(bx, nvars, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            cu_fab(i,j,k,n) = (2./3.) *( 0.5* cu_fab(i,j,k,n) + cup2_fab(i,j,k,n) - dt *
                                    (   AMREX_D_TERM(  (xflux_fab(i+1,j,k,n) - xflux_fab(i,j,k,n)) / dx[0],
                                                     + (yflux_fab(i,j+1,k,n) - yflux_fab(i,j,k,n)) / dx[1],
                                                     + (zflux_fab(i,j,k+1,n) - zflux_fab(i,j,k,n)) / dx[2])
                                                                                                            )
                                    + dt*source_fab(i,j,k,n) );

        });

        //        amrex::ParallelFor(bx, nvars, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        //        {
        //            cu_fab(i,j,k,n) = (2./3.) *( 0.5* cu_fab(i,j,k,n) + cup2_fab(i,j,k,n) - dt *
        //                                    (     (xflux_fab(i+1,j,k,n) - xflux_fab(i,j,k,n)) / dx[0])
        //
        //                                    + dt*source_fab(i,j,k,n) );
        //
        //        });

        //  what to do about tests

        //              if(cup(i,j,k,1) .lt. 0) then
        //                print *, "Aborting. Negative density at", i,j,k
        //                call exit(0)
        //              endif
        //              if(cup(i,j,k,5) .lt. 0) then
        //                print *, "Aborting. Negative energy at", i,j,k
        //                call exit(0)
        //              endif
        //


        amrex::ParallelFor(bx, 3, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            cu_fab(i,j,k,n+1) += 2./3.* dt * cup2_fab(i,j,k,0)*grav[n];
        });

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cu_fab(i,j,k,4) += 2./3. * dt * (  grav[0]*cup2_fab(i,j,k,1)
                                             + grav[1]*cup2_fab(i,j,k,2)
                                             + grav[2]*cup2_fab(i,j,k,3) );
        });

    }

    conservedToPrimitive(prim, cu);

    // Set BC: 1) fill boundary 2) physical
    cu.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());

    //doMembrane(cu,prim,flux,geom,dxp,dt);

    setBC(prim, cu);

}