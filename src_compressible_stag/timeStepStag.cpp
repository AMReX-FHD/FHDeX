#include "compressible_functions.H"
#include "compressible_functions_stag.H"

#include "common_functions.H"

#include "rng_functions.H"

void RK3stepStag(MultiFab& cu, 
                 std::array< MultiFab, AMREX_SPACEDIM >& cumom,
                 MultiFab& prim, std::array< MultiFab, AMREX_SPACEDIM >& facevel,
                 MultiFab& source,
                 MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
                 MultiFab& chi, MultiFab& D,
                 std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                 std::array< MultiFab, 2 >& edgeflux_x,
                 std::array< MultiFab, 2 >& edgeflux_y,
                 std::array< MultiFab, 2 >& edgeflux_z,
                 std::array< MultiFab, 3 >& cenflux,
                 std::array<MultiFab, AMREX_SPACEDIM>& stochFlux,
                 MultiFab& rancorn,
                 const amrex::Geometry geom, const amrex::Real* dxp, const amrex::Real dt)
{
    BL_PROFILE_VAR("RK3stepStag()",RK3stepStag);

    MultiFab cup (cu.boxArray(),cu.DistributionMap(),nvars,ngc);
    MultiFab cup2(cu.boxArray(),cu.DistributionMap(),nvars,ngc);
    cup.setVal(0.0,0,nvars,ngc);
    cup2.setVal(0.0,0,nvars,ngc);
    cup.setVal(rho0,0,1,ngc);
    cup2.setVal(rho0,0,1,ngc);

    std::array< MultiFab, AMREX_SPACEDIM > cupmom;
    std::array< MultiFab, AMREX_SPACEDIM > cup2mom;
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        cupmom[d].define(convert(cu.boxArray(),nodal_flag_dir[d]), cu.DistributionMap(), 1, 1);
        cupmom[d].setVal(0.);
        cup2mom[d].define(convert(cu.boxArray(),nodal_flag_dir[d]), cu.DistributionMap(), 1, 1);
        cup2mom[d].setVal(0.);
    }

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{AMREX_D_DECL(grav[0], grav[1], grav[2])};

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

    // fill random numbers (can skip density component 0)
    for(int d=0;d<AMREX_SPACEDIM;d++) {
    	for(int i=1;i<nvars;i++) {
    	    MultiFabFillRandom(stochFlux_A[d], i, 1.0, geom);
	    MultiFabFillRandom(stochFlux_B[d], i, 1.0, geom);
        }
    }

    MultiFabFillRandom(rancorn_A, 0, 1.0, geom);
    MultiFabFillRandom(rancorn_B, 0, 1.0, geom);

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

    calculateFluxStag(cu, cumom, prim, facevel, eta, zeta, kappa, chi, D, faceflux, edgeflux_x, edgeflux_y, edgeflux_z, cenflux, stochFlux, 
                      rancorn, geom, stoch_weights, dxp, dt);

    for ( MFIter mfi(cu,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.tilebox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<Real> & cu_fab = cu.array(mfi);
        const Array4<Real> & cup_fab = cup.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        AMREX_D_TERM(const Array4<Real>& momx = cumom[0].array(mfi);,
                     const Array4<Real>& momy = cumom[1].array(mfi);,
                     const Array4<Real>& momz = cumom[2].array(mfi););

        AMREX_D_TERM(const Array4<Real>& mompx = cupmom[0].array(mfi);,
                     const Array4<Real>& mompy = cupmom[1].array(mfi);,
                     const Array4<Real>& mompz = cupmom[2].array(mfi););

        AMREX_D_TERM(Array4<Real const> const& xflux_fab = faceflux[0].array(mfi);,
                     Array4<Real const> const& yflux_fab = faceflux[1].array(mfi);,
                     Array4<Real const> const& zflux_fab = faceflux[2].array(mfi););

        Array4<Real const> const& edgex_v = edgeflux_x[0].array(mfi);
        Array4<Real const> const& edgex_w = edgeflux_x[1].array(mfi);
        Array4<Real const> const& edgey_u = edgeflux_y[0].array(mfi);
        Array4<Real const> const& edgey_w = edgeflux_y[1].array(mfi);
        Array4<Real const> const& edgez_u = edgeflux_z[0].array(mfi);
        Array4<Real const> const& edgez_v = edgeflux_z[1].array(mfi);

        Array4<Real const> const& cenx_u = cenflux[0].array(mfi);
        Array4<Real const> const& ceny_v = cenflux[1].array(mfi);
        Array4<Real const> const& cenz_w = cenflux[2].array(mfi);

        amrex::ParallelFor(bx, nvars, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            cup_fab(i,j,k,n) = cu_fab(i,j,k,n) - dt *
                ( AMREX_D_TERM(  (xflux_fab(i+1,j,k,n) - xflux_fab(i,j,k,n)) / dx[0],
                               + (yflux_fab(i,j+1,k,n) - yflux_fab(i,j,k,n)) / dx[1],
                               + (zflux_fab(i,j,k+1,n) - zflux_fab(i,j,k,n)) / dx[2])
                                                                                       )
                + dt*source_fab(i,j,k,n);
        }); // [1:3 indices are not valuable -- momentum flux]

        // momentum flux
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            mompx(i,j,k) = momx(i,j,k) 
                    -dt*(cenx_u(i,j,k) - cenx_u(i-1,j,k))/dx[0]
                    -dt*(edgey_u(i,j+1,k) - edgey_u(i,j,k))/dx[1]
                    -dt*(edgez_u(i,j,k+1) - edgez_u(i,j,k))/dx[2]
                    +0.5*dt*grav_gpu[0]*(cu_fab(i-1,j,k,0)+cu_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            mompy(i,j,k) = momy(i,j,k)
                    -dt*(edgex_v(i+1,j,k) - edgex_v(i,j,k))/dx[0]
                    -dt*(ceny_v(i,j,k) - cenx_u(i,j-1,k))/dx[1]
                    -dt*(edgez_v(i,j,k+1) - edgez_v(i,j,k))/dx[2]
                    +0.5*dt*grav_gpu[1]*(cu_fab(i,j-1,k,0)+cu_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            mompz(i,j,k) = momz(i,j,k)
                    -dt*(edgex_w(i+1,j,k) - edgex_w(i,j,k))/dx[0]
                    -dt*(edgey_w(i,j+1,k) - edgey_w(i,j,k))/dx[1]
                    -dt*(cenz_w(i,j,k) - cenx_u(i,j,k-1))/dx[2]
                    +0.5*dt*grav_gpu[2]*(cu_fab(i,j,k-1,0)+cu_fab(i,j,k,0));
        });
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cup_fab(i,j,k,4) += 0.5 * dt * (  grav_gpu[0]*(momx(i+1,j,k)+momx(i,j,k))
                                            + grav_gpu[1]*(momy(i,j+1,k)+momy(i,j,k))
                                            + grav_gpu[2]*(momz(i,j,k+1)+momz(i,j,k)) );
        });
    }

    conservedToPrimitiveStag(prim, facevel, cup, cupmom);
    
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

    calculateFluxStag(cup, cupmom, prim, facevel, eta, zeta, kappa, chi, D, faceflux, edgeflux_x, edgeflux_y, edgeflux_z, cenflux, stochFlux, 
                      rancorn, geom, stoch_weights, dxp, dt);

    for ( MFIter mfi(cu,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.tilebox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<Real> & cu_fab = cu.array(mfi);
        const Array4<Real> & cup_fab = cup.array(mfi);
        const Array4<Real> & cup2_fab = cup2.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        AMREX_D_TERM(const Array4<Real>& momx = cumom[0].array(mfi);,
                     const Array4<Real>& momy = cumom[1].array(mfi);,
                     const Array4<Real>& momz = cumom[2].array(mfi););

        AMREX_D_TERM(const Array4<Real>& mompx = cupmom[0].array(mfi);,
                     const Array4<Real>& mompy = cupmom[1].array(mfi);,
                     const Array4<Real>& mompz = cupmom[2].array(mfi););

        AMREX_D_TERM(const Array4<Real>& momp2x = cup2mom[0].array(mfi);,
                     const Array4<Real>& momp2y = cup2mom[1].array(mfi);,
                     const Array4<Real>& momp2z = cup2mom[2].array(mfi););

        AMREX_D_TERM(Array4<Real const> const& xflux_fab = faceflux[0].array(mfi);,
                     Array4<Real const> const& yflux_fab = faceflux[1].array(mfi);,
                     Array4<Real const> const& zflux_fab = faceflux[2].array(mfi););

        Array4<Real const> const& edgex_v = edgeflux_x[0].array(mfi);
        Array4<Real const> const& edgex_w = edgeflux_x[1].array(mfi);
        Array4<Real const> const& edgey_u = edgeflux_y[0].array(mfi);
        Array4<Real const> const& edgey_w = edgeflux_y[1].array(mfi);
        Array4<Real const> const& edgez_u = edgeflux_z[0].array(mfi);
        Array4<Real const> const& edgez_v = edgeflux_z[1].array(mfi);

        Array4<Real const> const& cenx_u = cenflux[0].array(mfi);
        Array4<Real const> const& ceny_v = cenflux[1].array(mfi);
        Array4<Real const> const& cenz_w = cenflux[2].array(mfi);

        amrex::ParallelFor(bx, nvars, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            cup2_fab(i,j,k,n) = 0.25*( 3.0* cu_fab(i,j,k,n) + cup_fab(i,j,k,n) - dt *
                                       ( AMREX_D_TERM(  (xflux_fab(i+1,j,k,n) - xflux_fab(i,j,k,n)) / dx[0],
                                                      + (yflux_fab(i,j+1,k,n) - yflux_fab(i,j,k,n)) / dx[1],
                                                      + (zflux_fab(i,j,k+1,n) - zflux_fab(i,j,k,n)) / dx[2])
                                                                                                                )
                                       +dt*source_fab(i,j,k,n)  );
        }); // [1:3 indices are not valuable -- momentum flux]

        // momentum flux
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momp2x(i,j,k) = 0.25*3.0*momx(i,j,k) + 0.25*mompx(i,j,k)
                    -0.25*dt*(cenx_u(i,j,k) - cenx_u(i-1,j,k))/dx[0]
                    -0.25*dt*(edgey_u(i,j+1,k) - edgey_u(i,j,k))/dx[1]
                    -0.25*dt*(edgez_u(i,j,k+1) - edgez_u(i,j,k))/dx[2]
                    +0.5*0.25*dt*grav_gpu[0]*(cup_fab(i-1,j,k,0)+cup_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momp2y(i,j,k) = 0.25*3.0*momy(i,j,k) + 0.25*mompy(i,j,k)
                    -0.25*dt*(edgex_v(i+1,j,k) - edgex_v(i,j,k))/dx[0]
                    -0.25*dt*(ceny_v(i,j,k) - cenx_u(i,j-1,k))/dx[1]
                    -0.25*dt*(edgez_v(i,j,k+1) - edgez_v(i,j,k))/dx[2]
                    +0.5*0.25*dt*grav_gpu[1]*(cup_fab(i,j-1,k,0)+cup_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momp2z(i,j,k) = 0.25*3.0*momz(i,j,k) + 0.25*mompz(i,j,k)
                    -0.25*dt*(edgex_w(i+1,j,k) - edgex_w(i,j,k))/dx[0]
                    -0.25*dt*(edgey_w(i,j+1,k) - edgey_w(i,j,k))/dx[1]
                    -0.25*dt*(cenz_w(i,j,k) - cenx_u(i,j,k-1))/dx[2]
                    +0.5*0.25*dt*grav_gpu[2]*(cup_fab(i,j,k-1,0)+cup_fab(i,j,k,0));
        });
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cup2_fab(i,j,k,4) += 0.5 * 0.25 * dt * (  grav_gpu[0]*(momx(i+1,j,k)+momx(i,j,k))
                                                    + grav_gpu[1]*(momy(i,j+1,k)+momy(i,j,k))
                                                    + grav_gpu[2]*(momz(i,j,k+1)+momz(i,j,k)) );
        });
    }
        
    conservedToPrimitiveStag(prim, facevel, cup2, cup2mom);

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

    calculateFluxStag(cup2, cup2mom, prim, facevel, eta, zeta, kappa, chi, D, faceflux, edgeflux_x, edgeflux_y, edgeflux_z, cenflux, stochFlux,
                      rancorn, geom, stoch_weights, dxp, dt);

    for ( MFIter mfi(cu,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.tilebox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<Real> & cu_fab = cu.array(mfi);
        const Array4<Real> & cup2_fab = cup2.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        AMREX_D_TERM(const Array4<Real>& momx = cumom[0].array(mfi);,
                     const Array4<Real>& momy = cumom[1].array(mfi);,
                     const Array4<Real>& momz = cumom[2].array(mfi););

        AMREX_D_TERM(const Array4<Real>& momp2x = cup2mom[0].array(mfi);,
                     const Array4<Real>& momp2y = cup2mom[1].array(mfi);,
                     const Array4<Real>& momp2z = cup2mom[2].array(mfi););

        AMREX_D_TERM(Array4<Real const> const& xflux_fab = faceflux[0].array(mfi);,
                     Array4<Real const> const& yflux_fab = faceflux[1].array(mfi);,
                     Array4<Real const> const& zflux_fab = faceflux[2].array(mfi););

        Array4<Real const> const& edgex_v = edgeflux_x[0].array(mfi);
        Array4<Real const> const& edgex_w = edgeflux_x[1].array(mfi);
        Array4<Real const> const& edgey_u = edgeflux_y[0].array(mfi);
        Array4<Real const> const& edgey_w = edgeflux_y[1].array(mfi);
        Array4<Real const> const& edgez_u = edgeflux_z[0].array(mfi);
        Array4<Real const> const& edgez_v = edgeflux_z[1].array(mfi);

        Array4<Real const> const& cenx_u = cenflux[0].array(mfi);
        Array4<Real const> const& ceny_v = cenflux[1].array(mfi);
        Array4<Real const> const& cenz_w = cenflux[2].array(mfi);

        amrex::ParallelFor(bx, nvars, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            cu_fab(i,j,k,n) = (2./3.) *( 0.5* cu_fab(i,j,k,n) + cup2_fab(i,j,k,n) - dt *
                                    (   AMREX_D_TERM(  (xflux_fab(i+1,j,k,n) - xflux_fab(i,j,k,n)) / dx[0],
                                                     + (yflux_fab(i,j+1,k,n) - yflux_fab(i,j,k,n)) / dx[1],
                                                     + (zflux_fab(i,j,k+1,n) - zflux_fab(i,j,k,n)) / dx[2]) 
                                                                                                            )
                                    + dt*source_fab(i,j,k,n) );
            
        }); // [1:3 indices are not valuable -- momentum flux]

        // momentum flux
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momx(i,j,k) = (2./3.)*(0.5*momx(i,j,k) + momp2x(i,j,k))
                  -(2./3.)*dt*(cenx_u(i,j,k) - cenx_u(i-1,j,k))/dx[0]
                  -(2./3.)*dt*(edgey_u(i,j+1,k) - edgey_u(i,j,k))/dx[1]
                  -(2./3.)*dt*(edgez_u(i,j,k+1) - edgez_u(i,j,k))/dx[2]
                  +0.5*(2./3.)*dt*grav_gpu[0]*(cup2_fab(i-1,j,k,0)+cup2_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momy(i,j,k) = (2./3.)*(0.5*momy(i,j,k) + momp2y(i,j,k))
                  -(2./3.)*dt*(edgex_v(i+1,j,k) - edgex_v(i,j,k))/dx[0]
                  -(2./3.)*dt*(ceny_v(i,j,k) - cenx_u(i,j-1,k))/dx[1]
                  -(2./3.)*dt*(edgez_v(i,j,k+1) - edgez_v(i,j,k))/dx[2]
                  +0.5*(2/3.)*dt*grav_gpu[1]*(cup2_fab(i,j-1,k,0)+cup2_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momz(i,j,k) = (2./3.)*(0.5*momz(i,j,k) + momp2z(i,j,k))
                  -(2./3.)*dt*(edgex_w(i+1,j,k) - edgex_w(i,j,k))/dx[0]
                  -(2./3.)*dt*(edgey_w(i,j+1,k) - edgey_w(i,j,k))/dx[1]
                  -(2./3.)*dt*(cenz_w(i,j,k) - cenx_u(i,j,k-1))/dx[2]
                  +0.5*(2./3.)*dt*grav_gpu[2]*(cup2_fab(i,j,k-1,0)+cup2_fab(i,j,k,0));
        });
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cup2_fab(i,j,k,4) += 0.5 * (2./3.) * dt * (  grav_gpu[0]*(momx(i+1,j,k)+momx(i,j,k))
                                                    + grav_gpu[1]*(momy(i,j+1,k)+momy(i,j,k))
                                                    + grav_gpu[2]*(momz(i,j,k+1)+momz(i,j,k)) );
        });
    }

    conservedToPrimitiveStag(prim, facevel, cup, cupmom);

    // Set BC: 1) fill boundary 2) physical
    cu.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());

    setBC(prim, cu);
}

