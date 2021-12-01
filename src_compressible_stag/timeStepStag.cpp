#include "compressible_functions.H"
#include "compressible_functions_stag.H"

#include "common_functions.H"

#include "rng_functions.H"
#include <AMReX_VisMF.H>

void RK3stepStag(MultiFab& cu, 
                 std::array< MultiFab, AMREX_SPACEDIM >& cumom,
                 MultiFab& prim, std::array< MultiFab, AMREX_SPACEDIM >& vel,
                 MultiFab& source,
                 MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
                 MultiFab& chi, MultiFab& D,
                 std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                 std::array< MultiFab, 2 >& edgeflux_x,
                 std::array< MultiFab, 2 >& edgeflux_y,
                 std::array< MultiFab, 2 >& edgeflux_z,
                 std::array< MultiFab, AMREX_SPACEDIM>& cenflux,
                 const amrex::Geometry geom, const amrex::Real dt, const int step)
{
    BL_PROFILE_VAR("RK3stepStag()",RK3stepStag);

    MultiFab cup (cu.boxArray(),cu.DistributionMap(),nvars,ngc);
    MultiFab cup2(cu.boxArray(),cu.DistributionMap(),nvars,ngc);
    cup.setVal(0.0,0,nvars,ngc);
    cup2.setVal(0.0,0,nvars,ngc);
    //cup.setVal(rho0,0,1,ngc);
    //cup2.setVal(rho0,0,1,ngc);

    std::array< MultiFab, AMREX_SPACEDIM > cupmom;
    std::array< MultiFab, AMREX_SPACEDIM > cup2mom;
    AMREX_D_TERM(cupmom[0].define(convert(cu.boxArray(),nodal_flag_x), cu.DistributionMap(), 1, ngc);,
                 cupmom[1].define(convert(cu.boxArray(),nodal_flag_y), cu.DistributionMap(), 1, ngc);,
                 cupmom[2].define(convert(cu.boxArray(),nodal_flag_z), cu.DistributionMap(), 1, ngc););

    AMREX_D_TERM(cup2mom[0].define(convert(cu.boxArray(),nodal_flag_x), cu.DistributionMap(), 1, ngc);,
                 cup2mom[1].define(convert(cu.boxArray(),nodal_flag_y), cu.DistributionMap(), 1, ngc);,
                 cup2mom[2].define(convert(cu.boxArray(),nodal_flag_z), cu.DistributionMap(), 1, ngc););
    
    AMREX_D_TERM(cupmom[0].setVal(0.0);,
                 cupmom[1].setVal(0.0);,
                 cupmom[2].setVal(0.0););

    AMREX_D_TERM(cup2mom[0].setVal(0.0);,
                 cup2mom[1].setVal(0.0);,
                 cup2mom[2].setVal(0.0););

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    
    /////////////////////////////////////////////////////
    // Setup stochastic flux MultiFabs
    std::array< MultiFab, AMREX_SPACEDIM > stochface;
    AMREX_D_TERM(stochface[0].define(convert(cu.boxArray(),nodal_flag_x), cu.DistributionMap(), nvars, 0);,
                 stochface[1].define(convert(cu.boxArray(),nodal_flag_y), cu.DistributionMap(), nvars, 0);,
                 stochface[2].define(convert(cu.boxArray(),nodal_flag_z), cu.DistributionMap(), nvars, 0););
    
    std::array< MultiFab, 2 > stochedge_x; 
    std::array< MultiFab, 2 > stochedge_y; 
    std::array< MultiFab, 2 > stochedge_z; 

    stochedge_x[0].define(convert(cu.boxArray(),nodal_flag_xy), cu.DistributionMap(), 1, 0); 
    stochedge_x[1].define(convert(cu.boxArray(),nodal_flag_xz), cu.DistributionMap(), 1, 0);
             
    stochedge_y[0].define(convert(cu.boxArray(),nodal_flag_xy), cu.DistributionMap(), 1, 0);
    stochedge_y[1].define(convert(cu.boxArray(),nodal_flag_yz), cu.DistributionMap(), 1, 0);

    stochedge_z[0].define(convert(cu.boxArray(),nodal_flag_xz), cu.DistributionMap(), 1, 0);
    stochedge_z[1].define(convert(cu.boxArray(),nodal_flag_yz), cu.DistributionMap(), 1, 0);

    std::array< MultiFab, AMREX_SPACEDIM > stochcen;
    AMREX_D_TERM(stochcen[0].define(cu.boxArray(),cu.DistributionMap(),1,1);,
                 stochcen[1].define(cu.boxArray(),cu.DistributionMap(),1,1);,
                 stochcen[2].define(cu.boxArray(),cu.DistributionMap(),1,1););
    /////////////////////////////////////////////////////

    /////////////////////////////////////////////////////
    // Initialize white noise weighted fields
    // weights for stochastic fluxes; swgt2 changes each stage
    amrex::Vector< amrex::Real > stoch_weights;
    amrex::Real swgt1, swgt2;
    swgt1 = 1.0;

    // field "A"
    std::array< MultiFab, AMREX_SPACEDIM > stochface_A;
    AMREX_D_TERM(stochface_A[0].define(stochface[0].boxArray(), stochface[0].DistributionMap(), nvars, 0);,
                 stochface_A[1].define(stochface[1].boxArray(), stochface[1].DistributionMap(), nvars, 0);,
                 stochface_A[2].define(stochface[2].boxArray(), stochface[2].DistributionMap(), nvars, 0););

    AMREX_D_TERM(stochface_A[0].setVal(0.0);,
                 stochface_A[1].setVal(0.0);,
                 stochface_A[2].setVal(0.0););
    
    std::array< MultiFab, 2 > stochedge_x_A; 
    std::array< MultiFab, 2 > stochedge_y_A; 
    std::array< MultiFab, 2 > stochedge_z_A; 

    stochedge_x_A[0].define(stochedge_x[0].boxArray(), stochedge_x[0].DistributionMap(), 1, 0); 
    stochedge_x_A[1].define(stochedge_x[1].boxArray(), stochedge_x[1].DistributionMap(), 1, 0);
             
    stochedge_y_A[0].define(stochedge_y[0].boxArray(), stochedge_y[0].DistributionMap(), 1, 0);
    stochedge_y_A[1].define(stochedge_y[1].boxArray(), stochedge_y[1].DistributionMap(), 1, 0);

    stochedge_z_A[0].define(stochedge_z[0].boxArray(), stochedge_z[0].DistributionMap(), 1, 0);
    stochedge_z_A[1].define(stochedge_z[1].boxArray(), stochedge_z[1].DistributionMap(), 1, 0);

    stochedge_x_A[0].setVal(0.0); stochedge_x_A[1].setVal(0.0);
    stochedge_y_A[0].setVal(0.0); stochedge_y_A[1].setVal(0.0);
    stochedge_z_A[0].setVal(0.0); stochedge_z_A[1].setVal(0.0);

    std::array< MultiFab, AMREX_SPACEDIM > stochcen_A;
    AMREX_D_TERM(stochcen_A[0].define(stochcen[0].boxArray(),stochcen[0].DistributionMap(),1,1);,
                 stochcen_A[1].define(stochcen[1].boxArray(),stochcen[1].DistributionMap(),1,1);,
                 stochcen_A[2].define(stochcen[2].boxArray(),stochcen[2].DistributionMap(),1,1););

    AMREX_D_TERM(stochcen_A[0].setVal(0.0);,
                 stochcen_A[1].setVal(0.0);,
                 stochcen_A[2].setVal(0.0););

    // field "B"
    std::array< MultiFab, AMREX_SPACEDIM > stochface_B;
    AMREX_D_TERM(stochface_B[0].define(stochface[0].boxArray(), stochface[0].DistributionMap(), nvars, 0);,
                 stochface_B[1].define(stochface[1].boxArray(), stochface[1].DistributionMap(), nvars, 0);,
                 stochface_B[2].define(stochface[2].boxArray(), stochface[2].DistributionMap(), nvars, 0););

    AMREX_D_TERM(stochface_B[0].setVal(0.0);,
                 stochface_B[1].setVal(0.0);,
                 stochface_B[2].setVal(0.0););
    
    std::array< MultiFab, 2 > stochedge_x_B; 
    std::array< MultiFab, 2 > stochedge_y_B; 
    std::array< MultiFab, 2 > stochedge_z_B; 

    stochedge_x_B[0].define(stochedge_x[0].boxArray(), stochedge_x[0].DistributionMap(), 1, 0); 
    stochedge_x_B[1].define(stochedge_x[1].boxArray(), stochedge_x[1].DistributionMap(), 1, 0);
             
    stochedge_y_B[0].define(stochedge_y[0].boxArray(), stochedge_y[0].DistributionMap(), 1, 0);
    stochedge_y_B[1].define(stochedge_y[1].boxArray(), stochedge_y[1].DistributionMap(), 1, 0);

    stochedge_z_B[0].define(stochedge_z[0].boxArray(), stochedge_z[0].DistributionMap(), 1, 0);
    stochedge_z_B[1].define(stochedge_z[1].boxArray(), stochedge_z[1].DistributionMap(), 1, 0);

    stochedge_x_B[0].setVal(0.0); stochedge_x_B[1].setVal(0.0);
    stochedge_y_B[0].setVal(0.0); stochedge_y_B[1].setVal(0.0);
    stochedge_z_B[0].setVal(0.0); stochedge_z_B[1].setVal(0.0);

    std::array< MultiFab, AMREX_SPACEDIM > stochcen_B;
    AMREX_D_TERM(stochcen_B[0].define(stochcen[0].boxArray(),stochcen[0].DistributionMap(),1,1);,
                 stochcen_B[1].define(stochcen[1].boxArray(),stochcen[1].DistributionMap(),1,1);,
                 stochcen_B[2].define(stochcen[2].boxArray(),stochcen[2].DistributionMap(),1,1););

    AMREX_D_TERM(stochcen_B[0].setVal(0.0);,
                 stochcen_B[1].setVal(0.0);,
                 stochcen_B[2].setVal(0.0););

    // fill random numbers (can skip density component 0)
    if (do_1D) { // 1D need only for x- face 
        for(int i=1;i<nvars;i++) {
            MultiFabFillRandom(stochface_A[0], i, 1.0, geom);
            MultiFabFillRandom(stochface_B[0], i, 1.0, geom);
        }
    }
    else if (do_2D) { // 2D need only for x- and y- faces
        for(int i=1;i<nvars;i++) {
            MultiFabFillRandom(stochface_A[0], i, 1.0, geom);
            MultiFabFillRandom(stochface_B[0], i, 1.0, geom);
            MultiFabFillRandom(stochface_A[1], i, 1.0, geom);
            MultiFabFillRandom(stochface_B[1], i, 1.0, geom);
        }
    }
    else { // 3D
        for(int d=0;d<AMREX_SPACEDIM;d++) {
            for(int i=1;i<nvars;i++) {
                MultiFabFillRandom(stochface_A[d], i, 1.0, geom);
                MultiFabFillRandom(stochface_B[d], i, 1.0, geom);
            }
        }
    }

    if (do_1D) { // 1D no transverse shear fluxes
    }
    else if (do_2D) { // 2D only xy-shear
        MultiFabFillRandom(stochedge_x_A[0], 0, 1.0, geom);
        MultiFabFillRandom(stochedge_x_B[0], 0, 1.0, geom);
        MultiFabFillRandom(stochedge_y_A[0], 0, 1.0, geom);
        MultiFabFillRandom(stochedge_y_B[0], 0, 1.0, geom);
    }
    else { // 3D
        for (int i=0; i<2; i++) {
            MultiFabFillRandom(stochedge_x_A[i], 0, 1.0, geom);
            MultiFabFillRandom(stochedge_x_B[i], 0, 1.0, geom);
            MultiFabFillRandom(stochedge_y_A[i], 0, 1.0, geom);
            MultiFabFillRandom(stochedge_y_B[i], 0, 1.0, geom);
            MultiFabFillRandom(stochedge_z_A[i], 0, 1.0, geom);
            MultiFabFillRandom(stochedge_z_B[i], 0, 1.0, geom);
        }
    }

    if (do_1D) { // 1D no v_x and w_z stochastic terms
        MultiFabFillRandom(stochcen_A[0], 0, 1.0, geom, 1);
        MultiFabFillRandom(stochcen_B[0], 0, 1.0, geom, 1);
    }
    else if (do_2D) { // 2D simulation no w_z stochastic term
        MultiFabFillRandom(stochcen_A[0], 0, 2.0, geom, 1);
        MultiFabFillRandom(stochcen_B[0], 0, 2.0, geom, 1);
        MultiFabFillRandom(stochcen_A[1], 0, 2.0, geom, 1);
        MultiFabFillRandom(stochcen_B[1], 0, 2.0, geom, 1);
    }
    else { // 3D
        for (int i=0; i<3; i++) {
            MultiFabFillRandom(stochcen_A[i], 0, 2.0, geom, 1);
            MultiFabFillRandom(stochcen_B[i], 0, 2.0, geom, 1);
        }
    }

    /////////////////////////////////////////////////////

    /////////////////////////////////////////////////////
    // Perform weighting of white noise fields
    // Set stochastic weights
    swgt2 = ( 2.0*std::sqrt(2.0) + 1.0*std::sqrt(3.0) ) / 5.0;
    stoch_weights = {swgt1, swgt2};

    // apply weights (only energy and ns-1 species)
    AMREX_D_TERM(stochface[0].setVal(0.0);,
                 stochface[1].setVal(0.0);,
                 stochface[2].setVal(0.0););

    stochedge_x[0].setVal(0.0); stochedge_x[1].setVal(0.0);
    stochedge_y[0].setVal(0.0); stochedge_y[1].setVal(0.0);
    stochedge_z[0].setVal(0.0); stochedge_z[1].setVal(0.0);

    AMREX_D_TERM(stochcen[0].setVal(0.0);,
                 stochcen[1].setVal(0.0);,
                 stochcen[2].setVal(0.0););

    for (int d=0;d<AMREX_SPACEDIM;d++) {
	    MultiFab::LinComb(stochface[d], 
            stoch_weights[0], stochface_A[d], 1, 
            stoch_weights[1], stochface_B[d], 1,
            1, nvars-1, 0);
    }
    for (int i=0;i<2;i++) {
        MultiFab::LinComb(stochedge_x[i],
            stoch_weights[0], stochedge_x_A[i], 0,
            stoch_weights[1], stochedge_x_B[i], 0,
            0, 1, 0);
        MultiFab::LinComb(stochedge_y[i],
            stoch_weights[0], stochedge_y_A[i], 0,
            stoch_weights[1], stochedge_y_B[i], 0,
            0, 1, 0);
        MultiFab::LinComb(stochedge_z[i],
            stoch_weights[0], stochedge_z_A[i], 0,
            stoch_weights[1], stochedge_z_B[i], 0,
            0, 1, 0);
    }
    for (int i=0;i<3;i++) {
        MultiFab::LinComb(stochcen[i],
            stoch_weights[0], stochcen_A[i], 0,
            stoch_weights[1], stochcen_B[i], 0,
            0, 1, 1);
    }
    /////////////////////////////////////////////////////

    calculateTransportCoeffs(prim, eta, zeta, kappa, chi, D);

    calculateFluxStag(cu, cumom, prim, vel, eta, zeta, kappa, chi, D, 
        faceflux, edgeflux_x, edgeflux_y, edgeflux_z, cenflux, 
        stochface, stochedge_x, stochedge_y, stochedge_z, stochcen, 
        geom, stoch_weights,dt);

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
                    +0.5*dt*grav[0]*(cu_fab(i-1,j,k,0)+cu_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            mompy(i,j,k) = momy(i,j,k)
                    -dt*(edgex_v(i+1,j,k) - edgex_v(i,j,k))/dx[0]
                    -dt*(ceny_v(i,j,k) - ceny_v(i,j-1,k))/dx[1]
                    -dt*(edgez_v(i,j,k+1) - edgez_v(i,j,k))/dx[2]
                    +0.5*dt*grav[1]*(cu_fab(i,j-1,k,0)+cu_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            mompz(i,j,k) = momz(i,j,k)
                    -dt*(edgex_w(i+1,j,k) - edgex_w(i,j,k))/dx[0]
                    -dt*(edgey_w(i,j+1,k) - edgey_w(i,j,k))/dx[1]
                    -dt*(cenz_w(i,j,k) - cenz_w(i,j,k-1))/dx[2]
                    +0.5*dt*grav[2]*(cu_fab(i,j,k-1,0)+cu_fab(i,j,k,0));
        });
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cup_fab(i,j,k,4) += 0.5 * dt * (  grav[0]*(momx(i+1,j,k)+momx(i,j,k))
                                            + grav[1]*(momy(i,j+1,k)+momy(i,j,k))
                                            + grav[2]*(momz(i,j,k+1)+momz(i,j,k)) );
        });
    }

    // Set the correct momentum at the walls and ghost 
    for (int i=0; i<AMREX_SPACEDIM; i++) {
        BCMassTempPress(prim, cup, geom, i);
        BCMomNormal(cupmom[i], vel[i], cup, geom, i);
        BCMomTrans(cupmom[i], vel[i], geom, i);
    }

    // Fill boundaries for conserved variables
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        cupmom[d].FillBoundary(geom.periodicity());
    }
    cup.FillBoundary(geom.periodicity());

    // Conserved to primitive conversion (also writes momemtun at cell centers as averages of neighboring faces)
    conservedToPrimitiveStag(prim, vel, cup, cupmom);

    // Fill boundaries for primitive variables (also do for conserved, since cell-centered momentum is written above)
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        vel[d].FillBoundary(geom.periodicity());
    }
    prim.FillBoundary(geom.periodicity());
    cup.FillBoundary(geom.periodicity());

    // Correctly set momentum and velocity at the walls & temperature, pressure, density & mass/mole fractions in ghost cells
    setBCStag(prim, cup, cupmom, vel, geom);

    // Compute transport coefs after setting BCs
    calculateTransportCoeffs(prim, eta, zeta, kappa, chi, D);

    ///////////////////////////////////////////////////////////
    // Perform weighting of white noise fields

    // Set stochastic weights
    swgt2 = ( -4.0*std::sqrt(2.0) + 3.0*std::sqrt(3.0) ) / 5.0;
    stoch_weights = {swgt1, swgt2};

    // apply weights (only energy and ns-1 species)
    AMREX_D_TERM(stochface[0].setVal(0.0);,
                 stochface[1].setVal(0.0);,
                 stochface[2].setVal(0.0););

    stochedge_x[0].setVal(0.0); stochedge_x[1].setVal(0.0);
    stochedge_y[0].setVal(0.0); stochedge_y[1].setVal(0.0);
    stochedge_z[0].setVal(0.0); stochedge_z[1].setVal(0.0);

    AMREX_D_TERM(stochcen[0].setVal(0.0);,
                 stochcen[1].setVal(0.0);,
                 stochcen[2].setVal(0.0););

    for (int d=0;d<AMREX_SPACEDIM;d++) {
	    MultiFab::LinComb(stochface[d], 
            stoch_weights[0], stochface_A[d], 1, 
            stoch_weights[1], stochface_B[d], 1,
            1, nvars-1, 0);
    }
    for (int i=0;i<2;i++) {
        MultiFab::LinComb(stochedge_x[i],
            stoch_weights[0], stochedge_x_A[i], 0,
            stoch_weights[1], stochedge_x_B[i], 0,
            0, 1, 0);
        MultiFab::LinComb(stochedge_y[i],
            stoch_weights[0], stochedge_y_A[i], 0,
            stoch_weights[1], stochedge_y_B[i], 0,
            0, 1, 0);
        MultiFab::LinComb(stochedge_z[i],
            stoch_weights[0], stochedge_z_A[i], 0,
            stoch_weights[1], stochedge_z_B[i], 0,
            0, 1, 0);
    }
    for (int i=0;i<3;i++) {
        MultiFab::LinComb(stochcen[i],
            stoch_weights[0], stochcen_A[i], 0,
            stoch_weights[1], stochcen_B[i], 0,
            0, 1, 1);
    }
    ///////////////////////////////////////////////////////////

    calculateFluxStag(cup, cupmom, prim, vel, eta, zeta, kappa, chi, D, 
        faceflux, edgeflux_x, edgeflux_y, edgeflux_z, cenflux, 
        stochface, stochedge_x, stochedge_y, stochedge_z, stochcen, 
        geom, stoch_weights,dt);

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
                    +0.5*0.25*dt*grav[0]*(cup_fab(i-1,j,k,0)+cup_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momp2y(i,j,k) = 0.25*3.0*momy(i,j,k) + 0.25*mompy(i,j,k)
                    -0.25*dt*(edgex_v(i+1,j,k) - edgex_v(i,j,k))/dx[0]
                    -0.25*dt*(ceny_v(i,j,k) - ceny_v(i,j-1,k))/dx[1]
                    -0.25*dt*(edgez_v(i,j,k+1) - edgez_v(i,j,k))/dx[2]
                    +0.5*0.25*dt*grav[1]*(cup_fab(i,j-1,k,0)+cup_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momp2z(i,j,k) = 0.25*3.0*momz(i,j,k) + 0.25*mompz(i,j,k)
                    -0.25*dt*(edgex_w(i+1,j,k) - edgex_w(i,j,k))/dx[0]
                    -0.25*dt*(edgey_w(i,j+1,k) - edgey_w(i,j,k))/dx[1]
                    -0.25*dt*(cenz_w(i,j,k) - cenz_w(i,j,k-1))/dx[2]
                    +0.5*0.25*dt*grav[2]*(cup_fab(i,j,k-1,0)+cup_fab(i,j,k,0));
        });
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cup2_fab(i,j,k,4) += 0.5 * 0.25 * dt * (  grav[0]*(mompx(i+1,j,k)+mompx(i,j,k))
                                                    + grav[1]*(mompy(i,j+1,k)+mompy(i,j,k))
                                                    + grav[2]*(mompz(i,j,k+1)+mompz(i,j,k)) );
        });
    }
        
    // Set the correct momentum at the walls 
    for (int i=0; i<AMREX_SPACEDIM; i++) {
        BCMassTempPress(prim, cup2, geom, i);
        BCMomNormal(cup2mom[i], vel[i], cup2, geom, i);
        BCMomTrans(cup2mom[i], vel[i], geom, i);
    }

    // Fill  boundaries for conserved variables
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        cup2mom[d].FillBoundary(geom.periodicity());
    }
    cup2.FillBoundary(geom.periodicity()); 

    // Conserved to primitive conversion (also writes momemtun at cell centers as averages of neighboring faces)
    conservedToPrimitiveStag(prim, vel, cup2, cup2mom);

    // Fill  boundaries for primitive variables (also do for conserved, since cell-centered momentum is written above)
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        vel[d].FillBoundary(geom.periodicity());
    }
    prim.FillBoundary(geom.periodicity());
    cup2.FillBoundary(geom.periodicity());  

    // Correctly set momentum and velocity at the walls & temperature, pressure, density & mass/mole fractions in ghost cells
    setBCStag(prim, cup2, cup2mom, vel, geom);

    // Compute transport coefs after setting BCs
    calculateTransportCoeffs(prim, eta, zeta, kappa, chi, D);

    ///////////////////////////////////////////////////////////
    // Perform weighting of white noise fields

    // Set stochastic weights
    swgt2 = ( 1.0*std::sqrt(2.0) - 2.0*std::sqrt(3.0) ) / 10.0;
    stoch_weights = {swgt1, swgt2};

    // apply weights (only energy and ns-1 species)
    AMREX_D_TERM(stochface[0].setVal(0.0);,
                 stochface[1].setVal(0.0);,
                 stochface[2].setVal(0.0););

    stochedge_x[0].setVal(0.0); stochedge_x[1].setVal(0.0);
    stochedge_y[0].setVal(0.0); stochedge_y[1].setVal(0.0);
    stochedge_z[0].setVal(0.0); stochedge_z[1].setVal(0.0);

    AMREX_D_TERM(stochcen[0].setVal(0.0);,
                 stochcen[1].setVal(0.0);,
                 stochcen[2].setVal(0.0););

    for (int d=0;d<AMREX_SPACEDIM;d++) {
	    MultiFab::LinComb(stochface[d], 
            stoch_weights[0], stochface_A[d], 1, 
            stoch_weights[1], stochface_B[d], 1,
            1, nvars-1, 0);
    }
    for (int i=0;i<2;i++) {
        MultiFab::LinComb(stochedge_x[i],
            stoch_weights[0], stochedge_x_A[i], 0,
            stoch_weights[1], stochedge_x_B[i], 0,
            0, 1, 0);
        MultiFab::LinComb(stochedge_y[i],
            stoch_weights[0], stochedge_y_A[i], 0,
            stoch_weights[1], stochedge_y_B[i], 0,
            0, 1, 0);
        MultiFab::LinComb(stochedge_z[i],
            stoch_weights[0], stochedge_z_A[i], 0,
            stoch_weights[1], stochedge_z_B[i], 0,
            0, 1, 0);
    }
    for (int i=0;i<3;i++) {
        MultiFab::LinComb(stochcen[i],
            stoch_weights[0], stochcen_A[i], 0,
            stoch_weights[1], stochcen_B[i], 0,
            0, 1, 1);
    }
    ///////////////////////////////////////////////////////////

    calculateFluxStag(cup2, cup2mom, prim, vel, eta, zeta, kappa, chi, D, 
        faceflux, edgeflux_x, edgeflux_y, edgeflux_z, cenflux, 
        stochface, stochedge_x, stochedge_y, stochedge_z, stochcen, 
        geom, stoch_weights,dt);

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
                  +0.5*(2./3.)*dt*grav[0]*(cup2_fab(i-1,j,k,0)+cup2_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momy(i,j,k) = (2./3.)*(0.5*momy(i,j,k) + momp2y(i,j,k))
                  -(2./3.)*dt*(edgex_v(i+1,j,k) - edgex_v(i,j,k))/dx[0]
                  -(2./3.)*dt*(ceny_v(i,j,k) - ceny_v(i,j-1,k))/dx[1]
                  -(2./3.)*dt*(edgez_v(i,j,k+1) - edgez_v(i,j,k))/dx[2]
                  +0.5*(2/3.)*dt*grav[1]*(cup2_fab(i,j-1,k,0)+cup2_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momz(i,j,k) = (2./3.)*(0.5*momz(i,j,k) + momp2z(i,j,k))
                  -(2./3.)*dt*(edgex_w(i+1,j,k) - edgex_w(i,j,k))/dx[0]
                  -(2./3.)*dt*(edgey_w(i,j+1,k) - edgey_w(i,j,k))/dx[1]
                  -(2./3.)*dt*(cenz_w(i,j,k) - cenz_w(i,j,k-1))/dx[2]
                  +0.5*(2./3.)*dt*grav[2]*(cup2_fab(i,j,k-1,0)+cup2_fab(i,j,k,0));
        });
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cu_fab(i,j,k,4) += 0.5 * (2./3.) * dt * (  grav[0]*(momp2x(i+1,j,k)+momp2x(i,j,k))
                                                    + grav[1]*(momp2y(i,j+1,k)+momp2y(i,j,k))
                                                    + grav[2]*(momp2z(i,j,k+1)+momp2z(i,j,k)) );
        });
    }

    // Set the correct momentum at the walls 
    for (int i=0; i<AMREX_SPACEDIM; i++) {
        BCMassTempPress(prim, cu, geom, i);
        BCMomNormal(cumom[i], vel[i], cu, geom, i);
        BCMomTrans(cumom[i], vel[i], geom, i);
    }

    // Fill  boundaries for conserved variables
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        cumom[d].FillBoundary(geom.periodicity());
    }
    cu.FillBoundary(geom.periodicity());

    // Conserved to primitive conversion (also writes momemtun at cell centers as averages of neighboring faces)
    conservedToPrimitiveStag(prim, vel, cu, cumom);

    // Fill  boundaries for primitive variables (also do for conserved, since cell-centered momentum is written above)
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        vel[d].FillBoundary(geom.periodicity());
    }
    prim.FillBoundary(geom.periodicity());
    cu.FillBoundary(geom.periodicity()); 

    // Membrane setup
    if (membrane_cell >= 0) {
        doMembraneStag(cu,cumom,prim,vel,faceflux,geom,dt);
    }

    // Correctly set momentum and velocity at the walls & temperature, pressure, density & mass/mole fractions in ghost cells
    setBCStag(prim, cu, cumom, vel, geom);

}
