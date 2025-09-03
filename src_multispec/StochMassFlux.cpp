#include "rng_functions.H"
#include "gmres_functions.H"
#include "common_functions.H"
#include "multispec_functions.H"

#include "StochMassFlux.H"


// initialize n_rngs, geom
// build MultiFabs to hold random numbers
StochMassFlux::StochMassFlux(BoxArray ba_in, DistributionMapping dmap_in, Geometry geom_in,
                             int n_rngs_in) {

    BL_PROFILE_VAR("StochMassFlux()",StochMassFlux);

    n_rngs = n_rngs_in;
    geom = geom_in;

    stoch_W_fc.resize(n_rngs);

    // Here we store all the random number stages at all spatial locations
    for (int i=0; i<n_rngs; ++i) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            stoch_W_fc[i][d].define(convert(ba_in,nodal_flag_dir[d]),dmap_in,nspecies,0);
            stoch_W_fc[i][d].setVal(0.);
        }
    }

    // Temporary storage for linear combinations of random number stages
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        stoch_W_fc_weighted[d].define(convert(ba_in,nodal_flag_dir[d]),dmap_in,nspecies,0);
        stoch_W_fc_weighted[d].setVal(0.);
    }
}

// create weighted sum of stage RNGs and store in stoch_W_fc_weighted
void StochMassFlux::weightMassFlux(Vector< amrex::Real > weights) {

    BL_PROFILE_VAR("weightMassFlux()",weightMassFlux);

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        stoch_W_fc_weighted[d].setVal(0.);
    }

    // add weighted contribution of fluxes
    for (int i=0; i<n_rngs; ++i) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Saxpy(stoch_W_fc_weighted[d], weights[i], stoch_W_fc[i][d], 0, 0, nspecies, 0);
        }
    }
}

// fill stoch_W_fc with random numbers
void StochMassFlux::fillMassStochastic() {

    BL_PROFILE_VAR("fillMassStochastic()",fillMassStochastic);

    for (int i=0; i<n_rngs; ++i) {
        for (int n=0; n<AMREX_SPACEDIM; ++n) {
            for (int comp=0; comp<nspecies; comp++) {
                MultiFabFillRandom(stoch_W_fc[i][n],comp,1.0,geom);
            }
        }
    }
}

// scale random numbers that lie on physical boundaries appropriately
void StochMassFlux::StochMassFluxBC() {

    BL_PROFILE_VAR("StochMassFluxBC()",StochMassFluxBC);

    // lo-x domain boundary
    if (bc_mass_lo[0] == 1 || bc_mass_lo[0] == 2 || bc_mass_lo[0] == 4) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = reservoir   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_lo[0] == 1 || bc_mass_lo[0] == 4) ? 0. : sqrt(2.);

        // domain grown nodally based on stoch_W_fc_weighted[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), stoch_W_fc_weighted[0].ixType());

        // this is the x-lo domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xlo = amrex::bdryNode(dom_x, Orientation(0, Orientation::low));

        for (MFIter mfi(stoch_W_fc_weighted[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xlo;
            Array4<Real> const& flux = (stoch_W_fc_weighted[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    flux(i,j,k,n) *= factor;
                });
            }
        }
    }
    else if (bc_mass_lo[0] != -1) {
        Abort("StochMassFluxBC unsupported bc type");
    }

    // hi-x domain boundary
    if (bc_mass_hi[0] == 1 || bc_mass_hi[0] == 2 || bc_mass_hi[0] == 4) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = reservoir   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_hi[0] == 1 || bc_mass_hi[0] == 4) ? 0. : sqrt(2.);

        // domain grown nodally based on stoch_W_fc_weighted[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), stoch_W_fc_weighted[0].ixType());

        // this is the x-hi domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xhi = amrex::bdryNode(dom_x, Orientation(0, Orientation::high));

        for (MFIter mfi(stoch_W_fc_weighted[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xhi;
            Array4<Real> const& flux = (stoch_W_fc_weighted[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    flux(i,j,k,n) *= factor;
                });
            }
        }

    }
    else if (bc_mass_hi[0] != -1) {
        Abort("StochMassFluxBC unsupported bc type");
    }

    // lo-y domain boundary
    if (bc_mass_lo[1] == 1 || bc_mass_lo[1] == 2 || bc_mass_lo[1] == 4) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = reservoir   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_lo[1] == 1 || bc_mass_lo[1] == 4) ? 0. : sqrt(2.);

        // domain grown nodally based on stoch_W_fc_weighted[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), stoch_W_fc_weighted[1].ixType());

        // this is the y-lo domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_ylo = amrex::bdryNode(dom_y, Orientation(1, Orientation::low));

        for (MFIter mfi(stoch_W_fc_weighted[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_ylo;
            Array4<Real> const& flux = (stoch_W_fc_weighted[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    flux(i,j,k,n) *= factor;
                });
            }
        }
    }
    else if (bc_mass_lo[1] != -1) {
        Abort("StochMassFluxBC unsupported bc type");
    }

    // hi-y domain boundary
    if (bc_mass_hi[1] == 1 || bc_mass_hi[1] == 2 || bc_mass_hi[1] == 4) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = reservoir   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_hi[1] == 1 || bc_mass_hi[1] == 4) ? 0. : sqrt(2.);

        // domain grown nodally based on stoch_W_fc_weighted[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), stoch_W_fc_weighted[1].ixType());

        // this is the y-hi domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yhi = amrex::bdryNode(dom_y, Orientation(1, Orientation::high));

        for (MFIter mfi(stoch_W_fc_weighted[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yhi;
            Array4<Real> const& flux = (stoch_W_fc_weighted[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    flux(i,j,k,n) *= factor;
                });
            }
        }
    }
    else if (bc_mass_hi[1] != -1) {
        Abort("StochMassFluxBC unsupported bc type");
    }

#if (AMREX_SPACEDIM == 3)

    // lo-z domain boundary
    if (bc_mass_lo[2] == 1 || bc_mass_lo[2] == 2 || bc_mass_lo[2] == 4) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = reservoir   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_lo[2] == 1 || bc_mass_lo[2] == 4) ? 0. : sqrt(2.);

        // domain grown nodally based on stoch_W_fc_weighted[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), stoch_W_fc_weighted[2].ixType());

        // this is the z-lo domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zlo = amrex::bdryNode(dom_z, Orientation(2, Orientation::low));

        for (MFIter mfi(stoch_W_fc_weighted[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zlo;
            Array4<Real> const& flux = (stoch_W_fc_weighted[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    flux(i,j,k,n) *= factor;
                });
            }
        }
    }
    else if (bc_mass_lo[2] != -1) {
        Abort("StochMassFluxBC unsupported bc type");
    }

    // hi-z domain boundary
    if (bc_mass_hi[2] == 1 || bc_mass_hi[2] == 2 || bc_mass_hi[2] == 4) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = reservoir   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_hi[2] == 1 || bc_mass_hi[2] == 4) ? 0. : sqrt(2.);

        // domain grown nodally based on stoch_W_fc_weighted[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), stoch_W_fc_weighted[2].ixType());

        // this is the z-hi domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zhi = amrex::bdryNode(dom_z, Orientation(2, Orientation::high));

        for (MFIter mfi(stoch_W_fc_weighted[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zhi;
            Array4<Real> const& flux = (stoch_W_fc_weighted[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    flux(i,j,k,n) *= factor;
                });
            }
        }
    }
    else if (bc_mass_hi[2] != -1) {
        Abort("StochMassFluxBC unsupported bc type");
    }
#endif

}

void StochMassFlux::StochMassFluxDiv(const MultiFab& rho,
                                     const MultiFab& rhotot,
                                     const std::array<MultiFab, AMREX_SPACEDIM >& sqrtLonsager_fc,
                                     MultiFab& stoch_mass_fluxdiv,
                                     std::array<MultiFab, AMREX_SPACEDIM >& stoch_mass_flux,
                                     const amrex::Real& dt,
                                     const Vector< amrex::Real >& weights,
                                     int increment) {

    BL_PROFILE_VAR("StochMassFluxDiv()",StochMassFluxDiv);

    // take linear combination of mflux multifabs at each stage
    // store result in stoch_W_fc_weighted
    // note we do this here instead of in fillMassStochastic because we have weights now
    StochMassFlux::weightMassFlux(weights);

    // multiply noise stored in stoch_W_fc_weighted
    // on walls by 0 or on reservoirs by sqrt(2)
    StochMassFlux::StochMassFluxBC();

    // copy weighted stochastic noise stages into stoch_mass_flux
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(stoch_mass_flux[d], stoch_W_fc_weighted[d], 0, 0, nspecies, 0);
    }

    const Real* dx = geom.CellSize();
    Real dVol = (AMREX_SPACEDIM==2) ? dx[0]*dx[1]*cell_depth : dx[0]*dx[1]*dx[2];
    Real variance = sqrt(2.*k_B*variance_coef_mass/(dVol*dt));

    // compute variance X sqrtLonsager_fc X stoch_mass_flux X variance
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MatvecMul(stoch_mass_flux[d], sqrtLonsager_fc[d]);
        stoch_mass_flux[d].mult(variance);
    }

    // sync the fluxes at the boundaries
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        stoch_mass_flux[d].OverrideSync(geom.periodicity());
        stoch_mass_flux[d].FillBoundary(geom.periodicity());
    }

    // If there are walls with zero-flux boundary conditions
    if (is_nonisothermal == 1) {
        Abort("StochMassFlux: is_nonisothermal==1 not supported yet");
    }

    // correct fluxes to ensure mass conservation to roundoff
    if (correct_flux == 1 && nspecies > 1) {
        CorrectionFlux(rho,rhotot,stoch_mass_flux);
    }

    // compute divergence of stochastic flux
    ComputeDiv(stoch_mass_fluxdiv,stoch_mass_flux,0,0,nspecies,geom,increment);

}

