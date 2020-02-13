#include "rng_functions.H"
#include "gmres_functions.H"
#include "common_functions.H"
#include "hydro_functions_F.H"
#include "StochMassFlux.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>


// initialize n_rngs, geom
// build MultiFabs to hold random numbers
StochMassFlux::StochMassFlux(BoxArray ba_in, DistributionMapping dmap_in, Geometry geom_in,
                             int n_rngs_in) {

    BL_PROFILE_VAR("StochMassFlux::StochMassFlux()",StochMassFlux);

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


void StochMassFlux::weightMassFlux(Vector< amrex::Real > weights) {

    /*
    mflux_cc_weighted.setVal(0.0);
    for (int d=0; d<NUM_EDGE; ++d) {
        mflux_ed_weighted[d].setVal(0.0);
    }

    // add weighted contribution of fluxes
    for (int i=0; i<n_rngs; ++i) {
        MultiFab::Saxpy(mflux_cc_weighted, weights[i], mflux_cc[i], 0, 0, AMREX_SPACEDIM, std::max(1,filtering_width));
        for (int d=0; d<NUM_EDGE; ++d) {
            MultiFab::Saxpy(mflux_ed_weighted[d], weights[i], mflux_ed[i][d], 0, 0, ncomp_ed, filtering_width);
        }
    }
    */
}

// fill stoch_W_fc with random numbers
void StochMassFlux::fillMassStochastic() {

    BL_PROFILE_VAR("StochMassFlux::fillMassStochastic()",StochMassFlux);

    for (int i=0; i<n_rngs; ++i) {
        for (int n=0; n<AMREX_SPACEDIM; ++n) {
            for (int comp=0; comp<nspecies; comp++) {
                MultiFABFillRandom(stoch_W_fc[i][n],comp,1.0,geom);
            }
        }
    }    
}

void StochMassFlux::StochMassFluxBC() {

    // lo-x domain boundary
    if (bc_mass_lo[0] == 1 || bc_mass_lo[0] == 2) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = reservoir   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_lo[0] == 1) ? 0. : sqrt(2.);

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
                AMREX_HOST_DEVICE_FOR_4D(b, nspecies, i, j, k, n,
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
    if (bc_mass_hi[0] == 1 || bc_mass_hi[0] == 2) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = reservoir   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_hi[0] == 1) ? 0. : sqrt(2.);

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
                AMREX_HOST_DEVICE_FOR_4D(b, nspecies, i, j, k, n,
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
    if (bc_mass_lo[1] == 1 || bc_mass_lo[1] == 2) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = reservoir   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_lo[1] == 1) ? 0. : sqrt(2.);

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
                AMREX_HOST_DEVICE_FOR_4D(b, nspecies, i, j, k, n,
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
    if (bc_mass_hi[1] == 1 || bc_mass_hi[1] == 2) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = reservoir   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_hi[1] == 1) ? 0. : sqrt(2.);

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
                AMREX_HOST_DEVICE_FOR_4D(b, nspecies, i, j, k, n,
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
    if (bc_mass_lo[2] == 1 || bc_mass_lo[2] == 2) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = reservoir   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_lo[2] == 1) ? 0. : sqrt(2.);

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
                AMREX_HOST_DEVICE_FOR_4D(b, nspecies, i, j, k, n,
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
    if (bc_mass_hi[2] == 1 || bc_mass_hi[2] == 2) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = reservoir   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_hi[2] == 1) ? 0. : sqrt(2.);

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
                AMREX_HOST_DEVICE_FOR_4D(b, nspecies, i, j, k, n,
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

void StochMassFlux::StochMassFluxDiv(std::array< MultiFab, AMREX_SPACEDIM >& m_force,
                                     const int& increment,
                                     const MultiFab& eta_cc,
                                     const std::array< MultiFab, NUM_EDGE >& eta_ed,
                                     const MultiFab& temp_cc,
                                     const std::array< MultiFab, NUM_EDGE >& temp_ed,
                                     const Vector< amrex::Real >& weights,
                                     const amrex::Real& dt) {

    BL_PROFILE_VAR("StochMassFlux::StochMassFluxDiv()",StochMassFluxDiv);

    // Take linear combination of mflux multifabs at each stage
    StochMassFlux::weightMassFlux(weights);

    // multiply noise stored in mflux_ed_weighted
    // on walls by 0 or on reservoirs by sqrt(2)
    StochMassFlux::StochMassFluxBC();

    // calculate divergence and add to stoch_m_force
    Real dxinv = 1./(geom.CellSize()[0]);
/*
    // Loop over boxes
    for (MFIter mfi(mflux_cc_weighted); mfi.isValid(); ++mfi) {

        const Array4<Real const> & flux_cc = mflux_cc_weighted.array(mfi);
#if (AMREX_SPACEDIM == 2)
        const Array4<Real const> & flux_nd = mflux_ed_weighted[0].array(mfi);
#elif (AMREX_SPACEDIM == 3)
        const Array4<Real const> & flux_xy = mflux_ed_weighted[0].array(mfi);
        const Array4<Real const> & flux_xz = mflux_ed_weighted[1].array(mfi);
        const Array4<Real const> & flux_yz = mflux_ed_weighted[2].array(mfi);
#endif

        AMREX_D_TERM(const Array4<Real> & divx = m_force[0].array(mfi);,
                     const Array4<Real> & divy = m_force[1].array(mfi);,
                     const Array4<Real> & divz = m_force[2].array(mfi););

        AMREX_D_TERM(Box bx_x = mfi.validbox();,
                     Box bx_y = mfi.validbox();,
                     Box bx_z = mfi.validbox(););

        AMREX_D_TERM(bx_x.growHi(0);,
                     bx_y.growHi(1);,
                     bx_z.growHi(2););

#if (AMREX_SPACEDIM == 2)
        if (increment == 1) {
            amrex::ParallelFor(bx_x,bx_y,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   divx(i,j,k) += (flux_cc(i,j,k,0) - flux_cc(i-1,j,k,0) +
                                                   flux_nd(i,j+1,k,0) - flux_nd(i,j,k,0)) * dxinv;
                               },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   divy(i,j,k) += (flux_nd(i+1,j,k,1) - flux_nd(i,j,k,1) +
                                                   flux_cc(i,j,k,1) - flux_cc(i,j-1,k,1)) * dxinv;
                               });
	}
	else if (increment == 0) {
            amrex::ParallelFor(bx_x,bx_y,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   divx(i,j,k) = (flux_cc(i,j,k,0) - flux_cc(i-1,j,k,0) +
                                                  flux_nd(i,j+1,k,0) - flux_nd(i,j,k,0)) * dxinv;
                               },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   divy(i,j,k) = (flux_nd(i+1,j,k,1) - flux_nd(i,j,k,1) +
                                                  flux_cc(i,j,k,1) - flux_cc(i,j-1,k,1)) * dxinv;
                               });
	}
#elif (AMREX_SPACEDIM == 3)
        if (increment == 1) {
            amrex::ParallelFor(bx_x,bx_y,bx_z,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   divx(i,j,k) += (flux_cc(i,j,k,0) - flux_cc(i-1,j,k,0) +
                                                   flux_xy(i,j+1,k,0) - flux_xy(i,j,k,0) +
                                                   flux_xz(i,j,k+1,0) - flux_xz(i,j,k,0)) * dxinv;
                               },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   divy(i,j,k) += (flux_xy(i+1,j,k,1) - flux_xy(i,j,k,1) +
                                                   flux_cc(i,j,k,1) - flux_cc(i,j-1,k,1) +
                                                   flux_yz(i,j,k+1,0) - flux_yz(i,j,k,0)) * dxinv;
                               },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   divz(i,j,k) += (flux_xz(i+1,j,k,1) - flux_xz(i,j,k,1) +
                                                   flux_yz(i,j+1,k,1) - flux_yz(i,j,k,1) +
                                                   flux_cc(i,j,k,2) - flux_cc(i,j,k-1,2)) * dxinv;
                               });
	}
	else if (increment == 0) {
            amrex::ParallelFor(bx_x,bx_y,bx_z,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   divx(i,j,k) = (flux_cc(i,j,k,0) - flux_cc(i-1,j,k,0) +
                                                  flux_xy(i,j+1,k,0) - flux_xy(i,j,k,0) +
                                                  flux_xz(i,j,k+1,0) - flux_xz(i,j,k,0)) * dxinv;
                               },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   divy(i,j,k) = (flux_xy(i+1,j,k,1) - flux_xy(i,j,k,1) +
                                                  flux_cc(i,j,k,1) - flux_cc(i,j-1,k,1) +
                                                  flux_yz(i,j,k+1,0) - flux_yz(i,j,k,0)) * dxinv;
                               },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   divz(i,j,k) = (flux_xz(i+1,j,k,1) - flux_xz(i,j,k,1) +
                                                  flux_yz(i,j+1,k,1) - flux_yz(i,j,k,1) +
                                                  flux_cc(i,j,k,2) - flux_cc(i,j,k-1,2)) * dxinv;
                               });
	}
#endif
    }
*/
    
}

