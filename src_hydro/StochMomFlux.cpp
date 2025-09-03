#include "rng_functions.H"
#include "common_functions.H"
#include "hydro_functions.H"

#include "StochMomFlux.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>

// initialize n_rngs, geom
// build MultiFabs to hold random numbers
StochMomFlux::StochMomFlux(BoxArray ba_in, DistributionMapping dmap_in, Geometry geom_in,
                           int n_rngs_in) {

    BL_PROFILE_VAR("StochMomFlux()",StochMomFlux);

    if (filtering_width != 0) {
        Abort("StochMomFlux: filtering_width != 0 not fully implemented yet");
    }

    // number of random number stages
    n_rngs = n_rngs_in;

    // keep a local geometry object so we won't always have to pass one in
    geom = geom_in;

    // resize these to hold the number of RNG stages
    mflux_cc.resize(n_rngs);
    mflux_ed.resize(n_rngs);
    //filtering_width=1;
    // Here we store all the random number stages at all spatial locations
    for (int i=0; i<n_rngs; ++i) {
        mflux_cc[i].define(ba_in, dmap_in, AMREX_SPACEDIM, amrex::max(1,filtering_width));
        mflux_cc[i].setVal(0.);
#if (AMREX_SPACEDIM == 2)
        mflux_ed[i][0].define(convert(ba_in,nodal_flag), dmap_in, ncomp_ed, filtering_width);
#elif (AMREX_SPACEDIM == 3)
        mflux_ed[i][0].define(convert(ba_in,nodal_flag_xy), dmap_in, ncomp_ed, filtering_width);
        mflux_ed[i][1].define(convert(ba_in,nodal_flag_xz), dmap_in, ncomp_ed, filtering_width);
        mflux_ed[i][2].define(convert(ba_in,nodal_flag_yz), dmap_in, ncomp_ed, filtering_width);
#endif
        for (int d=0; d<NUM_EDGE; ++d) {
            mflux_ed[i][d].setVal(0.);
        }
    }

    // Temporary storage for linear combinations of random number stages
    mflux_cc_weighted.define(ba_in, dmap_in, AMREX_SPACEDIM, amrex::max(1,filtering_width));
    mflux_cc_weighted.setVal(0.);
#if (AMREX_SPACEDIM == 2)
    mflux_ed_weighted[0].define(convert(ba_in,nodal_flag), dmap_in, ncomp_ed, filtering_width);
#elif (AMREX_SPACEDIM == 3)
    mflux_ed_weighted[0].define(convert(ba_in,nodal_flag_xy), dmap_in, ncomp_ed, filtering_width);
    mflux_ed_weighted[1].define(convert(ba_in,nodal_flag_xz), dmap_in, ncomp_ed, filtering_width);
    mflux_ed_weighted[2].define(convert(ba_in,nodal_flag_yz), dmap_in, ncomp_ed, filtering_width);
    for (int d=0; d<NUM_EDGE; ++d) {
        mflux_ed_weighted[d].setVal(0.);
    }
    //filtering_width=0;

#endif
}

// fill mflux_cc and mflux_ed with random numbers
// we need these to hold Z = (W + W^T) / sqrt(2)
// For the symmetric case, the diagonal terms (stored at cell-centers)
// look like (W_11 + W_11) / sqrt(2)
// W_11 + W_11 has variance 4, and can be obtained by generating an RNG
// and multiplying by 2.  Then you have to divide by sqrt(2), so the net
// effect is an RNG multiplied by sqrt(2).  So we fill Z_11 with an
// RNG with variance 2.
// The off-diagnoal terms look like (W_12 + W_21) / sqrt(2)
// W_12 + W_21 has variance 2, and can be obtained by generating an RNG
// and multiplfying by sqrt(2).  Then you have to divide by sqrt(2), so the net
// effect is an RNG multiplied by 1.  So we will Z_12 with a unit-variance RNG
// Z is symmetric so we only store the lower-diagonal terms
void StochMomFlux::fillMomStochastic() {

    BL_PROFILE_VAR("fillMomStochastic()",StochMomFlux);

    for (int i=0; i<n_rngs; ++i) {

        switch(stoch_stress_form) {

        case 0: // Non-symmetric
            for (int n=0; n<AMREX_SPACEDIM; ++n) {
                MultiFabFillRandom(mflux_cc[i],n,1.0,geom);
            }

            for (int d=0; d<NUM_EDGE; ++d) {
                for (int n=0; n<ncomp_ed; ++n) {
                    MultiFabFillRandom(mflux_ed[i][d],n,1.0,geom);
                }
            }
            break;

        default: // Symmetric
            for (int n=0; n<AMREX_SPACEDIM; ++n) {
                MultiFabFillRandom(mflux_cc[i],n,2.0,geom);
            }

            for (int d=0; d<NUM_EDGE; ++d) {
                MultiFabFillRandom(mflux_ed[i][d],0,1.0,geom);
                MultiFab::Copy(mflux_ed[i][d], mflux_ed[i][d], 0, 1, ncomp_ed-1, 0);
            }
            break;
        }
    }
}


// create weighted sum of stage RNGs and store in mflux_cc_weighted and mflux_ed_weighted
void StochMomFlux::weightMomflux(Vector< amrex::Real > weights) {

    BL_PROFILE_VAR("weightMomFlux()",weightMomFlux);

    mflux_cc_weighted.setVal(0.0);
    for (int d=0; d<NUM_EDGE; ++d) {
        mflux_ed_weighted[d].setVal(0.0);
    }

    for (int i=0; i<n_rngs; ++i) {
        MultiFab::Saxpy(mflux_cc_weighted, weights[i], mflux_cc[i], 0, 0, AMREX_SPACEDIM, amrex::max(1,filtering_width));
        for (int d=0; d<NUM_EDGE; ++d) {
            MultiFab::Saxpy(mflux_ed_weighted[d], weights[i], mflux_ed[i][d], 0, 0, ncomp_ed, filtering_width);
        }
    }
}

// scale random numbers that lie on physical boundaries appropriately
void StochMomFlux::MomFluxBC() {

    BL_PROFILE_VAR("MomFluxBC()",MomFluxBC);

#if (AMREX_SPACEDIM == 2)

    // lo-x domain boundary
    if (bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) {

        // 1 = slip wall   : multiply fluxes on wall by 0
        // 2 = no-slip wall: multiply fluxes on wall by sqrt(2)
        Real factor = (bc_vel_lo[0] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on mflux_ed_weighted[0] nodality (xy)
        const Box& dom_nd = amrex::convert(geom.Domain(), mflux_ed_weighted[0].ixType());

        // this is the x-lo domain boundary box (xy nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_nd_xlo = amrex::bdryNode(dom_nd, Orientation(0, Orientation::low));

        for (MFIter mfi(mflux_ed_weighted[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_nd_xlo;
            Array4<Real> const& mflux_nd = (mflux_ed_weighted[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    mflux_nd(i,j,k,n) *= factor;
                });
            }
        }
    }
    else if (bc_vel_lo[0] != -1) {
        Abort("MomFluxBC unsupported bc type");
    }

    // hi-x domain boundary
    if (bc_vel_hi[0] == 1 || bc_vel_hi[0] == 2) {

        // 1 = slip wall   : multiply fluxes on wall by 0
        // 2 = no-slip wall: multiply fluxes on wall by sqrt(2)
        Real factor = (bc_vel_hi[0] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on mflux_ed_weighted[0] nodality (xy)
        const Box& dom_nd = amrex::convert(geom.Domain(), mflux_ed_weighted[0].ixType());

        // this is the x-hi domain boundary box (xy nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_nd_xhi = amrex::bdryNode(dom_nd, Orientation(0, Orientation::high));

        for (MFIter mfi(mflux_ed_weighted[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_nd_xhi;
            Array4<Real> const& mflux_nd = (mflux_ed_weighted[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    mflux_nd(i,j,k,n) *= factor;
                });
            }
        }
    }
    else if (bc_vel_hi[0] != -1) {
        Abort("MomFluxBC unsupported bc type");
    }

    // lo-y domain boundary
    if (bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) {

        // 1 = slip wall   : multiply fluxes on wall by 0
        // 2 = no-slip wall: multiply fluxes on wall by sqrt(2)
        Real factor = (bc_vel_lo[1] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on mflux_ed_weighted[0] nodality (xy)
        const Box& dom_nd = amrex::convert(geom.Domain(), mflux_ed_weighted[0].ixType());

        // this is the y-lo domain boundary box (xy nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_nd_ylo = amrex::bdryNode(dom_nd, Orientation(1, Orientation::low));

        for (MFIter mfi(mflux_ed_weighted[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_nd_ylo;
            Array4<Real> const& mflux_nd = (mflux_ed_weighted[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    mflux_nd(i,j,k,n) *= factor;
                });
            }
        }
    }
    else if (bc_vel_lo[1] != -1) {
        Abort("MomFluxBC unsupported bc type");
    }

    // hi-y domain boundary
    if (bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) {

        // 1 = slip wall   : multiply fluxes on wall by 0
        // 2 = no-slip wall: multiply fluxes on wall by sqrt(2)
        Real factor = (bc_vel_hi[1] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on mflux_ed_weighted[0] nodality (xy)
        const Box& dom_nd = amrex::convert(geom.Domain(), mflux_ed_weighted[0].ixType());

        // this is the y-hi domain boundary box (xy nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_nd_yhi = amrex::bdryNode(dom_nd, Orientation(1, Orientation::high));

        for (MFIter mfi(mflux_ed_weighted[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_nd_yhi;
            Array4<Real> const& mflux_nd = (mflux_ed_weighted[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    mflux_nd(i,j,k,n) *= factor;
                });
            }
        }
    }
    else if (bc_vel_hi[1] != -1) {
        Abort("MomFluxBC unsupported bc type");
    }

#elif (AMREX_SPACEDIM == 3)

    // lo-x domain boundary, y-facing fluxes
    // lo-x domain boundary, z-facing fluxes
    if (bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) {

        // 1 = slip wall   : multiply fluxes on wall by 0
        // 2 = no-slip wall: multiply fluxes on wall by sqrt(2)
        Real factor = (bc_vel_lo[0] == 1) ? 0. : sqrt(2.);

        ////////////////////////////////////////////////

        // domain grown nodally based on mflux_ed_weighted[0] nodality (xy)
        const Box& dom_xy = amrex::convert(geom.Domain(), mflux_ed_weighted[0].ixType());

        // this is the x-lo domain boundary box (xy nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xy_xlo = amrex::bdryNode(dom_xy, Orientation(0, Orientation::low));

        for (MFIter mfi(mflux_ed_weighted[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xy_xlo;
            Array4<Real> const& mflux_xy = (mflux_ed_weighted[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    mflux_xy(i,j,k,n) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////

        // domain grown nodally based on mflux_ed_weighted[1] nodality (xz)
        const Box& dom_xz = amrex::convert(geom.Domain(), mflux_ed_weighted[1].ixType());

        // this is the x-lo domain boundary box (xz nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xz_xlo = amrex::bdryNode(dom_xz, Orientation(0, Orientation::low));

        for (MFIter mfi(mflux_ed_weighted[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xz_xlo;
            Array4<Real> const& mflux_xz = (mflux_ed_weighted[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    mflux_xz(i,j,k,n) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
    }
    else if (bc_vel_lo[0] != -1) {
        Abort("MomFluxBC unsupported bc type");
    }

    // hi-x domain boundary, y-facing fluxes
    // hi-x domain boundary, z-facing fluxes
    if (bc_vel_hi[0] == 1 || bc_vel_hi[0] == 2) {

        // 1 = slip wall   : multiply fluxes on wall by 0
        // 2 = no-slip wall: multiply fluxes on wall by sqrt(2)
        Real factor = (bc_vel_hi[0] == 1) ? 0. : sqrt(2.);

        ////////////////////////////////////////////////

        // domain grown nodally based on mflux_ed_weighted[0] nodality (xy)
        const Box& dom_xy = amrex::convert(geom.Domain(), mflux_ed_weighted[0].ixType());

        // this is the x-hi domain boundary box (xy nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xy_xhi = amrex::bdryNode(dom_xy, Orientation(0, Orientation::high));

        for (MFIter mfi(mflux_ed_weighted[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xy_xhi;
            Array4<Real> const& mflux_xy = (mflux_ed_weighted[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    mflux_xy(i,j,k,n) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////

        // domain grown nodally based on mflux_ed_weighted[1] nodality (xz)
        const Box& dom_xz = amrex::convert(geom.Domain(), mflux_ed_weighted[1].ixType());

        // this is the x-hi domain boundary box (xz nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xz_xhi = amrex::bdryNode(dom_xz, Orientation(0, Orientation::high));

        for (MFIter mfi(mflux_ed_weighted[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xz_xhi;
            Array4<Real> const& mflux_xz = (mflux_ed_weighted[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    mflux_xz(i,j,k,n) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
    }
    else if (bc_vel_hi[0] != -1) {
        Abort("MomFluxBC unsupported bc type");
    }

    // lo-y domain boundary, x-facing fluxes
    // lo-y domain boundary, z-facing fluxes
    if (bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) {

        // 1 = slip wall   : multiply fluxes on wall by 0
        // 2 = no-slip wall: multiply fluxes on wall by sqrt(2)
        Real factor = (bc_vel_lo[1] == 1) ? 0. : sqrt(2.);

        ////////////////////////////////////////////////

        // domain grown nodally based on mflux_ed_weighted[0] nodality (xy)
        const Box& dom_xy = amrex::convert(geom.Domain(), mflux_ed_weighted[0].ixType());

        // this is the y-lo domain boundary box (xy nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xy_ylo = amrex::bdryNode(dom_xy, Orientation(1, Orientation::low));

        for (MFIter mfi(mflux_ed_weighted[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xy_ylo;
            Array4<Real> const& mflux_xy = (mflux_ed_weighted[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    mflux_xy(i,j,k,n) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////

        // domain grown nodally based on mflux_ed_weighted[2] nodality (yz)
        const Box& dom_yz = amrex::convert(geom.Domain(), mflux_ed_weighted[2].ixType());

        // this is the y-lo domain boundary box (yz nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yz_ylo = amrex::bdryNode(dom_yz, Orientation(1, Orientation::low));

        for (MFIter mfi(mflux_ed_weighted[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yz_ylo;
            Array4<Real> const& mflux_yz = (mflux_ed_weighted[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    mflux_yz(i,j,k,n) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
    }
    else if (bc_vel_lo[1] != -1) {
        Abort("MomFluxBC unsupported bc type");
    }

    // hi-y domain boundary, x-facing fluxes
    // hi-y domain boundary, z-facing fluxes
    if (bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) {

        // 1 = slip wall   : multiply fluxes on wall by 0
        // 2 = no-slip wall: multiply fluxes on wall by sqrt(2)
        Real factor = (bc_vel_hi[1] == 1) ? 0. : sqrt(2.);

        ////////////////////////////////////////////////

        // domain grown nodally based on mflux_ed_weighted[0] nodality (xy)
        const Box& dom_xy = amrex::convert(geom.Domain(), mflux_ed_weighted[0].ixType());

        // this is the y-hi domain boundary box (xy nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xy_yhi = amrex::bdryNode(dom_xy, Orientation(1, Orientation::high));

        for (MFIter mfi(mflux_ed_weighted[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xy_yhi;
            Array4<Real> const& mflux_xy = (mflux_ed_weighted[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    mflux_xy(i,j,k,n) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////

        // domain grown nodally based on mflux_ed_weighted[2] nodality (yz)
        const Box& dom_yz = amrex::convert(geom.Domain(), mflux_ed_weighted[2].ixType());

        // this is the y-hi domain boundary box (yz nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yz_yhi = amrex::bdryNode(dom_yz, Orientation(1, Orientation::high));

        for (MFIter mfi(mflux_ed_weighted[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yz_yhi;
            Array4<Real> const& mflux_yz = (mflux_ed_weighted[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    mflux_yz(i,j,k,n) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
    }
    else if (bc_vel_hi[1] != -1) {
        Abort("MomFluxBC unsupported bc type");
    }

    // lo-z domain boundary, x-facing fluxes
    // lo-z domain boundary, y-facing fluxes
    if (bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) {

        // 1 = slip wall   : multiply fluxes on wall by 0
        // 2 = no-slip wall: multiply fluxes on wall by sqrt(2)
        Real factor = (bc_vel_lo[2] == 1) ? 0. : sqrt(2.);

        ////////////////////////////////////////////////

        // domain grown nodally based on mflux_ed_weighted[1] nodality (xz)
        const Box& dom_xz = amrex::convert(geom.Domain(), mflux_ed_weighted[1].ixType());

        // this is the z-lo domain boundary box (xz nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xz_zlo = amrex::bdryNode(dom_xz, Orientation(2, Orientation::low));

        for (MFIter mfi(mflux_ed_weighted[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xz_zlo;
            Array4<Real> const& mflux_xz = (mflux_ed_weighted[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    mflux_xz(i,j,k,n) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////

        // domain grown nodally based on mflux_ed_weighted[2] nodality (yz)
        const Box& dom_yz = amrex::convert(geom.Domain(), mflux_ed_weighted[2].ixType());

        // this is the y-lo domain boundary box (yz nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yz_zlo = amrex::bdryNode(dom_yz, Orientation(2, Orientation::low));

        for (MFIter mfi(mflux_ed_weighted[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yz_zlo;
            Array4<Real> const& mflux_yz = (mflux_ed_weighted[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    mflux_yz(i,j,k,n) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
    }
    else if (bc_vel_lo[2] != -1) {
        Abort("MomFluxBC unsupported bc type");
    }

    // hi-z domain boundary, x-facing fluxes
    // hi-z domain boundary, y-facing fluxes
    if (bc_vel_hi[2] == 1 || bc_vel_hi[2] == 2) {

        // 1 = slip wall   : multiply fluxes on wall by 0
        // 2 = no-slip wall: multiply fluxes on wall by sqrt(2)
        Real factor = (bc_vel_hi[2] == 1) ? 0. : sqrt(2.);

        ////////////////////////////////////////////////

        // domain grown nodally based on mflux_ed_weighted[1] nodality (xz)
        const Box& dom_xz = amrex::convert(geom.Domain(), mflux_ed_weighted[1].ixType());

        // this is the z-hi domain boundary box (xz nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xz_zhi = amrex::bdryNode(dom_xz, Orientation(2, Orientation::high));

        for (MFIter mfi(mflux_ed_weighted[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xz_zhi;
            Array4<Real> const& mflux_xz = (mflux_ed_weighted[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    mflux_xz(i,j,k,n) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////

        // domain grown nodally based on mflux_ed_weighted[2] nodality (yz)
        const Box& dom_yz = amrex::convert(geom.Domain(), mflux_ed_weighted[2].ixType());

        // this is the y-hi domain boundary box (yz nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yz_zhi = amrex::bdryNode(dom_yz, Orientation(2, Orientation::high));

        for (MFIter mfi(mflux_ed_weighted[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yz_zhi;
            Array4<Real> const& mflux_yz = (mflux_ed_weighted[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    mflux_yz(i,j,k,n) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
    }
    else if (bc_vel_hi[2] != -1) {
        Abort("MomFluxBC unsupported bc type");
    }

#endif

}

// Multiply mflux_weighted by sqrt(eta*temperature)
void StochMomFlux::multbyVarSqrtEtaTemp(const MultiFab& eta_cc,
                                        const std::array< MultiFab, NUM_EDGE >& eta_ed,
                                        const MultiFab& temp_cc,
                                        const std::array< MultiFab, NUM_EDGE >& temp_ed,
                                        const amrex::Real& dt) {

    BL_PROFILE_VAR("multbyVarSqrtEtaTemp()",multbyVarSqrtEtaTemp);

    const Real* dx = geom.CellSize();
    Real dVol = (AMREX_SPACEDIM==2) ? dx[0]*dx[1]*cell_depth : dx[0]*dx[1]*dx[2];

    // Compute variance using computed differential volume
    Real variance = sqrt(variance_coef_mom*2.0*k_B/(dVol*dt));

    // Scale mflux_weighted by variance
    mflux_cc_weighted.mult(variance, filtering_width);
    for (int d=0; d<NUM_EDGE; d++) {
        mflux_ed_weighted[d].mult(variance, filtering_width);
    }

    // Loop over boxes
    for (MFIter mfi(mflux_cc_weighted); mfi.isValid(); ++mfi) {
        // Note: Make sure that multifab is cell-centered
        const Box& bx = mfi.growntilebox(1);

        const Array4<Real> & mflux_cc_fab = mflux_cc_weighted.array(mfi);
        const Array4<Real const> & eta_cc_fab = eta_cc.array(mfi);
        const Array4<Real const> &  temp_cc_fab = temp_cc.array(mfi);

        const Box & bx_xy = mfi.tilebox(nodal_flag_xy);
#if (AMREX_SPACEDIM == 3)
        const Box & bx_xz = mfi.tilebox(nodal_flag_xz);
        const Box & bx_yz = mfi.tilebox(nodal_flag_yz);
#endif

        const Array4<Real> & mflux_xy_fab = mflux_ed_weighted[0].array(mfi);
        const Array4<Real const> & eta_xy_fab = eta_ed[0].array(mfi);
        const Array4<Real const> & temp_xy_fab = temp_ed[0].array(mfi);

#if (AMREX_SPACEDIM == 3)
        const Array4<Real> & mflux_xz_fab = mflux_ed_weighted[1].array(mfi);
        const Array4<Real> & mflux_yz_fab = mflux_ed_weighted[2].array(mfi);
        const Array4<Real const> & eta_xz_fab = eta_ed[1].array(mfi);
        const Array4<Real const> & eta_yz_fab = eta_ed[2].array(mfi);
        const Array4<Real const> & temp_xz_fab = temp_ed[1].array(mfi);
        const Array4<Real const> & temp_yz_fab = temp_ed[2].array(mfi);
#endif

        amrex::ParallelFor(bx, AMREX_SPACEDIM, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            mflux_cc_fab(i,j,k,n) *= sqrt(eta_cc_fab(i,j,k)*temp_cc_fab(i,j,k));
        });

        amrex::ParallelFor(bx_xy, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            mflux_xy_fab(i,j,k,n) *= sqrt(eta_xy_fab(i,j,k)*temp_xy_fab(i,j,k));

        });

#if (AMREX_SPACEDIM == 3)

        amrex::ParallelFor(bx_xz, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            mflux_xz_fab(i,j,k,n) *= sqrt(eta_xz_fab(i,j,k)*temp_xz_fab(i,j,k));
        });

        amrex::ParallelFor(bx_yz, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            mflux_yz_fab(i,j,k,n) *= sqrt(eta_yz_fab(i,j,k)*temp_yz_fab(i,j,k));
        });

#endif
    }
}

// compute stochastic momentum flux divergence
void StochMomFlux::StochMomFluxDiv(std::array< MultiFab, AMREX_SPACEDIM >& m_force,
                                   const int& increment,
                                   const MultiFab& eta_cc,
                                   const std::array< MultiFab, NUM_EDGE >& eta_ed,
                                   const MultiFab& temp_cc,
                                   const std::array< MultiFab, NUM_EDGE >& temp_ed,
                                   const Vector< amrex::Real >& weights,
                                   const amrex::Real& dt) {

    BL_PROFILE_VAR("StochMomFluxDiv()",StochMomFluxDiv);

    // Take linear combination of mflux multifabs at each stage
    StochMomFlux::weightMomflux(weights);

    // Multiply weighted mflux (cc & edge) by sqrt(eta*temperature)
    StochMomFlux::multbyVarSqrtEtaTemp(eta_cc,eta_ed,temp_cc,temp_ed,dt);

    // multiply noise stored in mflux_ed_weighted
    // on walls by 0 (for slip) or sqrt(2) (for no-slip)
    StochMomFlux::MomFluxBC();

    // sync up random numbers at boundaries and ghost cells
    for (int d=0; d<NUM_EDGE; ++d) {
        mflux_ed_weighted[d].OverrideSync(geom.periodicity());
        mflux_ed_weighted[d].FillBoundary(geom.periodicity());
    }
    mflux_cc_weighted.FillBoundary(geom.periodicity());

    if (filtering_width > 0) {
        Abort("StochMomFlux: filtering_width != 0 not fully implemented yet");
        // need calls to filter_stoch_m_flux for mflux_ed and mflux_cc
        /*

         */
        mflux_cc_weighted.FillBoundary(geom.periodicity());
    }

    // calculate divergence and add to stoch_m_force
    Real dxinv = 1./(geom.CellSize()[0]);

    // if not incrementing, initialize data to zero
    if (increment == 0) {
        for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
            m_force[dir].setVal(0.,0,1,0);
        }
    }

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
        amrex::ParallelFor(bx_x,bx_y, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            divx(i,j,k) += (flux_cc(i,j,k,0) - flux_cc(i-1,j,k,0) +
                            flux_nd(i,j+1,k,0) - flux_nd(i,j,k,0)) * dxinv;
        },
                                      [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            divy(i,j,k) += (flux_nd(i+1,j,k,1) - flux_nd(i,j,k,1) +
                            flux_cc(i,j,k,1) - flux_cc(i,j-1,k,1)) * dxinv;
        });

#elif (AMREX_SPACEDIM == 3)
        amrex::ParallelFor(bx_x,bx_y,bx_z, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            divx(i,j,k) += (flux_cc(i,j,k,0) - flux_cc(i-1,j,k,0) +
                            flux_xy(i,j+1,k,0) - flux_xy(i,j,k,0) +
                            flux_xz(i,j,k+1,0) - flux_xz(i,j,k,0)) * dxinv;
        },
                                           [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            divy(i,j,k) += (flux_xy(i+1,j,k,1) - flux_xy(i,j,k,1) +
                            flux_cc(i,j,k,1) - flux_cc(i,j-1,k,1) +
                            flux_yz(i,j,k+1,0) - flux_yz(i,j,k,0)) * dxinv;
        },
                                           [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            divz(i,j,k) += (flux_xz(i+1,j,k,1) - flux_xz(i,j,k,1) +
                            flux_yz(i,j+1,k,1) - flux_yz(i,j,k,1) +
                            flux_cc(i,j,k,2) - flux_cc(i,j,k-1,2)) * dxinv;
        });
#endif
    }

    // m_force does not have ghost cells
    // set the value on physical boundaries to zero
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFabPhysBCDomainVel(m_force[d], geom, d);
    }
}

// compute stochastic momentum flux divergence
void StochMomFlux::StochMomFluxDivWideSplit(std::array< MultiFab, AMREX_SPACEDIM >& m_force,
                                   const int& increment,
                                   const MultiFab& eta_cc,
                                   const std::array< MultiFab, NUM_EDGE >& eta_ed,
                                   const MultiFab& temp_cc,
                                   const std::array< MultiFab, NUM_EDGE >& temp_ed,
                                   const Vector< amrex::Real >& weights,
                                   const amrex::Real& dt) {

    BL_PROFILE_VAR("StochMomFluxDiv()",StochMomFluxDiv);

    // Take linear combination of mflux multifabs at each stage
    StochMomFlux::weightMomflux(weights);

    // Multiply weighted mflux (cc & edge) by sqrt(eta*temperature)
    StochMomFlux::multbyVarSqrtEtaTemp(eta_cc,eta_ed,temp_cc,temp_ed,dt);

    // multiply noise stored in mflux_ed_weighted
    // on walls by 0 (for slip) or sqrt(2) (for no-slip)
    StochMomFlux::MomFluxBC();

    // sync up random numbers at boundaries and ghost cells
    for (int d=0; d<NUM_EDGE; ++d) {
        mflux_ed_weighted[d].OverrideSync(geom.periodicity());
        mflux_ed_weighted[d].FillBoundary(geom.periodicity());
    }
    mflux_cc_weighted.FillBoundary(geom.periodicity());

    if (filtering_width > 0) {
        Abort("StochMomFlux: filtering_width != 0 not fully implemented yet");
        // need calls to filter_stoch_m_flux for mflux_ed and mflux_cc
        /*

         */
        mflux_cc_weighted.FillBoundary(geom.periodicity());
    }

    // calculate divergence and add to stoch_m_force
    Real dxinv = 1./(geom.CellSize()[0]);

    int splitCell = 128;

    // if not incrementing, initialize data to zero
    if (increment == 0) {
        for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
            m_force[dir].setVal(0.,0,1,0);
        }
    }

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
        amrex::ParallelFor(bx_x,bx_y, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            divx(i,j,k) += (flux_cc(i,j,k,0) - flux_cc(i-1,j,k,0) +
                            flux_nd(i,j+1,k,0) - flux_nd(i,j,k,0)) * dxinv;
        },
                                      [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            divy(i,j,k) += (flux_nd(i+1,j,k,1) - flux_nd(i,j,k,1) +
                            flux_cc(i,j,k,1) - flux_cc(i,j-1,k,1)) * dxinv;
        });

#elif (AMREX_SPACEDIM == 3)
        amrex::ParallelFor(bx_x,bx_y,bx_z, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if(i!=splitCell  && j!=splitCell && k!=splitCell )
            //if(true)
            {
                divx(i,j,k) += (flux_cc(i,j,k,0) - flux_cc(i-1,j,k,0) +
                                flux_xy(i,j+1,k,0) - flux_xy(i,j,k,0) +
                                flux_xz(i,j,k+1,0) - flux_xz(i,j,k,0)) * dxinv;
            }else
            {
                divx(i,j,k) += (flux_cc(i+1,j,k,0) - flux_cc(i-2,j,k,0) +
                            flux_xy(i,j+2,k,0) - flux_xy(i,j-1,k,0) +
                            flux_xz(i,j,k+2,0) - flux_xz(i,j,k-1,0)) * dxinv * 0.5;
            }
        },
                                            [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if(true)
            {
                divy(i,j,k) += (flux_xy(i+1,j,k,1) - flux_xy(i,j,k,1) +
                                flux_cc(i,j,k,1) - flux_cc(i,j-1,k,1) +
                                flux_yz(i,j,k+1,0) - flux_yz(i,j,k,0)) * dxinv;
            }else
            {
                divy(i,j,k) += (flux_xy(i+2,j,k,1) - flux_xy(i-1,j,k,1) +
                            flux_cc(i,j+1,k,1) - flux_cc(i,j-2,k,1) +
                            flux_yz(i,j,k+2,0) - flux_yz(i,j,k-1,0)) * dxinv * 0.5;
            }
        },
                                           [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if(true )
            {
                divz(i,j,k) += (flux_xz(i+1,j,k,1) - flux_xz(i,j,k,1) +
                                flux_yz(i,j+1,k,1) - flux_yz(i,j,k,1) +
                                flux_cc(i,j,k,2) - flux_cc(i,j,k-1,2)) * dxinv;
            }else
            {
                divz(i,j,k) += (flux_xz(i+2,j,k,1) - flux_xz(i-1,j,k,1) +
                            flux_yz(i,j+2,k,1) - flux_yz(i,j-1,k,1) +
                            flux_cc(i,j,k+1,2) - flux_cc(i,j,k-2,2)) * dxinv * 0.5;

            }
        });
#endif
    }

    // m_force does not have ghost cells
    // set the value on physical boundaries to zero
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFabPhysBCDomainVel(m_force[d], geom, d);
    }
}

// compute stochastic momentum flux divergence
void StochMomFlux::StochMomFluxDivWide(std::array< MultiFab, AMREX_SPACEDIM >& m_force,
                                   const int& increment,
                                   const MultiFab& eta_cc,
                                   const std::array< MultiFab, NUM_EDGE >& eta_ed,
                                   const MultiFab& temp_cc,
                                   const std::array< MultiFab, NUM_EDGE >& temp_ed,
                                   const Vector< amrex::Real >& weights,
                                   const amrex::Real& dt) {

    BL_PROFILE_VAR("StochMomFluxDiv()",StochMomFluxDiv);

    // Take linear combination of mflux multifabs at each stage
    StochMomFlux::weightMomflux(weights);

    // Multiply weighted mflux (cc & edge) by sqrt(eta*temperature)
    StochMomFlux::multbyVarSqrtEtaTemp(eta_cc,eta_ed,temp_cc,temp_ed,dt);

    // multiply noise stored in mflux_ed_weighted
    // on walls by 0 (for slip) or sqrt(2) (for no-slip)
    StochMomFlux::MomFluxBC();

    // sync up random numbers at boundaries and ghost cells
    for (int d=0; d<NUM_EDGE; ++d) {
        mflux_ed_weighted[d].OverrideSync(geom.periodicity());
        mflux_ed_weighted[d].FillBoundary(geom.periodicity());
    }
    mflux_cc_weighted.FillBoundary(geom.periodicity());

    if (filtering_width > 0) {
        Abort("StochMomFlux: filtering_width != 0 not fully implemented yet");
        // need calls to filter_stoch_m_flux for mflux_ed and mflux_cc
        /*

         */
        mflux_cc_weighted.FillBoundary(geom.periodicity());
    }

    // calculate divergence and add to stoch_m_force
    Real dxinv = 1./(geom.CellSize()[0]);

    // if not incrementing, initialize data to zero
    if (increment == 0) {
        for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
            m_force[dir].setVal(0.,0,1,0);
        }
    }

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
        amrex::ParallelFor(bx_x,bx_y, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            divx(i,j,k) += (flux_cc(i,j,k,0) - flux_cc(i-1,j,k,0) +
                            flux_nd(i,j+1,k,0) - flux_nd(i,j,k,0)) * dxinv;
        },
                                      [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            divy(i,j,k) += (flux_nd(i+1,j,k,1) - flux_nd(i,j,k,1) +
                            flux_cc(i,j,k,1) - flux_cc(i,j-1,k,1)) * dxinv;
        });

#elif (AMREX_SPACEDIM == 3)
        amrex::ParallelFor(bx_x,bx_y,bx_z, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            divx(i,j,k) += (flux_cc(i+1,j,k,0) - flux_cc(i-2,j,k,0) +
                            flux_xy(i,j+2,k,0) - flux_xy(i,j-1,k,0) +
                            flux_xz(i,j,k+2,0) - flux_xz(i,j,k-1,0)) * dxinv * 0.5;
        },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k)

        {
            divy(i,j,k) += (flux_xy(i+2,j,k,1) - flux_xy(i-1,j,k,1) +
                            flux_cc(i,j+1,k,1) - flux_cc(i,j-2,k,1) +
                            flux_yz(i,j,k+2,0) - flux_yz(i,j,k-1,0)) * dxinv * 0.5;
        },
                                           [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            divz(i,j,k) += (flux_xz(i+2,j,k,1) - flux_xz(i-1,j,k,1) +
                            flux_yz(i,j+2,k,1) - flux_yz(i,j-1,k,1) +
                            flux_cc(i,j,k+1,2) - flux_cc(i,j,k-2,2)) * dxinv * 0.5;
        });
#endif
    }

    // m_force does not have ghost cells
    // set the value on physical boundaries to zero
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFabPhysBCDomainVel(m_force[d], geom, d);
    }
}

// compute stochastic momentum flux divergence
void StochMomFlux::StochMomFluxDivOrder3(std::array< MultiFab, AMREX_SPACEDIM >& m_force,
                                   const int& increment,
                                   const MultiFab& eta_cc,
                                   const std::array< MultiFab, NUM_EDGE >& eta_ed,
                                   const MultiFab& temp_cc,
                                   const std::array< MultiFab, NUM_EDGE >& temp_ed,
                                   const Vector< amrex::Real >& weights,
                                   const amrex::Real& dt) {

    BL_PROFILE_VAR("StochMomFluxDiv()",StochMomFluxDiv);

    // Take linear combination of mflux multifabs at each stage
    StochMomFlux::weightMomflux(weights);

    // Multiply weighted mflux (cc & edge) by sqrt(eta*temperature)
    StochMomFlux::multbyVarSqrtEtaTemp(eta_cc,eta_ed,temp_cc,temp_ed,dt);

    // multiply noise stored in mflux_ed_weighted
    // on walls by 0 (for slip) or sqrt(2) (for no-slip)
    StochMomFlux::MomFluxBC();

    // sync up random numbers at boundaries and ghost cells
    for (int d=0; d<NUM_EDGE; ++d) {
        mflux_ed_weighted[d].OverrideSync(geom.periodicity());
        mflux_ed_weighted[d].FillBoundary(geom.periodicity());
    }
    mflux_cc_weighted.FillBoundary(geom.periodicity());

    if (filtering_width > 0) {
        Abort("StochMomFlux: filtering_width != 0 not fully implemented yet");
        // need calls to filter_stoch_m_flux for mflux_ed and mflux_cc
        /*

         */
        mflux_cc_weighted.FillBoundary(geom.periodicity());
    }

    // calculate divergence and add to stoch_m_force
    Real dxinv = 1./(geom.CellSize()[0]);

    // if not incrementing, initialize data to zero
    if (increment == 0) {
        for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
            m_force[dir].setVal(0.,0,1,0);
        }
    }

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

        AMREX_D_TERM(const auto xlo = lbound(bx_x);,
                 const auto ylo = lbound(bx_y);,
                 const auto zlo = lbound(bx_z););

        AMREX_D_TERM(const auto xhi = ubound(bx_x);,
                 const auto yhi = ubound(bx_y);,
                 const auto zhi = ubound(bx_z););

        Real nineOver8 = 9.0/8.0;
        Real oneOver24 = 1.0/24.0;

        Real preFac = sqrt(1.0/(nineOver8*nineOver8 + oneOver24*oneOver24));
        preFac = 1;


#if (AMREX_SPACEDIM == 2)
        amrex::ParallelFor(bx_x,bx_y, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            divx(i,j,k) += (flux_cc(i,j,k,0) - flux_cc(i-1,j,k,0) +
                            flux_nd(i,j+1,k,0) - flux_nd(i,j,k,0)) * dxinv;
        },
                                      [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            divy(i,j,k) += (flux_nd(i+1,j,k,1) - flux_nd(i,j,k,1) +
                            flux_cc(i,j,k,1) - flux_cc(i,j-1,k,1)) * dxinv;
        });

#elif (AMREX_SPACEDIM == 3)
        amrex::ParallelFor(bx_x,bx_y,bx_z, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
//            if(i==xlo.x || i==xhi.x || j==xlo.y || j==xhi.y || k==xlo.z || k==xhi.z)
//            {
//                divx(i,j,k) += (flux_cc(i,j,k,0) - flux_cc(i-1,j,k,0) +
//                                flux_xy(i,j+1,k,0) - flux_xy(i,j,k,0) +
//                                flux_xz(i,j,k+1,0) - flux_xz(i,j,k,0)) * dxinv;
//            }else
//            {
                divx(i,j,k) += preFac*(nineOver8*flux_cc(i,j,k,0) - nineOver8*flux_cc(i-1,j,k,0) - oneOver24*flux_cc(i+1,j,k,0) + oneOver24*flux_cc(i-2,j,k,0) +
                     nineOver8*flux_xy(i,j+1,k,0) -  nineOver8*flux_xy(i,j,k,0) - oneOver24*flux_xy(i,j+2,k,0) + oneOver24*flux_xy(i,j-1,k,0) +
                    nineOver8*flux_xz(i,j,k+1,0) - nineOver8*flux_xz(i,j,k,0) - oneOver24*flux_xz(i,j,k+2,0) + oneOver24*flux_xz(i,j,k-1,0)) * dxinv;

//            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
//            if(i==ylo.x || i==yhi.x || j==ylo.y || j==yhi.y || k==ylo.z || k==yhi.z)
//            {
//                divy(i,j,k) += (flux_xy(i+1,j,k,1) - flux_xy(i,j,k,1) +
//                                flux_cc(i,j,k,1) - flux_cc(i,j-1,k,1) +
//                                flux_yz(i,j,k+1,0) - flux_yz(i,j,k,0)) * dxinv;
//            }else
//            {
                divy(i,j,k) += preFac*(nineOver8*flux_xy(i+1,j,k,1) - nineOver8*flux_xy(i,j,k,1) - oneOver24*flux_xy(i+2,j,k,1) + oneOver24*flux_xy(i-1,j,k,1) +
                                nineOver8*flux_cc(i,j,k,1) - nineOver8*flux_cc(i,j-1,k,1) - oneOver24*flux_cc(i,j+1,k,1) + oneOver24*flux_cc(i,j-2,k,1) +
                                nineOver8*flux_yz(i,j,k+1,0) - nineOver8*flux_yz(i,j,k,0) - oneOver24*flux_yz(i,j,k+2,0) + oneOver24*flux_yz(i,j,k-1,0)) * dxinv;
//            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
//            if(i==zlo.x || i==zhi.x || j==zlo.y || j==zhi.y || k==zlo.z || k==zhi.z)
//            {
//                divz(i,j,k) += (flux_xz(i+1,j,k,1) - flux_xz(i,j,k,1) +
//                                flux_yz(i,j+1,k,1) - flux_yz(i,j,k,1) +
//                                flux_cc(i,j,k,2) - flux_cc(i,j,k-1,2)) * dxinv;

//            }else
//            {
                divz(i,j,k) += preFac*(nineOver8*flux_xz(i+1,j,k,1) - nineOver8*flux_xz(i,j,k,1) - oneOver24*flux_xz(i+2,j,k,1) + oneOver24*flux_xz(i-1,j,k,1) +
                                nineOver8*flux_yz(i,j+1,k,1) - nineOver8*flux_yz(i,j,k,1) - oneOver24*flux_yz(i,j+2,k,1) + oneOver24*flux_yz(i,j-1,k,1) +
                                nineOver8*flux_cc(i,j,k,2) - nineOver8*flux_cc(i,j,k-1,2) - oneOver24*flux_cc(i,j,k+1,2) + oneOver24*flux_cc(i,j,k-2,2)) * dxinv;
//            }
        });
#endif
    }

    // m_force does not have ghost cells
    // set the value on physical boundaries to zero
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFabPhysBCDomainVel(m_force[d], geom, d);
    }
}


// compute stochastic momentum flux divergence
void StochMomFlux::StochMomFluxDivOrder3Split(std::array< MultiFab, AMREX_SPACEDIM >& m_force,
                                   const int& increment,
                                   const MultiFab& eta_cc,
                                   const std::array< MultiFab, NUM_EDGE >& eta_ed,
                                   const MultiFab& temp_cc,
                                   const std::array< MultiFab, NUM_EDGE >& temp_ed,
                                   const Vector< amrex::Real >& weights,
                                   const amrex::Real& dt) {

    BL_PROFILE_VAR("StochMomFluxDiv()",StochMomFluxDiv);

    // Take linear combination of mflux multifabs at each stage
    StochMomFlux::weightMomflux(weights);

    // Multiply weighted mflux (cc & edge) by sqrt(eta*temperature)
    StochMomFlux::multbyVarSqrtEtaTemp(eta_cc,eta_ed,temp_cc,temp_ed,dt);

    // multiply noise stored in mflux_ed_weighted
    // on walls by 0 (for slip) or sqrt(2) (for no-slip)
    StochMomFlux::MomFluxBC();

    // sync up random numbers at boundaries and ghost cells
    for (int d=0; d<NUM_EDGE; ++d) {
        mflux_ed_weighted[d].OverrideSync(geom.periodicity());
        mflux_ed_weighted[d].FillBoundary(geom.periodicity());
    }
    mflux_cc_weighted.FillBoundary(geom.periodicity());

    if (filtering_width > 0) {
        Abort("StochMomFlux: filtering_width != 0 not fully implemented yet");
        // need calls to filter_stoch_m_flux for mflux_ed and mflux_cc
        /*

         */
        mflux_cc_weighted.FillBoundary(geom.periodicity());
    }

    // calculate divergence and add to stoch_m_force
    Real dxinv = 1./(geom.CellSize()[0]);

    // if not incrementing, initialize data to zero
    if (increment == 0) {
        for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
            m_force[dir].setVal(0.,0,1,0);
        }
    }

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

        AMREX_D_TERM(const auto xlo = lbound(bx_x);,
                 const auto ylo = lbound(bx_y);,
                 const auto zlo = lbound(bx_z););

        AMREX_D_TERM(const auto xhi = ubound(bx_x);,
                 const auto yhi = ubound(bx_y);,
                 const auto zhi = ubound(bx_z););

        Real nineOver8 = 9.0/8.0;
        Real oneOver24 = 1.0/24.0;

        Real preFac = sqrt(1.0/(nineOver8*nineOver8 + oneOver24*oneOver24));

        int splitCell = 250;


#if (AMREX_SPACEDIM == 2)
        amrex::ParallelFor(bx_x,bx_y, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            divx(i,j,k) += (flux_cc(i,j,k,0) - flux_cc(i-1,j,k,0) +
                            flux_nd(i,j+1,k,0) - flux_nd(i,j,k,0)) * dxinv;
        },
                                      [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            divy(i,j,k) += (flux_nd(i+1,j,k,1) - flux_nd(i,j,k,1) +
                            flux_cc(i,j,k,1) - flux_cc(i,j-1,k,1)) * dxinv;
        });

#elif (AMREX_SPACEDIM == 3)
        amrex::ParallelFor(bx_x,bx_y,bx_z, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if(i<splitCell || i==xhi.x || j==xlo.y || j==xhi.y || k==xlo.z || k==xhi.z )
            {
                divx(i,j,k) += (flux_cc(i,j,k,0) - flux_cc(i-1,j,k,0) +
                                flux_xy(i,j+1,k,0) - flux_xy(i,j,k,0) +
                                flux_xz(i,j,k+1,0) - flux_xz(i,j,k,0)) * dxinv;
            }else
            {
                divx(i,j,k) += preFac*(nineOver8*flux_cc(i,j,k,0) - nineOver8*flux_cc(i-1,j,k,0) - oneOver24*flux_cc(i+1,j,k,0) + oneOver24*flux_cc(i-2,j,k,0) +
                     nineOver8*flux_xy(i,j+1,k,0) -  nineOver8*flux_xy(i,j,k,0) - oneOver24*flux_xy(i,j+2,k,0) + oneOver24*flux_xy(i,j-1,k,0) +
                    nineOver8*flux_xz(i,j,k+1,0) - nineOver8*flux_xz(i,j,k,0) - oneOver24*flux_xz(i,j,k+2,0) + oneOver24*flux_xz(i,j,k-1,0)) * dxinv;

            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if(i==ylo.x || i==yhi.x || j==ylo.y || j==yhi.y || k==ylo.z || k==yhi.z)
            {
                divy(i,j,k) += (flux_xy(i+1,j,k,1) - flux_xy(i,j,k,1) +
                                flux_cc(i,j,k,1) - flux_cc(i,j-1,k,1) +
                                flux_yz(i,j,k+1,0) - flux_yz(i,j,k,0)) * dxinv;
            }else
            {
                divy(i,j,k) += preFac*(nineOver8*flux_xy(i+1,j,k,1) - nineOver8*flux_xy(i,j,k,1) - oneOver24*flux_xy(i+2,j,k,1) + oneOver24*flux_xy(i-1,j,k,1) +
                                nineOver8*flux_cc(i,j,k,1) - nineOver8*flux_cc(i,j-1,k,1) - oneOver24*flux_cc(i,j+1,k,1) + oneOver24*flux_cc(i,j-2,k,1) +
                                nineOver8*flux_yz(i,j,k+1,0) - nineOver8*flux_yz(i,j,k,0) - oneOver24*flux_yz(i,j,k+2,0) + oneOver24*flux_yz(i,j,k-1,0)) * dxinv;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if(i==zlo.x || i==zhi.x || j==zlo.y || j==zhi.y || k==zlo.z || k==zhi.z)
            {
                divz(i,j,k) += (flux_xz(i+1,j,k,1) - flux_xz(i,j,k,1) +
                                flux_yz(i,j+1,k,1) - flux_yz(i,j,k,1) +
                                flux_cc(i,j,k,2) - flux_cc(i,j,k-1,2)) * dxinv;

            }else
            {
                divz(i,j,k) += preFac*(nineOver8*flux_xz(i+1,j,k,1) - nineOver8*flux_xz(i,j,k,1) - oneOver24*flux_xz(i+2,j,k,1) + oneOver24*flux_xz(i-1,j,k,1) +
                                nineOver8*flux_yz(i,j+1,k,1) - nineOver8*flux_yz(i,j,k,1) - oneOver24*flux_yz(i,j+2,k,1) + oneOver24*flux_yz(i,j-1,k,1) +
                                nineOver8*flux_cc(i,j,k,2) - nineOver8*flux_cc(i,j,k-1,2) - oneOver24*flux_cc(i,j,k+1,2) + oneOver24*flux_cc(i,j,k-2,2)) * dxinv;
            }
        });
#endif
    }

    // m_force does not have ghost cells
    // set the value on physical boundaries to zero
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFabPhysBCDomainVel(m_force[d], geom, d);
    }
}

// utility to write out random number MultiFabs to plotfiles
void StochMomFlux::writeMFs(std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv) {

    BL_PROFILE_VAR("writeMFs()",writeMFs);

    std::string plotfilename;
    std::string dimStr = "xyz";

    // Write out original fluxes
    for (int i=0; i<n_rngs; ++i){
        plotfilename = "a_mfluxcc_stage"+std::to_string(i);
        VisMF::Write(mflux_cc[i],plotfilename);

        for (int d=0; d<NUM_EDGE; ++d) {
            plotfilename = "a_mfluxnd_stage"+std::to_string(i)+"_";
            plotfilename += dimStr[d];
            VisMF::Write(mflux_ed[i][d],plotfilename);
        }
    }

    // Write out weighted fluxes
    plotfilename = "a_mfluxcc_weighted";
    VisMF::Write(mflux_cc_weighted,plotfilename);

    for (int d=0; d<NUM_EDGE; ++d) {
        plotfilename = "a_mfluxnd_weighted_";
        plotfilename += dimStr[d];
        VisMF::Write(mflux_ed_weighted[d],plotfilename);
    }

    // Write out fluxdiv
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        plotfilename = "a_mfluxdiv_";
        plotfilename += dimStr[d];
        VisMF::Write(mfluxdiv[d],plotfilename);
    }
}
