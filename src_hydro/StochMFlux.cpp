#include "rng_functions.H"
#include "gmres_functions.H"
#include "common_functions.H"
#include "hydro_functions_F.H"
#include "StochMFlux.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>


StochMFlux::StochMFlux(BoxArray ba_in, DistributionMapping dmap_in, Geometry geom_in,
		       int n_rngs_in) {

  BL_PROFILE_VAR("StochMFlux::StochMFlux()",StochMFlux);

  if (filtering_width != 0) {
    Abort("StochMFlux: filtering_width != 0 not fully implemented yet");
  }

  n_rngs = n_rngs_in;
  geom = geom_in;

  mflux_cc.resize(n_rngs);
  mflux_ed.resize(n_rngs);

  // Here we store all the random number stages at all spatial locations
  for (int i=0; i<n_rngs; ++i) {
      mflux_cc[i].define(ba_in, dmap_in, AMREX_SPACEDIM, std::max(1,filtering_width));
      mflux_cc[i].setVal(0.);
#if (AMREX_SPACEDIM == 2)
      mflux_ed[i][0].define(convert(ba_in,nodal_flag), dmap_in, ncomp_ed, filtering_width);
      mflux_ed[i][0].setVal(0.);
#elif (AMREX_SPACEDIM == 3)
      mflux_ed[i][0].define(convert(ba_in,nodal_flag_xy), dmap_in, ncomp_ed, filtering_width);
      mflux_ed[i][1].define(convert(ba_in,nodal_flag_xz), dmap_in, ncomp_ed, filtering_width);
      mflux_ed[i][2].define(convert(ba_in,nodal_flag_yz), dmap_in, ncomp_ed, filtering_width);

      for (int d=0; d<AMREX_SPACEDIM; ++d) mflux_ed[i][d].setVal(0.);
#endif
  }

  // Temporary storage for linear combinations of random number stages
  mflux_cc_weighted.define(ba_in, dmap_in, AMREX_SPACEDIM, std::max(1,filtering_width));
  mflux_cc_weighted.setVal(0.);
#if (AMREX_SPACEDIM == 2)
  mflux_ed_weighted[0].define(convert(ba_in,nodal_flag), dmap_in, ncomp_ed, filtering_width);
  mflux_ed_weighted[0].setVal(0.);
#elif (AMREX_SPACEDIM == 3)
  mflux_ed_weighted[0].define(convert(ba_in,nodal_flag_xy), dmap_in, ncomp_ed, filtering_width);
  mflux_ed_weighted[1].define(convert(ba_in,nodal_flag_xz), dmap_in, ncomp_ed, filtering_width);
  mflux_ed_weighted[2].define(convert(ba_in,nodal_flag_yz), dmap_in, ncomp_ed, filtering_width);

  for (int d=0; d<AMREX_SPACEDIM; ++d) mflux_ed_weighted[d].setVal(0.);
#endif
}


void StochMFlux::weightMflux(Vector< amrex::Real > weights) {

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
}

void StochMFlux::fillMStochastic() {

    BL_PROFILE_VAR("StochMFlux::fillMStochastic()",StochMFlux);

    for (int i=0; i<n_rngs; ++i) {

        switch(stoch_stress_form) {

        case 0: // Non-symmetric
            // Print() << "Non-symmetric \n";
            for (int n=0; n<AMREX_SPACEDIM; ++n) {
                MultiFABFillRandom(mflux_cc[i],n,1.0,geom);
            }

            for (int d=0; d<NUM_EDGE; ++d) {
                for (int n=0; n<ncomp_ed; ++n) {
                    MultiFABFillRandom(mflux_ed[i][d],n,1.0,geom);
                }
            }
            break;

        default: // Symmetric
            // Print() << "Symmetric \n";
            for (int n=0; n<AMREX_SPACEDIM; ++n) {
                MultiFABFillRandom(mflux_cc[i],n,2.0,geom);
            }

            for (int d=0; d<NUM_EDGE; ++d) {
                MultiFABFillRandom(mflux_ed[i][d],0,1.0,geom);
                MultiFab::Copy(mflux_ed[i][d], mflux_ed[i][d], 0, 1, ncomp_ed-1, 0);
            }
            break;
        }
    }
}

void StochMFlux::MfluxBC() {

#if (AMREX_SPACEDIM == 2)

    // lo-x domain boundary
    if (bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) {

        // 0 = slip wall   : multiply fluxes on wall by 0
        // 1 = no-slip wall: multiply fluxes on wall by sqrt(2)
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
                AMREX_HOST_DEVICE_FOR_4D(b, 2, i, j, k, n,
                {
                    mflux_nd(i,j,k,n) *= factor;
                });
            }
        }
    }
    else if (bc_vel_lo[0] != -1) {
        // DEBUGING deterministic boundary conditions
        // Abort("MfluxBC unsupported bc type");
    }

    // hi-x domain boundary
    if (bc_vel_hi[0] == 1 || bc_vel_hi[0] == 2) {

        // 0 = slip wall   : multiply fluxes on wall by 0
        // 1 = no-slip wall: multiply fluxes on wall by sqrt(2)
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
                AMREX_HOST_DEVICE_FOR_4D(b, 2, i, j, k, n,
                {
                    mflux_nd(i,j,k,n) *= factor;
                });
            }
        }
    }
    else if (bc_vel_hi[0] != -1) {
        // DEBUGING deterministic boundary conditions
        // Abort("MfluxBC unsupported bc type");
    }

    // lo-y domain boundary
    if (bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) {

        // 0 = slip wall   : multiply fluxes on wall by 0
        // 1 = no-slip wall: multiply fluxes on wall by sqrt(2)
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
                AMREX_HOST_DEVICE_FOR_4D(b, 2, i, j, k, n,
                {
                    mflux_nd(i,j,k,n) *= factor;
                });
            }
        }
    }
    else if (bc_vel_lo[1] != -1) {
        // DEBUGING deterministic boundary conditions
        // Abort("MfluxBC unsupported bc type");
    }

    // hi-y domain boundary
    if (bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) {

        // 0 = slip wall   : multiply fluxes on wall by 0
        // 1 = no-slip wall: multiply fluxes on wall by sqrt(2)
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
                AMREX_HOST_DEVICE_FOR_4D(b, 2, i, j, k, n,
                {
                    mflux_nd(i,j,k,n) *= factor;
                });
            }
        }
    }
    else if (bc_vel_hi[1] != -1) {
        // DEBUGING deterministic boundary conditions
        // Abort("MfluxBC unsupported bc type");
    }

#elif (AMREX_SPACEDIM == 3)

    // lo-x domain boundary, y-facing fluxes
    // lo-x domain boundary, z-facing fluxes
    if (bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) {

        // 0 = slip wall   : multiply fluxes on wall by 0
        // 1 = no-slip wall: multiply fluxes on wall by sqrt(2)
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
                AMREX_HOST_DEVICE_FOR_4D(b, 2, i, j, k, n,
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
                AMREX_HOST_DEVICE_FOR_4D(b, 2, i, j, k, n,
                {
                    mflux_xz(i,j,k,n) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
    }
    else if (bc_vel_lo[0] != -1) {
        // DEBUGING deterministic boundary conditions
        // Abort("MfluxBC unsupported bc type");
    }

    // hi-x domain boundary, y-facing fluxes
    // hi-x domain boundary, z-facing fluxes
    if (bc_vel_hi[0] == 1 || bc_vel_hi[0] == 2) {

        // 0 = slip wall   : multiply fluxes on wall by 0
        // 1 = no-slip wall: multiply fluxes on wall by sqrt(2)
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
                AMREX_HOST_DEVICE_FOR_4D(b, 2, i, j, k, n,
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
                AMREX_HOST_DEVICE_FOR_4D(b, 2, i, j, k, n,
                {
                    mflux_xz(i,j,k,n) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
    }
    else if (bc_vel_hi[0] != -1) {
        // DEBUGING deterministic boundary conditions
        // Abort("MfluxBC unsupported bc type");
    }

    // lo-y domain boundary, x-facing fluxes
    // lo-y domain boundary, z-facing fluxes
    if (bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) {

        // 0 = slip wall   : multiply fluxes on wall by 0
        // 1 = no-slip wall: multiply fluxes on wall by sqrt(2)
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
                AMREX_HOST_DEVICE_FOR_4D(b, 2, i, j, k, n,
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
                AMREX_HOST_DEVICE_FOR_4D(b, 2, i, j, k, n,
                {
                    mflux_yz(i,j,k,n) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
    }
    else if (bc_vel_lo[1] != -1) {
        // DEBUGING deterministic boundary conditions
        // Abort("MfluxBC unsupported bc type");
    }

    // hi-y domain boundary, x-facing fluxes
    // hi-y domain boundary, z-facing fluxes
    if (bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) {

        // 0 = slip wall   : multiply fluxes on wall by 0
        // 1 = no-slip wall: multiply fluxes on wall by sqrt(2)
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
                AMREX_HOST_DEVICE_FOR_4D(b, 2, i, j, k, n,
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
                AMREX_HOST_DEVICE_FOR_4D(b, 2, i, j, k, n,
                {
                    mflux_yz(i,j,k,n) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
    }
    else if (bc_vel_hi[1] != -1) {
        // DEBUGING deterministic boundary conditions
        // Abort("MfluxBC unsupported bc type");
    }

    // lo-z domain boundary, x-facing fluxes
    // lo-z domain boundary, y-facing fluxes
    if (bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) {

        // 0 = slip wall   : multiply fluxes on wall by 0
        // 1 = no-slip wall: multiply fluxes on wall by sqrt(2)
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
                AMREX_HOST_DEVICE_FOR_4D(b, 2, i, j, k, n,
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
                AMREX_HOST_DEVICE_FOR_4D(b, 2, i, j, k, n,
                {
                    mflux_yz(i,j,k,n) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
    }
    else if (bc_vel_lo[2] != -1) {
        // DEBUGING deterministic boundary conditions
        // Abort("MfluxBC unsupported bc type");
    }

    // hi-z domain boundary, x-facing fluxes
    // hi-z domain boundary, y-facing fluxes
    if (bc_vel_hi[2] == 1 || bc_vel_hi[2] == 2) {

        // 0 = slip wall   : multiply fluxes on wall by 0
        // 1 = no-slip wall: multiply fluxes on wall by sqrt(2)
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
                AMREX_HOST_DEVICE_FOR_4D(b, 2, i, j, k, n,
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
                AMREX_HOST_DEVICE_FOR_4D(b, 2, i, j, k, n,
                {
                    mflux_yz(i,j,k,n) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
    }
    else if (bc_vel_hi[2] != -1) {
        // DEBUGING deterministic boundary conditions
        // Abort("MfluxBC unsupported bc type");
    }

#endif

}


AMREX_GPU_HOST_DEVICE
inline
void mult_by_sqrt_eta_temp (const Box & tbx,
			    const Box& bx,
			    const Box& bx_xy,
#if (AMREX_SPACEDIM == 3)
			    const Box& bx_xz,
			    const Box& bx_yz,
#endif
			    const Array4<Real> mflux_cc,
			    const Array4<Real> mflux_xy,
#if (AMREX_SPACEDIM == 3)
			    const Array4<Real> mflux_xz,
			    const Array4<Real> mflux_yz,
#endif
			    const Array4<Real const> eta_cc,
			    const Array4<Real const> eta_xy,
#if (AMREX_SPACEDIM == 3)
			    const Array4<Real const> eta_xz,
			    const Array4<Real const> eta_yz,
#endif
			    const Array4<Real const> temp_cc,
			    const Array4<Real const> temp_xy
#if (AMREX_SPACEDIM == 3)
			    , const Array4<Real const> temp_xz,
			    const Array4<Real const> temp_yz
#endif
			    ) noexcept {

    // bx is the cell-centered box grown with 1 ghost cell
    // bx_xy, bx_xz, and bx_yz are the doubly nodal boxes

    // if running on the host: tlo is the minimal box contains the union of the
    // face-centered grid boxes

    // if running on the gpu: tlo is a box with a single point that comes from
    // the union of the face-centered grid boxes

    const auto tlo = lbound(tbx);
    const auto thi = ubound(tbx);

    // if running on the host, lo and hi are set to the lower/upper
    // bounds of the box of interest

    // if running on the gpu, lo and hi are set to the single point
    // defined by tlo, unless tlo is outside of the box of interest,
    // in which case they are set to values that make sure the loop
    // is not entered

    {
      const auto lo = amrex::elemwiseMax(tlo, lbound(bx));
      const auto hi = amrex::elemwiseMin(thi, ubound(bx));

      for (int n=0; n<AMREX_SPACEDIM; ++n) {
      for (int k=lo.z; k<=hi.z; ++k) {
      for (int j=lo.y; j<=hi.y; ++j) {
      AMREX_PRAGMA_SIMD
      for (int i=lo.x; i<=hi.x; ++i) {
	mflux_cc(i,j,k,n) *= sqrt(eta_cc(i,j,k)*temp_cc(i,j,k));
	}
	}
        }
	}
    }

    {
      const auto lo = amrex::elemwiseMax(tlo, lbound(bx_xy));
      const auto hi = amrex::elemwiseMin(thi, ubound(bx_xy));

      for (int n=0; n<2; ++n) {
      for (int k=lo.z; k<=hi.z; ++k) {
      for (int j=lo.y; j<=hi.y; ++j) {
      AMREX_PRAGMA_SIMD
      for (int i=lo.x; i<=hi.x; ++i) {
	mflux_xy(i,j,k,n) *= sqrt(eta_xy(i,j,k)*temp_xy(i,j,k));
      }
      }
      }
      }
    }

#if (AMREX_SPACEDIM == 3)

    {
      const auto lo = amrex::elemwiseMax(tlo, lbound(bx_xz));
      const auto hi = amrex::elemwiseMin(thi, ubound(bx_xz));

      for (int n=0; n<2; ++n) {
      for (int k=lo.z; k<=hi.z; ++k) {
      for (int j=lo.y; j<=hi.y; ++j) {
      AMREX_PRAGMA_SIMD
      for (int i=lo.x; i<=hi.x; ++i) {
	mflux_xz(i,j,k,n) *= sqrt(eta_xz(i,j,k)*temp_xz(i,j,k));
      }
      }
      }
      }
    }

    {
      const auto lo = amrex::elemwiseMax(tlo, lbound(bx_yz));
      const auto hi = amrex::elemwiseMin(thi, ubound(bx_yz));

      for (int n=0; n<2; ++n) {
      for (int k=lo.z; k<=hi.z; ++k) {
      for (int j=lo.y; j<=hi.y; ++j) {
      AMREX_PRAGMA_SIMD
      for (int i=lo.x; i<=hi.x; ++i) {
	mflux_yz(i,j,k,n) *= sqrt(eta_yz(i,j,k)*temp_yz(i,j,k));
      }
      }
      }
      }
    }

#endif

}

void StochMFlux::multbyVarSqrtEtaTemp(const MultiFab& eta_cc,
				      const std::array< MultiFab, NUM_EDGE >& eta_ed,
				      const MultiFab& temp_cc,
				      const std::array< MultiFab, NUM_EDGE >& temp_ed,
				      const amrex::Real& dt) {

  const Real* dx = geom.CellSize();

  Real dVol = dx[0]*dx[1];
  if (AMREX_SPACEDIM == 2) {
    dVol *= cell_depth;
  } else {
    if (AMREX_SPACEDIM == 3) {
      dVol *= dx[2];
    }
  }

  // Compute variance using computed differential volume
  Real variance = sqrt(variance_coef_mom*2.0*k_B/(dVol*dt));

  // Scale mflux_weighted by variance
  mflux_cc_weighted.mult(variance, filtering_width);
  for (int d=0; d<NUM_EDGE; d++) {
    mflux_ed_weighted[d].mult(variance, filtering_width);
  }

  // Multiply mflux_weighted by sqrt(eta*temperature)
  // Loop over boxes
  for (MFIter mfi(mflux_cc_weighted); mfi.isValid(); ++mfi) {
    // Note: Make sure that multifab is cell-centered
    const Box& bx = mfi.growntilebox(1);

    const Box & bx_xy = mfi.tilebox(nodal_flag_xy);
#if (AMREX_SPACEDIM == 3)
    const Box & bx_xz = mfi.tilebox(nodal_flag_xz);
    const Box & bx_yz = mfi.tilebox(nodal_flag_yz);
#endif

    const Array4<Real> & mflux_cc_fab = mflux_cc_weighted.array(mfi);
    const Array4<Real> & mflux_xy_fab = mflux_ed_weighted[0].array(mfi);
#if (AMREX_SPACEDIM == 3)
    const Array4<Real> & mflux_xz_fab = mflux_ed_weighted[1].array(mfi);
    const Array4<Real> & mflux_yz_fab = mflux_ed_weighted[2].array(mfi);
#endif

    const Array4<Real const> & eta_cc_fab = eta_cc.array(mfi);
    const Array4<Real const> & eta_xy_fab = eta_ed[0].array(mfi);
#if (AMREX_SPACEDIM == 3)
    const Array4<Real const> & eta_xz_fab = eta_ed[1].array(mfi);
    const Array4<Real const> & eta_yz_fab = eta_ed[2].array(mfi);
#endif

    const Array4<Real const> & temp_cc_fab = temp_cc.array(mfi);
    const Array4<Real const> & temp_xy_fab = temp_ed[0].array(mfi);
#if (AMREX_SPACEDIM == 3)
    const Array4<Real const> & temp_xz_fab = temp_ed[1].array(mfi);
    const Array4<Real const> & temp_yz_fab = temp_ed[2].array(mfi);
#endif

    AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
    {
      mult_by_sqrt_eta_temp(tbx, bx, bx_xy,
#if (AMREX_SPACEDIM == 3)
			    bx_xz, bx_yz,
#endif
			    mflux_cc_fab, mflux_xy_fab,
#if (AMREX_SPACEDIM == 3)
			    mflux_xz_fab, mflux_yz_fab,
#endif
			    eta_cc_fab, eta_xy_fab,
#if (AMREX_SPACEDIM == 3)
			    eta_xz_fab, eta_yz_fab,
#endif
			    temp_cc_fab, temp_xy_fab
#if (AMREX_SPACEDIM == 3)
			    , temp_xz_fab, temp_yz_fab
#endif
			    );
    });
  }
}

void StochMFlux::StochMFluxDiv(std::array< MultiFab, AMREX_SPACEDIM >& m_force,
                               const int& increment,
                               const MultiFab& eta_cc,
                               const std::array< MultiFab, NUM_EDGE >& eta_ed,
                               const MultiFab& temp_cc,
                               const std::array< MultiFab, NUM_EDGE >& temp_ed,
                               const Vector< amrex::Real >& weights,
                               const amrex::Real& dt) {

  BL_PROFILE_VAR("StochMFlux::StochMfluxDiv()",StochMfluxDiv);

  // Take linear combination of mflux multifabs at each stage
  StochMFlux::weightMflux(weights);

  // Multiply weighted mflux (cc & edge) by sqrt(eta*temperature)
  StochMFlux::multbyVarSqrtEtaTemp(eta_cc,eta_ed,temp_cc,temp_ed,dt);

  // multiply noise stored in mflux_ed_weighted
  // on walls by 0 (for slip) or sqrt(2) (for no-slip)
  StochMFlux::MfluxBC();

  // sync up random numbers at boundaries and ghost cells
  for (int d=0; d<NUM_EDGE; ++d) {
      mflux_ed_weighted[d].OverrideSync(geom.periodicity());
      mflux_ed_weighted[d].FillBoundary(geom.periodicity());
  }
  mflux_cc_weighted.FillBoundary(geom.periodicity());

  if (filtering_width > 0) {
      Abort("StochMFlux: filtering_width != 0 not fully implemented yet");
      // need calls to filter_stoch_m_flux for mflux_ed and mflux_cc
      /*

       */
      mflux_cc_weighted.FillBoundary(geom.periodicity());
  }

  // calculate divergence and add to stoch_m_force
  Real dxinv = 1./(geom.CellSize()[0]);

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

  // m_force does not have ghost cells
  // set the value on physical boundaries to zero
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    MultiFABPhysBCDomainVel(m_force[d], geom, d);
  }
}

void StochMFlux::addMfluctuations(std::array< MultiFab, AMREX_SPACEDIM >& umac,
				  const MultiFab& rhotot, const MultiFab& Temp,
				  const amrex::Real& variance) {

  std::array< MultiFab, AMREX_SPACEDIM > m_old;
  std::array< MultiFab, AMREX_SPACEDIM > rhotot_fc;
  std::array< MultiFab, AMREX_SPACEDIM > Temp_fc;

  for (int d=0; d<AMREX_SPACEDIM; d++) {
    m_old[d].define(     umac[d].boxArray(), umac[d].DistributionMap(), 1, 0);
    rhotot_fc[d].define( umac[d].boxArray(), umac[d].DistributionMap(), 1, 0);
    Temp_fc[d].define(   umac[d].boxArray(), umac[d].DistributionMap(), 1, 0);
  }

  // NOTE: these only operate on valid cells
  AverageCCToFace(rhotot, 0, rhotot_fc, 0, 1);
  AverageCCToFace(Temp,   0, Temp_fc,   0, 1);

  // Convert umac to momenta, rho*umac
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(     m_old[d], umac[d],      0, 0, 1, 0);
    MultiFab::Multiply( m_old[d], rhotot_fc[d], 0, 0, 1, 0);
  }

  addMfluctuations_stag(m_old, rhotot_fc, Temp_fc, variance);

  // Convert momenta to umac, (1/rho)*momentum
  for (int d=0; d<AMREX_SPACEDIM; d++) {
      MultiFab::Copy(   umac[d], m_old[d],     0, 0, 1, 0);
      MultiFab::Divide( umac[d], rhotot_fc[d], 0, 0, 1, 0);
  }
}

void StochMFlux::addMfluctuations_stag(std::array< MultiFab, AMREX_SPACEDIM >& m_old,
				       const std::array< MultiFab, AMREX_SPACEDIM >& rhotot_fc,
				       const std::array< MultiFab, AMREX_SPACEDIM >& Temp_fc,
				       const amrex::Real& variance) {

  const Real* dx = geom.CellSize();
  Real dVol = dx[0]*dx[1];
  if (AMREX_SPACEDIM == 2) {
    dVol *= cell_depth;
  } else {
    if (AMREX_SPACEDIM == 3) {
      dVol *= dx[2];
    }
  }

  // Initialize variances
  Real variance_mom = std::abs(variance)*k_B/dVol;

  std::array<MultiFab, AMREX_SPACEDIM> variance_mfab;
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    variance_mfab[d].define(m_old[d].boxArray(), m_old[d].DistributionMap(),1,0);
  }

  std::array< MultiFab, AMREX_SPACEDIM > mac_temp;
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    mac_temp[d].define(m_old[d].boxArray(), m_old[d].DistributionMap(),1,0);
  }

  // Fill momentum multifab with random numbers, scaled by equilibrium variances
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    // Set variance multifab to sqrt(rho*temp)
    MultiFab::Copy(     variance_mfab[d],rhotot_fc[d],0,0,1,0);
    MultiFab::Multiply( variance_mfab[d],Temp_fc[d],  0,0,1,0);
    SqrtMF(variance_mfab[d]);

    // Fill momentum with random numbers, scaled by sqrt(var*k_B/dV)
    MultiFABFillRandom(mac_temp[d],0,variance_mom,geom);

    // Scale random momenta further by factor of sqrt(rho*temp)
    MultiFab::Multiply(mac_temp[d],variance_mfab[d],0,0,1,0);

    MultiFab::Saxpy(m_old[d], 1.0, mac_temp[d],0,0,1,0);

    // For safety, although called by MultiFABFillRandom()
    m_old[d].OverrideSync(geom.periodicity());
    m_old[d].FillBoundary(geom.periodicity());
  }

  if (variance < 0.0) {
    // Ensure zero total momentum
    Vector<Real> av_mom;
    // take staggered sum & divide by number of cells
    SumStag(geom,m_old,0,av_mom,true);
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
      // subtract off average
      m_old[d].plus(-av_mom[d],1);
      m_old[d].OverrideSync(geom.periodicity());
      m_old[d].FillBoundary(geom.periodicity());
    }
  }

  for (int i=0; i<AMREX_SPACEDIM; i++) {
      m_old[i].FillBoundary(geom.periodicity());
      MultiFABPhysBCDomainVel(m_old[i], geom,i);
      MultiFABPhysBCMacVel(m_old[i], geom, i);
  }
}

void StochMFlux::writeMFs(std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv) {
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
