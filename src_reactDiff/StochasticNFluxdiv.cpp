#include "rng_functions.H"
#include "reactDiff_functions.H"

void StochasticNFluxdiv(MultiFab& n_in,
                        MultiFab& stoch_fluxdiv,
                        const std::array< MultiFab, AMREX_SPACEDIM >& diff_coef_face,
                        const Geometry& geom,
                        const Real& dt,
                        const Real& time,
                        int increment_div) {

    // single cell case set stochastic mass fluxdiv to zero
    // (or its increment if increment_in=T) and return
    long cell_count = (AMREX_SPACEDIM==2) ? n_cells[0]*n_cells[1] : n_cells[0]*n_cells[1]*n_cells[2];
    if (cell_count == 1 && increment_div==0) {
        stoch_fluxdiv.setVal(0.);
        return;
    }

    BoxArray ba = n_in.boxArray();
    DistributionMapping dmap = n_in.DistributionMap();

    std::array< MultiFab, AMREX_SPACEDIM > flux;
    AMREX_D_TERM(flux[0].define(convert(ba,nodal_flag_x), dmap, nspecies, 0);,
                 flux[1].define(convert(ba,nodal_flag_y), dmap, nspecies, 0);,
                 flux[2].define(convert(ba,nodal_flag_z), dmap, nspecies, 0););

    std::array< MultiFab, AMREX_SPACEDIM > rand;
    AMREX_D_TERM(rand[0].define(convert(ba,nodal_flag_x), dmap, nspecies, 0);,
                 rand[1].define(convert(ba,nodal_flag_y), dmap, nspecies, 0);,
                 rand[2].define(convert(ba,nodal_flag_z), dmap, nspecies, 0););

    const Real* dx = geom.CellSize();

    Real dv = (AMREX_SPACEDIM == 3) ? dx[0]*dx[1]*dx[2]*cell_depth : dx[0]*dx[1]*cell_depth;

    // average n_in to faces, store in flux
    for (MFIter mfi(n_in); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        const Array4<const Real>& n_arr = n_in.array(mfi);

        AMREX_D_TERM(const Array4<Real> & fluxx = flux[0].array(mfi);,
                     const Array4<Real> & fluxy = flux[1].array(mfi);,
                     const Array4<Real> & fluxz = flux[2].array(mfi););

        AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                     const Box & bx_y = mfi.nodaltilebox(1);,
                     const Box & bx_z = mfi.nodaltilebox(2););

        amrex::ParallelFor(bx_x, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            fluxx(i,j,k,n) = average_to_faces(n_arr(i-1,j,k,n),n_arr(i,j,k,n),dv);
        },
                           bx_y, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            fluxy(i,j,k,n) = average_to_faces(n_arr(i,j-1,k,n),n_arr(i,j,k,n),dv);
        }
#if (AMREX_SPACEDIM == 3)
                         , bx_z, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            fluxz(i,j,k,n) = average_to_faces(n_arr(i,j,k-1,n),n_arr(i,j,k,n),dv);
        }
#endif
        );
    }

    // generate random numbers
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        for (int n=0; n<nspecies; ++n) {
            MultiFabFillRandom(rand[i], n, 1., geom, 0);
        }
    }

    // assemble_stoch_n_fluxes
    for (MFIter mfi(n_in); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        const Array4<const Real>& n_arr = n_in.array(mfi);

        AMREX_D_TERM(const Array4<Real> & fluxx = flux[0].array(mfi);,
                     const Array4<Real> & fluxy = flux[1].array(mfi);,
                     const Array4<Real> & fluxz = flux[2].array(mfi););

        AMREX_D_TERM(const Array4<Real> & randx = rand[0].array(mfi);,
                     const Array4<Real> & randy = rand[1].array(mfi);,
                     const Array4<Real> & randz = rand[2].array(mfi););

        AMREX_D_TERM(const Array4<const Real> & coefx = diff_coef_face[0].array(mfi);,
                     const Array4<const Real> & coefy = diff_coef_face[1].array(mfi);,
                     const Array4<const Real> & coefz = diff_coef_face[2].array(mfi););

        AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                     const Box & bx_y = mfi.nodaltilebox(1);,
                     const Box & bx_z = mfi.nodaltilebox(2););

        amrex::ParallelFor(bx_x, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            fluxx(i,j,k,n) = std::sqrt(coefx(i,j,k,n)*fluxx(i,j,k,n)) * randx(i,j,k,n);
        },
                           bx_y, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            fluxy(i,j,k,n) = std::sqrt(coefy(i,j,k,n)*fluxy(i,j,k,n)) * randy(i,j,k,n);
        }
#if (AMREX_SPACEDIM == 3)
                         , bx_z, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            fluxz(i,j,k,n) = std::sqrt(coefz(i,j,k,n)*fluxz(i,j,k,n)) * randz(i,j,k,n);
        }
#endif
        );
    }

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_mass_lo[i] != -1 || bc_mass_hi[i] != -1) {
            Abort("StochasticNFluxdiv() - implement physical bc's for noise");
        }
    }

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        flux[i].mult(std::sqrt(2.*variance_coef_mass/(dv*dt)));
    }

    // compute flux divergence
    ComputeDiv(stoch_fluxdiv, flux, 0, 0, nspecies, geom, increment_div);

}

AMREX_GPU_HOST_DEVICE Real average_to_faces(const Real& value1,
                                            const Real& value2,
                                            const Real& dv) {

    if (avg_type == 1) { // Arithmetic with a C0-smoothed Heaviside

        if ( (value1 <= 0.) || (value2 <= 0.) ) {
            return 0.;
        } else {
            Real tmp1=std::min(dv*value1,1.);
            Real tmp2=std::min(dv*value2,1.);
            return (value1+value2)/2.*tmp1*tmp2;
        }

    } else if (avg_type == 2) { // Geometric

        Abort("average_to_faces: unimplemented avg_type");
        return 0;

    } else if (avg_type == 3) { // Harmonic
        // What we want here is the harmonic mean of max(value1,0) and max(value2,0)
        // Where we define the result to be zero if either one is zero
        // But numerically we want to avoid here division by zero

        Abort("average_to_faces: unimplemented avg_type");
        return 0;

    } else if (avg_type == 10) { // Arithmetic with (discontinuous) Heaviside

        Abort("average_to_faces: unimplemented avg_type");
        return 0;

    } else if (avg_type == 11) { // Arithmetic with C1-smoothed Heaviside

        Abort("average_to_faces: unimplemented avg_type");
        return 0;

    } else if (avg_type == 12) { // Arithmetic with C2-smoothed Heaviside

        Abort("average_to_faces: unimplemented avg_type");
        return 0;

    } else {

        Abort("average_to_faces: unimplemented avg_type");
        return 0;

    }

}
