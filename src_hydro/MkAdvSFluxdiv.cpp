#include "hydro_functions.H"

#include "common_functions.H"

using namespace amrex;

// takes cell-centered s_in
// fixme: not boundary-condition aware
void MkAdvSFluxdiv_cc(const std::array<MultiFab, AMREX_SPACEDIM>& umac_in,
    const MultiFab& s_in,
    MultiFab& s_update_in,
    const Geometry& geom,
    const int& scomp,
    const int& ncomp,
    const int& increment)
{

    BL_PROFILE_VAR("MkAdvSFluxdiv_cc()",MkAdvSFluxdiv_cc);

    Real dx = geom.CellSize(0);
    Real dxinv = 1./dx;

    // if not incrementing, initialize data to zero
    if (increment == 0) {
        s_update_in.setVal(0.,scomp,ncomp,0);
    }

    // Loop over boxes
    for (MFIter mfi(s_in); mfi.isValid(); ++mfi) {

        // Create cell-centered box from semi-nodal box
        const Box& bx = mfi.validbox();

        Array4<Real const> const& s = s_in.array(mfi);

        Array4<Real> const& s_update = s_update_in.array(mfi);

        AMREX_D_TERM(Array4<Real const> const& umac = umac_in[0].array(mfi);,
                     Array4<Real const> const& vmac = umac_in[1].array(mfi);,
                     Array4<Real const> const& wmac = umac_in[2].array(mfi););
        amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            s_update(i,j,k,scomp+n) -=
                + dxinv*( 0.5*(s(i+1,j,k,scomp+n)+s(i,j,k,scomp+n))*umac(i+1,j,k) - 0.5*(s(i,j,k,scomp+n)+s(i-1,j,k,scomp+n))*umac(i,j,k) )
                + dxinv*( 0.5*(s(i,j+1,k,scomp+n)+s(i,j,k,scomp+n))*vmac(i,j+1,k) - 0.5*(s(i,j,k,scomp+n)+s(i,j-1,k,scomp+n))*vmac(i,j,k) )
#if (AMREX_SPACEDIM == 3)
                + dxinv*( 0.5*(s(i,j,k+1,scomp+n)+s(i,j,k,scomp+n))*wmac(i,j,k+1) - 0.5*(s(i,j,k,scomp+n)+s(i,j,k-1,scomp+n))*wmac(i,j,k) )
#endif
                ;
        });
    }
}

// takes face-centered s_fc_in
void MkAdvSFluxdiv(const std::array<MultiFab, AMREX_SPACEDIM>& umac_in,
    const std::array<MultiFab, AMREX_SPACEDIM>& s_fc_in,
    MultiFab& s_update_in,
    const Geometry& geom,
    const int& scomp,
    const int& ncomp,
    const int& increment)
{

    BL_PROFILE_VAR("MkAdvSFluxdiv()",MkAdvSFluxdiv);

    Real dx = geom.CellSize(0);
    Real dxinv = 1./dx;


    // if not incrementing, initialize data to zero
    if (increment == 0) {
        s_update_in.setVal(0.,scomp,ncomp,0);
    }

    // Loop over boxes
    for (MFIter mfi(s_update_in); mfi.isValid(); ++mfi) {

        // Create cell-centered box from semi-nodal box
        const Box& bx = mfi.validbox();

        Array4<Real> const& s_update = s_update_in.array(mfi);

        AMREX_D_TERM(Array4<Real const> const& sx = s_fc_in[0].array(mfi);,
                     Array4<Real const> const& sy = s_fc_in[1].array(mfi);,
                     Array4<Real const> const& sz = s_fc_in[2].array(mfi););

        AMREX_D_TERM(Array4<Real const> const& umac = umac_in[0].array(mfi);,
                     Array4<Real const> const& vmac = umac_in[1].array(mfi);,
                     Array4<Real const> const& wmac = umac_in[2].array(mfi););

        amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            s_update(i,j,k,scomp+n) -=
                + dxinv*( sx(i+1,j,k,scomp+n)*umac(i+1,j,k) - sx(i,j,k,scomp+n)*umac(i,j,k) )
                + dxinv*( sy(i,j+1,k,scomp+n)*vmac(i,j+1,k) - sy(i,j,k,scomp+n)*vmac(i,j,k) )
#if (AMREX_SPACEDIM == 3)
                + dxinv*( sz(i,j,k+1,scomp+n)*wmac(i,j,k+1) - sz(i,j,k,scomp+n)*wmac(i,j,k) )
#endif
                ;
        });
    }
}