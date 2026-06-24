#include "hydro_functions.H"

#include "common_functions.H"


void MkAdvMFluxdiv(const std::array<MultiFab, AMREX_SPACEDIM>& umac_in,
    const std::array<MultiFab, AMREX_SPACEDIM>& m,
    std::array<MultiFab, AMREX_SPACEDIM>& m_update,
    const amrex::Real* dx,
    const int& increment)
{

    BL_PROFILE_VAR("MkAdvMFluxdiv()",MkAdvMFluxdiv);

    Real fourdxinv = 0.25/dx[0];

    // if not incrementing, initialize data to zero
    if (increment == 0) {
        for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
            m_update[dir].setVal(0.,0,1,0);
        }
    }

    // Loop over boxes
    for (MFIter mfi(umac_in[0],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        AMREX_D_TERM(const Array4<Real const> & umac = (umac_in[0]).array(mfi);,
                     const Array4<Real const> & vmac = (umac_in[1]).array(mfi);,
                     const Array4<Real const> & wmac = (umac_in[2]).array(mfi););

        AMREX_D_TERM(const Array4<Real const> & mx = m[0].array(mfi);,
                     const Array4<Real const> & my = m[1].array(mfi);,
                     const Array4<Real const> & mz = m[2].array(mfi););

        AMREX_D_TERM(const Array4<Real> & m_updatex = m_update[0].array(mfi);,
                     const Array4<Real> & m_updatey = m_update[1].array(mfi);,
                     const Array4<Real> & m_updatez = m_update[2].array(mfi););

        // since the MFIter is built on a nodal MultiFab we need to build the
        // nodal tileboxes for each direction in this way
        AMREX_D_TERM(Box bx_x = mfi.tilebox(nodal_flag_x);,
                     Box bx_y = mfi.tilebox(nodal_flag_y);,
                     Box bx_z = mfi.tilebox(nodal_flag_z););

#if (AMREX_SPACEDIM == 2)
        amrex::ParallelFor(bx_x,bx_y,[=]
            AMREX_GPU_DEVICE (int i, int j, int k) {
                Real fluxx_hi, fluxx_lo, fluxy_hi, fluxy_lo;
                fluxx_hi = (mx(i,j,k)+mx(i+1,j,k))*(umac(i,j,k)+umac(i+1,j,k));
                fluxx_lo = (mx(i-1,j,k)+mx(i,j,k))*(umac(i-1,j,k)+umac(i,j,k));
                fluxy_hi = (mx(i,j,k)+mx(i,j+1,k))*(vmac(i-1,j+1,k)+vmac(i,j+1,k));
                fluxy_lo = (mx(i,j-1,k)+mx(i,j,k))*(vmac(i-1,j,k)+vmac(i,j,k));
                m_updatex(i,j,k) -= ( fluxx_hi-fluxx_lo + fluxy_hi-fluxy_lo ) * fourdxinv;
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real fluxx_hi, fluxx_lo, fluxy_hi, fluxy_lo;
                fluxx_hi = (my(i,j,k)+my(i+1,j,k))*(umac(i+1,j-1,k)+umac(i+1,j,k));
                fluxx_lo = (my(i-1,j,k)+my(i,j,k))*(umac(i,j-1,k)+umac(i,j,k));
                fluxy_hi = (my(i,j,k)+my(i,j+1,k))*(vmac(i,j,k)+vmac(i,j+1,k));
                fluxy_lo = (my(i,j-1,k)+my(i,j,k))*(vmac(i,j-1,k)+vmac(i,j,k));
                m_updatey(i,j,k) -= ( fluxx_hi-fluxx_lo + fluxy_hi-fluxy_lo ) * fourdxinv;
        });

#elif (AMREX_SPACEDIM == 3)
        amrex::ParallelFor(bx_x,bx_y,bx_z,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real fluxx_hi, fluxx_lo, fluxy_hi, fluxy_lo, fluxz_hi, fluxz_lo;
                fluxx_hi = (mx(i,j,k)+mx(i+1,j,k))*(umac(i,j,k)+umac(i+1,j,k));
                fluxx_lo = (mx(i-1,j,k)+mx(i,j,k))*(umac(i-1,j,k)+umac(i,j,k));
                fluxy_hi = (mx(i,j,k)+mx(i,j+1,k))*(vmac(i-1,j+1,k)+vmac(i,j+1,k));
                fluxy_lo = (mx(i,j-1,k)+mx(i,j,k))*(vmac(i-1,j,k)+vmac(i,j,k));
                fluxz_hi = (mx(i,j,k)+mx(i,j,k+1))*(wmac(i-1,j,k+1)+wmac(i,j,k+1));
                fluxz_lo = (mx(i,j,k-1)+mx(i,j,k))*(wmac(i-1,j,k)+wmac(i,j,k));
                m_updatex(i,j,k) -= ( fluxx_hi-fluxx_lo + fluxy_hi-fluxy_lo + fluxz_hi-fluxz_lo ) * fourdxinv;
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real fluxx_hi, fluxx_lo, fluxy_hi, fluxy_lo, fluxz_hi, fluxz_lo;
                fluxx_hi = (my(i,j,k)+my(i+1,j,k))*(umac(i+1,j-1,k)+umac(i+1,j,k));
                fluxx_lo = (my(i-1,j,k)+my(i,j,k))*(umac(i,j-1,k)+umac(i,j,k));
                fluxy_hi = (my(i,j,k)+my(i,j+1,k))*(vmac(i,j,k)+vmac(i,j+1,k));
                fluxy_lo = (my(i,j-1,k)+my(i,j,k))*(vmac(i,j-1,k)+vmac(i,j,k));
                fluxz_hi = (my(i,j,k)+my(i,j,k+1))*(wmac(i,j-1,k+1)+wmac(i,j,k+1));
                fluxz_lo = (my(i,j,k-1)+my(i,j,k))*(wmac(i,j-1,k)+wmac(i,j,k));
                m_updatey(i,j,k) -= ( fluxx_hi-fluxx_lo + fluxy_hi-fluxy_lo + fluxz_hi-fluxz_lo ) * fourdxinv;
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real fluxx_hi, fluxx_lo, fluxy_hi, fluxy_lo, fluxz_hi, fluxz_lo;
                fluxx_hi = (mz(i,j,k)+mz(i+1,j,k))*(umac(i+1,j,k-1)+umac(i+1,j,k));
                fluxx_lo = (mz(i-1,j,k)+mz(i,j,k))*(umac(i,j,k-1)+umac(i,j,k));
                fluxy_hi = (mz(i,j,k)+mz(i,j+1,k))*(vmac(i,j+1,k-1)+vmac(i,j+1,k));
                fluxy_lo = (mz(i,j-1,k)+mz(i,j,k))*(vmac(i,j,k-1)+vmac(i,j,k));
                fluxz_hi = (mz(i,j,k)+mz(i,j,k+1))*(wmac(i,j,k)+wmac(i,j,k+1));
                fluxz_lo = (mz(i,j,k-1)+mz(i,j,k))*(wmac(i,j,k-1)+wmac(i,j,k));
                m_updatez(i,j,k) -= ( fluxx_hi-fluxx_lo + fluxy_hi-fluxy_lo + fluxz_hi-fluxz_lo ) * fourdxinv;
        });
#endif
    }
}