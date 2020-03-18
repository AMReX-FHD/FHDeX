#include "common_functions.H"

void AverageFaceToCC(const std::array<MultiFab, AMREX_SPACEDIM>& face,
                     MultiFab& cc, int cc_comp)
{

    BL_PROFILE_VAR("AverageFaceToCC()",AverageFaceToCC);

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(cc,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        AMREX_D_TERM(Array4<Real const> const& facex_fab = face[0].array(mfi);,
                     Array4<Real const> const& facey_fab = face[1].array(mfi);,
                     Array4<Real const> const& facez_fab = face[2].array(mfi););

        Array4<Real> const& cc_fab = cc.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            AMREX_D_TERM(cc_fab(i,j,k,cc_comp  ) = 0.5*(facex_fab(i+1,j,k) + facex_fab(i,j,k));,
                         cc_fab(i,j,k,cc_comp+1) = 0.5*(facey_fab(i,j+1,k) + facey_fab(i,j,k));,
                         cc_fab(i,j,k,cc_comp+2) = 0.5*(facez_fab(i,j,k+1) + facez_fab(i,j,k)););
        });
    }
}

void AverageCCToFace(const MultiFab& cc, int cc_comp,
                     std::array<MultiFab, AMREX_SPACEDIM>& face, int face_comp,
                     int ncomp)
{

    BL_PROFILE_VAR("AverageCCToFace()",AverageCCToFace);

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(cc,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Array4<Real const> & cc_fab = cc.array(mfi);

        AMREX_D_TERM(const Array4<Real> & facex_fab = face[0].array(mfi);,
                     const Array4<Real> & facey_fab = face[1].array(mfi);,
                     const Array4<Real> & facez_fab = face[2].array(mfi););

        AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                     const Box & bx_y = mfi.nodaltilebox(1);,
                     const Box & bx_z = mfi.nodaltilebox(2););

        amrex::ParallelFor(bx_x, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            facex_fab(i,j,k,face_comp+n) = 0.5*(cc_fab(i,j,k,cc_comp+n)+cc_fab(i-1,j,k,cc_comp+n));
        },
                           bx_y, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            facey_fab(i,j,k,face_comp+n) = 0.5*(cc_fab(i,j,k,cc_comp+n)+cc_fab(i,j-1,k,cc_comp+n));
        }
#if (AMREX_SPACEDIM == 3)
        ,
                           bx_z, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            facez_fab(i,j,k,face_comp+n) = 0.5*(cc_fab(i,j,k,cc_comp+n)+cc_fab(i,j,k-1,cc_comp+n));
        }
#endif
        );
    }

}


void ShiftFaceToCC(const MultiFab& face, int face_comp,
                     MultiFab& cc, int cc_comp,
                     int ncomp)
{

    BL_PROFILE_VAR("ShiftFaceToCC()",ShiftFaceToCC);

    if (!face.is_nodal(0) && !face.is_nodal(1) && !face.is_nodal(2)) {
        Abort("ShiftFaceToCC requires a face-centered MultiFab");
    }

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(cc,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        Array4<Real const> const& face_fab = face.array(mfi);

        Array4<Real> const& cc_fab = cc.array(mfi);

        amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            cc_fab(i,j,k,cc_comp+n) = face_fab(i,j,k,face_comp+n);
        });
    }
}
