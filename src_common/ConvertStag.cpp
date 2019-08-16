#include "common_functions.H"
#include "common_functions_F.H"

void AverageFaceToCC(const std::array<MultiFab, AMREX_SPACEDIM>& face,
                     MultiFab& cc, int cc_comp)
{

    BL_PROFILE_VAR("AverageFaceToCC()",AverageFaceToCC);

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(cc); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        AMREX_D_TERM(Array4<Real const> const& facex_fab = face[0].array(mfi);,
                     Array4<Real const> const& facey_fab = face[1].array(mfi);,
                     Array4<Real const> const& facez_fab = face[2].array(mfi););

        Array4<Real> const& cc_fab = cc.array(mfi);

        AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
        {
            AMREX_D_TERM(cc_fab(i,j,k,cc_comp  ) = 0.5*(facex_fab(i+1,j,k) + facex_fab(i,j,k));,
                         cc_fab(i,j,k,cc_comp+1) = 0.5*(facey_fab(i,j+1,k) + facey_fab(i,j,k));,
                         cc_fab(i,j,k,cc_comp+2) = 0.5*(facez_fab(i,j,k+1) + facez_fab(i,j,k)););
        });
    }
}


AMREX_GPU_HOST_DEVICE
inline
void avg_cc_to_fc (const Box & tbx,
                   const Box & xbx,
#if (AMREX_SPACEDIM == 2)
                   const Box & ybx,
#endif
#if (AMREX_SPACEDIM == 3)
                   const Box & zbx,
#endif
                   const Array4<Real> & fx,
#if (AMREX_SPACEDIM == 2)
                   const Array4<Real> & fy,
#endif
#if (AMREX_SPACEDIM == 3)
                   const Array4<Real> & fz,
#endif
                   const Array4<Real const> & cc,
                   int fcomp, int ccomp, int ncomp) noexcept
{
    // xbx, ybx, and zbx are the face-centered boxes

    // if running on the host: tlo is the minimal box contains the union of the
    // face-centered grid boxes

    // if running on the gpu: tlo is a box with a single point that comes from
    // the union of the face-centered grid boxes

    const auto tlo = lbound(tbx);
    const auto thi = ubound(tbx);

    // if running on the host, x/y/zlo and x/y/zhi are set to the lower/upper
    // bounds of x/y/zbx

    // if running on the gpu, x/y/zlo and x/y/zhi are set to the single point
    // defined by tlo, unless tlo is outside of the union of the face-centered
    // grid boxes, in which case they are set to values that make sure the loop
    // is not entered

    AMREX_D_TERM(const auto xlo = amrex::elemwiseMax(tlo, lbound(xbx));,
                 const auto ylo = amrex::elemwiseMax(tlo, lbound(ybx));,
                 const auto zlo = amrex::elemwiseMax(tlo, lbound(zbx)););

    AMREX_D_TERM(const auto xhi = amrex::elemwiseMin(thi, ubound(xbx));,
                 const auto yhi = amrex::elemwiseMin(thi, ubound(ybx));,
                 const auto zhi = amrex::elemwiseMin(thi, ubound(zbx)););

    for (int n = 0; n < ncomp; ++n) {
    for (int k = xlo.z; k <= xhi.z; ++k) {
    for (int j = xlo.y; j <= xhi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = xlo.x; i <= xhi.x; ++i) {
        fx(i,j,k,fcomp+n) = 0.5*(cc(i-1,j,k,ccomp+n) + cc(i,j,k,ccomp+n));
    }
    }
    }
    }

#if (AMREX_SPACEDIM == 2)
    for (int n = 0; n < ncomp; ++n) {
    for (int k = ylo.z; k <= yhi.z; ++k) {
    for (int j = ylo.y; j <= yhi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = ylo.x; i <= yhi.x; ++i) {
        fy(i,j,k,fcomp+n) = 0.5*(cc(i,j-1,k,ccomp+n) + cc(i,j,k,ccomp+n));
    }
    }
    }
    }
#endif

#if (AMREX_SPACEDIM == 3)
    for (int n = 0; n < ncomp; ++n) {
    for (int k = zlo.z; k <= zhi.z; ++k) {
    for (int j = zlo.y; j <= zhi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = zlo.x; i <= zhi.x; ++i) {
        fz(i,j,k,fcomp+n) = 0.5*(cc(i,j,k-1,ccomp+n) + cc(i,j,k,ccomp+n));
    }
    }
    }
    }
#endif

}

void AverageCCToFace(const MultiFab& cc, int cc_comp,
                     std::array<MultiFab, AMREX_SPACEDIM>& face, int face_comp,
                     int ncomp)
{

    BL_PROFILE_VAR("AverageCCToFace()",AverageCCToFace);

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(cc); mfi.isValid(); ++mfi) {

        const Box & bx = mfi.validbox();

        const Array4<Real const> & cc_fab = cc.array(mfi);

        AMREX_D_TERM(const Array4<Real> & facex_fab = face[0].array(mfi);,
                     const Array4<Real> & facey_fab = face[1].array(mfi);,
                     const Array4<Real> & facez_fab = face[2].array(mfi););

        AMREX_D_TERM(const Box & bx_x = mfi.nodaltilebox(0);,
                     const Box & bx_y = mfi.nodaltilebox(1);,
                     const Box & bx_z = mfi.nodaltilebox(2););

        const Box& index_bounds = amrex::getIndexBounds(AMREX_D_DECL(bx_x,bx_y,bx_z));

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA(index_bounds, tbx,
        {
            avg_cc_to_fc(tbx, AMREX_D_DECL(bx_x,bx_y,bx_z),
                         AMREX_D_DECL(facex_fab,facey_fab,facez_fab), cc_fab,
                         face_comp, cc_comp, ncomp);

        });

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
    for (MFIter mfi(cc); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        Array4<Real const> const& face_fab = face.array(mfi);

        Array4<Real> const& cc_fab = cc.array(mfi);

        AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, n,
        {
            cc_fab(i,j,k,cc_comp+n) = face_fab(i,j,k,face_comp+n);
        });
    }
}
