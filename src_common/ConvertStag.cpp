#include "common_functions.H"
#include "common_functions_F.H"

void AverageFaceToCC(const MultiFab& face, int face_comp,
                     MultiFab& cc, int cc_comp,
                     int ncomp)
{

    BL_PROFILE_VAR("AverageFaceToCC()",AverageFaceToCC);

    int av_dim;  // along which dimension to do the average

    if (face.is_nodal(0)) {
        av_dim = 0;
    }
    else if (face.is_nodal(1)) {
        av_dim = 1;
    }
    else if (face.is_nodal(2)) {
        av_dim = 2;
    }
    else {
        Abort("AverageFaceToCC requires a face-centered MultiFab");
    }
    
    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(cc); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        const auto& face_fab = (&face)->array(mfi);
        const auto& cc_fab = (&cc)->array(mfi);

        if (av_dim == 0) {
            AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, n,
            {
                cc_fab(i,j,k,cc_comp+n) = 0.5*(face_fab(i+1,j,k,face_comp+n) + face_fab(i,j,k,face_comp+n));
            });
        }
        else  if (av_dim == 1) {
            AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, n,
            {
                cc_fab(i,j,k,cc_comp+n) = 0.5*(face_fab(i,j+1,k,face_comp+n) + face_fab(i,j,k,face_comp+n));
            });
        }
        else {
            AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, n,
            {
                cc_fab(i,j,k,cc_comp+n) = 0.5*(face_fab(i,j,k+1,face_comp+n) + face_fab(i,j,k,face_comp+n));
            });
        }
    }
}


void AverageCCToFace(const MultiFab& cc, int cc_comp,
                     std::array<MultiFab, AMREX_SPACEDIM>& face, int face_comp,
                     int ncomp)
{

    BL_PROFILE_VAR("AverageCCToFace()",AverageCCToFace);

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(cc); mfi.isValid(); ++mfi) {

        const Box& validBox = mfi.validbox();

        average_cc_to_face(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                           BL_TO_FORTRAN_FAB(cc[mfi]),
                           BL_TO_FORTRAN_FAB(face[0][mfi]),
                           BL_TO_FORTRAN_FAB(face[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                           BL_TO_FORTRAN_FAB(face[2][mfi]),
#endif
                           &cc_comp, &face_comp, &ncomp);

    }

}

void ShiftFaceToCC(const MultiFab& face, int face_comp,
                     MultiFab& cc, int cc_comp,
                     int ncomp)
{

    BL_PROFILE_VAR("ShiftFaceToCC()",ShiftFaceToCC);

    int av_dim;  // along which dimension to do the shift

    if (face.is_nodal(0)) {
        av_dim = 0;
    }
    else if (face.is_nodal(1)) {
        av_dim = 1;
    }
    else if (face.is_nodal(2)) {
        av_dim = 2;
    }

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(cc); mfi.isValid(); ++mfi) {

        const Box& validBox = mfi.validbox();

        shift_face_to_cc(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                         BL_TO_FORTRAN_FAB(face[mfi]),
                         BL_TO_FORTRAN_FAB(cc[mfi]),
                         &face_comp, &cc_comp, &ncomp, &av_dim);
    }
}
