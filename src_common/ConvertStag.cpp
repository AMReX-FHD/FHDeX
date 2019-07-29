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

        const Box& bx = mfi.validbox();        
        
        const auto& cc_fab = (&cc)->array(mfi);
        const auto& facex_fab = (&face[0]) -> array(mfi);
        const auto& facey_fab = (&face[1]) -> array(mfi);
        const Box& bx_x = amrex::growHi(bx,0,1);
        const Box& bx_y = amrex::growHi(bx,1,1);
#if (AMREX_SPACEDIM == 3)        
        const auto& facez_fab = (&face[2]) -> array(mfi);
        const Box& bx_z = amrex::growHi(bx,2,1);
#endif

        AMREX_HOST_DEVICE_FOR_4D(bx_x, ncomp, i, j, k, n,
        {
            facex_fab(i,j,k,face_comp+n) = 0.5*(cc_fab(i,j,k,cc_comp+n) + cc_fab(i-1,j,k,cc_comp+n));
        });

        AMREX_HOST_DEVICE_FOR_4D(bx_y, ncomp, i, j, k, n,
        {
            facey_fab(i,j,k,face_comp+n) = 0.5*(cc_fab(i,j,k,cc_comp+n) + cc_fab(i,j-1,k,cc_comp+n));
        });

#if (AMREX_SPACEDIM == 3)
        AMREX_HOST_DEVICE_FOR_4D(bx_z, ncomp, i, j, k, n,
        {
            facez_fab(i,j,k,face_comp+n) = 0.5*(cc_fab(i,j,k,cc_comp+n) + cc_fab(i,j,k-1,cc_comp+n));
        });
#endif

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
    else {
        Abort("ShiftFaceToCC requires a face-centered MultiFab");
    }

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(cc); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        const auto& face_fab = (&face)->array(mfi);
        const auto& cc_fab = (&cc)->array(mfi);

        AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, n,
        {
            cc_fab(i,j,k,cc_comp+n) = face_fab(i,j,k,face_comp+n);
        });
    }
}
