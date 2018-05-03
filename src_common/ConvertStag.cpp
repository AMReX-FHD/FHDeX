#include "common_functions.H"
#include "common_functions_F.H"

void AverageFaceToCC(MultiFab& face, int face_comp,
                     MultiFab& cc, int cc_comp,
                     int ncomp) 
{
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


    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for ( MFIter mfi(cc); mfi.isValid(); ++mfi) {

        const Box& validBox = mfi.validbox();

        average_face_to_cc(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                           BL_TO_FORTRAN_FAB(face[mfi]),
                           BL_TO_FORTRAN_FAB(cc[mfi]),
                           &face_comp, &cc_comp, &ncomp, &av_dim);


    }

}
