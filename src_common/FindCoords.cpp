#include "common_functions.H"
#include "common_functions_F.H"


void FindFaceCoords(std::array<MultiFab, AMREX_SPACEDIM> & RealFaceCoords,
                    const Geometry & geom)
{
    BL_PROFILE_VAR("FindFaceCoords()", FindFaceCoords);   

    const RealBox & realDomain = geom.ProbDomain();
    const Real * dx            = geom.CellSize();

    for (MFIter mfi(RealFaceCoords[0]); mfi.isValid(); ++mfi) {

        find_face_coords(ZFILL(realDomain.lo()), ZFILL(realDomain.hi()),
                         BL_TO_FORTRAN_3D(RealFaceCoords[0][mfi]),
                         BL_TO_FORTRAN_3D(RealFaceCoords[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                         BL_TO_FORTRAN_3D(RealFaceCoords[2][mfi]),
#endif
                         ZFILL(dx) );

    }

}

void FindCenterCoords(MultiFab & RealCenterCoords, const Geometry & geom)
{
    BL_PROFILE_VAR("FindCenterCoords()", FindCenterCoords);  

    const RealBox & realDomain = geom.ProbDomain();
    const Real * dx            = geom.CellSize();

    for (MFIter mfi(RealCenterCoords); mfi.isValid(); ++mfi) {
        
        find_center_coords(ZFILL(realDomain.lo()), ZFILL(realDomain.hi()),
                           BL_TO_FORTRAN_3D(RealCenterCoords[mfi]),
                           ZFILL(dx) );

    }

}
