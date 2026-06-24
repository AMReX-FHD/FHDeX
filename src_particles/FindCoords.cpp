#include "ib_functions.H"

void FindFaceCoords(std::array<MultiFab, AMREX_SPACEDIM> & RealFaceCoords,
                    const Geometry & geom)
{
    BL_PROFILE_VAR("FindFaceCoords()", FindFaceCoords);

    const GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for (MFIter mfi(RealFaceCoords[0]); mfi.isValid(); ++mfi) {

        AMREX_D_TERM(const Array4<Real> & xface = (RealFaceCoords[0]).array(mfi);,
                     const Array4<Real> & yface = (RealFaceCoords[1]).array(mfi);,
                     const Array4<Real> & zface = (RealFaceCoords[2]).array(mfi););

        // since the MFIter is built on a nodal MultiFab we need to build the
        // nodal tileboxes for each direction in this way
        AMREX_D_TERM(Box bx_x = mfi.tilebox(nodal_flag_x,RealFaceCoords[0].nGrowVect());,
                     Box bx_y = mfi.tilebox(nodal_flag_y,RealFaceCoords[1].nGrowVect());,
                     Box bx_z = mfi.tilebox(nodal_flag_z,RealFaceCoords[2].nGrowVect()););

#if (AMREX_SPACEDIM == 2)
            amrex::ParallelFor(bx_x,bx_y,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   xface(i,j,k,0) = i*dx[0]+prob_lo[0];
                                   xface(i,j,k,1) = (j+0.5)*dx[1]+prob_lo[1];
                               },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   yface(i,j,k,0) = (i+0.5)*dx[0]+prob_lo[0];
                                   yface(i,j,k,1) = (j)*dx[1]+prob_lo[1];
                               });

#elif (AMREX_SPACEDIM == 3)
            amrex::ParallelFor(bx_x,bx_y,bx_z,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   xface(i,j,k,0) = i*dx[0]+prob_lo[0];
                                   xface(i,j,k,1) = (j+0.5)*dx[1]+prob_lo[1];
                                   xface(i,j,k,2) = (k+0.5)*dx[2]+prob_lo[2];
                               },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   yface(i,j,k,0) = (i+0.5)*dx[0]+prob_lo[0];
                                   yface(i,j,k,1) = (j)*dx[1]+prob_lo[1];
                                   yface(i,j,k,2) = (k+0.5)*dx[2]+prob_lo[2];
                               },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   zface(i,j,k,0) = (i+0.5)*dx[0]+prob_lo[0];
                                   zface(i,j,k,1) = (j+0.5)*dx[1]+prob_lo[1];
                                   zface(i,j,k,2) = (k)*dx[2]+prob_lo[2];
                               });

#endif
    }
}

void FindCenterCoords(MultiFab & RealCenterCoords, const Geometry & geom)
{
    BL_PROFILE_VAR("FindCenterCoords()", FindCenterCoords);

    const GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for (MFIter mfi(RealCenterCoords); mfi.isValid(); ++mfi) {

        const Array4<Real> & centers = RealCenterCoords.array(mfi);

        Box bx = mfi.growntilebox(RealCenterCoords.nGrowVect());

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            centers(i,j,k,0) = (i+0.5)*dx[0]+prob_lo[0];
            centers(i,j,k,1) = (j+0.5)*dx[1]+prob_lo[1];
#if (AMREX_SPACEDIM == 3)
            centers(i,j,k,2) = (k+0.5)*dx[2]+prob_lo[2];
#endif
        });

    }

}
