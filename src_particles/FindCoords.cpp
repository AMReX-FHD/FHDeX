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

void SimpleShear(std::array<MultiFab, AMREX_SPACEDIM> & externalV, const Geometry & geom, const std::array<MultiFab, AMREX_SPACEDIM> & coords, const Real pin_y)
{
    BL_PROFILE_VAR("SimpleShear()", SimpleShear);

    const GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    Real prob_lo_x = prob_lo[0];
    Real prob_lo_y = prob_lo[1];
    Real prob_lo_z = prob_lo[2];


    if (pin_y != prob_lo[1]) {
       prob_lo_y = pin_y;
       //AllPrint() << "(global in) Replace prob_lo_y with " << prob_lo_y << std::endl;
    }
    
    //AllPrint() << "(global out) Replace prob_lo_y with " << pin_y << std::endl;

    for (MFIter mfi(externalV[0]); mfi.isValid(); ++mfi) {

        AMREX_D_TERM(const Array4<Real> & xface = (externalV[0]).array(mfi);,
                     const Array4<Real> & yface = (externalV[1]).array(mfi);,
                     const Array4<Real> & zface = (externalV[2]).array(mfi););

        // since the MFIter is built on a nodal MultiFab we need to build the
        // nodal tileboxes for each direction in this way
        AMREX_D_TERM(Box bx_x = mfi.tilebox(nodal_flag_x,externalV[0].nGrowVect());,
                     Box bx_y = mfi.tilebox(nodal_flag_y,externalV[1].nGrowVect());,
                     Box bx_z = mfi.tilebox(nodal_flag_z,externalV[2].nGrowVect()););
	
        Array4<const Real> const& coords_x = coords[0].array(mfi);
        Array4<const Real> const& coords_y = coords[1].array(mfi);
        Array4<const Real> const& coords_z = coords[2].array(mfi);

#if (AMREX_SPACEDIM == 2)
            amrex::ParallelFor(bx_x,bx_y,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
			           if (coords_x(i,j,k,1) >= prob_lo_y) {
                                   	xface(i,j,k) = (wallspeed_y_hi[0]-wallspeed_y_lo[0])/(prob_hi[1]-prob_lo_y)*(coords_x(i,j,k,1)-prob_lo_y)+wallspeed_y_lo[0];
				   }
				   else {
				   	xface(i,j,k) = 0.0;
				   }
                               },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
			           if (coords_y(i,j,k,0) >= prob_lo_x) {
                                   	yface(i,j,k) = (wallspeed_x_hi[1]-wallspeed_x_lo[1])/(prob_hi[0]-prob_lo_x)*(coords_y(i,j,k,0)-prob_lo_x)+wallspeed_x_lo[1];
				   }
				   else {
				   	yface(i,j,k) = 0.0;
				   }
                               });

#elif (AMREX_SPACEDIM == 3)
            amrex::ParallelFor(bx_x,bx_y,bx_z,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   //AllPrint() << "Replace prob_lo_y with " << prob_lo_y << std::endl;
			           if (coords_x(i,j,k,1) >= prob_lo_y && coords_x(i,j,k,2) >= prob_lo_z) {
                                        //AllPrint() << "xface if 1" << std::endl;
                                   	xface(i,j,k) = (wallspeed_y_hi[0]-wallspeed_y_lo[0])/(prob_hi[1]-prob_lo_y)*(coords_x(i,j,k,1)-prob_lo_y)+wallspeed_y_lo[0] + (wallspeed_z_hi[0]-wallspeed_z_lo[0])/(prob_hi[2]-prob_lo_z)*(coords_x(i,j,k,2)-prob_lo_z)+wallspeed_z_lo[0];
	   			   }
				   else if (coords_x(i,j,k,2) >= prob_lo_z && coords_x(i,j,k,1) < prob_lo_y) {
                                        //AllPrint() << "xface if 2" << std::endl;
                                        xface(i,j,k) = (wallspeed_z_hi[0]-wallspeed_z_lo[0])/(prob_hi[2]-prob_lo_z)*(coords_x(i,j,k,2)-prob_lo_z)+wallspeed_z_lo[0];
				   }
				   else if (coords_x(i,j,k,2) < prob_lo_z && coords_x(i,j,k,1) >= prob_lo_y) {
                                        //AllPrint() << "xface if 3" << std::endl;
                                        xface(i,j,k) = (wallspeed_y_hi[0]-wallspeed_y_lo[0])/(prob_hi[1]-prob_lo_y)*(coords_x(i,j,k,1)-prob_lo_y)+wallspeed_y_lo[0];
				   }
				   else {
                                        //AllPrint() << "xface if 4" << std::endl;
				   	xface(i,j,k) = 0.0;
				   }
                                   //AllPrint() << "coords_x(" << i << "," << j << "," << k << ",0) is: " << coords_x(i,j,k,0) << std::endl;
                                   //AllPrint() << "coords_x(" << i << "," << j << "," << k << ",1) is: " << coords_x(i,j,k,1) << std::endl;
                                   //AllPrint() << "coords_x(" << i << "," << j << "," << k << ",2) is: " << coords_x(i,j,k,2) << std::endl;
                                   //AllPrint() << "umac_x(" << i << "," << j << "," << k << ") is: " << xface(i,j,k) << std::endl;
                               },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
			           if (coords_y(i,j,k,0) >= prob_lo_x && coords_y(i,j,k,2) >= prob_lo_z) {
                                   	yface(i,j,k) = (wallspeed_x_hi[1]-wallspeed_x_lo[1])/(prob_hi[0]-prob_lo_x)*(coords_y(i,j,k,0)-prob_lo_x)+wallspeed_x_lo[1] + (wallspeed_z_hi[1]-wallspeed_z_lo[1])/(prob_hi[2]-prob_lo_z)*(coords_y(i,j,k,2)-prob_lo_z)+wallspeed_z_lo[1]; 
			           }
				   else if (coords_y(i,j,k,2) >= prob_lo_z && coords_y(i,j,k,0) < prob_lo_x) {
					yface(i,j,k) = (wallspeed_z_hi[1]-wallspeed_z_lo[1])/(prob_hi[2]-prob_lo_z)*(coords_y(i,j,k,2)-prob_lo_z)+wallspeed_z_lo[1];
				   }
				   else if (coords_y(i,j,k,2) < prob_lo_z && coords_y(i,j,k,0) >= prob_lo_x) {
				        yface(i,j,k) = (wallspeed_x_hi[1]-wallspeed_x_lo[1])/(prob_hi[0]-prob_lo_x)*(coords_y(i,j,k,0)-prob_lo_x)+wallspeed_x_lo[1];
				   }
				   else {
				   	yface(i,j,k) = 0.0;
				   }
                               },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
			           if (coords_z(i,j,k,0) >= prob_lo_x && coords_z(i,j,k,1) >= prob_lo_y) {
                                   	zface(i,j,k) = (wallspeed_x_hi[2]-wallspeed_x_lo[2])/(prob_hi[0]-prob_lo_x)*(coords_z(i,j,k,0)-prob_lo_x)+wallspeed_x_lo[2] + (wallspeed_y_hi[2]-wallspeed_y_lo[2])/(prob_hi[1]-prob_lo_y)*(coords_z(i,j,k,1)-prob_lo_y)+wallspeed_y_lo[2];
				   }
				   else if (coords_z(i,j,k,1) >= prob_lo_y && coords_z(i,j,k,0) < prob_lo_x) {
				       	zface(i,j,k) = (wallspeed_y_hi[2]-wallspeed_y_lo[2])/(prob_hi[1]-prob_lo_y)*(coords_z(i,j,k,1)-prob_lo_y)+wallspeed_y_lo[2];
				   }
				   else if (coords_z(i,j,k,1) < prob_lo_y && coords_z(i,j,k,0) >= prob_lo_x) {
				        zface(i,j,k) = (wallspeed_x_hi[2]-wallspeed_x_lo[2])/(prob_hi[0]-prob_lo_x)*(coords_z(i,j,k,0)-prob_lo_x)+wallspeed_x_lo[2];
				   }
				   else {
				   	zface(i,j,k) = 0.0;
				   }
                               });

#endif
    }
 
}

void UniformVel(std::array<MultiFab, AMREX_SPACEDIM> & externalV, const Geometry & geom)
{
    BL_PROFILE_VAR("SimpleShear()", SimpleShear);

    const GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for (MFIter mfi(externalV[0]); mfi.isValid(); ++mfi) {

        AMREX_D_TERM(const Array4<Real> & xface = (externalV[0]).array(mfi);,
                     const Array4<Real> & yface = (externalV[1]).array(mfi);,
                     const Array4<Real> & zface = (externalV[2]).array(mfi););

        // since the MFIter is built on a nodal MultiFab we need to build the
        // nodal tileboxes for each direction in this way
        AMREX_D_TERM(Box bx_x = mfi.tilebox(nodal_flag_x,externalV[0].nGrowVect());,
                     Box bx_y = mfi.tilebox(nodal_flag_y,externalV[1].nGrowVect());,
                     Box bx_z = mfi.tilebox(nodal_flag_z,externalV[2].nGrowVect()););
	
#if (AMREX_SPACEDIM == 2)
            amrex::ParallelFor(bx_x,bx_y,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   	xface(i,j,k) = wallspeed_y_hi[0]-wallspeed_y_lo[0];
                               },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   	yface(i,j,k) = wallspeed_x_hi[1]-wallspeed_x_lo[1];
                               });

#elif (AMREX_SPACEDIM == 3)
            amrex::ParallelFor(bx_x,bx_y,bx_z,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   	xface(i,j,k) = wallspeed_y_hi[0]-wallspeed_y_lo[0] + wallspeed_z_hi[0]-wallspeed_z_lo[0];
                               },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   	yface(i,j,k) = wallspeed_x_hi[1]-wallspeed_x_lo[1] + wallspeed_z_hi[1]-wallspeed_z_lo[1]; 
                               },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                   	zface(i,j,k) = wallspeed_x_hi[2]-wallspeed_x_lo[2] + wallspeed_y_hi[2]-wallspeed_y_lo[2];
                               });

#endif
    }
 
}
