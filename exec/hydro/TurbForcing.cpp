#include "TurbForcing.H"

// initialize n_rngs, geom
// build MultiFabs to hold random numbers
TurbForcing::TurbForcing(BoxArray ba_in, DistributionMapping dmap_in, Geometry geom_in,
                         const Real& a_in, const Real& b_in)
{

    BL_PROFILE_VAR("TurbForcing()",TurbForcing);

    forcing_a = a_in;
    forcing_b = b_in;
    
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        sines  [i].define(convert(ba_in,nodal_flag_dir[i]), dmap_in, 22, 0);
        cosines[i].define(convert(ba_in,nodal_flag_dir[i]), dmap_in, 22, 0);
    }

    Real pi = 3.1415926535897932;
    Real L = prob_hi[0] - prob_lo[0];

    if (prob_hi[0] - prob_lo[0] !=  prob_hi[1] - prob_lo[1]) {
        Abort("TurbForce requires square domain for now");
    }
#if (AMREX_SPACEDIM == 3)
    if (prob_hi[0] - prob_lo[0] !=  prob_hi[2] - prob_lo[2]) {
        Abort("TurbForce requires square domain for now");
    }
#endif

    const GpuArray<Real,AMREX_SPACEDIM> dx = geom_in.CellSizeArray();
    
    // Loop over boxes
    for (MFIter mfi(sines[0],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        AMREX_D_TERM(const Array4<Real> & sin_x = sines[0].array(mfi);,
                     const Array4<Real> & sin_y = sines[1].array(mfi);,
                     const Array4<Real> & sin_z = sines[2].array(mfi););
        
        AMREX_D_TERM(const Array4<Real> & cos_x = cosines[0].array(mfi);,
                     const Array4<Real> & cos_y = cosines[1].array(mfi);,
                     const Array4<Real> & cos_z = cosines[2].array(mfi););

        // since the MFIter is built on a nodal MultiFab we need to build the
        // nodal tileboxes for each direction in this way
        AMREX_D_TERM(Box bx_x = mfi.tilebox(nodal_flag_x);,
                     Box bx_y = mfi.tilebox(nodal_flag_y);,
                     Box bx_z = mfi.tilebox(nodal_flag_z););
    
#if (AMREX_SPACEDIM == 2)
        amrex::ParallelFor(bx_x, bx_y, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real x = prob_lo[0] + i*dx[0];
                Real y = prob_lo[1] + (j+0.5)*dx[1];
                for (int d=0; d<22; ++d) {
                    sin_x(i,j,k,d) = sin(2.*pi*(kx[d]*x + ky[d]*y) / L);
                    cos_x(i,j,k,d) = cos(2.*pi*(kx[d]*x + ky[d]*y) / L);
                }
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real x = prob_lo[0] + (i+0.5)*dx[0];
                Real y = prob_lo[1] + j*dx[1];
                for (int d=0; d<22; ++d) {
                    sin_y(i,j,k,d) = sin(2.*pi*(kx[d]*x + ky[d]*y) / L);
                    cos_y(i,j,k,d) = cos(2.*pi*(kx[d]*x + ky[d]*y) / L);
                }
            });
#elif (AMREX_SPACEDIM ==3)
        amrex::ParallelFor(bx_x, bx_y, bx_z, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real x = prob_lo[0] + i*dx[0];
                Real y = prob_lo[1] + (j+0.5)*dx[1];
                Real z = prob_lo[2] + (k+0.5)*dx[2];
                for (int d=0; d<22; ++d) {
                    sin_x(i,j,k,d) = sin(2.*pi*(kx[d]*x + ky[d]*y + kz[d]*z) / L);
                    cos_x(i,j,k,d) = cos(2.*pi*(kx[d]*x + ky[d]*y + kz[d]*z) / L);
                }
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real x = prob_lo[0] + (i+0.5)*dx[0];
                Real y = prob_lo[1] + j*dx[1];
                Real z = prob_lo[2] + (k+0.5)*dx[2];
                for (int d=0; d<22; ++d) {
                    sin_y(i,j,k,d) = sin(2.*pi*(kx[d]*x + ky[d]*y + kz[d]*z) / L);
                    cos_y(i,j,k,d) = cos(2.*pi*(kx[d]*x + ky[d]*y + kz[d]*z) / L);
                }
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real x = prob_lo[0] + (i+0.5)*dx[0];
                Real y = prob_lo[1] + (j+0.5)*dx[1];
                Real z = prob_lo[2] + k*dx[2];
                for (int d=0; d<22; ++d) {
                    sin_z(i,j,k,d) = sin(2.*pi*(kx[d]*x + ky[d]*y + kz[d]*z) / L);
                    cos_z(i,j,k,d) = cos(2.*pi*(kx[d]*x + ky[d]*y + kz[d]*z) / L);
                }
            });
#endif
    }
    
}

void TurbForcing::AddTurbForcing(std::array< MultiFab, AMREX_SPACEDIM >& gmres_rhs_u,
                                 const Real& dt,
                                 const int& reset_rng)
{    
    // constants for OU process
    const Real a = 0.1;
    const Real b = 0.1;

    if (reset_rng == 1) {
        Vector<Real> rngs_tmp;
        rngs_tmp.resize(132);

        // compute random numbers on IOProcessor
        if (ParallelDescriptor::IOProcessor()) {
            for (int i=0; i<132; ++i) {
                rngs_tmp[i] = amrex::RandomNormal(0.,1.);
            }
        }

        // broadcast random numbers fo all processors
        amrex::BroadcastArray(rngs_tmp,
                              ParallelDescriptor::MyProc(),
                              ParallelDescriptor::IOProcessorNumber(),
                              ParallelDescriptor::Communicator());

        // copy random numbers into GpuArray
        for (int i=0; i<132; ++i) {
            rngs[i] = rngs_tmp[i];
        }
    }

    // Loop over boxes
    for (MFIter mfi(sines[0],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        AMREX_D_TERM(const Array4<Real> & sin_x = sines[0].array(mfi);,
                     const Array4<Real> & sin_y = sines[1].array(mfi);,
                     const Array4<Real> & sin_z = sines[2].array(mfi););
        
        AMREX_D_TERM(const Array4<Real> & cos_x = cosines[0].array(mfi);,
                     const Array4<Real> & cos_y = cosines[1].array(mfi);,
                     const Array4<Real> & cos_z = cosines[2].array(mfi););
        
        AMREX_D_TERM(const Array4<Real> & rhs_x = gmres_rhs_u[0].array(mfi);,
                     const Array4<Real> & rhs_y = gmres_rhs_u[1].array(mfi);,
                     const Array4<Real> & rhs_z = gmres_rhs_u[2].array(mfi););

        // since the MFIter is built on a nodal MultiFab we need to build the
        // nodal tileboxes for each direction in this way
        AMREX_D_TERM(Box bx_x = mfi.tilebox(nodal_flag_x);,
                     Box bx_y = mfi.tilebox(nodal_flag_y);,
                     Box bx_z = mfi.tilebox(nodal_flag_z););
    
#if (AMREX_SPACEDIM == 2)
        Warning("2D AddTurbForcing not defined yet");
        amrex::ParallelFor(bx_x, bx_y, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            }
            );
#elif (AMREX_SPACEDIM ==3)
        amrex::ParallelFor(bx_x, bx_y, bx_z, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                for (int d=0; d<22; ++d) {
                    rhs_x(i,j,k) += (-a*dt + b*rngs[6*d+0]) * cos_x(i,j,k,d);
                    rhs_x(i,j,k) += (-a*dt + b*rngs[6*d+1]) * sin_x(i,j,k,d);
                }
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                for (int d=0; d<22; ++d) {
                    rhs_y(i,j,k) += (-a*dt + b*rngs[6*d+2]) * cos_y(i,j,k,d);
                    rhs_y(i,j,k) += (-a*dt + b*rngs[6*d+3]) * sin_y(i,j,k,d);
                }
                
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                for (int d=0; d<22; ++d) {
                    rhs_z(i,j,k) += (-a*dt + b*rngs[6*d+4]) * cos_z(i,j,k,d);
                    rhs_z(i,j,k) += (-a*dt + b*rngs[6*d+5]) * sin_z(i,j,k,d);
                }                
            });
#endif
    }    
}
