#include "TurbForcing.H"

// initialize n_rngs, geom
// build MultiFabs to hold random numbers
TurbForcing::TurbForcing()
{}

void TurbForcing::define(BoxArray ba_in, DistributionMapping dmap_in,
                    const Real& a_in, const Real& b_in)
{
    BL_PROFILE_VAR("TurbForcing::define()",TurbForcingDefine);

    for (int i=0; i<132; ++i) {
        forcing_U[i] = 0.;
    }

    forcing_a = a_in;
    forcing_b = b_in;

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        sines  [i].define(convert(ba_in,nodal_flag_dir[i]), dmap_in, 22, 0);
        cosines[i].define(convert(ba_in,nodal_flag_dir[i]), dmap_in, 22, 0);
    }
}

void TurbForcing::Initialize(const Geometry& geom_in) {

    GpuArray<int,22> kx{1, 0, 0, 1, 1, 0, 1, 2, 0, 0, 2, 2, 1, 0, 1, 0, 2, 1, 1, 2, 2, 0};
    GpuArray<int,22> ky{0, 1, 0, 1, 0, 1, 1, 0, 2, 0, 1, 0, 2, 2, 0, 1, 1, 2, 1, 2, 0, 2};
    GpuArray<int,22> kz{0, 0, 1, 0, 1, 1, 1, 0, 0, 2, 0, 1, 0, 1, 2, 2, 1, 1, 2, 0, 2, 2};

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

    GpuArray<Real,AMREX_SPACEDIM> prob_lo_gpu;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        prob_lo_gpu[d] = prob_lo[d];
    }

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
                Real x = prob_lo_gpu[0] + i*dx[0];
                Real y = prob_lo_gpu[1] + (j+0.5)*dx[1];
                for (int d=0; d<22; ++d) {
                    sin_x(i,j,k,d) = std::sin(2.*pi*(kx[d]*x + ky[d]*y) / L);
                    cos_x(i,j,k,d) = std::cos(2.*pi*(kx[d]*x + ky[d]*y) / L);
                }
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real x = prob_lo_gpu[0] + (i+0.5)*dx[0];
                Real y = prob_lo_gpu[1] + j*dx[1];
                for (int d=0; d<22; ++d) {
                    sin_y(i,j,k,d) = std::sin(2.*pi*(kx[d]*x + ky[d]*y) / L);
                    cos_y(i,j,k,d) = std::cos(2.*pi*(kx[d]*x + ky[d]*y) / L);
                }
            });
#elif (AMREX_SPACEDIM ==3)
        amrex::ParallelFor(bx_x, bx_y, bx_z, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real x = prob_lo_gpu[0] + i*dx[0];
                Real y = prob_lo_gpu[1] + (j+0.5)*dx[1];
                Real z = prob_lo_gpu[2] + (k+0.5)*dx[2];
                for (int d=0; d<22; ++d) {
                    sin_x(i,j,k,d) = std::sin(2.*pi*(kx[d]*x + ky[d]*y + kz[d]*z) / L);
                    cos_x(i,j,k,d) = std::cos(2.*pi*(kx[d]*x + ky[d]*y + kz[d]*z) / L);
                }
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real x = prob_lo_gpu[0] + (i+0.5)*dx[0];
                Real y = prob_lo_gpu[1] + j*dx[1];
                Real z = prob_lo_gpu[2] + (k+0.5)*dx[2];
                for (int d=0; d<22; ++d) {
                    sin_y(i,j,k,d) = std::sin(2.*pi*(kx[d]*x + ky[d]*y + kz[d]*z) / L);
                    cos_y(i,j,k,d) = std::cos(2.*pi*(kx[d]*x + ky[d]*y + kz[d]*z) / L);
                }
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real x = prob_lo_gpu[0] + (i+0.5)*dx[0];
                Real y = prob_lo_gpu[1] + (j+0.5)*dx[1];
                Real z = prob_lo_gpu[2] + k*dx[2];
                for (int d=0; d<22; ++d) {
                    sin_z(i,j,k,d) = std::sin(2.*pi*(kx[d]*x + ky[d]*y + kz[d]*z) / L);
                    cos_z(i,j,k,d) = std::cos(2.*pi*(kx[d]*x + ky[d]*y + kz[d]*z) / L);
                }
            });
#endif
    }

}

void TurbForcing::AddTurbForcing(std::array< MultiFab, AMREX_SPACEDIM >& gmres_rhs_u,
                                 const Real& dt,
                                 const int& update_U)
{

    Real sqrtdt = std::sqrt(dt);

    // update U = U - a*dt + b*sqrt(dt)*Z
    if (update_U == 1) {

        Vector<Real> rngs(132);

        if (ParallelDescriptor::IOProcessor()) {
            // compute random numbers on IOProcessor
            for (int i=0; i<132; ++i) {
                rngs[i] = amrex::RandomNormal(0.,1.);
            }
        }

        // broadcast random numbers to all processors
        amrex::BroadcastArray(rngs,
                              ParallelDescriptor::MyProc(),
                              ParallelDescriptor::IOProcessorNumber(),
                              ParallelDescriptor::Communicator());

        // update forcing_U
        for (int i=0; i<132; ++i) {
            forcing_U[i] += -forcing_a*forcing_U[i]*dt + forcing_b*sqrtdt*rngs[i];
        }
    }

    GpuArray<Real,132> forcing_U_gpu;
    for (int i=0; i<132; ++i) {
        forcing_U_gpu[i] = forcing_U[i];
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
                    rhs_x(i,j,k) += forcing_U_gpu[d] * cos_x(i,j,k,d);
                    rhs_x(i,j,k) += forcing_U_gpu[d+22] * sin_x(i,j,k,d);
                }
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                for (int d=0; d<22; ++d) {
                    rhs_y(i,j,k) += forcing_U_gpu[d+44] * cos_y(i,j,k,d);
                    rhs_y(i,j,k) += forcing_U_gpu[d+66] * sin_y(i,j,k,d);
                }

            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                for (int d=0; d<22; ++d) {
                    rhs_z(i,j,k) += forcing_U_gpu[d+88] * cos_z(i,j,k,d);
                    rhs_z(i,j,k) += forcing_U_gpu[d+110] * sin_z(i,j,k,d);
                }
            });
#endif
    }
}

Real TurbForcing::getU(const int& i) {
    return forcing_U[i];
}

void TurbForcing::setU(const int& i, Real x) {
    forcing_U[i] = x;
    return;
}
