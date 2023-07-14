#include "TurbForcingComp.H"

// initialize n_rngs, geom
// build MultiFabs to hold random numbers
TurbForcingComp::TurbForcingComp()
{}

void TurbForcingComp::define(BoxArray ba_in, DistributionMapping dmap_in,
                    const Real& a_in, const Real& b_in, const Real& c_in, const Real& d_in, const Real& alpha_in)
{
    BL_PROFILE_VAR("TurbForcingComp::define()",TurbForcingCompDefine);

    for (int i=0; i<132; ++i) {
        forcing.forcing_S[i] = 0.;
        forcing.forcing_C[i] = 0.;
    }
    
    forcing_a = a_in;
    forcing_b = b_in;
    forcing_c = c_in;
    forcing_d = d_in;
    alpha     = alpha_in;
    
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        sines  [i].define(convert(ba_in,nodal_flag_dir[i]), dmap_in, 22, 0);
        cosines[i].define(convert(ba_in,nodal_flag_dir[i]), dmap_in, 22, 0);
    }
}

void TurbForcingComp::Initialize(const Geometry& geom_in) {

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
    
    KVEC kvec_gpu;
    // copy the object to the device and pass a device pointer
    KVEC* kvec_gpu_ptr = (KVEC*)The_Arena()->alloc(sizeof(KVEC));
    Gpu::htod_memcpy_async(kvec_gpu_ptr, &kvec_gpu, sizeof(KVEC));
    Gpu::streamSynchronize();
    
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
                    sin_x(i,j,k,d) = std::sin(2.*pi*(kvec.kx[d]*x + kvec.ky[d]*y) / L);
                    cos_x(i,j,k,d) = std::cos(2.*pi*(kvec.kx[d]*x + kvec.ky[d]*y) / L);
                }
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real x = prob_lo_gpu[0] + (i+0.5)*dx[0];
                Real y = prob_lo_gpu[1] + j*dx[1];
                for (int d=0; d<22; ++d) {
                    sin_y(i,j,k,d) = std::sin(2.*pi*(kvec.kx[d]*x + kvec.ky[d]*y) / L);
                    cos_y(i,j,k,d) = std::cos(2.*pi*(kvec.kx[d]*x + kvec.ky[d]*y) / L);
                }
            });
#elif (AMREX_SPACEDIM ==3)
        amrex::ParallelFor(bx_x, bx_y, bx_z, [prob_lo_gpu,kvec_gpu_ptr,L,sin_x,cos_x,dx] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real pi = 3.1415926535897932;
                Real x = prob_lo_gpu[0] + i*dx[0];
                Real y = prob_lo_gpu[1] + (j+0.5)*dx[1];
                Real z = prob_lo_gpu[2] + (k+0.5)*dx[2];
                for (int d=0; d<22; ++d) {
                    sin_x(i,j,k,d) = std::sin(2.*pi*(kvec_gpu_ptr->kx[d]*x + kvec_gpu_ptr->ky[d]*y + kvec_gpu_ptr->kz[d]*z) / L);
                    cos_x(i,j,k,d) = std::cos(2.*pi*(kvec_gpu_ptr->kx[d]*x + kvec_gpu_ptr->ky[d]*y + kvec_gpu_ptr->kz[d]*z) / L);
                }
            },
            [prob_lo_gpu,kvec_gpu_ptr,L,sin_y,cos_y,dx] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real pi = 3.1415926535897932;
                Real x = prob_lo_gpu[0] + (i+0.5)*dx[0];
                Real y = prob_lo_gpu[1] + j*dx[1];
                Real z = prob_lo_gpu[2] + (k+0.5)*dx[2];
                for (int d=0; d<22; ++d) {
                    sin_y(i,j,k,d) = std::sin(2.*pi*(kvec_gpu_ptr->kx[d]*x + kvec_gpu_ptr->ky[d]*y + kvec_gpu_ptr->kz[d]*z) / L);
                    cos_y(i,j,k,d) = std::cos(2.*pi*(kvec_gpu_ptr->kx[d]*x + kvec_gpu_ptr->ky[d]*y + kvec_gpu_ptr->kz[d]*z) / L);
                }
            },
            [prob_lo_gpu,kvec_gpu_ptr,L,sin_z,cos_z,dx] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real pi = 3.1415926535897932;
                Real x = prob_lo_gpu[0] + (i+0.5)*dx[0];
                Real y = prob_lo_gpu[1] + (j+0.5)*dx[1];
                Real z = prob_lo_gpu[2] + k*dx[2];
                for (int d=0; d<22; ++d) {
                    sin_z(i,j,k,d) = std::sin(2.*pi*(kvec_gpu_ptr->kx[d]*x + kvec_gpu_ptr->ky[d]*y + kvec_gpu_ptr->kz[d]*z) / L);
                    cos_z(i,j,k,d) = std::cos(2.*pi*(kvec_gpu_ptr->kx[d]*x + kvec_gpu_ptr->ky[d]*y + kvec_gpu_ptr->kz[d]*z) / L);
                }
            });
#endif
    }
    
}

void TurbForcingComp::CalcTurbForcingComp(std::array< MultiFab, AMREX_SPACEDIM >& vel_f,
                                 const Real& dt,
                                 const int& update)
{

    Real sqrtdt = std::sqrt(dt);

    // update U = U - a*dt + b*sqrt(dt)*Z
    if (update == 1) {
        
        Vector<Real> rngs_s(132); // solenoidal
        Vector<Real> rngs_c(132); // comopressional

        if (ParallelDescriptor::IOProcessor()) {
            // compute random numbers on IOProcessor
            for (int i=0; i<132; ++i) {
                rngs_s[i] = amrex::RandomNormal(0.,1.);
                rngs_c[i] = amrex::RandomNormal(0.,1.);
            }
        }

        // broadcast random numbers to all processors
        amrex::BroadcastArray(rngs_s,
                              ParallelDescriptor::MyProc(),
                              ParallelDescriptor::IOProcessorNumber(),
                              ParallelDescriptor::Communicator());
        amrex::BroadcastArray(rngs_c,
                              ParallelDescriptor::MyProc(),
                              ParallelDescriptor::IOProcessorNumber(),
                              ParallelDescriptor::Communicator());

        // update forcing (OU)
        for (int i=0; i<132; ++i) {
            forcing.forcing_S[i] += forcing_a*forcing.forcing_S[i]*dt + forcing_b*sqrtdt*rngs_s[i];
            forcing.forcing_C[i] += forcing_c*forcing.forcing_C[i]*dt + forcing_d*sqrtdt*rngs_c[i];
        }

    }

    
    KVEC kvec_gpu;
    // copy the object to the device and pass a device pointer
    KVEC* kvec_gpu_ptr = (KVEC*)The_Arena()->alloc(sizeof(KVEC));
    Gpu::htod_memcpy_async(kvec_gpu_ptr, &kvec_gpu, sizeof(KVEC));
    Gpu::streamSynchronize();
    
    FORCING forcing_gpu;    
    // copy the object to the device and pass a device pointer
    FORCING* forcing_gpu_ptr = (FORCING*)The_Arena()->alloc(sizeof(FORCING));
    Gpu::htod_memcpy_async(forcing_gpu_ptr, &forcing_gpu, sizeof(FORCING));
    Gpu::streamSynchronize();
    
    for (int i=0; i<132; ++i) {
        forcing_gpu_ptr->forcing_S[i] = forcing.forcing_S[i];
        forcing_gpu_ptr->forcing_C[i] = forcing.forcing_C[i];
    }

    Real alpha_gpu = alpha;

    // Loop over boxes
    for (MFIter mfi(sines[0],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        AMREX_D_TERM(const Array4<Real> & sin_x = sines[0].array(mfi);,
                     const Array4<Real> & sin_y = sines[1].array(mfi);,
                     const Array4<Real> & sin_z = sines[2].array(mfi););
        
        AMREX_D_TERM(const Array4<Real> & cos_x = cosines[0].array(mfi);,
                     const Array4<Real> & cos_y = cosines[1].array(mfi);,
                     const Array4<Real> & cos_z = cosines[2].array(mfi););
        
        AMREX_D_TERM(const Array4<Real> & vel_x = vel_f[0].array(mfi);,
                     const Array4<Real> & vel_y = vel_f[1].array(mfi);,
                     const Array4<Real> & vel_z = vel_f[2].array(mfi););

        // since the MFIter is built on a nodal MultiFab we need to build the
        // nodal tileboxes for each direction in this way
        AMREX_D_TERM(Box bx_x = mfi.tilebox(nodal_flag_x);,
                     Box bx_y = mfi.tilebox(nodal_flag_y);,
                     Box bx_z = mfi.tilebox(nodal_flag_z););
    
#if (AMREX_SPACEDIM == 2)
        Warning("2D CalcTurbForcingComp not defined yet");
        amrex::ParallelFor(bx_x, bx_y, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            }
            );
#elif (AMREX_SPACEDIM ==3)
        amrex::ParallelFor(bx_x, bx_y, bx_z, [vel_x,kvec_gpu_ptr,forcing_gpu_ptr,alpha_gpu,sin_x,cos_x] AMREX_GPU_DEVICE (int i, int j, int k) {
                for (int d=0; d<22; ++d) {
                    Real kx = kvec_gpu_ptr->kx[d];
                    Real ky = kvec_gpu_ptr->ky[d];
                    Real kz = kvec_gpu_ptr->kz[d];
                    Real kk = kx*kx + ky*ky + kz*kz;
                    vel_x(i,j,k) = alpha_gpu*cos_x(i,j,k,d)*
                        (forcing_gpu_ptr->forcing_S[d]*(1.0-kx*kx/kk)    - forcing_gpu_ptr->forcing_S[d+22]*(kx*ky/kk)     - forcing_gpu_ptr->forcing_S[d+44]*(kx*kz/kk)); // solenoidal
                    vel_x(i,j,k) += alpha_gpu*sin_x(i,j,k,d)*
                        (forcing_gpu_ptr->forcing_S[d+66]*(1.0-kx*kx/kk) - forcing_gpu_ptr->forcing_S[d+88]*(kx*ky/kk)     - forcing_gpu_ptr->forcing_S[d+110]*(kx*kz/kk)); // solenoidal
                    vel_x(i,j,k) += (1.0-alpha_gpu)*cos_x(i,j,k,d)*
                        (forcing_gpu_ptr->forcing_C[d]*(kx*kx/kk)        + forcing_gpu_ptr->forcing_C[d+22]*(kx*ky/kk)     + forcing_gpu_ptr->forcing_C[d+44]*(kx*kz/kk)); // compressional
                    vel_x(i,j,k) += (1.0-alpha_gpu)*sin_x(i,j,k,d)*
                        (forcing_gpu_ptr->forcing_C[d+66]*(kx*kx/kk)     + forcing_gpu_ptr->forcing_C[d+88]*(kx*ky/kk)     + forcing_gpu_ptr->forcing_C[d+110]*(kx*kz/kk)); // compressional
                }
            },
            [vel_y,kvec_gpu_ptr,forcing_gpu_ptr,alpha_gpu,sin_y,cos_y] AMREX_GPU_DEVICE (int i, int j, int k) {
                for (int d=0; d<22; ++d) {
                    Real kx = kvec_gpu_ptr->kx[d];
                    Real ky = kvec_gpu_ptr->ky[d];
                    Real kz = kvec_gpu_ptr->kz[d];
                    Real kk = kx*kx + ky*ky + kz*kz;
                    vel_y(i,j,k) = alpha_gpu*cos_y(i,j,k,d)*
                        (-forcing_gpu_ptr->forcing_S[d]*(kx*ky/kk)       + forcing_gpu_ptr->forcing_S[d+22]*(1.0-ky*ky/kk) - forcing_gpu_ptr->forcing_S[d+44]*(ky*kz/kk)); // solenoidal
                    vel_y(i,j,k) += alpha_gpu*sin_y(i,j,k,d)*
                        (-forcing_gpu_ptr->forcing_S[d+66]*(kx*ky/kk)    + forcing_gpu_ptr->forcing_S[d+88]*(1.0-ky*ky/kk) - forcing_gpu_ptr->forcing_S[d+110]*(ky*kz/kk)); // solenoidal
                    vel_y(i,j,k) += (1.0-alpha_gpu)*cos_y(i,j,k,d)*
                        (forcing_gpu_ptr->forcing_C[d]*(kx*ky/kk)        + forcing_gpu_ptr->forcing_C[d+22]*(ky*ky/kk)     + forcing_gpu_ptr->forcing_C[d+44]*(ky*kz/kk)); // compressional
                    vel_y(i,j,k) += (1.0-alpha_gpu)*sin_y(i,j,k,d)*
                        (forcing_gpu_ptr->forcing_C[d+66]*(kx*ky/kk)     + forcing_gpu_ptr->forcing_C[d+88]*(ky*ky/kk)     + forcing_gpu_ptr->forcing_C[d+110]*(ky*kz/kk)); // compressional
                }
                
            },
            [vel_z,kvec_gpu_ptr,forcing_gpu_ptr,alpha_gpu,sin_z,cos_z] AMREX_GPU_DEVICE (int i, int j, int k) {
                for (int d=0; d<22; ++d) {
                    Real kx = kvec_gpu_ptr->kx[d];
                    Real ky = kvec_gpu_ptr->ky[d];
                    Real kz = kvec_gpu_ptr->kz[d];
                    Real kk = kx*kx + ky*ky + kz*kz;
                    vel_z(i,j,k) = alpha_gpu*cos_z(i,j,k,d)*
                        (-forcing_gpu_ptr->forcing_S[d]*(kx*kz/kk)       - forcing_gpu_ptr->forcing_S[d+22]*(ky*kz/kk)     + forcing_gpu_ptr->forcing_S[d+44]*(1.0-kz*kz/kk)); // solenoidal
                    vel_z(i,j,k) += alpha_gpu*sin_z(i,j,k,d)*
                        (-forcing_gpu_ptr->forcing_S[d+66]*(kx*kz/kk)    - forcing_gpu_ptr->forcing_S[d+88]*(ky*kz/kk)     + forcing_gpu_ptr->forcing_S[d+110]*(1.0-kz*kz/kk)); // solenoidal
                    vel_z(i,j,k) += (1.0-alpha_gpu)*cos_z(i,j,k,d)*
                        (forcing_gpu_ptr->forcing_C[d]*(kx*kz/kk)        + forcing_gpu_ptr->forcing_C[d+22]*(ky*kz/kk)     + forcing_gpu_ptr->forcing_C[d+44]*(kz*kz/kk)    ); // compressional
                    vel_z(i,j,k) += (1.0-alpha_gpu)*sin_z(i,j,k,d)*
                        (forcing_gpu_ptr->forcing_C[d+66]*(kx*kz/kk)     + forcing_gpu_ptr->forcing_C[d+88]*(ky*kz/kk)     + forcing_gpu_ptr->forcing_C[d+110]*(kz*kz/kk)    ); // compressional
                }                
            });
#endif
    }    
}

std::tuple<amrex::Real, Real> TurbForcingComp::getU(const int& i) {
    return {forcing.forcing_S[i], forcing.forcing_C[i]};
}

void TurbForcingComp::setU(const int& i, Real fs, Real fc) {
    forcing.forcing_S[i] = fs;
    forcing.forcing_C[i] = fc;
}
