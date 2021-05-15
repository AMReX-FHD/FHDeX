#include "hydro_test_functions.H"

using namespace amrex;

void InitVel(std::array< MultiFab, AMREX_SPACEDIM >& umac,
             const Geometry& geom) {

    Real zshft = (AMREX_SPACEDIM == 2) ? 0. : 0.5;

    GpuArray<Real,AMREX_SPACEDIM> reallo = geom.ProbLoArray();
    GpuArray<Real,AMREX_SPACEDIM> realhi = geom.ProbHiArray();
    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    GpuArray<Real,AMREX_SPACEDIM> center;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        center[d] = 0.5*(realhi[d]-reallo[d]);
    }

    // IC parameters
    Real L_hlf = center[0];
    // k1 & k2 determine steepness of velocity profile:
    Real k1 = 0.1*L_hlf;
    // k1 = 1d-6*L_hlf
    Real k2 = k1;
    Real k1_inv = 1./k1;
    Real k2_inv = 1./k2;

    // Vortex:
    // [r_a r_b] defines radial bounds of velocity bump:
    Real r_a = 0.35*L_hlf;
    Real r_b = L_hlf - r_a;

    // Kelvin-Helmholtz:
    Real pi = 3.1415926535897932;
    Real freq = 02.*pi/L_hlf;
    Real amp = 2.0e-3*L_hlf;
    // amp = 2.0d-1*L_hlf
    Real width1 = L_hlf/2.;

    if (prob_type == 0) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            umac[d].setVal(0.);
        }
        return;
    }
    
    for ( MFIter mfi(umac[0]); mfi.isValid(); ++mfi ) {

        AMREX_D_TERM(const Array4<Real> & u = (umac[0]).array(mfi);,
                     const Array4<Real> & v = (umac[1]).array(mfi);,
                     const Array4<Real> & w = (umac[2]).array(mfi););

        // since the MFIter is built on a nodal MultiFab we need to build the
        // nodal tileboxes for each direction in this way
        AMREX_D_TERM(Box bx_x = mfi.tilebox(nodal_flag_x);,
                     Box bx_y = mfi.tilebox(nodal_flag_y);,
                     Box bx_z = mfi.tilebox(nodal_flag_z););

        
        amrex::ParallelFor(bx_x,
                           bx_y,
#if (AMREX_SPACEDIM == 3)
                           bx_z,
#endif
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            GpuArray<Real,AMREX_SPACEDIM> itVec;

            AMREX_D_TERM(itVec[0] = i*dx[0];,
                         itVec[1] = (j+0.5)*dx[1];,
                         itVec[2] = (k+zshft)*dx[2];);

            GpuArray<Real,AMREX_SPACEDIM> relpos;
            Real rad2 = 0.;
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                relpos[d] = reallo[d] + itVec[d] - center[d];
            }

            // only sum first two velocities
            for (int d=0; d<2; ++d) {
                rad2 += relpos[d]*relpos[d];
            }

            Real rad = std::sqrt(rad2);
            
            // note prob_type == 0 handled above
            if (prob_type == 1) {

                // Multiply velocity magnitude by sin(theta)
                u(i,j,k) = 0.25*(1.+std::tanh(k1_inv*(rad-r_a)))*(1.+std::tanh(k2_inv*(r_b-rad)))
                    *(relpos[1]/rad);

            } else if (prob_type == 2) {

                Real perturb = amp*sin(freq*relpos[0]);
                u(i,j,k) = 0.25*(1.+std::tanh(k1_inv*(relpos[1] - (-width1/2.+perturb))))
                    *(1.+std::tanh(k2_inv*((width1/2.+perturb) - relpos[1])));

            } else if (prob_type == 3) {
                
                u(i,j,k) = 0.25*(1.+std::tanh(k1_inv*(relpos[1] - (-width1/2.))))
                    *(1.+std::tanh(k2_inv*((width1/2.) - relpos[1])));

            }
        },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            GpuArray<Real,AMREX_SPACEDIM> itVec;

            AMREX_D_TERM(itVec[0] = (i+0.5)*dx[0];,
                         itVec[1] = j*dx[1];,
                         itVec[2] = (k+zshft)*dx[2];);

            GpuArray<Real,AMREX_SPACEDIM> relpos;
            Real rad2 = 0.;
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                relpos[d] = reallo[d] + itVec[d] - center[d];
            }

            // only sum first two velocities
            for (int d=0; d<2; ++d) {
                rad2 += relpos[d]*relpos[d];
            }

            Real rad = std::sqrt(rad2);

            // note prob_type == 0 handled above
            if (prob_type == 1) {

                // Multiply velocity magnitude by -cos(theta)
                v(i,j,k) = 0.25*(1.+std::tanh(k1_inv*(rad-r_a)))*(1.+std::tanh(k2_inv*(r_b-rad)))
                    *(-relpos[0]/rad);
                
            } else if (prob_type == 2) {
                
                Real perturb = amp*sin(freq*relpos[0]);
                Real slope = amp*freq*cos(freq*relpos[0]);
                Real fun_ptrb = 0.25*(1.*tanh(k1_inv*(relpos[1] - (-width1/2.+perturb))))
                    *(1.+tanh(k2_inv*((width1/2.+perturb) - relpos[1])));
                v(i,j,k) = slope*fun_ptrb;
            } else if (prob_type == 3) {
                v(i,j,k) = 0.;
            }
        }
#if (AMREX_SPACEDIM == 3)
                         , [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            w(i,j,k) = 0.;
        }
#endif
        );
    }
}

void InitTracer(MultiFab& tracer,
                const Geometry& geom) {

    GpuArray<Real,AMREX_SPACEDIM> reallo = geom.ProbLoArray();
    GpuArray<Real,AMREX_SPACEDIM> realhi = geom.ProbHiArray();
    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    GpuArray<Real,AMREX_SPACEDIM> center;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        center[d] = 0.5*(realhi[d]-reallo[d]);
    }
            
    Real L_hlf = center[0];
    
    Real pi = 3.1415926535897932;

    // k1 & k2 determine steepness of profile:
    Real k1 = 1.e-2*L_hlf;
    Real k2 = k1;
    Real k1_inv = 1./k1;
    Real k2_inv = 1./k2;

    // Vortex:
    // [r_a r_b] defines radial bounds of bump:
    Real r_a = 0.5*L_hlf;
    // r_b = L_hlf - r_a

    // Stream:
    // freq = 3.d0*pi/L_hlf
    // amp = 2.0d-1*L_hlf
    Real width1 = L_hlf/2.;
            
    for ( MFIter mfi(tracer); mfi.isValid(); ++mfi ) {

        const Array4<Real> & phic = tracer.array(mfi);

        Box bx = mfi.tilebox();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            GpuArray<Real,AMREX_SPACEDIM> itVec;

            AMREX_D_TERM(itVec[0] = i*dx[0];,
                         itVec[1] = j*dx[1];,
                         itVec[2] = k*dx[2];);

            GpuArray<Real,AMREX_SPACEDIM> relpos;
            Real rad2 = 0.;
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                relpos[d] = reallo[d] + itVec[d] - center[d];
                rad2 += relpos[d]*relpos[d];
            }

            Real rad = std::sqrt(rad2);
            
            // Circle
            // phic(i,j,k) = 0.5d0*(1d0+tanh(k2_inv*(r_a-rad)))

            // Stream:
            // perturb = amp*sin(freq*relpos(1))
            // slope = amp*freq*cos(freq*relpos(1))
            Real perturb = 0.; 

            phic(i,j,k) = 0.25*(1.+std::tanh(k1_inv*(relpos[1] - (-width1/2.+perturb)))) 
                *(1.+std::tanh(k2_inv*((width1/2.+perturb) - relpos[1])));
        });
     }
}
