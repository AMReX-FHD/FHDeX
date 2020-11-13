#include "compressible_functions.H"

void setBC(MultiFab& prim_in, MultiFab& cons_in)
{
    BL_PROFILE_VAR("setBC()",setBC);

    int ng_c = cons_in.nGrow();
    int ng_p = prim_in.nGrow();
    if (ng_c != ng_p) {
        Abort("setBC: prim and cons need the same number of ghost cells");
    }

    int nprimvars_gpu = nprimvars;

    GpuArray<Real,MAX_SPECIES> bc_Yk_x_lo_gpu;
    GpuArray<Real,MAX_SPECIES> bc_Yk_x_hi_gpu;
    GpuArray<Real,MAX_SPECIES> bc_Yk_y_lo_gpu;
    GpuArray<Real,MAX_SPECIES> bc_Yk_y_hi_gpu;

    for (int n=0; n<nspecies; ++n) {
        bc_Yk_x_lo_gpu[n] = bc_Yk[n*LOHI*AMREX_SPACEDIM];
        bc_Yk_x_hi_gpu[n] = bc_Yk[AMREX_SPACEDIM + n*LOHI*AMREX_SPACEDIM];
        bc_Yk_y_lo_gpu[n] = bc_Yk[1 + n*LOHI*AMREX_SPACEDIM];
        bc_Yk_y_hi_gpu[n] = bc_Yk[1 + AMREX_SPACEDIM + n*LOHI*AMREX_SPACEDIM];
    }

    GpuArray<Real,MAX_SPECIES> bc_Xk_x_lo_gpu;
    GpuArray<Real,MAX_SPECIES> bc_Xk_x_hi_gpu;
    GpuArray<Real,MAX_SPECIES> bc_Xk_y_lo_gpu;
    GpuArray<Real,MAX_SPECIES> bc_Xk_y_hi_gpu;

    for (int n=0; n<nspecies; ++n) {
        bc_Xk_x_lo_gpu[n] = bc_Xk[n*LOHI*AMREX_SPACEDIM];
        bc_Xk_x_hi_gpu[n] = bc_Xk[AMREX_SPACEDIM + n*LOHI*AMREX_SPACEDIM];
        bc_Xk_y_lo_gpu[n] = bc_Xk[1 + n*LOHI*AMREX_SPACEDIM];
        bc_Xk_y_hi_gpu[n] = bc_Xk[1 + AMREX_SPACEDIM + n*LOHI*AMREX_SPACEDIM];
    }
    
    GpuArray<Real,MAX_SPECIES> molmass_gpu;
    for (int n=0; n<nspecies; ++n) {
        molmass_gpu[n] = molmass[n];
    }
    
    GpuArray<Real,MAX_SPECIES> hcv_gpu;
    for (int n=0; n<nspecies; ++n) {
        hcv_gpu[n] = hcv[n];
    }
    
    Real Runiv_gpu = Runiv;

    int nspecies_gpu = nspecies;
    
    // Loop over boxes
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        const Box& gbx = mfi.growntilebox(ng_c);
                
        const Array4<Real> cons = cons_in.array(mfi);
        const Array4<Real> prim = prim_in.array(mfi);

        // LOWER X
        
        // mass fractions, wall
        if (bc_mass_lo[0] == 1 && algorithm_type == 2) {
            amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < 0) {
                    for (int n=5; n<nprimvars; ++n) {
                        prim(i,j,k,n) = prim(-i-1,j,k,n);
                    }
                }
            });
        }

        // mass fracations, reservoir
        if (bc_mass_lo[0] == 2 && algorithm_type == 2) {
            amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < 0) {
                    for (int n=0; n<nspecies_gpu; ++n) {
                        prim(i,j,k,6+n)          = 2.*bc_Yk_x_lo_gpu[n] - prim(-i-1,j,k,6+n);
                        prim(i,j,k,6+nspecies+n) = 2.*bc_Xk_x_lo_gpu[n] - prim(-i-1,j,k,6+nspecies+n);
                    }
                }
            });
        }

        // temperature and pressure, adiabatic
        if (bc_therm_lo[0] == 1) {
            amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < 0) {
                    for (int n=4; n<6; ++n) {
                        prim(i,j,k,n) = prim(-i-1,j,k,n);
                    }
                }
            });
        } else if (bc_therm_lo[0] == 2) { // isothermal
            amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < 0) {
                    prim(i,j,k,4) = -prim(-i-1,j,k,4) + 2.*t_lo[0];
                    prim(i,j,k,5) = prim(-i-1,j,k,5);
                }
            });
        }

        // momentum, velocity, rho, rhoY, rhoE
        if (bc_vel_lo[0] == 1) { // slip
            amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < 0) {

                    cons(i,j,k,1) = -cons(-i-1,j,k,1);
                    cons(i,j,k,2) = cons(-i-1,j,k,2);
                    cons(i,j,k,3) = cons(-i-1,j,k,3);
                        
                    prim(i,j,k,1) = -prim(-i-1,j,k,1);
                    prim(i,j,k,2) = prim(-i-1,j,k,2);
                    prim(i,j,k,3) = prim(-i-1,j,k,3);

                    // thermal & species (+pressure) BCs must be enforced first
                    GpuArray<Real,MAX_SPECIES> fracvec_gpu;
                    for (int n=0; n<nspecies_gpu; ++n) {
                        fracvec_gpu[n] = prim(i,j,k,6+n);
                    }
                    Real temp = prim(i,j,k,4);
                    Real pt = prim(i,j,k,5);
                    Real rho;
                    Real intenergy;
                    
                    GetDensity(pt,rho,temp,fracvec_gpu,molmass_gpu,Runiv_gpu,nspecies_gpu);
                    GetEnergy(intenergy,fracvec_gpu,temp,hcv_gpu,nspecies_gpu);

                    // total density depends on temperature
                    prim(i,j,k,0) = rho;
                    cons(i,j,k,0) = rho;
                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies_gpu; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }

                    // must be last BC enforced: depends on rho, vel, & temp
                    cons(i,j,k,4) = rho*intenergy + 0.5*rho*(prim(i,j,k,1)*prim(i,j,k,1) + 
                                                             prim(i,j,k,2)*prim(i,j,k,2) +
                                                             prim(i,j,k,3)*prim(i,j,k,3));
                    
                }
            });
        }
        
        
        const Box& vbx = mfi.tilebox();
        
        set_bc(ARLIM_3D(vbx.loVect()), ARLIM_3D(vbx.hiVect()),
                       cons_in[mfi].dataPtr(),
                       prim_in[mfi].dataPtr());

    }
}

