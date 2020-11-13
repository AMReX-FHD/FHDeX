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
                    for (int n=0; n<nspecies; ++n) {
                        prim(i,j,k,6+n)          = 2.*bc_Yk_x_lo_gpu[n] - prim(-i-1,j,k,6+n);
                        prim(i,j,k,6+nspecies+n) = 2.*bc_Xk_x_lo_gpu[n] - prim(-i-1,j,k,6+nspecies+n);
                    }
                }
            });
        }

        // temperature and pressure, adiabatic
        if (bc_therm_lo[0] == 1) {


        }


        
        const Box& vbx = mfi.tilebox();
        
        set_bc(ARLIM_3D(vbx.loVect()), ARLIM_3D(vbx.hiVect()),
                       cons_in[mfi].dataPtr(),
                       prim_in[mfi].dataPtr());

    }
}

