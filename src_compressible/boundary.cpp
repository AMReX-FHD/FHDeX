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
    
    // Loop over boxes
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        const Box& gbx = mfi.growntilebox(ng_c);
                
        const Array4<Real> cons = cons_in.array(mfi);
        const Array4<Real> prim = prim_in.array(mfi);

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
        
        const Box& vbx = mfi.tilebox();
        
        set_bc(ARLIM_3D(vbx.loVect()), ARLIM_3D(vbx.hiVect()),
                       cons_in[mfi].dataPtr(),
                       prim_in[mfi].dataPtr());

    }
}

