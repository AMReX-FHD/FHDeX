#include "reactDiff_functions.H"
#include "chemistry_functions.H"

void AdvanceTimestep(MultiFab& n_old,
                     MultiFab& n_new,
                     const Real& dt,
                     const Real& time,
                     const Geometry& geom) {

    if (temporal_integrator > 0 && reactDiff_reaction_type != 0) {
        if (reaction_type == 2) {
            Abort("SSA (reaction_type==2) requires reactDiff_reaction_type=0 for split schemes");
        }
    }

    // external source term for diffusion/reaction solvers for inhomogeneous bc algorithm
    MultiFab Rn_steady(n_old.boxArray(), n_old.DistributionMap(), nspecies, 0);
    
    if (temporal_integrator < 0) {
        // unsplit schemes



    } else {

        if (inhomogeneous_bc_fix) {
            Abort("inhomogeneous_bc_fix not implemented yet");
        } else {
            Rn_steady.setVal(0.);
        }

        if (temporal_integrator == 0) {
            // D + R
            AdvanceDiffusion(n_old,n_new,Rn_steady,dt,time,geom);
            MultiFab::Copy(n_old,n_new,0,0,nspecies,1);
            AdvanceReaction(n_old,n_new,Rn_steady,dt,time,geom);

        } else if (temporal_integrator == 1) {
            // (1/2)R + D + (1/2)R
/*
          call advance_reaction (mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower,Rn_steady)
          ! swap n_new/n_old to avoid calling copy()
          call advance_diffusion(mla,n_new,n_old,dx,dt      ,the_bc_tower,Rn_steady)  
          call advance_reaction (mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower,Rn_steady)
*/

        } else if (temporal_integrator == 2) {
            // (1/2)D + R + (1/2)D
/*
          call advance_diffusion(mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower,Rn_steady)
          ! swap n_new/n_old to avoid calling copy()
          call advance_reaction (mla,n_new,n_old,dx,dt      ,the_bc_tower,Rn_steady)
          call advance_diffusion(mla,n_old,n_new,dx,0.5d0*dt,the_bc_tower,Rn_steady)
*/

        } else {
            Abort("AdvanceTimestep(): invalid temporal_integrator");
        }

    }
    
}
