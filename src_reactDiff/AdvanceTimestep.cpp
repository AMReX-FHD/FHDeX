#include "reactDiff_functions.H"
#include "chemistry_functions.H"

void AdvanceTimestep(MultiFab& n_old,
                     MultiFab& n_new,
                     const Real& dt,
                     const Real& time,
                     const Geometry& geom) {

    if (temporal_integrator >= 0 && reactDiff_reaction_type != 0) {
        if (reaction_type == 2) {
            Abort("SSA (reaction_type==2) requires reactDiff_reaction_type=0 for split schemes");
        }
    }

    // external source term for diffusion/reaction solvers for inhomogeneous bc algorithm
    MultiFab Rn_steady(n_old.boxArray(), n_old.DistributionMap(), nspecies, 0);

    if (temporal_integrator < 0) {

        Rn_steady.setVal(0.);

        // unsplit schemes
        AdvanceReactionDiffusion(n_old,n_new,Rn_steady,dt,time,geom);

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
            AdvanceReaction(n_old,n_new,Rn_steady,0.5*dt,time,geom);
            // swap n_new/n_old to avoid calling copy()
            AdvanceDiffusion(n_new,n_old,Rn_steady,dt,time,geom);
            AdvanceReaction(n_old,n_new,Rn_steady,0.5*dt,time,geom);

        } else if (temporal_integrator == 2) {
            // (1/2)D + R + (1/2)D
            AdvanceDiffusion(n_old,n_new,Rn_steady,0.5*dt,time,geom);
            // swap n_new/n_old to avoid calling copy()
            AdvanceReaction(n_new,n_old,Rn_steady,dt,time,geom);
            AdvanceDiffusion(n_old,n_new,Rn_steady,0.5*dt,time,geom);

        } else {
            Abort("AdvanceTimestep(): invalid temporal_integrator");
        }

    }

}
