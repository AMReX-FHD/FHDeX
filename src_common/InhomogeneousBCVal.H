AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE

Real InhomogeneousBCVal ( int comp,
                          Real x,
                          Real y,
                          Real z,
                          Real time)
{
    if (comp == RHO_BC_COMP) {

        // only the correct thing for homogeneous Neumann
        // need to fix to support Dirichlet or inhomogeneous conditions
        return 0;
        
    } else if (comp >= SPEC_BC_COMP && comp <= SPEC_END_BC_COMP) {

        // only the correct thing for homogeneous Neumann
        // need to fix to support Dirichlet or inhomogeneous conditions
        return 0;
        
    } else if (comp == TEMP_BC_COMP) {
        
        if (x == prob_lo[0]) {
            return t_lo[0];
        } else if (x == prob_hi[0]) {
            return t_hi[0];
        } else if (y == prob_lo[1]) {
            return t_lo[1];
        } else if (y == prob_hi[1]) {
            return t_hi[1];
        } else if (z == prob_lo[2]) {
            return t_lo[2];
        } else if (z == prob_hi[2]) {
            return t_hi[2];
        }  

    } else if (comp == EPOT_BC_COMP) {

        if (x == prob_lo[0]) {
            return potential_lo[0];
        } else if (x == prob_hi[0]) {
            return potential_hi[0];
        } else if (y == prob_lo[1]) {
            return potential_lo[1];
        } else if (y == prob_hi[1]) {
            return potential_hi[1];
        } else if (z == prob_lo[2]) {
            return potential_lo[2];
        } else if (z == prob_hi[2]) {
            return potential_hi[2];
        }

    } else {

        printf("comp %d",comp);
        Abort("InhomogeneousBCVal: invalid comp");

    }

    
    return 0.;
};
