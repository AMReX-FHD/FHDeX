AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE

Real InhomogeneousBCVal ( int comp,
                          Real x,
                          Real y,
                          Real z,
                          Real /*time*/)
{
    if (comp == PRES_BC_COMP) {

        if (x == prob_lo[0]) {
            return p_lo[0];
        } else if (x == prob_hi[0]) {
            return p_hi[0];
        } else if (y == prob_lo[1]) {
            return p_lo[1];
        } else if (y == prob_hi[1]) {
            return p_hi[1];
        } else if (z == prob_lo[2]) {
            return p_lo[2];
        } else if (z == prob_hi[2]) {
            return p_hi[2];
        }

    } else if (comp == RHO_BC_COMP) {

        if (x == prob_lo[0]) {
            return bc_rhotot_x_lo;
        } else if (x == prob_hi[0]) {
            return bc_rhotot_x_hi;
        } else if (y == prob_lo[1]) {
            return bc_rhotot_y_lo;
        } else if (y == prob_hi[1]) {
            return bc_rhotot_y_hi;
        } else if (z == prob_lo[2]) {
            return bc_rhotot_z_lo;
        } else if (z == prob_hi[2]) {
            return bc_rhotot_z_hi;
        }

    } else if (comp == TRAN_BC_COMP) {

        if (x == prob_lo[0]) {
            return 0.;
        } else if (x == prob_hi[0]) {
            return 0.;
        } else if (y == prob_lo[1]) {
            return 0.;
        } else if (y == prob_hi[1]) {
            return 0.;
        } else if (z == prob_lo[2]) {
            return 0.;
        } else if (z == prob_hi[2]) {
            return 0.;
        }

    } else if (comp >= SPEC_BC_COMP && comp <= SPEC_END_BC_COMP) {

        if (x == prob_lo[0]) {
            return bc_Yk_x_lo[comp-SPEC_BC_COMP];
        } else if (x == prob_hi[0]) {
            return bc_Yk_x_hi[comp-SPEC_BC_COMP];
        } else if (y == prob_lo[1]) {
            return bc_Yk_y_lo[comp-SPEC_BC_COMP];
        } else if (y == prob_hi[1]) {
            return bc_Yk_y_hi[comp-SPEC_BC_COMP];
        } else if (z == prob_lo[2]) {
            return bc_Yk_z_lo[comp-SPEC_BC_COMP];
        } else if (z == prob_hi[2]) {
            return bc_Yk_z_hi[comp-SPEC_BC_COMP];
        }

    } else if (comp >= MOLFRAC_BC_COMP && comp <= MOLFRAC_END_BC_COMP) {

        if (x == prob_lo[0]) {
            return bc_Xk_x_lo[comp-MOLFRAC_BC_COMP];
        } else if (x == prob_hi[0]) {
            return bc_Xk_x_hi[comp-MOLFRAC_BC_COMP];
        } else if (y == prob_lo[1]) {
            return bc_Xk_y_lo[comp-MOLFRAC_BC_COMP];
        } else if (y == prob_hi[1]) {
            return bc_Xk_y_hi[comp-MOLFRAC_BC_COMP];
        } else if (z == prob_lo[2]) {
            return bc_Xk_z_lo[comp-MOLFRAC_BC_COMP];
        } else if (z == prob_hi[2]) {
            return bc_Xk_z_hi[comp-MOLFRAC_BC_COMP];
        }

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
}
