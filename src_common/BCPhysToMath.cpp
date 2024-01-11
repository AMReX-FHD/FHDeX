#include "common_functions.H"

// this is loosely analogous to what "define_bc_tower" used to do in pure fortran codes
// it converts "physical" boundary descriptions to "mathematical" descriptions used
// to fill ghost cells in PhysBC routines

void BCPhysToMath(int bccomp, amrex::Vector<int>& bc_lo, amrex::Vector<int>& bc_hi) {
    
    BL_PROFILE_VAR("BCPhysToMath()",BCPhysToMath);

    // set to interior/periodic by default; overwrite below
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        bc_lo[i] = bc_hi[i] = INT_DIR;
    }
    
    if (bccomp == PRES_BC_COMP) { // PRESSURE
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            if (bc_vel_lo[i] == 1 || bc_vel_lo[i] == 2) {
                // wall -> first-order extrapolation
                bc_lo[i] = FOEXTRAP;
            } else if (bc_vel_lo[i] == -2) {
                // pressure inflow
                bc_lo[i] = EXT_DIR;
            }
            if (bc_vel_hi[i] == 1 || bc_vel_hi[i] == 2) {
                // wall -> first-order extrapolation
                bc_hi[i] = FOEXTRAP;
            } else if (bc_vel_hi[i] == -2) {
                // pressure inflow
                bc_hi[i] = EXT_DIR;
            }
        }
    }
    else if (bccomp == RHO_BC_COMP ) { // density
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            if (bc_mass_lo[i] == 1) {
                // wall -> first-order extrapolation
                bc_lo[i] = FOEXTRAP;
            }
            else if (bc_mass_lo[i] == 2) {
                // reservoir -> dirichlet
                bc_lo[i] = EXT_DIR;
            }
            else if (bc_mass_lo[i] == 4) {
                // reservoir -> dirichlet
                bc_lo[i] = FOEXTRAP;
            }
            if (bc_mass_hi[i] == 1) {
                // wall -> first-order extrapolation
                bc_hi[i] = FOEXTRAP;
            }
            else if (bc_mass_hi[i] == 2) {
                // reservoir -> dirichlet
                bc_hi[i] = EXT_DIR;
            }
            else if (bc_mass_hi[i] == 4) {
                // reservoir -> dirichlet
                bc_hi[i] = FOEXTRAP;
            }
        }
    }
    else if (bccomp == SPEC_BC_COMP ||
             bccomp == MOLFRAC_BC_COMP) { // density, species, or mole fraction
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            if (bc_mass_lo[i] == 1) {
                // wall -> first-order extrapolation
                bc_lo[i] = FOEXTRAP;
            }
            else if (bc_mass_lo[i] == 2) {
                // reservoir -> dirichlet
                bc_lo[i] = EXT_DIR;
            }
            else if (bc_mass_lo[i] == 4) {
                // reservoir -> dirichlet
                bc_lo[i] = SPEC_CONTACT_BC;
            }
            if (bc_mass_hi[i] == 1) {
                // wall -> first-order extrapolation
                bc_hi[i] = FOEXTRAP;
            }
            else if (bc_mass_hi[i] == 2) {
                // reservoir -> dirichlet
                bc_hi[i] = EXT_DIR;
            }
            else if (bc_mass_hi[i] == 4) {
                // reservoir -> dirichlet
                bc_hi[i] = SPEC_CONTACT_BC;
            }
        }
    }
    else if (bccomp == TEMP_BC_COMP) { // TEMPERATURE
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            if (bc_therm_lo[i] == 1) {
                // adiabtic -> first-order extrapolation
                bc_lo[i] = FOEXTRAP;
            }
            else if (bc_therm_lo[i] == 2) {
                // isothermal -> dirichlet
                bc_lo[i] = EXT_DIR;
            }
            if (bc_therm_hi[i] == 1) {
                // adiabatic -> first-order extrapolation
                bc_hi[i] = FOEXTRAP;
            }
            else if (bc_therm_hi[i] == 2) {
                // isothermal -> dirichlet
                bc_hi[i] = EXT_DIR;
            }
        }
    }
    else if (bccomp == EPOT_BC_COMP) { // ELECTRIC POTENTIAL
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            if (bc_es_lo[i] == 1) {
                // dirichlet
                bc_lo[i] = EXT_DIR;
            }
            else if (bc_es_lo[i] == 2) {
                // neumann
                bc_lo[i] = FOEXTRAP;
            }
            if (bc_es_hi[i] == 1) {
                // dirichlet
                bc_hi[i] = EXT_DIR;
            }
            else if (bc_es_hi[i] == 2) {
                // neumann
                bc_hi[i] = FOEXTRAP;
            }
        }
    }
}
   
