#include "common_functions.H"

// this is loosely analogous to what "define_bc_tower" used to do in pure fortran codes
// it converts "physical" boundary descriptions to "mathematical" descriptions used
// to fill ghost cells in PhysBC routines

void BCPhysToMath(int varType, amrex::Vector<int>& bc_lo, amrex::Vector<int>& bc_hi) {
    
    BL_PROFILE_VAR("BCPhysToMath()",BCPhysToMath);

    // varType -1: density
    // varType  0: pressure
    // varType  1: species
    // varType  2: temperature
    // varType  3: electric potential

    // set to interior/periodic by default; overwrite below
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        bc_lo[i] = bc_hi[i] = INT_DIR;
    }
    
    if (varType == 0) { // PRESSURE
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            if (bc_vel_lo[i] == 1 || bc_vel_lo[i] == 2) {
                // wall -> first-order extrapolation
                bc_lo[i] = FOEXTRAP;
            }
            if (bc_vel_hi[i] == 1 || bc_vel_hi[i] == 2) {
                // wall -> first-order extrapolation
                bc_hi[i] = FOEXTRAP;
            }
        }
    }
    else if (varType == -1 || varType == 1) { // density or species
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            if (bc_mass_lo[i] == 1) {
                // wall -> first-order extrapolation
                bc_lo[i] = FOEXTRAP;
            }
            else if (bc_mass_lo[i] == 2) {
                // reservoir -> dirichlet
                bc_lo[i] = EXT_DIR;
            }
            if (bc_mass_hi[i] == 1) {
                // wall -> first-order extrapolation
                bc_hi[i] = FOEXTRAP;
            }
            else if (bc_mass_hi[i] == 2) {
                // reservoir -> dirichlet
                bc_hi[i] = EXT_DIR;
            }
        }
    }
    else if (varType == 2) { // TEMPERATURE
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
    else if (varType == 3) { // ELECTRIC POTENTIAL
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
   
