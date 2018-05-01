#include "common_functions.H"
#include "common_functions_F.H"
#include "common_params.H"

using namespace common;

void set_common_params() {

    prob_lo.resize(AMREX_SPACEDIM);
    prob_hi.resize(AMREX_SPACEDIM);
    n_cells.resize(AMREX_SPACEDIM);
    max_grid_size.resize(AMREX_SPACEDIM);
    grav.resize(AMREX_SPACEDIM);
    molmass.resize(MAX_SPECIES);
    rhobar.resize(MAX_SPECIES);
    u_init.resize(2);
    bc_lo.resize(AMREX_SPACEDIM);
    bc_hi.resize(AMREX_SPACEDIM);
    wallspeed_lo.resize((AMREX_SPACEDIM-1)*AMREX_SPACEDIM);
    wallspeed_hi.resize((AMREX_SPACEDIM-1)*AMREX_SPACEDIM);
    density_weights.resize(MAX_SPECIES);
    shift_cc_to_boundary.resize(AMREX_SPACEDIM*LOHI);
    
    copy_common_params_to_c(prob_lo.dataPtr(),prob_hi.dataPtr());


}
