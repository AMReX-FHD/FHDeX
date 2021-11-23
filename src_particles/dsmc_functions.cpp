#include "dsmc_functions.H"


#include "AMReX_ParmParse.H"

int         dsmc::sim_type;
amrex::Real dsmc::hbar;
amrex::Real dsmc::debye_group_velocity;

void InitializeDsmcNamespace() {


    sim_type = 0;

    // use the viscosity-based BFBt Schur complement (from Georg Stadler)
    hbar = 1.0546e-27;

    // weighting of pressure when computing norms and inner products
    debye_group_velocity = 1;

    ParmParse pp;

    // pp.query searches for optional parameters
    // pp.get aborts if the parameter is not found
    // pp.getarr and queryarr("string",inputs,start_indx,count); can be used for arrays

    pp.query("sim_type",sim_type);
    pp.query("hbar",hbar);
    pp.query("debye_group_velocity",debye_group_velocity);
}
