#include "MFsurfchem_functions.H"
#include "AMReX_ParmParse.H"

AMREX_GPU_MANAGED int MFsurfchem::ads_spec;
AMREX_GPU_MANAGED amrex::Real MFsurfchem::surfcov0;

AMREX_GPU_MANAGED amrex::Real MFsurfchem::ads_rate_const;
AMREX_GPU_MANAGED amrex::Real MFsurfchem::des_rate;

void InitializeMFSurfchemNamespace()
{
    // extract inputs parameters
    ParmParse pp;

    ads_spec = -1;
    // get the species number that undergoes adsoprtion/desorption 
    pp.query("ads_spec",ads_spec);

    // if ads_spec is set to -1 or not defined in the inputs file, quit the routine
    if (ads_spec==-1) return;

    surfcov0 = 0.;
    // get initial surface coverage
    pp.query("surfcov0",surfcov0);

    ads_rate_const = 0.;
    // get adsorption rate const 
    pp.query("ads_rate_const",ads_rate_const);

    des_rate = 0.;
    // get desoprtion rate
    pp.query("des_rate",des_rate);

    return;
}
