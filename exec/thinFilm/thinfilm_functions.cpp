#include "thinfilm_functions.H"

#include "AMReX_ParmParse.H"

AMREX_GPU_MANAGED amrex::Real thinfilm::thinfilm_h0;
AMREX_GPU_MANAGED amrex::Real thinfilm::thinfilm_gamma;

void InitializeThinfilmNamespace() {

    BL_PROFILE_VAR("InitializeThinfilmNamespace()",InitializeThinfilmNameSpace);
  
    ParmParse pp;
    
    pp.get("thinfilm_h0",thinfilm_h0);
    pp.get("thinfilm_gamma",thinfilm_gamma);
    
}
