#include "thinfilm_functions.H"

#include "AMReX_ParmParse.H"

AMREX_GPU_MANAGED amrex::Real thinfilm::thinfilm_h0;
AMREX_GPU_MANAGED amrex::Real thinfilm::thinfilm_gamma;
AMREX_GPU_MANAGED amrex::Real thinfilm::thinfilm_pertamp;

AMREX_GPU_MANAGED int thinfilm::thinfilm_icorr;
AMREX_GPU_MANAGED int thinfilm::thinfilm_jcorr;

void InitializeThinfilmNamespace() {

    BL_PROFILE_VAR("InitializeThinfilmNamespace()",InitializeThinfilmNameSpace);

    // defaults
    thinfilm_icorr = 0;
    thinfilm_jcorr = 0;
    thinfilm_pertamp = 0;
    
    ParmParse pp;
    
    pp.get("thinfilm_h0",thinfilm_h0);
    pp.get("thinfilm_gamma",thinfilm_gamma);
    pp.get("thinfilm_pertamp",thinfilm_pertamp);

    pp.query("thinfilm_icorr",thinfilm_icorr);
    pp.query("thinfilm_icorr",thinfilm_jcorr);
    
}
