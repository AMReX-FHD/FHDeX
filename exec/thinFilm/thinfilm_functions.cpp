#include "thinfilm_functions.H"

#include "AMReX_ParmParse.H"

AMREX_GPU_MANAGED amrex::Real thinfilm::thinfilm_h0;
AMREX_GPU_MANAGED amrex::Real thinfilm::thinfilm_gamma;
AMREX_GPU_MANAGED amrex::Real thinfilm::thinfilm_pertamp;
AMREX_GPU_MANAGED amrex::Real thinfilm::thinfilm_hamaker;

AMREX_GPU_MANAGED int thinfilm::thinfilm_icorr;
AMREX_GPU_MANAGED int thinfilm::thinfilm_jcorr;

AMREX_GPU_MANAGED int thinfilm::do_fft_diag;

void InitializeThinfilmNamespace() {

    BL_PROFILE_VAR("InitializeThinfilmNamespace()",InitializeThinfilmNameSpace);

    // defaults
    thinfilm_icorr = 0;
    thinfilm_jcorr = 0;
    thinfilm_pertamp = 0.;
    thinfilm_hamaker = 0.;

    do_fft_diag = 1;

    ParmParse pp;

    pp.get("thinfilm_h0",thinfilm_h0);
    pp.get("thinfilm_gamma",thinfilm_gamma);

    pp.query("thinfilm_pertamp",thinfilm_pertamp);
    pp.query("thinfilm_hamaker",thinfilm_hamaker);

    pp.query("thinfilm_icorr",thinfilm_icorr);
    pp.query("thinfilm_jcorr",thinfilm_jcorr);

    pp.query("do_fft_diag",do_fft_diag);

}

#ifdef AMREX_USE_CUDA
std::string cufftErrorToString (const cufftResult& err)
{
    switch (err) {
    case CUFFT_SUCCESS:  return "CUFFT_SUCCESS";
    case CUFFT_INVALID_PLAN: return "CUFFT_INVALID_PLAN";
    case CUFFT_ALLOC_FAILED: return "CUFFT_ALLOC_FAILED";
    case CUFFT_INVALID_TYPE: return "CUFFT_INVALID_TYPE";
    case CUFFT_INVALID_VALUE: return "CUFFT_INVALID_VALUE";
    case CUFFT_INTERNAL_ERROR: return "CUFFT_INTERNAL_ERROR";
    case CUFFT_EXEC_FAILED: return "CUFFT_EXEC_FAILED";
    case CUFFT_SETUP_FAILED: return "CUFFT_SETUP_FAILED";
    case CUFFT_INVALID_SIZE: return "CUFFT_INVALID_SIZE";
    case CUFFT_UNALIGNED_DATA: return "CUFFT_UNALIGNED_DATA";
    default: return std::to_string(err) + " (unknown error code)";
    }
}
#endif
