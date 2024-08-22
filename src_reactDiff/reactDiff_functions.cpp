#include "reactDiff_functions.H"
#include "AMReX_ParmParse.H"

// only used for split schemes (temporal_integrator>=0)
// 0=explicit trapezoidal predictor/corrector
// 1=Crank-Nicolson semi-implicit
// 2=explicit midpoint
// 3=multinomial diffusion
// 4=forward Euler  

AMREX_GPU_MANAGED int reactDiff::diffusion_type;

void InitializeReactDiffNamespace()
{
    // extract inputs parameters
    ParmParse pp;

    diffusion_type = 0;
    pp.query("diffusion_type",diffusion_type);
    
    return;
}
