#include "reactDiff_functions.H"
#include "AMReX_ParmParse.H"

// 0=D + R (first-order splitting)
// 1=(1/2)R + D + (1/2)R (Strang option 1)
// 2=(1/2)D + R + (1/2)D (Strang option 2)
// -1=unsplit forward Euler
// -2=unsplit explicit midpoint
// -3=unsplit multinomial diffusion
// -4=unsplit implicit midpoint
AMREX_GPU_MANAGED int reactDiff::temporal_integrator;

// only used for split schemes (temporal_integrator>=0)
// 0=explicit trapezoidal predictor/corrector
// 1=Crank-Nicolson semi-implicit
// 2=explicit midpoint
// 3=multinomial diffusion
// 4=forward Euler
AMREX_GPU_MANAGED int reactDiff::reactDiff_diffusion_type;

// only used for split schemes (temporal_integrator>=0)
// 0=first-order (deterministic, tau leaping, CLE, or SSA)
// 1=second-order (determinisitc, tau leaping, or CLE only)
AMREX_GPU_MANAGED int reactDiff::reactDiff_reaction_type;

// only used for midpoint diffusion schemes (split as well as unsplit)
// corrector formulation of noise
// 1 = K(nold) * W1 + K(nold)         * W2
// 2 = K(nold) * W1 + K(npred)        * W2
// 3 = K(nold) * W1 + K(2*npred-nold) * W2
AMREX_GPU_MANAGED int reactDiff::midpoint_stoch_flux_type;

// how to compute n on faces for stochastic weighting
// 1=arithmetic (with C0-Heaviside), 2=geometric, 3=harmonic
// 10=arithmetic average with discontinuous Heaviside function
// 11=arithmetic average with C1-smoothed Heaviside function
// 12=arithmetic average with C2-smoothed Heaviside function
AMREX_GPU_MANAGED int reactDiff::avg_type;

// use the Einkemmer boundary condition fix (split schemes only)
AMREX_GPU_MANAGED int reactDiff::inhomogeneous_bc_fix;

// volume multiplier (dv = product(dx(1:MAX_SPACEDIM))*volume_factor)
// only really intended for 3D since in 2D one can control the cell depth
AMREX_GPU_MANAGED amrex::Real reactDiff::volume_factor;

// initial values to be used in init_n.f90
AMREX_GPU_MANAGED Array2D<amrex::Real, 0, 2 ,0, MAX_SPECIES> reactDiff::n_init_in;

// initialize from model file
AMREX_GPU_MANAGED int reactDiff::model_file_init;

// initialize with all number of molecules strictly integer
AMREX_GPU_MANAGED int reactDiff::integer_populations;

// Fickian diffusion coeffs
AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> reactDiff::D_Fick;

// diffusion boundary stencil order
AMREX_GPU_MANAGED int reactDiff::diffusion_stencil_order;

// implicit diffusion solve verbosity
AMREX_GPU_MANAGED int reactDiff::diffusion_verbose;

// implicit diffusion solve bottom solver verbosity
AMREX_GPU_MANAGED int reactDiff::diffusion_bottom_verbose;

// relative eps for implicit diffusion solve
AMREX_GPU_MANAGED amrex::Real reactDiff::implicit_diffusion_rel_eps;

// absolute eps for implicit diffusion solve
AMREX_GPU_MANAGED amrex::Real reactDiff::implicit_diffusion_abs_eps;

void InitializeReactDiffNamespace()
{
    // extract inputs parameters
    ParmParse pp;

    int temp_max = std::max(MAX_SPECIES,MAX_REACTION);

    amrex::Vector<amrex::Real> temp    (temp_max,0.);
    amrex::Vector<int>         temp_int(temp_max,0 );

    temporal_integrator = 0;
    pp.query("temporal_integrator",temporal_integrator);

    reactDiff_diffusion_type = 0;
    pp.query("reactDiff_diffusion_type",reactDiff_diffusion_type);

    reactDiff_reaction_type = 0;
    pp.query("reactDiff_reaction_type",reactDiff_reaction_type);

    midpoint_stoch_flux_type = 1;
    pp.query("midpoint_stoch_flux_type",midpoint_stoch_flux_type);

    avg_type = 1;
    pp.query("avg_type",avg_type);

    inhomogeneous_bc_fix = 0;
    pp.query("inhomogeneous_bc_fix",inhomogeneous_bc_fix);

    volume_factor = 1.;
    pp.query("volume_factor",volume_factor);

    if (pp.queryarr("n_init_in_1",temp)) {
        for (int i=0; i<nspecies; ++i) {
            n_init_in(0,i) = temp[i];
        }
    }
    if (pp.queryarr("n_init_in_2",temp)) {
        for (int i=0; i<nspecies; ++i) {
            n_init_in(1,i) = temp[i];
        }
    }

    model_file_init = 0;
    pp.query("model_file_init",model_file_init);

    integer_populations = 0;
    pp.query("integer_populations",integer_populations);

    if (pp.queryarr("D_Fick",temp)) {
        for (int i=0; i<nspecies; ++i) {
            D_Fick[i] = temp[i];
        }
    }

    diffusion_stencil_order = 1;
    pp.query("diffusion_stencil_order",diffusion_stencil_order);

    diffusion_verbose = 0;
    pp.query("diffusion_verbose",diffusion_verbose);

    diffusion_bottom_verbose = 0;
    pp.query("diffusion_bottom_verbose",diffusion_bottom_verbose);

    implicit_diffusion_rel_eps = 1.e-10;
    pp.query("implicit_diffusion_rel_eps",implicit_diffusion_rel_eps);

    implicit_diffusion_abs_eps = -1.;
    pp.query("implicit_diffusion_abs_eps",implicit_diffusion_abs_eps);

    return;
}
