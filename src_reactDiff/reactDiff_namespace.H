namespace reactDiff {

    // 0=D + R (first-order splitting)
    // 1=(1/2)R + D + (1/2)R (Strang option 1)
    // 2=(1/2)D + R + (1/2)D (Strang option 2)
    // -1=unsplit forward Euler
    // -2=unsplit explicit midpoint 
    // -3=unsplit multinomial diffusion
    // -4=unsplit implicit midpoint
    extern AMREX_GPU_MANAGED int temporal_integrator;

    // only used for split schemes (temporal_integrator>=0)
    // 0=explicit trapezoidal predictor/corrector
    // 1=Crank-Nicolson semi-implicit
    // 2=explicit midpoint
    // 3=multinomial diffusion
    // 4=forward Euler  
    extern AMREX_GPU_MANAGED int reactDiff_diffusion_type;

    // only used for split schemes (temporal_integrator>=0)
    // 0=first-order (deterministic, tau leaping, CLE, or SSA)
    // 1=second-order (determinisitc, tau leaping, or CLE only)
    extern AMREX_GPU_MANAGED int reactDiff_reaction_type;

    // only used for midpoint diffusion schemes (split as well as unsplit)
    // corrector formulation of noise
    // 1 = K(nold) * W1 + K(nold)         * W2
    // 2 = K(nold) * W1 + K(npred)        * W2
    // 3 = K(nold) * W1 + K(2*npred-nold) * W2
    extern AMREX_GPU_MANAGED int midpoint_stoch_flux_type;

    // how to compute n on faces for stochastic weighting
    // 1=arithmetic (with C0-Heaviside), 2=geometric, 3=harmonic
    // 10=arithmetic average with discontinuous Heaviside function
    // 11=arithmetic average with C1-smoothed Heaviside function
    // 12=arithmetic average with C2-smoothed Heaviside function
    extern AMREX_GPU_MANAGED int avg_type;

    // use the Einkemmer boundary condition fix (split schemes only)
    extern AMREX_GPU_MANAGED int inhomogeneous_bc_fix;

    // volume multiplier (dv = product(dx(1:MAX_SPACEDIM))*volume_factor)
    // only really intended for 3D since in 2D one can control the cell depth
    extern AMREX_GPU_MANAGED amrex::Real volume_factor;

    // initial values to be used in init_n.f90
    extern AMREX_GPU_MANAGED Array2D<amrex::Real, 0, 2 ,0, MAX_SPECIES> n_init_in;

    // initialize from model file
    extern AMREX_GPU_MANAGED int model_file_init;

    // initialize with all number of molecules strictly integer
    extern AMREX_GPU_MANAGED int integer_populations;

    // Fickian diffusion coeffs
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> D_Fick;

    // diffusion boundary stencil order
    extern AMREX_GPU_MANAGED int diffusion_stencil_order;

    // implicit diffusion solve verbosity
    extern AMREX_GPU_MANAGED int diffusion_verbose;

    // implicit diffusion solve bottom solver verbosity
    extern AMREX_GPU_MANAGED int diffusion_bottom_verbose;

    // relative eps for implicit diffusion solve
    extern AMREX_GPU_MANAGED amrex::Real implicit_diffusion_rel_eps;

    // absolute eps for implicit diffusion solve
    extern AMREX_GPU_MANAGED amrex::Real implicit_diffusion_abs_eps;

}
