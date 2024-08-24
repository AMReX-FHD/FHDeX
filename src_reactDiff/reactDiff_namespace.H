namespace reactDiff {

    // see reactDiff_functions.cpp for descriptions and default values
    extern AMREX_GPU_MANAGED int temporal_integrator;
    extern AMREX_GPU_MANAGED int reactDiff_diffusion_type;
    extern AMREX_GPU_MANAGED int reactDiff_reaction_type;
    extern AMREX_GPU_MANAGED int midpoint_stoch_flux_type;
    extern AMREX_GPU_MANAGED int avg_type;
    extern AMREX_GPU_MANAGED int inhomogeneous_bc_fix;
    extern AMREX_GPU_MANAGED amrex::Real volume_factor;
    extern AMREX_GPU_MANAGED Array2D<amrex::Real, 0, 2 ,0, MAX_SPECIES> n_init_in;
    extern AMREX_GPU_MANAGED int model_file_init;
    extern AMREX_GPU_MANAGED int integer_populations;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> D_Fick;
    extern AMREX_GPU_MANAGED int diffusion_stencil_order;
    extern AMREX_GPU_MANAGED int diffusion_verbose;
    extern AMREX_GPU_MANAGED int diffusion_bottom_verbose;
    extern AMREX_GPU_MANAGED amrex::Real implicit_diffusion_rel_eps;
    extern AMREX_GPU_MANAGED amrex::Real implicit_diffusion_abs_eps;
}
