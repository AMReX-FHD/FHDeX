namespace MFsurfchem {

    extern AMREX_GPU_MANAGED int n_ads_spec;
    extern AMREX_GPU_MANAGED int ads_wall_dir;

    extern AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_SPECIES> surfcov0;
    extern AMREX_GPU_MANAGED amrex::Real surf_site_num_dens;

    extern AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_SPECIES> ads_rate_const;
    extern AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_SPECIES> des_rate;

    extern AMREX_GPU_MANAGED int stoch_surfcov0;
    extern AMREX_GPU_MANAGED int stoch_MFsurfchem;

    extern AMREX_GPU_MANAGED amrex::Real k_beta;
    extern AMREX_GPU_MANAGED amrex::Real e_beta;

    extern AMREX_GPU_MANAGED int splitting_MFsurfchem;
    extern AMREX_GPU_MANAGED int conversion_MFsurfchem;

}
