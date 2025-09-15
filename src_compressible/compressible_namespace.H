namespace compressible {

    extern AMREX_GPU_MANAGED int transport_type;
    extern AMREX_GPU_MANAGED int membrane_cell;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> transmission;
    extern AMREX_GPU_MANAGED int do_1D;
    extern AMREX_GPU_MANAGED int do_2D;
    extern AMREX_GPU_MANAGED int all_correl;
    extern AMREX_GPU_MANAGED int nspec_surfcov;
    extern AMREX_GPU_MANAGED int turbRestartRun;
    extern AMREX_GPU_MANAGED bool do_reservoir;
    extern AMREX_GPU_MANAGED amrex::Real zeta_ratio;
    extern AMREX_GPU_MANAGED amrex::Real p_init;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> Xk_init;
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM + 2 + MAX_SPECIES> cons_SF_pairA;
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM + 2 + MAX_SPECIES> cons_SF_pairB;
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM + 3 + 2*MAX_SPECIES> prim_SF_pairA;
    extern AMREX_GPU_MANAGED amrex::GpuArray<int, AMREX_SPACEDIM + 3 + 2*MAX_SPECIES> prim_SF_pairB;
    extern AMREX_GPU_MANAGED int do_SF_pair_cons;
    extern AMREX_GPU_MANAGED int do_SF_pair_prim;
}

