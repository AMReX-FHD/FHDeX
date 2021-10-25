namespace compressible {

    extern AMREX_GPU_MANAGED int transport_type;
    extern AMREX_GPU_MANAGED int membrane_cell;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, MAX_SPECIES> transmission;
    extern AMREX_GPU_MANAGED int do_1D;
    extern AMREX_GPU_MANAGED int do_2D;
    extern AMREX_GPU_MANAGED int all_correl;

}

