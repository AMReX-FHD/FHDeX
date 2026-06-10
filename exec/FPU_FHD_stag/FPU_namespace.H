namespace FPU {

    extern AMREX_GPU_MANAGED int enable_fluctuations;
    extern AMREX_GPU_MANAGED int nonlinear_fhd;
    extern AMREX_GPU_MANAGED int diag_int;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, 3*3> A;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, 3*3> D;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, 3>   B;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, 3*3> R;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, 3*3> H1;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, 3*3> H2;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, 3*3> H3;
}
