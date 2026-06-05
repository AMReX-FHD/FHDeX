namespace FPU {

    extern AMREX_GPU_MANAGED int enable_fluctuations;
    extern int diag_int;
    extern AMREX_GPU_MANAGED amrex::Real r0;
    extern AMREX_GPU_MANAGED amrex::Real p0;
    extern AMREX_GPU_MANAGED amrex::Real e0;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, 3*3> A;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, 3*3> D;
    extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real, 3>   B;
}
