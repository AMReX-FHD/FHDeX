namespace chemistry {
#define MAX_REACTION 5

    extern AMREX_GPU_MANAGED int nreaction;

    extern AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_REACTION> rate_const; 
    extern AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_REACTION> alpha_param; 
    extern AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_REACTION> beta_param; 
    extern AMREX_GPU_MANAGED amrex::Real T0_chem;

    extern AMREX_GPU_MANAGED Array2D<int,0, MAX_REACTION,0, MAX_SPECIES> stoich_coeffs_R; 
    extern AMREX_GPU_MANAGED Array2D<int,0, MAX_REACTION,0, MAX_SPECIES> stoich_coeffs_P; 
    extern AMREX_GPU_MANAGED Array2D<int,0, MAX_REACTION,0, MAX_SPECIES> stoich_coeffs_PR; 

    extern AMREX_GPU_MANAGED int reaction_type;
}
