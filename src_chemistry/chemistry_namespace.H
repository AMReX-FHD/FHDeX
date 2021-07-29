namespace chemistry {
#define MAX_REACTION 5

    extern AMREX_GPU_MANAGED int nreaction;

    extern AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_REACTION> rate_const; 
    extern AMREX_GPU_MANAGED Array2D<int,0, MAX_SPECIES,0, MAX_REACTION> stoich_coeffs_R; 
    extern AMREX_GPU_MANAGED Array2D<int,0, MAX_SPECIES,0, MAX_REACTION> stoich_coeffs_P; 
    extern AMREX_GPU_MANAGED Array2D<int,0, MAX_SPECIES,0, MAX_REACTION> stoich_coeffs_PR; 

    extern AMREX_GPU_MANAGED int reaction_type;
}
