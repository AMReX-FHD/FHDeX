AMREX_GPU_MANAGED int chemistry::nreaction;

AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_REACTION> chemistry::rate_const;

AMREX_GPU_MANAGED Array2D<int,0, MAX_REACTION,0, MAX_SPECIES> chemistry::stoich_coeffs_R; 
AMREX_GPU_MANAGED Array2D<int,0, MAX_REACTION,0, MAX_SPECIES> chemistry::stoich_coeffs_P;
AMREX_GPU_MANAGED Array2D<int,0, MAX_REACTION,0, MAX_SPECIES> chemistry::stoich_coeffs_PR;

AMREX_GPU_MANAGED int chemistry::reaction_type;
