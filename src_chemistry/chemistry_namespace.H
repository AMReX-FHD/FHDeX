namespace chemistry {

    extern AMREX_GPU_MANAGED int nreaction;

    // from the fortran code, stoich_coeffs_R = stoichiometric_factors(spec,1,reac)
    // from the fortran code, stoich_coeffs_P = stoichiometric_factors(spec,2,reac)
    // stoich_coeffs_PR = stoich_coeffs_P - stoich_coeffs_R
    extern AMREX_GPU_MANAGED Array2D<int,0, MAX_REACTION,0, MAX_SPECIES> stoich_coeffs_R;
    extern AMREX_GPU_MANAGED Array2D<int,0, MAX_REACTION,0, MAX_SPECIES> stoich_coeffs_P;
    extern AMREX_GPU_MANAGED Array2D<int,0, MAX_REACTION,0, MAX_SPECIES> stoich_coeffs_PR;

    // reaction rate constant for each reaction (assuming Law of Mass Action holds)
    // using rate_multiplier, reaction rates can be changed by the same factor
    // if include_discrete_LMA_correction, n^2 and n^3 in rate expressions become
    // n*(n-1/dv) and n*(n-1/dv)*(n-2/dv).
    extern AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_REACTION> rate_const;
    extern AMREX_GPU_MANAGED amrex::Real rate_multiplier;
    extern AMREX_GPU_MANAGED int include_discrete_LMA_correction;

    // if n is positive, exclude species n (=solvent) when computing reaction rates
    // in this case, the concentration of the solvent is assumed to be constant,
    // which should be reflected on rate constants.
    // if 0, no species is excluded
    // e.g. U + S -> 2U, if exclude_solvent_comput_rates=0, rate=k*n_U*n_S
    //                   if exclude_solvent_comput_rates=2, rate=k_new*n_U where k_new=k*n_S
    extern AMREX_GPU_MANAGED int exclude_solvent_comput_rates;

    // from the fortran code this was use_Poisson_rng (0=CLE; 1=tau leaping; -1=deterministic; 2=SSA)
    // here it's being used as reaction_type (0=deterministic; 1=CLE; 2=SSA; 3=tau leap)
    extern AMREX_GPU_MANAGED int reaction_type;

    // use mole fraction based LMA
    extern AMREX_GPU_MANAGED int use_mole_frac_LMA;

    // specific to compressible codes
    extern AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_REACTION> alpha_param;
    extern AMREX_GPU_MANAGED GpuArray<amrex::Real, MAX_REACTION> beta_param;
    extern AMREX_GPU_MANAGED amrex::Real T0_chem;
}
