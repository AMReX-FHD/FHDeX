namespace chemistry {

#define MAX_SPECIES 5
#define MAX_REACTION 5

    extern int nspecies;
    extern int nreaction;

    extern std::array<amrex::Real, MAX_REACTION> rate_const; 
    
    extern std::array<std::array<int, MAX_SPECIES>, MAX_REACTION> stoich_coeffs_R; 
    extern std::array<std::array<int, MAX_SPECIES>, MAX_REACTION> stoich_coeffs_P; 
    extern std::array<std::array<int, MAX_SPECIES>, MAX_REACTION> stoich_coeffs_PR; 

    extern int reaction_type;
}
