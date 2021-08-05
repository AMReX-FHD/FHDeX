int chemistry::nspecies;
int chemistry::nreaction;

std::array<amrex::Real, MAX_REACTION> chemistry::rate_const; 

std::array<std::array<int, MAX_SPECIES>, MAX_REACTION> chemistry::stoich_coeffs_R; 
std::array<std::array<int, MAX_SPECIES>, MAX_REACTION> chemistry::stoich_coeffs_P; 
std::array<std::array<int, MAX_SPECIES>, MAX_REACTION> chemistry::stoich_coeffs_PR; 

int chemistry::reaction_type;
