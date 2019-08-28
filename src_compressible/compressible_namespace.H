namespace compressible {

  // Boundary conditions

  extern amrex::Vector<int> bc_mass_lo;   // mass boundary conditions
  extern amrex::Vector<int> bc_mass_hi;
  extern amrex::Vector<int> bc_therm_lo;  // temperature boundary conditions
  extern amrex::Vector<int> bc_therm_hi;
  extern amrex::Vector<int> bc_vel_lo;    // velocity field boundary conditions
  extern amrex::Vector<int> bc_vel_hi;
  extern amrex::Vector<amrex::Real> bc_Yk;   // wall mass fractions
  extern amrex::Vector<amrex::Real> bc_Xk;   // wall mole fractions

}
