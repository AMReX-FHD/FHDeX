namespace compressible {

  // Boundary conditions

  extern amrex::Vector<int> mass_bc_lo;   // mass boundary conditions
  extern amrex::Vector<int> mass_bc_hi;
  extern amrex::Vector<int> therm_bc_lo;  // temperature boundary conditions
  extern amrex::Vector<int> therm_bc_hi;
  extern amrex::Vector<int> vel_bc_lo;    // velocity field boundary conditions
  extern amrex::Vector<int> vel_bc_hi;
  extern amrex::Vector<amrex::Real> Yk_bc;   // wall mass fractions
  extern amrex::Vector<amrex::Real> Xk_bc;   // wall mole fractions

}
