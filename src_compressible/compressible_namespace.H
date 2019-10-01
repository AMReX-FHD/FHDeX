namespace compressible {

    // Boundary conditions
    extern amrex::Vector<amrex::Real> bc_Yk;   // wall mass fractions
    extern amrex::Vector<amrex::Real> bc_Xk;   // wall mole fractions

    extern int plot_means; // write out means to plotfile
    extern int plot_vars;  // write out variances to plotfile
    
}
