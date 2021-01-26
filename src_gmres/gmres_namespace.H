namespace gmres {

    // preconditioner type
    // 1 = projection preconditioner
    //-1 = projection preconditioner with expensive pressure update
    // 2 = lower triangular preconditioner
    //-2 = lower triangular preconditioner with negative sign
    // 3 = upper triangular preconditioner
    //-3 = upper triangular preconditioner with negative sign
    // 4 = Block diagonal preconditioner
    //-4 = Block diagonal preconditioner with negative sign
    extern int         precon_type;

    // use the viscosity-based BFBt Schur complement (from Georg Stadler)
    extern int         visc_schur_approx;

    // weighting of pressure when computing norms and inner products
    extern amrex::Real p_norm_weight;

    // scale theta_alpha, beta, gamma, and b_u by this, and then scale x_p by the inverse
    extern amrex::Real scale_factor;

    // MAC projection solver parameters:
    extern int         mg_verbose;            // multigrid verbosity
    extern int         cg_verbose;            // BiCGStab (mg_bottom_solver=1) verbosity
    extern int         mg_max_vcycles;        // maximum number of V-cycles
    extern int         mg_minwidth;           // length of box at coarsest multigrid level
    extern int         mg_bottom_solver;      // bottom solver type
    // 0 = smooths only, controlled by mg_nsmooths_bottom
    // 1 = BiCGStab with agglomeration
    extern int         mg_nsmooths_down;      // number of smooths at each level on the way down
    extern int         mg_nsmooths_up;        // number of smooths at each level on the way up
    extern int         mg_nsmooths_bottom;    // number of smooths at the bottom (only if mg_bottom_solver=0)
    extern int         mg_max_bottom_nlevels; // for mg_bottom_solver=4, number of additional levels of multigrid
    extern amrex::Real mg_rel_tol;            // rel_tol for Poisson solve
    extern amrex::Real mg_abs_tol;            // abs_tol for Poisson solve

    // Staggered multigrid solver parameters
    extern int         stag_mg_verbosity;          // verbosity
    extern int         stag_mg_max_vcycles;        // max number of v-cycles
    extern int         stag_mg_minwidth;           // length of box at coarsest multigrid level
    extern int         stag_mg_bottom_solver;      // bottom solver type
    // 0 = smooths only, controlled by mg_nsmooths_bottom
    extern int         stag_mg_nsmooths_down;      // number of smooths at each level on the way down
    extern int         stag_mg_nsmooths_up;        // number of smooths at each level on the way up
    extern int         stag_mg_nsmooths_bottom;    // number of smooths at the bottom
    extern int         stag_mg_max_bottom_nlevels; // for stag_mg_bottom_solver=4, number of additional levels of multigrid
    extern amrex::Real stag_mg_omega;              // weighted-jacobi omega coefficient
    extern int         stag_mg_smoother;           // 0 = jacobi; 1 = 2*dm-color Gauss-Seidel
    extern amrex::Real stag_mg_rel_tol;            // relative tolerance stopping criteria

    // GMRES solver parameters
    extern amrex::Real gmres_rel_tol;         // relative tolerance stopping criteria
    extern amrex::Real gmres_abs_tol;         // absolute tolerance stopping criteria
    extern int         gmres_verbose;         // gmres verbosity; if greater than 1, more residuals will be printed out
    extern int         gmres_max_outer;       // max number of outer iterations
    extern int         gmres_max_inner;       // max number of inner iterations, or restart number
    extern int         gmres_max_iter;        // max number of gmres iterations
    extern int         gmres_min_iter;        // min number of gmres iterations

    extern int         gmres_spatial_order;   // spatial order of viscous and gradient operators in matrix "A"
}

