#include "gmres_functions.H"

#include "AMReX_ParmParse.H"

int         gmres::precon_type;
int         gmres::visc_schur_approx;
amrex::Real gmres::p_norm_weight;
amrex::Real gmres::scale_factor;
int         gmres::mg_verbose;
int         gmres::cg_verbose;
int         gmres::mg_max_vcycles;
int         gmres::mg_minwidth;
int         gmres::mg_bottom_solver;
int         gmres::mg_nsmooths_down;
int         gmres::mg_nsmooths_up;
int         gmres::mg_nsmooths_bottom;
int         gmres::mg_max_bottom_nlevels;
amrex::Real gmres::mg_rel_tol;
amrex::Real gmres::mg_abs_tol;
int         gmres::stag_mg_verbosity;
int         gmres::stag_mg_max_vcycles;
int         gmres::stag_mg_minwidth;
int         gmres::stag_mg_bottom_solver;
int         gmres::stag_mg_nsmooths_down;
int         gmres::stag_mg_nsmooths_up;
int         gmres::stag_mg_nsmooths_bottom;
int         gmres::stag_mg_max_bottom_nlevels;
amrex::Real gmres::stag_mg_omega;
int         gmres::stag_mg_smoother;
amrex::Real gmres::stag_mg_rel_tol;
amrex::Real gmres::gmres_rel_tol;
amrex::Real gmres::gmres_abs_tol;
int         gmres::gmres_verbose;
int         gmres::gmres_max_outer;
int         gmres::gmres_max_inner;
int         gmres::gmres_max_iter;
int         gmres::gmres_min_iter;
int         gmres::gmres_spatial_order;

void InitializeGmresNamespace() {

    // specify default values first, then read in values from inputs file

    // preconditioner type
    // 1 = projection preconditioner
    //-1 = projection preconditioner with expensive pressure update
    // 2 = lower triangular preconditioner
    //-2 = lower triangular preconditioner with negative sign
    // 3 = upper triangular preconditioner
    //-3 = upper triangular preconditioner with negative sign
    // 4 = Block diagonal preconditioner
    //-4 = Block diagonal preconditioner with negative sign
    precon_type = 1;

    // use the viscosity-based BFBt Schur complement (from Georg Stadler)
    visc_schur_approx = 0;

    // weighting of pressure when computing norms and inner products
    p_norm_weight = 1.;

    // scale theta_alpha, beta, gamma, and b_u by this, and then scale x_p by the inverse
    scale_factor = 1.;

    // MAC projection solver parameters:
    mg_verbose = 0;             // multigrid verbosity
    cg_verbose = 0;             // BiCGStab (mg_bottom_solver=1) verbosity
    mg_max_vcycles = 1;         // maximum number of V-cycles
    mg_minwidth = 2;            // length of box at coarsest multigrid level
    mg_bottom_solver = 0;       // bottom solver type
    // 0 = smooths only, controlled by mg_nsmooths_bottom
    // 1 = BiCGStab with agglomeration
    mg_nsmooths_down = 2;       // number of smooths at each level on the way down
    mg_nsmooths_up = 2;         // number of smooths at each level on the way up
    mg_nsmooths_bottom = 8;     // number of smooths at the bottom (only if mg_bottom_solver=0)
    mg_max_bottom_nlevels = 10; // for mg_bottom_solver 4, number of additional levels of multigrid
    mg_rel_tol = 1.e-9;         // rel_tol for Poisson solve
    mg_abs_tol = 1.e-14;        // abs_tol for Poisson solve

    // Staggered multigrid solver parameters
    stag_mg_verbosity = 0;           // verbosity
    stag_mg_max_vcycles = 1;         // max number of v-cycles
    stag_mg_minwidth = 2;            // length of box at coarsest multigrid level
    stag_mg_bottom_solver = 0;       // bottom solver type
    // 0 = smooths only, controlled by mg_nsmooths_bottom
    stag_mg_nsmooths_down = 2;       // number of smooths at each level on the way down
    stag_mg_nsmooths_up = 2;         // number of smooths at each level on the way up
    stag_mg_nsmooths_bottom = 8;     // number of smooths at the bottom
    stag_mg_max_bottom_nlevels = 10; // for stag_mg_bottom_solver 4, number of additional levels of multigrid
    stag_mg_omega = 1.;              // weightee-jacobi omega coefficient
    stag_mg_smoother = 1;            // 0 = jacobi; 1 = 2*dm-color Gauss-Seidel
    stag_mg_rel_tol = 1.e-9;         // relative tolerance stopping criteria

    // GMRES solver parameters
    gmres_rel_tol = 1.e-9;     // relative tolerance stopping criteria
    gmres_abs_tol = 0.;        // absolute tolerance stopping criteria
    gmres_verbose = 1;         // gmres verbosity; if greater than 1, more residuals will be printed out
    gmres_max_outer = 20;      // max number of outer iterations
    gmres_max_inner = 5;       // max number of inner iterations, or restart number
    gmres_max_iter = 100;      // max number of gmres iterations
    gmres_min_iter = 1;        // min number of gmres iterations

    gmres_spatial_order = 2;   // spatial order of viscous and gradient operators in matrix "A"

    ParmParse pp;

    // pp.query searches for optional parameters
    // pp.get aborts if the parameter is not found
    // pp.getarr and queryarr("string",inputs,start_indx,count); can be used for arrays

    pp.query("precon_type",precon_type);
    pp.query("visc_schur_approx",visc_schur_approx);
    pp.query("p_norm_weight",p_norm_weight);
    pp.query("scale_factor",scale_factor);
    pp.query("mg_verbose",mg_verbose);
    pp.query("cg_verbose",cg_verbose);
    pp.query("mg_max_vcycles",mg_max_vcycles);
    pp.query("mg_minwidth",mg_minwidth);
    pp.query("mg_bottom_solver",mg_bottom_solver);
    pp.query("mg_nsmooths_down",mg_nsmooths_down);
    pp.query("mg_nsmooths_up",mg_nsmooths_up);
    pp.query("mg_nsmooths_bottom",mg_nsmooths_bottom);
    pp.query("mg_max_bottom_nlevels",mg_max_bottom_nlevels);
    pp.query("mg_rel_tol",mg_rel_tol);
    pp.query("mg_abs_tol",mg_abs_tol);
    pp.query("stag_mg_verbosity",stag_mg_verbosity);
    pp.query("stag_mg_max_vcycles",stag_mg_max_vcycles);
    pp.query("stag_mg_minwidth",stag_mg_minwidth);
    pp.query("stag_mg_bottom_solver",stag_mg_bottom_solver);
    pp.query("stag_mg_nsmooths_down",stag_mg_nsmooths_down);
    pp.query("stag_mg_nsmooths_up",stag_mg_nsmooths_up);
    pp.query("stag_mg_nsmooths_bottom",stag_mg_nsmooths_bottom);
    pp.query("stag_mg_max_bottom_nlevels",stag_mg_max_bottom_nlevels);
    pp.query("stag_mg_omega",stag_mg_omega);
    pp.query("stag_mg_smoother",stag_mg_smoother);
    pp.query("stag_mg_rel_tol",stag_mg_rel_tol);
    pp.query("gmres_rel_tol",gmres_rel_tol);
    pp.query("gmres_abs_tol",gmres_abs_tol);
    pp.query("gmres_verbose",gmres_verbose);
    pp.query("gmres_max_outer",gmres_max_outer);
    pp.query("gmres_max_inner",gmres_max_inner);
    pp.query("gmres_max_iter",gmres_max_iter);
    pp.query("gmres_min_iter",gmres_min_iter);
    pp.query("gmres_spatial_order",gmres_spatial_order);

}
