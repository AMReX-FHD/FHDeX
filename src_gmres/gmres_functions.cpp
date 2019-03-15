#include "gmres_functions.H"
#include "gmres_functions_F.H"
#include "gmres_namespace.H"

using namespace gmres;

void InitializeGmresNamespace() {

    initialize_gmres_namespace(&precon_type, &visc_schur_approx, &p_norm_weight,
                               &scale_factor,
                               &mg_verbose, &cg_verbose, &mg_max_vcycles, &mg_minwidth,
                               &mg_bottom_solver, &mg_nsmooths_down, &mg_nsmooths_up,
                               &mg_nsmooths_bottom, &mg_max_bottom_nlevels, &mg_rel_tol,
                               &mg_abs_tol, &stag_mg_verbosity, &stag_mg_max_vcycles,
                               &stag_mg_minwidth, &stag_mg_bottom_solver,
                               &stag_mg_nsmooths_down, &stag_mg_nsmooths_up,
                               &stag_mg_nsmooths_bottom, &stag_mg_max_bottom_nlevels,
                               &stag_mg_omega, &stag_mg_smoother, &stag_mg_rel_tol,
                               &gmres_rel_tol, &gmres_abs_tol, &gmres_verbose,
                               &gmres_max_outer, &gmres_max_inner, &gmres_max_iter,
                               &gmres_min_iter, &gmres_spatial_order);

}
