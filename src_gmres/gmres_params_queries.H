// read in common parameters
{
    // ParmParse is way of reading inputs from the inputs file
    ParmParse pp;

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
