module gmres_params_module

  implicit none

  ! Begin the declarations of the ParmParse parameters

  integer,          save :: precon_type
  integer,          save :: visc_schur_approx
  double precision, save :: p_norm_weight
  double precision, save :: scale_factor
  integer,          save :: mg_verbose
  integer,          save :: cg_verbose
  integer,          save :: mg_max_vcycles
  integer,          save :: mg_minwidth
  integer,          save :: mg_bottom_solver
  integer,          save :: mg_nsmooths_down
  integer,          save :: mg_nsmooths_up
  integer,          save :: mg_nsmooths_bottom
  integer,          save :: mg_max_bottom_nlevels
  double precision, save :: mg_rel_tol
  double precision, save :: mg_abs_tol
  integer,          save :: stag_mg_verbosity
  integer,          save :: stag_mg_max_vcycles
  integer,          save :: stag_mg_minwidth
  integer,          save :: stag_mg_bottom_solver
  integer,          save :: stag_mg_nsmooths_down
  integer,          save :: stag_mg_nsmooths_up
  integer,          save :: stag_mg_nsmooths_bottom
  integer,          save :: stag_mg_max_bottom_nlevels
  double precision, save :: stag_mg_omega
  integer,          save :: stag_mg_smoother
  double precision, save :: stag_mg_rel_tol
  double precision, save :: gmres_rel_tol
  double precision, save :: gmres_abs_tol
  integer,          save :: gmres_verbose
  integer,          save :: gmres_max_outer
  integer,          save :: gmres_max_inner
  integer,          save :: gmres_max_iter
  integer,          save :: gmres_min_iter
  integer,          save :: gmres_spatial_order

  ! End the declarations of the ParmParse parameters

contains

  subroutine read_gmres_params() bind(C, name="read_gmres_params")

    use amrex_parmparse_module, only: amrex_parmparse_build, &
                                      amrex_parmparse_destroy, &
                                      amrex_parmparse

    type (amrex_parmparse) :: pp

    ! allocate arrays (if needed)

    ! default values
    precon_type = 1
    visc_schur_approx = 0
    p_norm_weight = 1.d0
    scale_factor = 1.d0
    mg_verbose = 0
    cg_verbose = 0
    mg_max_vcycles = 1
    mg_minwidth = 2
    mg_bottom_solver = 4
    mg_nsmooths_down = 2
    mg_nsmooths_up = 2
    mg_nsmooths_bottom = 8
    mg_max_bottom_nlevels = 10
    mg_rel_tol = 1.d-9
    mg_abs_tol = 1.d-14
    stag_mg_verbosity = 0
    stag_mg_max_vcycles = 1
    stag_mg_minwidth = 2
    stag_mg_bottom_solver = 4
    stag_mg_nsmooths_down = 2
    stag_mg_nsmooths_up = 2
    stag_mg_nsmooths_bottom = 8
    stag_mg_max_bottom_nlevels = 10
    stag_mg_omega = 1.d0
    stag_mg_smoother = 1
    stag_mg_rel_tol = 1.d-9
    gmres_rel_tol = 1.d-9
    gmres_abs_tol = 0.d0
    gmres_verbose = 1
    gmres_max_outer = 20
    gmres_max_inner = 5
    gmres_max_iter = 100
    gmres_min_iter = 1
    gmres_spatial_order = 2

    ! read in from inputs file
    call amrex_parmparse_build(pp)

    call pp%query("precon_type",precon_type);
    call pp%query("visc_schur_approx",visc_schur_approx);
    call pp%query("p_norm_weight",p_norm_weight);
    call pp%query("scale_factor",scale_factor);
    call pp%query("mg_verbose",mg_verbose);
    call pp%query("cg_verbose",cg_verbose);             
    call pp%query("mg_max_vcycles",mg_max_vcycles);
    call pp%query("mg_minwidth",mg_minwidth);
    call pp%query("mg_bottom_solver",mg_bottom_solver);
    call pp%query("mg_nsmooths_down",mg_nsmooths_down);
    call pp%query("mg_nsmooths_up",mg_nsmooths_up);
    call pp%query("mg_nsmooths_bottom",mg_nsmooths_bottom);
    call pp%query("mg_max_bottom_nlevels",mg_max_bottom_nlevels);
    call pp%query("mg_rel_tol",mg_rel_tol);
    call pp%query("mg_abs_tol",mg_abs_tol);
    call pp%query("stag_mg_verbosity",stag_mg_verbosity);
    call pp%query("stag_mg_max_vcycles",stag_mg_max_vcycles);
    call pp%query("stag_mg_minwidth",stag_mg_minwidth);
    call pp%query("stag_mg_bottom_solver",stag_mg_bottom_solver);
    call pp%query("stag_mg_nsmooths_down",stag_mg_nsmooths_down);
    call pp%query("stag_mg_nsmooths_up",stag_mg_nsmooths_up);
    call pp%query("stag_mg_nsmooths_bottom",stag_mg_nsmooths_bottom);
    call pp%query("stag_mg_max_bottom_nlevels",stag_mg_max_bottom_nlevels);
    call pp%query("stag_mg_omega",stag_mg_omega);
    call pp%query("stag_mg_smoother",stag_mg_smoother);
    call pp%query("stag_mg_rel_tol",stag_mg_rel_tol);
    call pp%query("gmres_rel_tol",gmres_rel_tol);
    call pp%query("gmres_abs_tol",gmres_abs_tol);
    call pp%query("gmres_verbose",gmres_verbose);
    call pp%query("gmres_max_outer",gmres_max_outer);
    call pp%query("gmres_max_inner",gmres_max_inner);
    call pp%query("gmres_max_iter",gmres_max_iter);
    call pp%query("gmres_min_iter",gmres_min_iter);
    call pp%query("gmres_spatial_order",gmres_spatial_order);

    call amrex_parmparse_destroy(pp)

  end subroutine read_gmres_params

  subroutine finalize_gmres_params() bind(C, name="finalize_gmres_params")

    ! deallocate parameters to quiet valgrind

  end subroutine finalize_gmres_params

end module gmres_params_module
