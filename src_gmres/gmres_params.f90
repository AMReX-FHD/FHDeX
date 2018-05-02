module gmres_params_module

  use iso_c_binding, only: c_char
  use amrex_string_module, only: amrex_string_c_to_f

  implicit none

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
  
  namelist /probin_gmres/ precon_type
  namelist /probin_gmres/ visc_schur_approx
  namelist /probin_gmres/ p_norm_weight
  namelist /probin_gmres/ scale_factor
  namelist /probin_gmres/ mg_verbose
  namelist /probin_gmres/ cg_verbose
  namelist /probin_gmres/ mg_max_vcycles
  namelist /probin_gmres/ mg_minwidth
  namelist /probin_gmres/ mg_bottom_solver
  namelist /probin_gmres/ mg_nsmooths_down
  namelist /probin_gmres/ mg_nsmooths_up
  namelist /probin_gmres/ mg_nsmooths_bottom
  namelist /probin_gmres/ mg_max_bottom_nlevels
  namelist /probin_gmres/ mg_rel_tol
  namelist /probin_gmres/ mg_abs_tol
  namelist /probin_gmres/ stag_mg_verbosity
  namelist /probin_gmres/ stag_mg_max_vcycles
  namelist /probin_gmres/ stag_mg_minwidth
  namelist /probin_gmres/ stag_mg_bottom_solver
  namelist /probin_gmres/ stag_mg_nsmooths_down
  namelist /probin_gmres/ stag_mg_nsmooths_up
  namelist /probin_gmres/ stag_mg_nsmooths_bottom
  namelist /probin_gmres/ stag_mg_max_bottom_nlevels
  namelist /probin_gmres/ stag_mg_omega
  namelist /probin_gmres/ stag_mg_smoother
  namelist /probin_gmres/ stag_mg_rel_tol
  namelist /probin_gmres/ gmres_rel_tol
  namelist /probin_gmres/ gmres_abs_tol
  namelist /probin_gmres/ gmres_verbose
  namelist /probin_gmres/ gmres_max_outer
  namelist /probin_gmres/ gmres_max_inner
  namelist /probin_gmres/ gmres_max_iter
  namelist /probin_gmres/ gmres_min_iter
  namelist /probin_gmres/ gmres_spatial_order

contains

  ! read in fortran namelist into gmres_params_module
  subroutine read_gmres_params(inputs_file,length) bind(C, name="read_gmres_params")

    integer               , value         :: length
    character(kind=c_char), intent(in   ) :: inputs_file(length)

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

    ! read in probin_gmres namelist
    open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
    read(unit=100, nml=probin_gmres)
    close(unit=100)

  end subroutine read_gmres_params

  ! copy contents of gmres_params_module to C++ gmres namespace
  subroutine copy_gmres_params_to_c() &
                                     bind(C, name="copy_gmres_params_to_c")

  end subroutine copy_gmres_params_to_c 

end module gmres_params_module
