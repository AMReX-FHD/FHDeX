module gmres_namelist_module

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
  
  ! preconditioner type
  ! 1 = projection preconditioner
  !-1 = projection preconditioner with expensive pressure update
  ! 2 = lower triangular preconditioner
  !-2 = lower triangular preconditioner with negative sign
  ! 3 = upper triangular preconditioner
  !-3 = upper triangular preconditioner with negative sign
  ! 4 = Block diagonal preconditioner
  !-4 = Block diagonal preconditioner with negative sign
  namelist /gmres/ precon_type

  ! use the viscosity-based BFBt Schur complement (from Georg Stadler)
  namelist /gmres/ visc_schur_approx

  ! weighting of pressure when computing norms and inner products
  namelist /gmres/ p_norm_weight

  ! scale theta_alpha, beta, gamma, and b_u by this, and then scale x_p by the inverse
  namelist /gmres/ scale_factor

  ! MAC projection solver parameters:
  namelist /gmres/ mg_verbose            ! multigrid verbosity
  namelist /gmres/ cg_verbose            ! BiCGStab (mg_bottom_solver=1) verbosity
  namelist /gmres/ mg_max_vcycles        ! maximum number of V-cycles
  namelist /gmres/ mg_minwidth           ! length of box at coarsest multigrid level
  namelist /gmres/ mg_bottom_solver      ! bottom solver type
  ! 0 = smooths only, controlled by mg_nsmooths_bottom
  ! 1 = BiCGStab with agglomeration
  namelist /gmres/ mg_nsmooths_down      ! number of smooths at each level on the way down
  namelist /gmres/ mg_nsmooths_up        ! number of smooths at each level on the way up
  namelist /gmres/ mg_nsmooths_bottom    ! number of smooths at the bottom (only if mg_bottom_solver=0)
  namelist /gmres/ mg_max_bottom_nlevels ! for mg_bottom_solver=4, number of additional levels of multigrid
  namelist /gmres/ mg_rel_tol            ! rel_tol for Poisson solve
  namelist /gmres/ mg_abs_tol            ! abs_tol for Poisson solve

  ! Staggered multigrid solver parameters
  namelist /gmres/ stag_mg_verbosity          ! verbosity
  namelist /gmres/ stag_mg_max_vcycles        ! max number of v-cycles
  namelist /gmres/ stag_mg_minwidth           ! length of box at coarsest multigrid level
  namelist /gmres/ stag_mg_bottom_solver      ! bottom solver type
  ! 0 = smooths only, controlled by mg_nsmooths_bottom
  namelist /gmres/ stag_mg_nsmooths_down      ! number of smooths at each level on the way down
  namelist /gmres/ stag_mg_nsmooths_up        ! number of smooths at each level on the way up
  namelist /gmres/ stag_mg_nsmooths_bottom    ! number of smooths at the bottom
  namelist /gmres/ stag_mg_max_bottom_nlevels ! for stag_mg_bottom_solver=4, number of additional levels of multigrid
  namelist /gmres/ stag_mg_omega              ! weighted-jacobi omega coefficient
  namelist /gmres/ stag_mg_smoother           ! 0 = jacobi; 1 = 2*dm-color Gauss-Seidel
  namelist /gmres/ stag_mg_rel_tol            ! relative tolerance stopping criteria

  ! GMRES solver parameters
  namelist /gmres/ gmres_rel_tol         ! relative tolerance stopping criteria
  namelist /gmres/ gmres_abs_tol         ! absolute tolerance stopping criteria
  namelist /gmres/ gmres_verbose         ! gmres verbosity; if greater than 1, more residuals will be printed out
  namelist /gmres/ gmres_max_outer       ! max number of outer iterations
  namelist /gmres/ gmres_max_inner       ! max number of inner iterations, or restart number
  namelist /gmres/ gmres_max_iter        ! max number of gmres iterations
  namelist /gmres/ gmres_min_iter        ! min number of gmres iterations

  namelist /gmres/ gmres_spatial_order   ! spatial order of viscous and gradient operators in matrix "A"

contains

  ! read in fortran namelist into gmres_namelist_module
  subroutine read_gmres_namelist(inputs_file,length) bind(C, name="read_gmres_namelist")

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
    mg_bottom_solver = 0
    mg_nsmooths_down = 2
    mg_nsmooths_up = 2
    mg_nsmooths_bottom = 8
    mg_max_bottom_nlevels = 10
    mg_rel_tol = 1.d-9
    mg_abs_tol = 1.d-14
    stag_mg_verbosity = 0
    stag_mg_max_vcycles = 1
    stag_mg_minwidth = 2
    stag_mg_bottom_solver = 0
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

    ! read in gmres namelist
    open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
    read(unit=100, nml=gmres)
    close(unit=100)

  end subroutine read_gmres_namelist

  ! copy contents of gmres_namelist_module to C++ gmres namespace
  subroutine initialize_gmres_namespace(precon_type_in, visc_schur_approx_in, &
                                        p_norm_weight_in, &
                                        scale_factor_in, mg_verbose_in, cg_verbose_in, &
                                        mg_max_vcycles_in, mg_minwidth_in, &
                                        mg_bottom_solver_in, &
                                        mg_nsmooths_down_in, mg_nsmooths_up_in, &
                                        mg_nsmooths_bottom_in, mg_max_bottom_nlevels_in, &
                                        mg_rel_tol_in, mg_abs_tol_in, stag_mg_verbosity_in, &
                                        stag_mg_max_vcycles_in, stag_mg_minwidth_in, &
                                        stag_mg_bottom_solver_in, stag_mg_nsmooths_down_in, &
                                        stag_mg_nsmooths_up_in, stag_mg_nsmooths_bottom_in, &
                                        stag_mg_max_bottom_nlevels_in, stag_mg_omega_in, &
                                        stag_mg_smoother_in, stag_mg_rel_tol_in, &
                                        gmres_rel_tol_in, gmres_abs_tol_in, &
                                        gmres_verbose_in, &
                                        gmres_max_outer_in, gmres_max_inner_in, &
                                        gmres_max_iter_in, gmres_min_iter_in, &
                                        gmres_spatial_order_in) &
                                        bind(C, name="initialize_gmres_namespace")
    
    integer,          intent(inout) :: precon_type_in
    integer,          intent(inout) :: visc_schur_approx_in
    double precision, intent(inout) :: p_norm_weight_in
    double precision, intent(inout) :: scale_factor_in
    integer,          intent(inout) :: mg_verbose_in
    integer,          intent(inout) :: cg_verbose_in
    integer,          intent(inout) :: mg_max_vcycles_in
    integer,          intent(inout) :: mg_minwidth_in
    integer,          intent(inout) :: mg_bottom_solver_in
    integer,          intent(inout) :: mg_nsmooths_down_in
    integer,          intent(inout) :: mg_nsmooths_up_in
    integer,          intent(inout) :: mg_nsmooths_bottom_in
    integer,          intent(inout) :: mg_max_bottom_nlevels_in
    double precision, intent(inout) :: mg_rel_tol_in
    double precision, intent(inout) :: mg_abs_tol_in
    integer,          intent(inout) :: stag_mg_verbosity_in
    integer,          intent(inout) :: stag_mg_max_vcycles_in
    integer,          intent(inout) :: stag_mg_minwidth_in
    integer,          intent(inout) :: stag_mg_bottom_solver_in
    integer,          intent(inout) :: stag_mg_nsmooths_down_in
    integer,          intent(inout) :: stag_mg_nsmooths_up_in
    integer,          intent(inout) :: stag_mg_nsmooths_bottom_in
    integer,          intent(inout) :: stag_mg_max_bottom_nlevels_in
    double precision, intent(inout) :: stag_mg_omega_in
    integer,          intent(inout) :: stag_mg_smoother_in
    double precision, intent(inout) :: stag_mg_rel_tol_in
    double precision, intent(inout) :: gmres_rel_tol_in
    double precision, intent(inout) :: gmres_abs_tol_in
    integer,          intent(inout) :: gmres_verbose_in
    integer,          intent(inout) :: gmres_max_outer_in
    integer,          intent(inout) :: gmres_max_inner_in
    integer,          intent(inout) :: gmres_max_iter_in
    integer,          intent(inout) :: gmres_min_iter_in
    integer,          intent(inout) :: gmres_spatial_order_in

    precon_type_in = precon_type
    visc_schur_approx_in = visc_schur_approx
    p_norm_weight_in = p_norm_weight
    scale_factor_in = scale_factor
    mg_verbose_in = mg_verbose
    cg_verbose_in = cg_verbose
    mg_max_vcycles_in = mg_max_vcycles
    mg_minwidth_in = mg_minwidth
    mg_bottom_solver_in = mg_bottom_solver
    mg_nsmooths_down_in = mg_nsmooths_down
    mg_nsmooths_up_in = mg_nsmooths_up
    mg_nsmooths_bottom_in = mg_nsmooths_bottom
    mg_max_bottom_nlevels_in = mg_max_bottom_nlevels
    mg_rel_tol_in = mg_rel_tol
    mg_abs_tol_in = mg_abs_tol
    stag_mg_verbosity_in = stag_mg_verbosity
    stag_mg_max_vcycles_in = stag_mg_max_vcycles
    stag_mg_minwidth_in = stag_mg_minwidth
    stag_mg_bottom_solver_in = stag_mg_bottom_solver
    stag_mg_nsmooths_down_in = stag_mg_nsmooths_down
    stag_mg_nsmooths_up_in = stag_mg_nsmooths_up
    stag_mg_nsmooths_bottom_in = stag_mg_nsmooths_bottom
    stag_mg_max_bottom_nlevels_in = stag_mg_max_bottom_nlevels
    stag_mg_omega_in = stag_mg_omega
    stag_mg_smoother_in = stag_mg_smoother
    stag_mg_rel_tol_in = stag_mg_rel_tol
    gmres_rel_tol_in = gmres_rel_tol
    gmres_abs_tol_in = gmres_abs_tol
    gmres_verbose_in = gmres_verbose
    gmres_max_outer_in = gmres_max_outer
    gmres_max_inner_in = gmres_max_inner
    gmres_max_iter_in = gmres_max_iter
    gmres_min_iter_in = gmres_min_iter
    gmres_spatial_order_in = gmres_spatial_order

  end subroutine initialize_gmres_namespace

end module gmres_namelist_module
