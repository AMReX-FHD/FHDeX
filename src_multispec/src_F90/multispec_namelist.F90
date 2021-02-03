module multispec_namelist_module

  use iso_c_binding, only: c_char
  use amrex_string_module, only: amrex_string_c_to_f, amrex_string_f_to_c
  use common_namelist_module, only: MAX_SPECIES

  implicit none

  integer, parameter :: MAX_ELEMENT=MAX_SPECIES*(MAX_SPECIES-1)/2

  integer,            save :: inverse_type
  integer,            save :: temp_type
  integer,            save :: chi_iterations
  double precision,   save :: start_time 
  double precision,   save :: Dbar(MAX_ELEMENT)
  double precision,   save :: Dtherm(MAX_SPECIES)
  double precision,   save :: H_offdiag(MAX_ELEMENT)
  double precision,   save :: H_diag(MAX_SPECIES)
  double precision,   save :: fraction_tolerance
  integer,            save :: correct_flux
  integer,            save :: print_error_norms
  integer,            save :: is_nonisothermal
  integer,            save :: is_ideal_mixture
  integer,            save :: use_lapack
  integer,            save :: use_multiphase
  double precision,   save :: kc_tension
  double precision,   save :: alpha_gex
  integer,            save :: n_gex
  double precision,   save :: c_init_1(MAX_SPECIES)
  double precision,   save :: c_init_2(MAX_SPECIES)
  
  integer,            save :: midpoint_stoch_mass_flux_type
  integer,            save :: avg_type
  integer,            save :: mixture_type

  ! charged
  integer,            save :: use_charged_fluid, print_debye_len
  double precision,   save :: dielectric_const
  integer,            save :: dielectric_type
  double precision,   save :: charge_per_mass(MAX_SPECIES)
  double precision,   save :: theta_pot
  integer,            save :: num_pot_iters
  double precision,   save :: dpdt_factor
  double precision,   save :: relxn_param_charge
  integer,            save :: E_ext_type
  double precision,   save :: E_ext_value(1:3)
  integer,            save :: electroneutral
  integer,            save :: induced_charge_eo 
  integer,            save :: zero_eps_on_wall_type
  integer,            save :: zero_charge_on_wall_type
  double precision,   save :: zero_eps_on_wall_left_end, zero_eps_on_wall_right_start
  integer,            save :: bc_function_type
  double precision,   save :: L_pos, L_trans, L_zero

  ! Physical properties:
  !----------------------
  namelist /multispec/ fraction_tolerance ! For roundoff errors in mass and mole fractions
  namelist /multispec/ start_time
  namelist /multispec/ inverse_type       ! Only for LAPACK:  1=inverse, 2=pseudo inverse
  namelist /multispec/ correct_flux       ! Manually ensure mass is conserved to roundoff 
  namelist /multispec/ print_error_norms   
  namelist /multispec/ is_ideal_mixture   ! If T assume Gamma=I (H=0) and simplify
  namelist /multispec/ is_nonisothermal   ! If T Soret effect will be included
  namelist /multispec/ use_lapack         ! Use LAPACK or iterative method for diffusion matrix (recommend False)
  namelist /multispec/ use_multiphase     ! for RTIL
  namelist /multispec/ kc_tension         ! for RTIL
  namelist /multispec/ alpha_gex          ! for RTIL
  namelist /multispec/ n_gex              ! for RTIL
  namelist /multispec/ chi_iterations     ! number of iterations used in Dbar2chi_iterative

  ! Initial and boundary conditions 
  !----------------------

  namelist /multispec/ temp_type  ! for initializing temperature
  namelist /multispec/ c_init_1   ! initial values for c
  namelist /multispec/ c_init_2
  
  ! Thermodynamic and transport properties:
  !----------------------

  ! These are lower-triangules of symmetric matrices represented as vectors
  ! Number of elements is (nspecies*(nspecies-1)/2)
  ! The values are read row by row starting from top going down (this allows easy addition/deletion of new species/rows)
  ! So D_12; D_13, D_23; D_14, D_24, D_34; ...
  namelist /multispec/ Dbar       ! Maxwell-Stefan diffusion constant  
  namelist /multispec/ Dtherm     ! thermo-diffusion coefficients, only differences among elements matter
  namelist /multispec/ H_offdiag
  namelist /multispec/ H_diag     ! Diagonal of H=d^2F/dx^2, these are vectors of length nspecies

  ! Algorithm control
  !----------------------
  namelist /multispec/ midpoint_stoch_mass_flux_type  ! 1 = Strato
                                                      ! 2 = Ito

  namelist /multispec/ avg_type   ! how to compute stochastc_mass_fluxdiv
                                  ! 1=arithmetic (with C0-Heaviside), 2=geometric, 3=harmonic
                                  ! 10=arithmetic average with discontinuous Heaviside function
                                  ! 11=arithmetic average with C1-smoothed Heaviside function
                                  ! 12=arithmetic average with C2-smoothed Heaviside function

  namelist /multispec/ mixture_type ! Model for how transport and thermodynamic coefficients depend on composition
                                    ! See compute_mixture_properties.f90 for values supported at present
                                    ! The default mixture_type=0 means no dependence on composition

  ! for charged fluid
  namelist /multispec/ use_charged_fluid
  namelist /multispec/ print_debye_len
  namelist /multispec/ dielectric_const
  namelist /multispec/ dielectric_type
  namelist /multispec/ charge_per_mass
  namelist /multispec/ bc_function_type   ! 0 = constant 
                                          ! 1 = cubic, see description below

  namelist /multispec/ L_pos              ! length of part of boundary where there is positive charge flux, if cubic function is imposed
  namelist /multispec/ L_trans            ! length of transition part of boundary, where the value varies like a cubic
  namelist /multispec/ L_zero             ! length of part of boundary where there is zero charge flux, if cubic function is imposed

  namelist /multispec/ theta_pot          ! for implicit algorithm_type=3, controls
                                               ! temporal discretization for potential term
  namelist /multispec/ num_pot_iters
  namelist /multispec/ dpdt_factor
  namelist /multispec/ E_ext_type         ! if 1, sets an external E field to E_ext_value
  namelist /multispec/ E_ext_value        ! spacedim-vector specifying external E field

  namelist /multispec/ electroneutral               ! use electroneutral diffusion fluxes
  namelist /multispec/ induced_charge_eo            ! are we simulating ICEO?  
  namelist /multispec/ relxn_param_charge           ! Used to prevent slow buildup of charge for electroneutral, keep at 1.0
  namelist /multispec/ zero_eps_on_wall_type        ! set eps=0 on certain Dirichlet walls
                                                    ! if we want homogeneous Neumann bc's on 
                                                    ! phi for part of a Dirichlet wall
  namelist /multispec/ zero_charge_on_wall_type     ! set sigma=0 on certain Neumann walls
                                                    ! if we want homogeneous Neumann bc's on 
                                                    ! phi for part of a Dirichlet wall
  namelist /multispec/ zero_eps_on_wall_left_end    ! eg if set to 0.25, then eps will be set to 0 on the wall from 0*Lx --> 0.25*Lx
  namelist /multispec/ zero_eps_on_wall_right_start ! eg if set to 0.75, then eps will be set to 0 on the wall from 0.75*Lx --> 1.*Lx

contains

  ! read in fortran namelist into multispec_params_module
  subroutine read_multispec_namelist(inputs_file,length) bind(C, name="read_multispec_namelist")

    integer               , value         :: length
    character(kind=c_char), intent(in   ) :: inputs_file(length)

    ! default values
    fraction_tolerance = 1.d-14  ! must be larger than machine eps of else the W=(1,0) case fails)
    start_time         = 0.0d0 
    inverse_type       = 1
    correct_flux       = 1
    print_error_norms  = 1
    is_ideal_mixture   = 1
    is_nonisothermal   = 0
    use_lapack         = 0
    use_multiphase     = 0
    kc_tension         = 0.d0
    alpha_gex          = 0.d0
    n_gex              = 1
    chi_iterations     = 10
    temp_type          = 0
    c_init_1(:)        = 1.0d0
    c_init_2(:)        = 1.0d0
    Dbar(:)            = 1.0d0
    Dtherm(:)          = 0.0d0
    H_offdiag(:)       = 0.0d0
    H_diag(:)          = 0.0d0
    midpoint_stoch_mass_flux_type = 1
    avg_type           = 1
    mixture_type       = 0

    ! charged
    use_charged_fluid  = 0
    print_debye_len    = 0
    dielectric_const   = 1.d0
    dielectric_type    = 0      ! 0 = assumes constant epsilon
                                ! 1 = (1+c1)*dielectric_const
                                ! see fluid_charge.f90:compute_permittivity()
    charge_per_mass(:)     = 0.d0
    bc_function_type       = 0 
    L_pos                  = 0.d0
    L_trans                = 0.d0
    L_zero                 = 0.d0
    theta_pot              = 0.5d0
    num_pot_iters          = 2
    dpdt_factor            = 0.d0
    E_ext_type             = 0
    E_ext_value(:)         = 0.d0

    electroneutral = 0
    induced_charge_eo = 0
    relxn_param_charge = 1.0d0
    zero_eps_on_wall_type = 0
    zero_charge_on_wall_type = 0
    zero_eps_on_wall_left_end = 0.25d0
    zero_eps_on_wall_right_start = 0.75d0
    
    ! read in multispec namelist
    open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
    read(unit=100, nml=multispec)
    close(unit=100)

  end subroutine read_multispec_namelist

  ! copy contents of multispec_params_module to C++ multispec namespace
  subroutine initialize_multispec_namespace( inverse_type_in, temp_type_in, chi_iterations_in, start_time_in, &
                                             Dbar_in, Dtherm_in, H_offdiag_in, H_diag_in, &
                                             fraction_tolerance_in, correct_flux_in, print_error_norms_in, &
                                             is_nonisothermal_in, is_ideal_mixture_in, &
                                             use_lapack_in, use_multiphase_in, &
                                             kc_tension_in, alpha_gex_in, n_gex_in, &
                                             c_init_1_in, c_init_2_in, &
                                             midpoint_stoch_mass_flux_type_in, &
                                             avg_type_in, mixture_type_in, &
                                             use_charged_fluid_in, print_debye_len_in, dielectric_const_in, &
                                             dielectric_type_in, charge_per_mass_in, &
                                             theta_pot_in, num_pot_iters_in, dpdt_factor_in, &
                                             relxn_param_charge_in, E_ext_type_in, E_ext_value_in, &
                                             electroneutral_in, induced_charge_eo_in, &
                                             zero_eps_on_wall_type_in, zero_charge_on_wall_type_in, &
                                             zero_eps_on_wall_left_end_in, zero_eps_on_wall_right_start_in, &
                                             bc_function_type_in, L_pos_in, L_trans_in, L_zero_in) &
                                             bind(C, name="initialize_multispec_namespace")

    integer,            intent(inout) :: inverse_type_in
    integer,            intent(inout) :: temp_type_in
    integer,            intent(inout) :: chi_iterations_in
    double precision,   intent(inout) :: start_time_in 
    double precision,   intent(inout) :: Dbar_in(MAX_ELEMENT)
    double precision,   intent(inout) :: Dtherm_in(MAX_SPECIES)
    double precision,   intent(inout) :: H_offdiag_in(MAX_ELEMENT)
    double precision,   intent(inout) :: H_diag_in(MAX_SPECIES)
    double precision,   intent(inout) :: fraction_tolerance_in
    integer,            intent(inout) :: correct_flux_in
    integer,            intent(inout) :: print_error_norms_in
    integer,            intent(inout) :: is_nonisothermal_in
    integer,            intent(inout) :: is_ideal_mixture_in
    integer,            intent(inout) :: use_lapack_in
    integer,            intent(inout) :: use_multiphase_in
    double precision,   intent(inout) :: kc_tension_in
    double precision,   intent(inout) :: alpha_gex_in
    integer,            intent(inout) :: n_gex_in
    double precision,   intent(inout) :: c_init_1_in(MAX_SPECIES)
    double precision,   intent(inout) :: c_init_2_in(MAX_SPECIES)

    integer,            intent(inout) :: midpoint_stoch_mass_flux_type_in
    integer,            intent(inout) :: avg_type_in
    integer,            intent(inout) :: mixture_type_in
    
    integer,            intent(inout) :: use_charged_fluid_in
    integer,            intent(inout) :: print_debye_len_in
    double precision,   intent(inout) :: dielectric_const_in
    integer,            intent(inout) :: dielectric_type_in
    double precision,   intent(inout) :: charge_per_mass_in(MAX_SPECIES)
    double precision,   intent(inout) :: theta_pot_in
    integer,            intent(inout) :: num_pot_iters_in
    double precision,   intent(inout) :: dpdt_factor_in
    double precision,   intent(inout) :: relxn_param_charge_in
    integer,            intent(inout) :: E_ext_type_in
    double precision,   intent(inout) :: E_ext_value_in(1:3)
    integer,            intent(inout) :: electroneutral_in
    integer,            intent(inout) :: induced_charge_eo_in
    integer,            intent(inout) :: zero_eps_on_wall_type_in
    integer,            intent(inout) :: zero_charge_on_wall_type_in
    double precision,   intent(inout) :: zero_eps_on_wall_left_end_in
    double precision,   intent(inout) :: zero_eps_on_wall_right_start_in
    integer,            intent(inout) :: bc_function_type_in
    double precision,   intent(inout) :: L_pos_in
    double precision,   intent(inout) :: L_trans_in
    double precision,   intent(inout) :: L_zero_in


    inverse_type_in = inverse_type
    temp_type_in = temp_type
    chi_iterations_in = chi_iterations
    start_time_in = start_time
    Dbar_in = Dbar
    Dtherm_in = Dtherm
    H_offdiag_in = H_offdiag
    H_diag_in = H_diag
    fraction_tolerance_in = fraction_tolerance
    correct_flux_in = correct_flux
    print_error_norms_in = print_error_norms
    is_nonisothermal_in = is_nonisothermal
    is_ideal_mixture_in = is_ideal_mixture
    use_lapack_in = use_lapack
    use_multiphase_in = use_multiphase
    kc_tension_in = kc_tension
    alpha_gex_in = alpha_gex
    n_gex_in = n_gex
    c_init_1_in = c_init_1
    c_init_2_in = c_init_2
    midpoint_stoch_mass_flux_type_in = midpoint_stoch_mass_flux_type
    avg_type_in = avg_type
    mixture_type_in = mixture_type

    use_charged_fluid_in = use_charged_fluid
    print_debye_len_in = print_debye_len
    dielectric_const_in = dielectric_const
    dielectric_type_in = dielectric_type
    charge_per_mass_in = charge_per_mass
    theta_pot_in = theta_pot
    num_pot_iters_in = num_pot_iters
    dpdt_factor_in = dpdt_factor
    relxn_param_charge_in = relxn_param_charge
    E_ext_type_in = E_ext_type
    E_ext_value_in = E_ext_value
    electroneutral_in = electroneutral
    induced_charge_eo_in = induced_charge_eo
    zero_eps_on_wall_type_in = zero_eps_on_wall_type
    zero_charge_on_wall_type_in = zero_charge_on_wall_type
    zero_eps_on_wall_left_end_in = zero_eps_on_wall_left_end
    zero_eps_on_wall_right_start_in = zero_eps_on_wall_right_start
    bc_function_type_in = bc_function_type
    L_pos_in = L_pos
    L_trans_in = L_trans
    L_zero_in = L_zero
    
  end subroutine initialize_multispec_namespace

end module multispec_namelist_module
