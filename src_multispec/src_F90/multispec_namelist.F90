module multispec_namelist_module

  use iso_c_binding, only: c_char
  use amrex_string_module, only: amrex_string_c_to_f, amrex_string_f_to_c
  use common_namelist_module, only: MAX_SPECIES

  implicit none

  integer, parameter :: max_element=MAX_SPECIES*(MAX_SPECIES-1)/2

  integer,            save :: inverse_type
  integer,            save :: temp_type
  integer,            save :: chi_iterations
  double precision,   save :: start_time 
  double precision,   save :: Dbar(max_element)
  double precision,   save :: Dtherm(MAX_SPECIES)
  double precision,   save :: H_offdiag(max_element)
  double precision,   save :: H_diag(MAX_SPECIES)
  double precision,   save :: fraction_tolerance
  integer,            save :: correct_flux
  integer,            save :: print_error_norms
  integer,            save :: is_nonisothermal
  integer,            save :: is_ideal_mixture
  integer,            save :: use_lapack
  double precision,   save :: c_init(2,MAX_SPECIES)
  double precision,   save :: c_bc(AMREX_SPACEDIM,2,MAX_SPECIES)
  
  integer,            save :: midpoint_stoch_mass_flux_type
  integer,            save :: avg_type
  integer,            save :: mixture_type

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
  namelist /multispec/ chi_iterations     ! number of iterations used in Dbar2chi_iterative

  ! Initial and boundary conditions 
  !----------------------

  namelist /multispec/ temp_type  ! for initializing temperature
  namelist /multispec/ c_init     ! initial values for c
  namelist /multispec/ c_bc       ! c_i boundary conditions (dir,lohi,species)
  
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
    chi_iterations     = 10
    temp_type          = 0
    c_init(:,:)        = 1.0d0
    c_bc(:,:,:)        = 0.d0
    Dbar(:)            = 1.0d0
    Dtherm(:)          = 0.0d0
    H_offdiag(:)       = 0.0d0
    H_diag(:)          = 0.0d0
    midpoint_stoch_mass_flux_type = 1
    avg_type           = 1
    mixture_type       = 0

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
                                             use_lapack_in, c_init_in, c_bc_in, &
                                             midpoint_stoch_mass_flux_type_in, &
                                             avg_type_in, mixture_type_in) &
                                             bind(C, name="initialize_multispec_namespace")

    integer,                intent(inout) :: inverse_type_in
    integer,                intent(inout) :: temp_type_in
    integer,                intent(inout) :: chi_iterations_in
    double precision,       intent(inout) :: start_time_in 
    double precision,       intent(inout) :: Dbar_in(max_element)
    double precision,       intent(inout) :: Dtherm_in(MAX_SPECIES)
    double precision,       intent(inout) :: H_offdiag_in(max_element)
    double precision,       intent(inout) :: H_diag_in(MAX_SPECIES)
    double precision,       intent(inout) :: fraction_tolerance_in
    integer,                intent(inout) :: correct_flux_in
    integer,                intent(inout) :: print_error_norms_in
    integer,                intent(inout) :: is_nonisothermal_in
    integer,                intent(inout) :: is_ideal_mixture_in
    integer,                intent(inout) :: use_lapack_in
    double precision,       intent(inout) :: c_init_in(2,MAX_SPECIES)
    double precision,       intent(inout) :: c_bc_in(AMREX_SPACEDIM,2,MAX_SPECIES)

    integer,                intent(inout) :: midpoint_stoch_mass_flux_type_in
    integer,                intent(inout) :: avg_type_in
    integer,                intent(inout) :: mixture_type_in

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
    c_init_in = c_init
    c_bc_in = c_bc
    midpoint_stoch_mass_flux_type_in = midpoint_stoch_mass_flux_type
    avg_type_in = avg_type
    mixture_type_in = mixture_type

  end subroutine initialize_multispec_namespace

end module multispec_namelist_module
