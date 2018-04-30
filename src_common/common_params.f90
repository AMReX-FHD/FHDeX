module common_params_module

  implicit none

  integer, parameter :: MAX_SPACEDIM = 3
  integer, parameter :: MAX_SPECIES = 6
  integer, parameter :: LOHI = 2

  double precision,   save :: prob_lo(MAX_SPACEDIM)
  double precision,   save :: prob_hi(MAX_SPACEDIM)
  integer,            save :: n_cells(MAX_SPACEDIM)
  integer,            save :: max_grid_size(MAX_SPACEDIM)
  double precision,   save :: fixed_dt
  double precision,   save :: cfl
  integer,            save :: max_step
  integer,            save :: plot_int
  character(len=128), save :: plot_base_name
  integer,            save :: chk_int
  character(len=128), save :: chk_base_name
  integer,            save :: prob_type
  integer,            save :: restart
  integer,            save :: print_int
  integer,            save :: project_eos_int
  double precision,   save :: grav(MAX_SPACEDIM)
  integer,            save :: nspecies
  double precision,   save :: molmass(MAX_SPECIES)
  double precision,   save :: rhobar(MAX_SPECIES)
  double precision,   save :: rho0
  double precision,   save :: variance_coef_mom
  double precision,   save :: variance_coef_mass
  double precision,   save :: k_B
  double precision,   save :: Runiv
  integer,            save :: algorithm_type
  integer,            save :: barodiffusion_type
  logical,            save :: use_bl_rng
  integer,            save :: seed
  integer,            save :: seed_momentum
  integer,            save :: seed_diffusion
  integer,            save :: seed_reaction
  integer,            save :: seed_init_mass
  integer,            save :: seed_init_momentum
  double precision,   save :: visc_coef
  integer,            save :: visc_type
  integer,            save :: advection_type
  integer,            save :: filtering_width
  integer,            save :: stoch_stress_form
  double precision,   save :: u_init(2)
  double precision,   save :: perturb_width
  double precision,   save :: smoothing_width
  double precision,   save :: initial_variance_mom
  double precision,   save :: initial_variance_mass
  integer,            save :: bc_lo(MAX_SPACEDIM)
  integer,            save :: bc_hi(MAX_SPACEDIM)
  double precision,   save :: wallspeed_lo(MAX_SPACEDIM-1,MAX_SPACEDIM)
  double precision,   save :: wallspeed_hi(MAX_SPACEDIM-1,MAX_SPACEDIM)
  integer,            save :: histogram_unit
  double precision,   save :: density_weights(MAX_SPECIES)
  integer,            save :: shift_cc_to_boundary(MAX_SPACEDIM,LOHI)
  
  namelist /probin_common/ prob_lo
  namelist /probin_common/ prob_hi
  namelist /probin_common/ n_cells
  namelist /probin_common/ max_grid_size
  namelist /probin_common/ fixed_dt
  namelist /probin_common/ cfl
  namelist /probin_common/ max_step
  namelist /probin_common/ plot_int
  namelist /probin_common/ plot_base_name
  namelist /probin_common/ chk_int
  namelist /probin_common/ chk_base_name
  namelist /probin_common/ prob_type
  namelist /probin_common/ restart
  namelist /probin_common/ print_int
  namelist /probin_common/ project_eos_int
  namelist /probin_common/ grav
  namelist /probin_common/ nspecies
  namelist /probin_common/ molmass
  namelist /probin_common/ rhobar
  namelist /probin_common/ rho0
  namelist /probin_common/ variance_coef_mom
  namelist /probin_common/ variance_coef_mass
  namelist /probin_common/ k_B
  namelist /probin_common/ Runiv
  namelist /probin_common/ algorithm_type
  namelist /probin_common/ barodiffusion_type
  namelist /probin_common/ use_bl_rng
  namelist /probin_common/ seed
  namelist /probin_common/ seed_momentum
  namelist /probin_common/ seed_diffusion
  namelist /probin_common/ seed_reaction
  namelist /probin_common/ seed_init_mass
  namelist /probin_common/ seed_init_momentum
  namelist /probin_common/ visc_coef
  namelist /probin_common/ visc_type
  namelist /probin_common/ advection_type
  namelist /probin_common/ filtering_width
  namelist /probin_common/ stoch_stress_form
  namelist /probin_common/ u_init
  namelist /probin_common/ perturb_width
  namelist /probin_common/ smoothing_width
  namelist /probin_common/ initial_variance_mom
  namelist /probin_common/ initial_variance_mass
  namelist /probin_common/ bc_lo
  namelist /probin_common/ bc_hi
  namelist /probin_common/ wallspeed_lo
  namelist /probin_common/ wallspeed_hi
  namelist /probin_common/ histogram_unit
  namelist /probin_common/ density_weights
  namelist /probin_common/ shift_cc_to_boundary

contains

  subroutine read_common_params(probin_file,length) bind(C, name="read_common_params")

    use iso_c_binding, only: c_char, c_null_char
    use amrex_string_module, only: amrex_string_c_to_f

    integer               , value         :: length
    character(kind=c_char), intent(in   ) :: probin_file(length)

    integer :: i
    character(len=length) :: probin_file_f

    ! default values
    prob_lo(:) = 0.d0
    prob_hi(:) = 1.d0
    n_cells(:) = 1
    max_grid_size(:) = 1
    fixed_dt = 1.
    cfl = 0.5
    max_step = 1
    plot_int = 0
    plot_base_name = "plt"
    chk_int = 0
    chk_base_name = "chk"
    prob_type = 1
    restart = -1
    print_int = 0
    project_eos_int = -1
    grav(:) = 0.d0
    nspecies = 2
    molmass(:) = 1.d0
    rhobar(:) = 1.d0
    rho0 = 1.
    variance_coef_mom = 1.
    variance_coef_mass = 1.
    k_B = 1.
    Runiv = 8.314462175d7
    algorithm_type = 0
    barodiffusion_type = 0
    use_bl_rng = .false.
    seed = 1
    seed_momentum = 1
    seed_diffusion = 1
    seed_reaction = 1
    seed_init_mass = 1
    seed_init_momentum = 1
    visc_coef = 1.
    visc_type = 1
    advection_type = 0
    filtering_width = 0
    stoch_stress_form = 1
    u_init(:) = 0.d0
    perturb_width = 0.
    smoothing_width = 1.
    initial_variance_mom = 0.
    initial_variance_mass = 0.
    bc_lo(:) = 0
    bc_hi(:) = 0
    wallspeed_lo(:,:) = 0
    wallspeed_hi(:,:) = 0
    histogram_unit = 0
    density_weights(:) = 0.d0
    shift_cc_to_boundary(:,:) = 0

    ! convert c string to a fortran string
    do i=1,length
       if (probin_file(i) == c_null_char) exit
       probin_file_f(i:i) = probin_file(i)
    end do

    open(unit=100, file=probin_file_f, status='old', action='read')
    read(unit=100, nml=probin_common)
    close(unit=100)

  end subroutine read_common_params

end module common_params_module
