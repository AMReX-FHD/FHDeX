module compressible_namelist_module

  use iso_c_binding, only: c_char
  use amrex_string_module, only: amrex_string_c_to_f, amrex_string_f_to_c

  implicit none

  integer, parameter :: MAX_SPECIES = 10
  integer, parameter :: LOHI = 2

  integer,            save :: bc_mass_lo(AMREX_SPACEDIM)
  integer,            save :: bc_mass_hi(AMREX_SPACEDIM)
  integer,            save :: bc_therm_lo(AMREX_SPACEDIM)
  integer,            save :: bc_therm_hi(AMREX_SPACEDIM)
  integer,            save :: bc_vel_lo(AMREX_SPACEDIM)
  integer,            save :: bc_vel_hi(AMREX_SPACEDIM)
  double precision,   save :: bc_Yk(AMREX_SPACEDIM,LOHI,MAX_SPECIES)
  double precision,   save :: bc_Xk(AMREX_SPACEDIM,LOHI,MAX_SPECIES)

  ! Boundary conditions
  namelist /compressible/ bc_mass_lo  ! mass boundary conditions
  namelist /compressible/ bc_mass_hi
  namelist /compressible/ bc_therm_lo ! temperature boundary conditions
  namelist /compressible/ bc_therm_hi
  namelist /compressible/ bc_vel_lo   ! velocity field boundary conditions
  namelist /compressible/ bc_vel_hi
  namelist /compressible/ bc_Yk       ! mass fraction wall boundary value
  namelist /compressible/ bc_Xk       ! mole fraction wall boundary value

contains

  ! read in fortran namelist into compressible_params_module
  subroutine read_compressible_namelist(inputs_file,length) bind(C, name="read_compressible_namelist")

    integer               , value         :: length
    character(kind=c_char), intent(in   ) :: inputs_file(length)

    ! default values
    bc_mass_lo(:) = 0
    bc_mass_hi(:) = 0
    bc_therm_lo(:) = 0
    bc_therm_hi(:) = 0
    bc_vel_lo(:) = 0
    bc_vel_hi(:) = 0
    bc_Yk(:,:,:) = 0.d0
    bc_Xk(:,:,:) = 0.d0

    ! read in compressible namelist
    open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
    read(unit=100, nml=compressible)
    close(unit=100)

  end subroutine read_compressible_namelist

  ! copy contents of compressible_params_module to C++ compressible namespace
  subroutine initialize_compressible_namespace( bc_mass_lo_in, bc_mass_hi_in, &
                                                bc_therm_lo_in, bc_therm_hi_in, &
                                                bc_vel_lo_in, bc_vel_hi_in, &
                                                bc_Yk_in, bc_Xk_in) &
                                                bind(C, name="initialize_compressible_namespace")

    double precision,       intent(inout) :: bc_mass_lo_in(AMREX_SPACEDIM), bc_mass_hi_in(AMREX_SPACEDIM), bc_therm_lo_in(AMREX_SPACEDIM), bc_therm_hi_in(AMREX_SPACEDIM), bc_vel_lo_in(AMREX_SPACEDIM), bc_vel_hi_in(AMREX_SPACEDIM),bc_Yk_in(AMREX_SPACEDIM,LOHI,MAX_SPECIES), bc_Xk_in(AMREX_SPACEDIM,LOHI,MAX_SPECIES)

    bc_mass_lo_in = bc_mass_lo
    bc_mass_hi_in = bc_mass_hi
    bc_therm_lo_in = bc_therm_lo
    bc_therm_hi_in = bc_therm_hi
    bc_vel_lo_in = bc_vel_lo
    bc_vel_hi_in = bc_vel_hi
    bc_Yk_in = bc_Yk
    bc_Xk_in = bc_Xk

  end subroutine initialize_compressible_namespace

end module compressible_namelist_module
