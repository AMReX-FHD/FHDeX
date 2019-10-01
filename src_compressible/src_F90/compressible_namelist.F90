module compressible_namelist_module

  use iso_c_binding, only: c_char
  use amrex_string_module, only: amrex_string_c_to_f, amrex_string_f_to_c

  implicit none

  integer, parameter :: MAX_SPECIES = 10
  integer, parameter :: LOHI = 2

  double precision,   save :: bc_Yk(AMREX_SPACEDIM,LOHI,MAX_SPECIES)
  double precision,   save :: bc_Xk(AMREX_SPACEDIM,LOHI,MAX_SPECIES)
  integer         ,   save :: plot_means
  integer         ,   save :: plot_vars

  ! Boundary conditions
  namelist /compressible/ bc_Yk       ! mass fraction wall boundary value
  namelist /compressible/ bc_Xk       ! mole fraction wall boundary value
  namelist /compressible/ plot_means  ! write out means to plotfile
  namelist /compressible/ plot_vars   ! write out variances to plotfile

contains

  ! read in fortran namelist into compressible_params_module
  subroutine read_compressible_namelist(inputs_file,length) bind(C, name="read_compressible_namelist")

    integer               , value         :: length
    character(kind=c_char), intent(in   ) :: inputs_file(length)

    ! default values
    bc_Yk(:,:,:) = 0.d0
    bc_Xk(:,:,:) = 0.d0

    plot_means = 0
    plot_vars = 0

    ! read in compressible namelist
    open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
    read(unit=100, nml=compressible)
    close(unit=100)

  end subroutine read_compressible_namelist

  ! copy contents of compressible_params_module to C++ compressible namespace
  subroutine initialize_compressible_namespace(bc_Yk_in, bc_Xk_in, plot_means_in, plot_vars_in) &
                                                bind(C, name="initialize_compressible_namespace")

    double precision, intent(inout) :: bc_Yk_in(AMREX_SPACEDIM,LOHI,MAX_SPECIES)
    double precision, intent(inout) :: bc_Xk_in(AMREX_SPACEDIM,LOHI,MAX_SPECIES)
    integer         , intent(inout) :: plot_means_in
    integer         , intent(inout) :: plot_vars_in
    
    bc_Yk_in = bc_Yk
    bc_Xk_in = bc_Xk

    plot_means_in = plot_means
    plot_vars_in = plot_vars

  end subroutine initialize_compressible_namespace

end module compressible_namelist_module
