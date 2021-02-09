module compressible_namelist_module

  use iso_c_binding, only: c_char
  use amrex_string_module, only: amrex_string_c_to_f, amrex_string_f_to_c
  use common_namelist_module, only: MAX_SPECIES
  
  implicit none

  integer, parameter :: LOHI = 2

  integer         ,   save :: plot_means
  integer         ,   save :: plot_vars
  integer         ,   save :: plot_covars
  integer         ,   save :: plot_cross

  ! Boundary conditions
  namelist /compressible/ plot_means  ! write out means to plotfile
  namelist /compressible/ plot_vars   ! write out variances to plotfile
  namelist /compressible/ plot_covars   ! write out covariances to plotfile
  namelist /compressible/ plot_cross   ! write out spatial correlation to textfile

contains

  ! read in fortran namelist into compressible_params_module
  subroutine read_compressible_namelist(inputs_file,length) bind(C, name="read_compressible_namelist")

    integer               , value         :: length
    character(kind=c_char), intent(in   ) :: inputs_file(length)

    ! default values
    plot_means = 0
    plot_vars = 0
    plot_covars = 0
    plot_cross = 0

    ! read in compressible namelist
    open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
    read(unit=100, nml=compressible)
    close(unit=100)

  end subroutine read_compressible_namelist

  ! copy contents of compressible_params_module to C++ compressible namespace
  subroutine initialize_compressible_namespace(plot_means_in, plot_vars_in, plot_covars_in, plot_cross_in) &
                                               bind(C, name="initialize_compressible_namespace")

    integer         , intent(inout) :: plot_means_in
    integer         , intent(inout) :: plot_vars_in
    integer         , intent(inout) :: plot_covars_in
    integer         , intent(inout) :: plot_cross_in
    
    plot_means_in = plot_means
    plot_vars_in = plot_vars
    plot_covars_in = plot_covars
    plot_cross_in = plot_cross

  end subroutine initialize_compressible_namespace

end module compressible_namelist_module
