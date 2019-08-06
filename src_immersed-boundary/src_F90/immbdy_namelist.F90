module immbdy_namelist_module

    use iso_c_binding, only: c_char
    use amrex_string_module, only: amrex_string_c_to_f, amrex_string_f_to_c

    implicit none

    ! Parameter controlling the minimum number of ghost cells in the immersed boundary
    ! marker container => gris can't be smaller thant this
    integer, parameter :: IBMC_MIN_NGHOST = 4


    integer,              save :: n_immbdy
    integer, allocatable, save :: n_marker


end module immbdy_namelist_module
