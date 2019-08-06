module immbdy_namelist_module

    use iso_c_binding, only: c_char
    use amrex_string_module, only: amrex_string_c_to_f, amrex_string_f_to_c

    implicit none

    ! Parameter controlling the minimum number of ghost cells in the immersed boundary
    ! marker container => gris can't be smaller thant this
    integer, parameter :: IBMC_MIN_NGHOST = 4


    integer,              save :: n_immbdy
    integer, allocatable, save :: n_marker(:) ! n_immbdy large

    namelist /immbdy/     n_immbdy
    namelist /ib_markers/ n_marker


contains

    subroutine read_immbdy_namelist(inputs_file, length) &
            bind(C, name="read_immbdy_namelist")

        integer,                intent(in), value :: length
        character(kind=c_char), intent(in)        :: inputs_file(length)


        ! default values
        n_immbdy = 0


        ! read in immbdy namelist
        open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
        read(unit=100, nml=immbdy)
        close(unit=100)


        allocate(n_marker(1:n_immbdy))

        ! default values
        n_marker(1:n_immbdy) = 0

        ! read in immbdy namelist
        open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
        read(unit=100, nml=ib_markers)
        close(unit=100)


        write(*,*) n_marker

    end subroutine read_immbdy_namelist



end module immbdy_namelist_module
