module immbdy_namelist_module

    use iso_c_binding,       only: c_char
    use amrex_string_module, only: amrex_string_c_to_f, amrex_string_f_to_c
    use amrex_fort_module,   only: amrex_real

    implicit none

    ! Parameter controlling the minimum number of ghost cells in the immersed boundary
    ! marker container => gris can't be smaller thant this
    integer, parameter :: IBMC_MIN_NGHOST = 4


    integer,                       save :: n_immbdy
    logical,                       save :: contains_flagellum
    integer,          allocatable, save :: n_marker(:)
    real(amrex_real), allocatable, save :: offset_0(:, :)


    namelist /immbdy/ n_immbdy
    namelist /immbdy/ contains_flagellum

    namelist /ib_flagellum/ n_marker
    namelist /ib_flagellum/ offset_0

contains

    subroutine read_immbdy_namelist(inputs_file, length) &
            bind(C, name="read_immbdy_namelist")

        integer,                intent(in), value :: length
        character(kind=c_char), intent(in)        :: inputs_file(length)

        ! default values
        n_immbdy           = 0
        contains_flagellum = .false.

        ! read in immbdy namelist
        open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
        read(unit=100, nml=immbdy)
        close(unit=100)


        if (contains_flagellum) then

            allocate(n_marker(n_immbdy))
            allocate(offset_0(n_immbdy, AMREX_SPACEDIM))

            ! default values
            n_marker(:)    = 0
            offset_0(:, :) = 0

            ! read in immbdy namelist
            open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
            read(unit=100, nml=ib_flagellum)
            close(unit=100)
        end if

    end subroutine read_immbdy_namelist



end module immbdy_namelist_module
