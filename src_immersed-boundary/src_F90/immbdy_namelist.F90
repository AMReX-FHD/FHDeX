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
    real(amrex_real), allocatable, save :: amplitude(:)
    real(amrex_real), allocatable, save :: frequency(:)
    real(amrex_real), allocatable, save :: length(:)
    real(amrex_real), allocatable, save :: wavelength(:)



    namelist /immbdy/ n_immbdy
    namelist /immbdy/ contains_flagellum

    namelist /ib_flagellum/ n_marker
    namelist /ib_flagellum/ offset_0
    namelist /ib_flagellum/ amplitude
    namelist /ib_flagellum/ frequency
    namelist /ib_flagellum/ length
    namelist /ib_flagellum/ wavelength

contains

    subroutine read_immbdy_namelist(inputs_file, f_length) &
            bind(C, name="read_immbdy_namelist")

        integer,                intent(in), value :: f_length
        character(kind=c_char), intent(in)        :: inputs_file(f_length)

        ! default values
        n_immbdy           = 0
        contains_flagellum = .false.

        ! read in immbdy namelist
        open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
        read(unit=100, nml=immbdy)
        close(unit=100)


        if (contains_flagellum) then

            allocate(n_marker(n_immbdy))
            allocate(offset_0(AMREX_SPACEDIM, n_immbdy))
            allocate(amplitude(n_immbdy))
            allocate(frequency(n_immbdy))
            allocate(length(n_immbdy))
            allocate(wavelength(n_immbdy))

            ! default values
            n_marker(:)    = 0
            offset_0(:, :) = 0
            amplitude(:)   = 0
            frequency(:)   = 0
            length(:)      = 0
            wavelength(:)  = 0

            ! read in immbdy namelist
            open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
            read(unit=100, nml=ib_flagellum)
            close(unit=100)

        end if

    end subroutine read_immbdy_namelist


    subroutine initialize_immbdy_namespace (n_immbdy_in, contains_flagellum_in) &
            bind(C, name="initialize_immbdy_namespace")

        integer, intent(inout) :: n_immbdy_in
        integer, intent(inout) :: contains_flagellum_in

        n_immbdy_in = n_immbdy

        contains_flagellum_in = 0
        if (contains_flagellum) contains_flagellum_in = 1

    end subroutine initialize_immbdy_namespace


    subroutine initialize_ib_flagellum_namespace (n_immbdy, n_marker_in, offset_0_in,    &
            &                                     amplitude_in, frequency_in, length_in, &
            &                                     wavelength_in                        ) &
            bind(C, name="initialize_ib_flagellum_namespace")

        integer, value, intent(in) :: n_immbdy

        integer,          intent(inout) :: n_marker_in(n_immbdy)
        real(amrex_real), intent(inout) :: offset_0_in(AMREX_SPACEDIM, n_immbdy)
        real(amrex_real), intent(inout) :: amplitude_in(n_immbdy)
        real(amrex_real), intent(inout) :: frequency_in(n_immbdy)
        real(amrex_real), intent(inout) :: length_in(n_immbdy)
        real(amrex_real), intent(inout) :: wavelength_in(n_immbdy)

        n_marker_in   = n_marker;
        offset_0_in   = offset_0;
        amplitude_in  = amplitude;
        frequency_in  = frequency;
        length_in     = length;
        wavelength_in = wavelength;

    end subroutine initialize_ib_flagellum_namespace

end module immbdy_namelist_module
