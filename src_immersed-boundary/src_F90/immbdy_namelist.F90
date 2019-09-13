module immbdy_namelist_module

    use iso_c_binding,       only: c_int, c_char, c_null_char
    use amrex_string_module, only: amrex_string_c_to_f, amrex_string_f_to_c
    use amrex_fort_module,   only: amrex_real

    implicit none

    ! Parameter controlling the minimum number of ghost cells in the immersed boundary
    ! marker container => gris can't be smaller thant this
    integer, parameter :: IBMC_MIN_NGHOST = 4
    ! Parameter controlling the string buffer size.
    integer, parameter :: STRING_BUFF_LEN = 150


    integer,                       save :: n_immbdy
    logical,                       save :: contains_flagellum
    logical,                       save :: contains_fourier

    integer,          allocatable, save :: n_marker(:)
    real(amrex_real), allocatable, save :: offset_0(:, :)
    real(amrex_real), allocatable, save :: amplitude(:)
    real(amrex_real), allocatable, save :: frequency(:)
    real(amrex_real), allocatable, save :: length(:)
    real(amrex_real), allocatable, save :: wavelength(:)
    real(amrex_real), allocatable, save :: k_spring(:)
    real(amrex_real), allocatable, save :: k_driving(:)

    integer,                       save :: fourier_coef_len
    real(amrex_real), allocatable, save :: a_coef(:,:,:)
    real(amrex_real), allocatable, save :: b_coef(:,:,:)
    character(:),     allocatable, save :: chlamy_flagellum_datafile


    ! most general immersed boundary information:
    namelist /immbdy/ n_immbdy
    namelist /immbdy/ contains_flagellum
    namelist /immbdy/ contains_fourier

    ! information relevant to generic flagellum
    namelist /ib_flagellum/ n_marker
    namelist /ib_flagellum/ offset_0
    namelist /ib_flagellum/ amplitude
    namelist /ib_flagellum/ frequency
    namelist /ib_flagellum/ length
    namelist /ib_flagellum/ wavelength
    namelist /ib_flagellum/ k_spring
    namelist /ib_flagellum/ k_driving

    ! general information about chlamy (fourier) data
    namelist /ib_flagellum/ fourier_coef_len
    namelist /ib_flagellum/ chlamy_flagellum_datafile

    ! fourier series analysis based on experimental data on chlamy
    namelist /chlamy_flagellum/ a_coef
    namelist /chlamy_flagellum/ b_coef



contains


    subroutine read_immbdy_namelist(inputs_file, f_length) &
            bind(C, name="read_immbdy_namelist")

        integer,                intent(in), value :: f_length
        character(kind=c_char), intent(in)        :: inputs_file(f_length)

        ! default values
        n_immbdy           = 0
        contains_flagellum = .false.
        contains_fourier   = .false.

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
            allocate(k_spring(n_immbdy))
            allocate(k_driving(n_immbdy))

            ! default values
            n_marker(:)    = 0
            offset_0(:, :) = 0
            amplitude(:)   = 0
            frequency(:)   = 0
            length(:)      = 0
            wavelength(:)  = 0
            k_spring(:)    = 0
            k_driving(:)   = 0

            ! ensure that the filename string buffer is large enough
            allocate(character(len=STRING_BUFF_LEN)::chlamy_flagellum_datafile)

            ! default file name for any additional chlamy data
            fourier_coef_len = 0
            chlamy_flagellum_datafile(1:4) = "none"
            chlamy_flagellum_datafile(5:) = " " ! delete any junk

            ! read in immbdy namelist
            open(unit=100, file=amrex_string_c_to_f(inputs_file), status='old', action='read')
            read(unit=100, nml=ib_flagellum)
            close(unit=100)
        end if

        if (contains_flagellum .and. contains_fourier) then

            ! fourier_coef_len has already been initialized before, but (just in
            ! case the `fourier_coef_len` field is overwritten in the chlamy data
            ! file, we'll load it again here

            open(unit=100, file=chlamy_flagellum_datafile, status='old', action='read')
            read(unit=100, nml=ib_flagellum)
            close(unit=100)

            ! now we can allocate the correct array sizes ahead of loading the
            ! coefficient tables. Also note that the last marker does not store
            ! any fourier values (as it does not have a next neighbor wrt which
            ! it can bend)

            allocate(a_coef(fourier_coef_len, maxval(n_marker)-1, n_immbdy))
            allocate(b_coef(fourier_coef_len, maxval(n_marker)-1, n_immbdy))

            open(unit=100, file=chlamy_flagellum_datafile, status='old', action='read')
            read(unit=100, nml=chlamy_flagellum)
            close(unit=100)
        end if

    end subroutine read_immbdy_namelist


    subroutine chlamy_flagellum_datafile_len(len_buffer) &
            bind(C, name="chlamy_flagellum_datafile_len")

        integer(c_int), intent(out) :: len_buffer

        len_buffer = len_trim(chlamy_flagellum_datafile)

    end subroutine chlamy_flagellum_datafile_len


    subroutine flagellum_max_markers(max_n_markers) &
            bind(C, name="flagellum_max_markers")

        integer(c_int), intent(out) :: max_n_markers

        max_n_markers = maxval(n_marker)

    end subroutine flagellum_max_markers


    subroutine initialize_immbdy_namespace (n_immbdy_in, contains_flagellum_in, &
            &                               contains_fourier_in               ) &
            bind(C, name="initialize_immbdy_namespace")

        integer(c_int), intent(inout) :: n_immbdy_in
        integer(c_int), intent(inout) :: contains_flagellum_in
        integer(c_int), intent(inout) :: contains_fourier_in

        n_immbdy_in = n_immbdy

        contains_flagellum_in = 0
        if (contains_flagellum) contains_flagellum_in = 1

        contains_fourier_in = 0
        if (contains_fourier) contains_fourier_in = 1

    end subroutine initialize_immbdy_namespace


    subroutine initialize_ib_flagellum_namespace (n_immbdy, n_marker_in, offset_0_in,       &
            &                                     amplitude_in, frequency_in, length_in,    &
            &                                     wavelength_in, k_spring_in, k_driving_in, &
            &                                     fourier_coef_len_in,                      &
            &                                     chlamy_flagellum_datafile_in             )&
            bind(C, name="initialize_ib_flagellum_namespace")

        integer :: i, n

        integer(c_int), value, intent(in) :: n_immbdy

        integer(c_int),   intent(inout) :: n_marker_in(n_immbdy)
        real(amrex_real), intent(inout) :: offset_0_in(AMREX_SPACEDIM, n_immbdy)
        real(amrex_real), intent(inout) :: amplitude_in(n_immbdy)
        real(amrex_real), intent(inout) :: frequency_in(n_immbdy)
        real(amrex_real), intent(inout) :: length_in(n_immbdy)
        real(amrex_real), intent(inout) :: wavelength_in(n_immbdy)
        real(amrex_real), intent(inout) :: k_spring_in(n_immbdy)
        real(amrex_real), intent(inout) :: k_driving_in(n_immbdy)

        integer(c_int), intent(inout) :: fourier_coef_len_in
        character(kind=c_char), intent(inout) :: chlamy_flagellum_datafile_in(*)


        ! copy general immersed-boundary flagellum parameters
        !------------------------------------------------------------------------

        n_marker_in   = n_marker;
        offset_0_in   = offset_0;
        amplitude_in  = amplitude;
        frequency_in  = frequency;
        length_in     = length;
        wavelength_in = wavelength;
        k_spring_in   = k_spring;
        k_driving_in  = k_driving;


        ! extra flagellum data relating to fourier expansion:
        !------------------------------------------------------------------------

        ! number of fourier modes
        fourier_coef_len_in = fourier_coef_len

        ! name of data file containing fourier modes (note: use len_trim => the
        ! _in pointer is only allocated to len_trim+1)
        n = len_trim(chlamy_flagellum_datafile)
        do i = 1, n
            ! you need to "substring" into a fortran string to access a single character:
           chlamy_flagellum_datafile_in(i) = chlamy_flagellum_datafile(i:i);
        end do
        chlamy_flagellum_datafile_in(n+1) = c_null_char

    end subroutine initialize_ib_flagellum_namespace


    subroutine copy_ib_fourier_data(i_immbdy, i_marker,   &
            &                       a_coef_in, b_coef_in) &
            bind(C, name="copy_ib_fourier_data")

        integer(c_int), intent(in), value :: i_immbdy, i_marker

        real(amrex_real), intent(inout) :: a_coef_in(fourier_coef_len)
        real(amrex_real), intent(inout) :: b_coef_in(fourier_coef_len)

        a_coef_in(:) = a_coef(:, i_marker, i_immbdy)
        b_coef_in(:) = b_coef(:, i_marker, i_immbdy)

    end subroutine copy_ib_fourier_data

end module immbdy_namelist_module
