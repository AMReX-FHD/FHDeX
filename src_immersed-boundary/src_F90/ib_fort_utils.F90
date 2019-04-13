module ib_fort_utils
    use amrex_fort_module, only: amrex_real, amrex_particle_real
    use iso_c_binding,     only: c_int

    implicit none

    public particle_info_t, particle_t

    type, bind(C) :: particle_info_t
        real(amrex_real)   :: pos(3)    !< Position
        real(amrex_real)   :: vel(3)    !< Velocity
        integer(c_int)     :: ind(3)    !< Index in grid
        real(amrex_real)   :: radius    !< Particle radius
        integer(c_int)     :: id        !< Unique index
        integer(c_int)     :: cpu       !< CPU at time of initialization
    end type particle_info_t

    ! Fortran analoque to the particle data structure. This _must_ have the same order as:
    ! ( pos(1), pos(2), pos(3), {IBP_realData}, id, cpu, {IBP_intData} )
    ! ... as defined in the IBParticleContainer.
    type, bind(C) :: particle_t
        real(amrex_real)   :: pos(3)    !< Position
        real(amrex_real)   :: radius    !< Particle radius
        real(amrex_real)   :: volume
        real(amrex_real)   :: mass
        real(amrex_real)   :: density
        real(amrex_real)   :: oneOverI
        real(amrex_real)   :: velx
        real(amrex_real)   :: vely
        real(amrex_real)   :: velz
        real(amrex_real)   :: omegax
        real(amrex_real)   :: omegay
        real(amrex_real)   :: omegaz
        real(amrex_real)   :: dragx
        real(amrex_real)   :: dragy
        real(amrex_real)   :: dragz
        integer(c_int)     :: id        !< Unique index
        integer(c_int)     :: cpu       !< CPU at time of initialization
        integer(c_int)     :: phase
        integer(c_int)     :: state
    end type particle_t

contains

    subroutine test_interface(info, np) &
            bind(C, name="test_interface")

        integer,               intent(in)         :: np
        type(particle_info_t), intent(in), target :: info(np)

        integer                        :: i
        type(particle_info_t), pointer :: p

        do i = 1, np
            p => info(i)

            write(*,*) "Particle Info recieved:"
            write(*,*) "   pos = ", p%pos
            write(*,*) "   vel = ", p%vel
            write(*,*) "   ind = ", p%ind
            write(*,*) "   rad = ", p%radius, ", id = ", p%id, ", cpu = ", p%cpu
        end do

    end subroutine test_interface


    pure subroutine tag_interface_ib(iface, iflo,  ifhi,  &
                                     phi,   philo, phihi, &
                                     tag,   taglo, taghi )&
                    bind(C, name="tag_interface_ib")

        integer(c_int), dimension(3), intent(in   ) :: iflo, ifhi, philo, phihi, taglo, taghi

        integer(c_int),   intent(  out) :: iface(iflo(1):ifhi(1), iflo(2):ifhi(2), iflo(3):ifhi(3))
        real(amrex_real), intent(in   ) :: phi(philo(1):phihi(1), philo(2):phihi(2), philo(3):phihi(3))
        integer(c_int),   intent(in   ) :: tag(taglo(1):taghi(1), taglo(2):taghi(2), taglo(3):taghi(3), 3)

        integer :: i, j, k

        do k = iflo(3), ifhi(3)
            do j = iflo(2), ifhi(2)
                do i = iflo(1), ifhi(1)
                   iface(i, j, k) = cell_contains_interface( &
                        i, j, k,           &
                        phi, philo, phihi, &
                        tag, taglo, taghi  &
                        )
                end do
             end do
        end do

    contains

        pure function cell_contains_interface(i,   j,     k,     &
                                              phi, philo, phihi, &
                                              tag, taglo, taghi  )

            ! ** output type
            integer :: cell_contains_interface

            ! ** input types
            integer(c_int), dimension(3), intent(in) :: philo, phihi, taglo, taghi

            real(amrex_real), intent(in) :: phi(philo(1):phihi(1), philo(2):phihi(2), philo(3):phihi(3))
            integer(c_int),   intent(in) :: tag(taglo(1):taghi(1), taglo(2):taghi(2), taglo(3):taghi(3), 3)
            integer,          intent(in) :: i, j, k


            ! ** declare local variables
            ! ii, jj, kk : loop variables itterating over neighbour stencil
            ! klo ... ihi: boundaries of stencil which will be checked for valid cells
            !              a cell is valid if phi <= 0
            integer :: ii, jj, kk, klo, khi, jlo, jhi, ilo, ihi
            logical :: negative, positive


            !------------------------------------------------------------------!
            ! build neighbour stencil to cover all 8 corners of the cell       !
            ! note: stencil could be out-of-bounds => bounds-checking          !
            !------------------------------------------------------------------!

            klo = k
            khi = k + 1
            if ( khi .gt. phihi(3) ) then
                khi = phihi(3)
            end if

            jlo = j
            jhi = j + 1
            if ( jhi .gt. phihi(2) ) then
                jhi = phihi(2)
            end if

            ilo = i
            ihi = i + 1
            if ( ihi .gt. phihi(1) ) then
                ihi = phihi(1)
            end if


            !-----------------------------------------------------------------------------------------------!
            ! check members of neighbour stencil:                                                           !
            !       cell is "valid" whenever at least one cell in the neighbour stencil has a level-set phi !
            !       less than, or equal to, 0                                                               !
            !-----------------------------------------------------------------------------------------------!

            positive = .false.
            negative = .false.

            ! note that we're not escaping the do-loop because we need to check _all_ corners
            do kk = klo, khi
                do jj = jlo, jhi
                    do ii = ilo, ihi
                        if ( ( phi(ii, jj, kk) .lt. 0 ) .and. tag(ii, jj, kk, 1) .gt. 0 ) then
                            negative = .true.
                        end if

                        if ( ( phi(ii, jj, kk) .ge. 0 ) .and. tag(ii, jj, kk, 1) .ge. 0 ) then
                            positive = .true.
                        end if
                    end do
                end do
            end do


            if (negative .and. positive) then
                cell_contains_interface = 1
            else if ( (.not. negative) .and. positive ) then
                cell_contains_interface = 2
            else
                cell_contains_interface = 0
            end if


        end function cell_contains_interface
    end subroutine tag_interface_ib


    subroutine fill_levelset_ib(lo,        hi,             &
                                part_info, np,             &
                                phi,       philo, phihi,   &
                                tag,       taglo, taghi,   &
                                vel,       vello, velhi,   &
                                dx                       ) &
               bind(C, name="fill_levelset_ib")

        integer(c_int),   dimension(3), intent(in   ) :: lo, hi, philo, phihi, taglo, taghi, vello, velhi
        integer(c_int),                 intent(in   ) :: np
        type(particle_info_t),          intent(in   ) :: part_info(np)
        real(amrex_real), dimension(3), intent(in   ) :: dx

        real(amrex_real), intent(  out) :: phi(philo(1):phihi(1), philo(2):phihi(2), philo(3):phihi(3))
        integer(c_int),   intent(  out) :: tag(taglo(1):taghi(1), taglo(2):taghi(2), taglo(3):taghi(3), 3)
        real(amrex_real), intent(  out) :: vel(vello(1):velhi(1), vello(2):velhi(2), vello(3):velhi(3), 3)


        integer                        :: i, j, k
        real(amrex_real), dimension(3) :: pos, orig, vel_min
        real(amrex_real)               :: signed_dist
        integer                        :: id_loc, cpu_loc, ind_loc

        if ( np .le. 0 ) return

        ! Level-set is the nodal MultiFab associated with the cell-centered particle MultiFab
        ! => the effective origin of the level-set MultiFab is actually -0.5*dx
        orig(:) = (/ 0., 0., 0. /) ! 0.5*dx(:)

        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    pos = (/ i*dx(1), j*dx(2), k*dx(3) /) - orig
                    call union_particle( signed_dist, id_loc, cpu_loc, ind_loc, vel_min, &
                                         pos, part_info, np                              )

                    phi(i, j, k)    = signed_dist
                    tag(i, j, k, :) = (/ id_loc, cpu_loc, ind_loc /)
                    vel(i, j, k, :) = vel_min(:)
                end do
            end do
        end do

    end subroutine fill_levelset_ib


    subroutine union_particle(min_dist, id_min, cpu_min, ind_min, vel_min, &
                              x, part_info, np                             )

        real(amrex_real),      intent(  out)         :: min_dist, vel_min(3)
        integer(c_int),        intent(  out)         :: id_min, cpu_min, ind_min
        integer(c_int),        intent(in   )         :: np
        type(particle_info_t), intent(in   ), target :: part_info(np)
        real(amrex_real),      intent(in   )         :: x(3)

        integer                        :: i
        real(amrex_real)               :: dist
        type(particle_info_t), pointer :: p

        min_dist = -huge(min_dist)

        do i = 1, np
            p    => part_info(i)
            dist =  dist_sphere(x, p%pos, p%radius)

            if ( dist > min_dist ) then
                min_dist   = dist
                id_min     = p%id
                cpu_min    = p%cpu
                ind_min    = i
                vel_min(:) = p%vel(:)
            end if
        end do


    end subroutine union_particle



    pure function dist_sphere(x, pos, r)

      ! ** output type
      real(amrex_real) :: dist_sphere

      ! ** input types
      real(amrex_real), dimension(3), intent(in   ) :: x, pos
      real(amrex_real),               intent(in   ) :: r

      dist_sphere = r - sqrt( sum( (pos(:) - x(:))**2 ) )

    end function dist_sphere



    pure function interp_lagrange(r, h)

      ! ** output type
      real(amrex_real) :: interp_lagrange

      ! ** input types
      real(amrex_real), intent(in   ) :: r, h


      ! ** internal parameters
      real(amrex_real)            :: r_scale
      real(amrex_real), parameter :: inv3 = 1./3.
      real(amrex_real), parameter :: inv6 = 1./6.

      r_scale = abs(r) / h

      if ( r_scale .le. 0.5 ) then
         interp_lagrange = inv3*( 1 + sqrt( 1-3*r_scale**2 ) )
      else if (r_scale .le. 1.5) then
         interp_lagrange = inv6*( 5 - 3*r_scale - sqrt( 1 - 3*(1-r_scale)**2 ) )
      else
         interp_lagrange = 0
      end if


    end function interp_lagrange



    pure function effective_tag(i,   j,     k,    dir, &
                                tag, taglo, taghi      )
      ! ** output type
      integer :: effective_tag

      ! ** input types
      integer,               intent(in) :: i, j, k, dir
      integer, dimension(3), intent(in) :: taglo, taghi
      integer,               intent(in) :: tag(taglo(1):taghi(1), taglo(2):taghi(2), taglo(3):taghi(3))


      effective_tag = 0


      ! Cell already tagged => no need to check neighbours
      if ( tag(i, j, k) .eq. 1 ) then
         effective_tag = 1
         return
      end if

      ! If dir .eq. 0 : check cells in all directions
      if ( dir .eq. 0) then
         if ( ( tag(i - 1, j, k) .ge. 1 ) .or. ( tag(i + 1, j, k) .ge. 1 ) .or. &
              ( tag(i, j - 1, k) .ge. 1 ) .or. ( tag(i, j + 1, k) .ge. 1 ) .or. &
              ( tag(i, j, k + 1) .ge. 1 ) .or. ( tag(i, j, k + 1) .ge. 1 )) then
            effective_tag = 1
            return
         end if
      end if

      ! Check neighboring cells in the specific directions
      if ( dir .eq. 1) then
         if ( ( tag(i - 1, j, k) .ge. 1 ) .or. ( tag(i + 1, j, k) .ge. 1 ) ) then
            effective_tag = 1
            return
         end if
      else if (dir .eq. 2) then
         if ( ( tag(i, j - 1, k) .ge. 1 ) .or. ( tag(i, j + 1, k) .ge. 1 ) ) then
            effective_tag = 1
            return
         end if
      else if (dir .eq. 3) then
         if ( ( tag(i, j, k - 1) .ge. 1 ) .or. ( tag(i, j, k + 1) .ge. 1 ) ) then
            effective_tag = 1
            return
         end if
      end if

      if ( tag(i, j, k) .eq. 2 ) then
         effective_tag = 2
      end if

    end function effective_tag



    subroutine interpolate_ib_staggered(lo,  hi,           &
                                        u_d, udlo, udhi,   &
                                        v_d, vdlo, vdhi,   &
                                        w_d, wdlo, wdhi,   &
                                        u_s, uslo, ushi,   &
                                        v_s, vslo, vshi,   &
                                        w_s, wslo, wshi,   &
                                        et,  etlo, ethi,   &
                                        phi, philo, phihi, &
                                        tag, taglo, taghi, &
                                        vel, vello, velhi )&
               bind(C, name="interpolate_ib_staggered")

        integer(c_int),   dimension(3), intent(in   ) :: lo, hi, &
            udlo, udhi, vdlo, vdhi, wdlo, wdhi, uslo, ushi, vslo, vshi, wslo, wshi, etlo, ethi
        real(amrex_real),               intent(  out) :: u_d(udlo(1):udhi(1), udlo(2):udhi(2), udlo(3):udhi(3))
        real(amrex_real),               intent(  out) :: v_d(vdlo(1):vdhi(1), vdlo(2):vdhi(2), vdlo(3):vdhi(3))
        real(amrex_real),               intent(  out) :: w_d(wdlo(1):wdhi(1), wdlo(2):wdhi(2), wdlo(3):wdhi(3))
        real(amrex_real),               intent(in   ) :: u_s(uslo(1):ushi(1), uslo(2):ushi(2), uslo(3):ushi(3))
        real(amrex_real),               intent(in   ) :: v_s(vslo(1):vshi(1), vslo(2):vshi(2), vslo(3):vshi(3))
        real(amrex_real),               intent(in   ) :: w_s(wslo(1):wshi(1), wslo(2):wshi(2), wslo(3):wshi(3))
        integer(c_int),                 intent(  out) :: et(etlo(1):ethi(1), etlo(2):ethi(2), etlo(3):ethi(3))

        integer(c_int), dimension(3), intent(in   ) :: philo, phihi, taglo, taghi, vello, velhi
        real(amrex_real), intent(in   ) :: phi(philo(1):phihi(1), philo(2):phihi(2), philo(3):phihi(3))
        integer(c_int),   intent(in   ) :: tag(taglo(1):taghi(1), taglo(2):taghi(2), taglo(3):taghi(3))
        real(amrex_real), intent(in   ) :: vel(vello(1):velhi(1), vello(2):velhi(2), vello(3):velhi(3), 3)


        integer :: i, j, k, eff_tag

        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    et(i, j, k) = 0

                    ! eff_tag = effective_tag(i, j, k, 0, tag, taglo, taghi)
                    ! Tag only interface cells
                    if (tag(i, j, k) .eq. 1) then
                        eff_tag = 1
                    else
                        eff_tag = 0
                    end if

                    ! x-components
                    !eff_tag = effective_tag(i, j, k, 1, tag, taglo, taghi)
                    if ( eff_tag .ge. 1 ) then
                        u_d(i, j, k) = vel(i, j, k, 1)
                        et(i, j, k)  = 1
                    else
                        u_d(i, j, k) = 0. ! u_s(i, j, k)
                    end if

                    ! y-components
                    !eff_tag = effective_tag(i, j, k, 2, tag, taglo, taghi)
                    if ( eff_tag .ge. 1 ) then
                        v_d(i, j, k) = vel(i, j, k, 2)
                        et(i, j, k)  = 1
                    else
                        v_d(i, j, k) = 0. ! v_s(i, j, k)
                    end if

                    ! z-components
                    !eff_tag = effective_tag(i, j, k, 3, tag, taglo, taghi)
                    if ( eff_tag .ge. 1 ) then
                        w_d(i, j, k) = vel(i, j, k, 3)
                        et(i, j, k)  = 1
                    else
                        w_d(i, j, k) = 0. ! w_s(i, j, k)
                    end if
                end do
            end do
        end do

    end subroutine interpolate_ib_staggered



    subroutine interpolate_ib_cc(lo,  hi,           &
                                 u_d, udlo, udhi,   &
                                 u_s, uslo, ushi,   &
                                 et,  etlo, ethi,   &
                                 phi, philo, phihi, &
                                 tag, taglo, taghi, &
                                 vel, vello, velhi )&
               bind(C, name="interpolate_ib_cc")

        integer(c_int),   dimension(3), intent(in   ) :: lo, hi, udlo, udhi, uslo, ushi, etlo, ethi
        real(amrex_real),               intent(  out) :: u_d(udlo(1):udhi(1), udlo(2):udhi(2), udlo(3):udhi(3), 3)
        real(amrex_real),               intent(in   ) :: u_s(uslo(1):ushi(1), uslo(2):ushi(2), uslo(3):ushi(3), 3)
        integer(c_int),                 intent(  out) :: et(etlo(1):ethi(1), etlo(2):ethi(2), etlo(3):ethi(3))

        integer(c_int), dimension(3), intent(in   ) :: philo, phihi, taglo, taghi, vello, velhi
        real(amrex_real), intent(in   ) :: phi(philo(1):phihi(1), philo(2):phihi(2), philo(3):phihi(3))
        integer(c_int),   intent(in   ) :: tag(taglo(1):taghi(1), taglo(2):taghi(2), taglo(3):taghi(3))
        real(amrex_real), intent(in   ) :: vel(vello(1):velhi(1), vello(2):velhi(2), vello(3):velhi(3), 3)


        integer :: i, j, k, eff_tag

        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    et(i, j, k) = 0

                    eff_tag = effective_tag(i, j, k, 0, tag, taglo, taghi)

                    if ( eff_tag .eq. 1 ) then
                       u_d(i, j, k, :) = vel(i, j, k, :)
                       et(i, j, k)     = 1
                    else
                       u_d(i, j, k, :) = u_s(i, j, k, :) ! Alternatively: 0
                    end if

                end do
            end do
        end do

    end subroutine interpolate_ib_cc



    subroutine fill_fgds_ib(lo,  hi,         &
                            f_u, fulo, fuhi, &
                            f_v, fvlo, fvhi, &
                            f_w, fwlo, fwhi, &
                            et,  etlo, ethi )&
               bind(C, name="fill_fgds_ib")

        integer(c_int),   dimension(3), intent(in   ) :: lo, hi, fulo, fuhi, fvlo, fvhi, fwlo, fwhi, etlo, ethi
        real(amrex_real),               intent(  out) :: f_u(fulo(1):fuhi(1), fulo(2):fuhi(2), fulo(3):fuhi(3))
        real(amrex_real),               intent(  out) :: f_v(fvlo(1):fvhi(1), fvlo(2):fvhi(2), fvlo(3):fvhi(3))
        real(amrex_real),               intent(  out) :: f_w(fwlo(1):fwhi(1), fwlo(2):fwhi(2), fwlo(3):fwhi(3))
        integer(c_int),                 intent(  out) :: et(etlo(1):ethi(1), etlo(2):ethi(2), etlo(3):ethi(3))


        integer          :: i, j, k

        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    if (  et(i, j, k)  .ge. 1 ) then
                        f_u(i, j, k) = 1.
                        f_v(i, j, k) = 1.
                        f_w(i, j, k) = 1.
                    else
                        f_u(i, j, k) = 0.
                        f_v(i, j, k) = 0.
                        f_w(i, j, k) = 0.
                    end if
                end do
            end do
        end do

    end subroutine fill_fgds_ib


    subroutine fill_force_ib_staggered(lo,  hi,         &
                                       f_u, fulo, fuhi, &
                                       f_v, fvlo, fvhi, &
                                       f_w, fwlo, fwhi, &
                                       u_g, uglo, ughi, &
                                       v_g, vglo, vghi, &
                                       w_g, wglo, wghi, &
                                       u_d, udlo, udhi, &
                                       v_d, vdlo, vdhi, &
                                       w_d, wdlo, wdhi, &
                                       et,  etlo, ethi, &
                                       dt              )&
               bind(C, name="fill_force_ib_staggered")

        integer(c_int),   dimension(3), intent(in   ) :: lo, hi, &
            fulo, fuhi, fvlo, fvhi, fwlo, fwhi, uglo, ughi, vglo, vghi, wglo, wghi, &
            udlo, udhi, vdlo, vdhi, wdlo, wdhi, etlo, ethi
        real(amrex_real),               intent(  out) :: f_u(fulo(1):fuhi(1), fulo(2):fuhi(2), fulo(3):fuhi(3))
        real(amrex_real),               intent(  out) :: f_v(fvlo(1):fvhi(1), fvlo(2):fvhi(2), fvlo(3):fvhi(3))
        real(amrex_real),               intent(  out) :: f_w(fwlo(1):fwhi(1), fwlo(2):fwhi(2), fwlo(3):fwhi(3))
        real(amrex_real),               intent(in   ) :: u_g(uglo(1):ughi(1), uglo(2):ughi(2), uglo(3):ughi(3))
        real(amrex_real),               intent(in   ) :: v_g(vglo(1):vghi(1), vglo(2):vghi(2), vglo(3):vghi(3))
        real(amrex_real),               intent(in   ) :: w_g(wglo(1):wghi(1), wglo(2):wghi(2), wglo(3):wghi(3))
        real(amrex_real),               intent(in   ) :: u_d(udlo(1):udhi(1), udlo(2):udhi(2), udlo(3):udhi(3))
        real(amrex_real),               intent(in   ) :: v_d(vdlo(1):vdhi(1), vdlo(2):vdhi(2), vdlo(3):vdhi(3))
        real(amrex_real),               intent(in   ) :: w_d(wdlo(1):wdhi(1), wdlo(2):wdhi(2), wdlo(3):wdhi(3))
        integer(c_int),                 intent(in   ) :: et(etlo(1):ethi(1), etlo(2):ethi(2), etlo(3):ethi(3))
        real(amrex_real),               intent(in   ) :: dt


        integer          :: i, j, k

        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    if ( et(i, j, k) .ge. 1 ) then
                        f_u(i, j, k) = (u_d(i, j, k) - u_g(i, j, k)) / dt
                        f_v(i, j, k) = (v_d(i, j, k) - v_g(i, j, k)) / dt
                        f_w(i, j, k) = (w_d(i, j, k) - w_g(i, j, k)) / dt
                    else
                        f_u(i, j, k) = 0.
                        f_v(i, j, k) = 0.
                        f_w(i, j, k) = 0.
                    end if
                end do
            end do
        end do

    end subroutine fill_force_ib_staggered


    subroutine fill_force_ib_cc(lo,  hi,         &
                                f_u, fulo, fuhi, &
                                u_g, uglo, ughi, &
                                u_d, udlo, udhi, &
                                et,  etlo, ethi, &
                                dt              )&
               bind(C, name="fill_force_ib_cc")

        integer(c_int), dimension(3), intent(in   ) :: lo, hi, fulo, fuhi, uglo, ughi, udlo, udhi, etlo, ethi
        real(amrex_real), intent(  out) :: f_u(fulo(1):fuhi(1), fulo(2):fuhi(2), fulo(3):fuhi(3), 3)
        real(amrex_real), intent(in   ) :: u_g(uglo(1):ughi(1), uglo(2):ughi(2), uglo(3):ughi(3), 3)
        real(amrex_real), intent(in   ) :: u_d(udlo(1):udhi(1), udlo(2):udhi(2), udlo(3):udhi(3), 3)
        integer(c_int),   intent(in   ) :: et(etlo(1):ethi(1), etlo(2):ethi(2), etlo(3):ethi(3))
        real(amrex_real), intent(in   ) :: dt


        integer          :: i, j, k

        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    if ( et(i, j, k) .eq. 1 ) then
                        f_u(i, j, k, :) = (u_d(i, j, k, :) - u_g(i, j, k, :)) / dt
                    else
                        f_u(i, j, k, :) = 0.
                    end if
                end do
            end do
        end do

    end subroutine fill_force_ib_cc


end module ib_fort_utils
