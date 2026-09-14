#include <AMReX_Config.H>


    module ib_fort_utils

    use amrex_fort_module,      only: amrex_real, amrex_particle_real
    use iso_c_binding,          only: c_int

    implicit none


    public particle_info_t, particle_t, marker_info_t, marker_t


    type, bind(C) :: particle_info_t
        ! Same types an order as struct IBP_info in IBParticleInfo.H
        real(amrex_real)   :: pos(3)    !< Position
        real(amrex_real)   :: vel(3)    !< Velocity
        integer(c_int)     :: ind(3)    !< Index in grid
        real(amrex_real)   :: ori(3)    !< Orientation
        real(amrex_real)   :: radius    !< Particle radius
        integer(c_int)     :: id        !< Unique index
        integer(c_int)     :: cpu       !< CPU at time of initialization
        integer(c_int)     :: real      !< Particle is real (not a neighbor)
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


    type, bind(C) :: marker_info_t
       ! Same types an order as struct IBM_info in IBParticleInfo.H
       real(amrex_real)   :: pos(3)    !< Position
       real(amrex_real)   :: vel(3)    !< Velocity
       real(amrex_real)   :: force(3)  !< Force
       integer(c_int)     :: id        !< Unique index
       integer(c_int)     :: cpu       !< CPU at time of initialization
       integer(c_int)     :: real      !< Particle is real (not a neighbor)
    end type marker_info_t


    ! Fortran analoque to the marker data structure. This _must_ have the same order as:
    ! ( pos(1), pos(2), pos(3), {IBM_realData}, id, cpu, {IBM_intData} )
    ! ... as defined in the IBMarkerContainer.
    type, bind(C) :: marker_t
        real(amrex_real)   :: pos(3)        !< Position
        real(amrex_real)   :: velx          !< Velocity x
        real(amrex_real)   :: vely          !< Velocity y
        real(amrex_real)   :: velz          !< Velocity z
        real(amrex_real)   :: forcex        !< Force x
        real(amrex_real)   :: forcey        !< Force y
        real(amrex_real)   :: forcez        !< Force z
        real(amrex_real)   :: pred_posx     !< Predictor position x
        real(amrex_real)   :: pred_posy     !< Predictor position y
        real(amrex_real)   :: pred_posz     !< Predictor position z
        real(amrex_real)   :: pred_velx     !< Predictor position x
        real(amrex_real)   :: pred_vely     !< Predictor position y
        real(amrex_real)   :: pred_velz     !< Predictor position z
        real(amrex_real)   :: pred_forcex   !< Predictor force x
        real(amrex_real)   :: pred_forcey   !< Predictor force y
        real(amrex_real)   :: pred_forcez   !< Predictor force z
        integer(c_int)     :: id            !< Unique index
        integer(c_int)     :: cpu           !< CPU at time of initialization
        integer(c_int)     :: id_0          !< ID of neighbor 0
        integer(c_int)     :: cpu_0         !< CPU (at time of init) of neighbor 0
        integer(c_int)     :: id_1          !< ID of neighbor 1
        integer(c_int)     :: cpu_1         !< CPU (at time of init) of neighbor 1
    end type marker_t

contains


    subroutine test_interface(info, np) &
            bind(C, name="test_interface")

        integer,               intent(in)         :: np
        type(particle_info_t), intent(in), target :: info(np)

        integer                        :: i
        type(particle_info_t), pointer :: p

!        do i = 1, np
!            p => info(i)
!
!            write(*,*) "Particle Info recieved:"
!            write(*,*) "   pos = ", p%pos
!            write(*,*) "   vel = ", p%vel
!            write(*,*) "   ind = ", p%ind
!            write(*,*) "   rad = ", p%radius, "id = ", p%id, "cpu = ", p%cpu
!            write(*,*) "  real = ", p%real
!        end do

    end subroutine test_interface


!    pure
 subroutine tag_interface_ib(iface, iflo,  ifhi,  &
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

                        if ( ( phi(ii, jj, kk) .ge. 0 ) .and. tag(ii, jj, kk, 1) .gt. 0 ) then
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

    subroutine tag_catalyst_interface(lo,        hi,           &
                                    part_info,   np,           &
                                    iface, iflo, ifhi,         &
                                    ctag,       ctaglo, ctaghi,&
                                    dx                       ) &
               bind(C, name="tag_catalyst_interface")

        !________________________________________________________________________
        ! ** work region
        integer(c_int),   dimension(3), intent(in   ) :: lo, hi

        ! ** IN:  particle info
        integer(c_int),                 intent(in   ) :: np
        type(particle_info_t),          intent(in   ), target :: part_info(np)
        ! type(particle_info_t),          intent(in   ), target :: part_info
        ! ** IN:  spatial discretization
        real(amrex_real), dimension(3), intent(in   ) :: dx
        integer(c_int), dimension(3), intent(in   ) :: iflo, ifhi
        ! ** IN:  location of the interface

        integer(c_int),   intent( in) :: iface(iflo(1):ifhi(1), iflo(2):ifhi(2), iflo(3):ifhi(3))

        ! ** OUT: (nodal) location of the catalyst on the colloid
        integer(c_int),   dimension(3), intent(in   ) :: ctaglo, ctaghi
        integer(c_int),   intent(  out) :: ctag(ctaglo(1):ctaghi(1), &
            &                                  ctaglo(2):ctaghi(2), &
            &                                  ctaglo(3):ctaghi(3))


        !________________________________________________________________________
        ! ** Internal variables:
        ! i, j, k => cell-centered indices
        ! pos     => (cell centered) position of the cell (i, j, k)
        ! pos1    => position of the center of the particle
        ! vect    => vector from pos1 to pos
        ! ori     => orientation of particle
        ! dot     => dot product of the vect and ori
        integer                        :: i, j, k
        real(amrex_real) :: dot
        real(amrex_real), dimension(3) ::pos, pos1, cent, vect, ori1
        integer                        :: m
        type(particle_info_t), pointer :: p
        do m = 1, np
            p => part_info(m)
           cent =p%pos
           ori1 = p%ori
        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    ! print *, " center ",cent
                    pos = (/ (i+0.5)*dx(1), (j+0.5)*dx(2), (k+0.5)*dx(3) /)
                    ! print *," ori ", ori1
                    vect=cent-pos !part_info%pos
                    dot=ori1(1)*vect(1)+ori1(2)*vect(2)+ori1(3)*vect(3)
                    ! if we are on the interface and on the " bottom half " of
                    ! the particle (with respect to the orientation  ie dot <=0)
                    ! then there is catalyst present in this cell, otherwise
                    ! there isn't
                    !if (iface(i,j,k)==1) then
                    !print *, iface(i,j,k)!" dot ", dot, "pos", pos(3), "vect", vect(3)
                    !end if

                    if ((dot<=0.) .and. (iface(i,j,k)==1)) then
                    ctag(i, j, k) = 1

                    else
                    ctag(i,j,k)=0
                    end if
                end do
            end do
        end do

        end do
    end subroutine tag_catalyst_interface


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


    subroutine fill_levelset_sphere(lo,        hi,             &
                                    part_info,                 &
                                    phi,       philo, phihi,   &
                                    dx                       ) &
               bind(C, name="fill_levelset_sphere")

        !________________________________________________________________________
        ! ** work region
        integer(c_int),   dimension(3), intent(in   ) :: lo, hi

        ! ** IN:  particle info
        type(particle_info_t),          intent(in   ) :: part_info
        ! ** IN:  spatial discretization
        real(amrex_real), dimension(3), intent(in   ) :: dx

        ! ** OUT: (nodal) level-set (signed distance from particle surface)
        integer(c_int),   dimension(3), intent(in   ) :: philo, phihi
        real(amrex_real), intent(  out) :: phi(philo(1):phihi(1), &
            &                                  philo(2):phihi(2), &
            &                                  philo(3):phihi(3))


        !________________________________________________________________________
        ! ** Internal variables:
        ! i, j, k => cell-centered indices
        ! pos     => (nodal) position of the cell (i, j, k)
        integer                        :: i, j, k
        real(amrex_real), dimension(3) :: pos


        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    pos = (/ i*dx(1), j*dx(2), k*dx(3) /)
                    phi(i, j, k) = dist_sphere(pos, part_info%pos, part_info%radius)
                end do
            end do
        end do

    end subroutine fill_levelset_sphere


    pure function dist_sphere(x, pos, r)

      ! ** output type
      real(amrex_real) :: dist_sphere

      ! ** input types
      real(amrex_real), dimension(3), intent(in   ) :: x, pos
      real(amrex_real),               intent(in   ) :: r

      dist_sphere = r - sqrt( sum( (pos(:) - x(:))**2 ) )

    end function dist_sphere



    pure function interp_lagrange(r, h)

        ! The 3-point kernel function based on the paper:
        ! >*An Adaptive Version of the Immersed Boundary Method*
        ! >Alexandre M Roma, Charles S Peskin, Marsha J Berger, *Journal of Computational Physics* **153**, 509 (1999)

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



    pure function kernel_6p(r_in)

        ! The 6-point kernel function, based on the paper:
        !
        ! >*A Gaussian-like immersed-boundary kernel with three continuous derivatives and improved translational invariance*
        ! >Yuanxun Bao, Jason Kaye, Charles S. Peskin, *Journal of Computational Physics* **316**, 139 (2016)
        ! >https://dx.doi.org/10.1016/j.jcp.2016.04.024
        !
        ! Note also: https://github.com/stochasticHydroTools/IBMethod/blob/master/IBKernels/Kernels.c because the paper above has
        ! mistakes (but the repo's OK)

        ! ** output type
        real(amrex_real) :: kernel_6p

        ! ** input types
        real(amrex_real), intent(in) :: r_in

        ! ** internal parameters
        real(amrex_real), parameter :: K = 59d0/60 - sqrt(29d0)/20
        real(amrex_real), parameter :: sgn = sign(1.d0, 3d0/2-K)

        ! ** pre-computed ratios
        real(amrex_real), parameter :: inv16 = 1./16.
        real(amrex_real), parameter :: inv8  = 1./8.
        real(amrex_real), parameter :: inv12 = 1./12.
        real(amrex_real), parameter :: inv4  = 1./4.
        real(amrex_real), parameter :: inv6  = 1./6.
        real(amrex_real), parameter :: rat58 = 5./8.

        ! ** internal variables
        real(amrex_real) :: r

        ! ** initialize r
        r = r_in


        ! ** compute kernel function
        if (r .le. -3) then
            kernel_6p = 0.
        else if (r .le. -2) then
            r = r + 3
            kernel_6p = phi1(r)
        else if (r .le. -1) then
            r = r + 2
            kernel_6p = -3*phi1(r) - inv16 + inv8*(K+r**2) + inv12*((3*K-1)*r+r**3)
        else if (r .le.  0) then
            r = r + 1
            kernel_6p = 2*phi1(r) + inv4 + inv6*((4-3*K)*r-r**3)
        else if (r .le.  1) then
            kernel_6p = 2*phi1(r) + rat58 - inv4*(K+r**2)
        else if (r .le.  2) then
            r = r - 1
            kernel_6p = -3*phi1(r) + inv4 - inv6*((4-3*K)*r-r**3)
        else if (r .le.  3) then
            r = r - 2
            kernel_6p = phi1(r) - inv16 + inv8*(K+r**2) - inv12*((3*K-1)*r+r**3)
        else
            kernel_6p = 0
        end if


        ! ** sub-functions: beta, gamma, phi1
        contains

            pure function beta(r)

                ! ** output type
                real(amrex_real) :: beta

                ! ** input types
                real(amrex_real), intent(in) :: r

                ! ** pre-computed ratios
                real(amrex_real), parameter :: a = 9d0/4
                real(amrex_real), parameter :: b = 3d0/2
                real(amrex_real), parameter :: c = 22d0/3
                real(amrex_real), parameter :: d = 7d0/3

                ! NOTE: mistake in the paper: b*(K+r**2)*r -> b*(K+r**2)
                ! beta = a - b*(K+r**2)*r + (c-7*K)*r - d*r**3
                beta = a - b*(K+r**2) + (c-7*K)*r - d*r**3


            end function beta


            pure function gamma(r)

                ! ** output type
                real(amrex_real) :: gamma

                ! ** input types
                real(amrex_real), intent(in) :: r

                ! ** pre-computed ratios
                real(amrex_real), parameter :: a = 11d0/32
                real(amrex_real), parameter :: b = 3d0/32
                real(amrex_real), parameter :: c = 1d0/72
                real(amrex_real), parameter :: d = 1d0/18

                gamma = - a*r**2 + b*(2*K+r**2)*r**2 &
                    &   + c*((3*K-1)*r+r**3)**2 + d*((4-3*K)*r-r**3)**2
            end function gamma

            pure function phi1(r)

                ! ** output type
                real(amrex_real) :: phi1

                ! ** input types
                real(amrex_real), intent(in) :: r

                ! ** parameters
                real(amrex_real), parameter :: alpha = 28

                ! ** pre-computed ratios
                real(amrex_real), parameter :: inv_alpha = 1./(2.*alpha)

                phi1 = inv_alpha * ( -beta(r) + sgn * sqrt(beta(r)**2 - 4*alpha*gamma(r)) )

            end function phi1


    end function kernel_6p

    pure function kernel_3p(r_in)

        ! The 3-point kernel function, based on Google?

        ! ** output type
        real(amrex_real) :: kernel_3p

        ! ** input types
        real(amrex_real), intent(in) :: r_in


        ! ** internal variables
        real(amrex_real) :: r
        real(amrex_real) :: r1
        real(amrex_real) :: r2

        ! ** initialize r
        r = r_in
        r1 = abs(r_in)
        r2 = r1*r1

        if(r1 .le. 0.5) then

          kernel_3p = (1 + sqrt(1-3*r2))/3

        elseif(r1 .le. 1.5) then

          kernel_3p = (5 - 3*r1 - sqrt(-3*(1-r1)*(1-r1)+1)) / 6.

        else

          kernel_3p = 0

        endif

    end function kernel_3p

    subroutine spread_kernel(lo,       hi,               &
            &                mf_x,     mfx_lo,   mfx_hi, &
            &                mf_y,     mfy_lo,   mfy_hi, &
            &                mf_z,     mfz_lo,   mfz_hi, &
            &                weight_x, wfx_lo,   wfx_hi, &
            &                weight_y, wfy_lo,   wfy_hi, &
            &                weight_z, wfz_lo,   wfz_hi, &
            &                coords_x, cx_lo,    cx_hi,  &
            &                coords_y, cy_lo,    cy_hi,  &
            &                coords_z, cz_lo,    cz_hi,  &
            &                pos, v_spread, dx, ghost, pkernel_fluid_in)


        !________________________________________________________________________
        ! ** work region
        integer(c_int), dimension(3), intent(in   ) :: lo, hi
        integer(c_int), intent(in   ) :: ghost
        integer(c_int), intent(in   ) :: pkernel_fluid_in(1)

        ! ** OUT: vector quantity `v_spread` is spread to (staggered) MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: mfx_lo, mfx_hi
        real(amrex_real), intent(inout) :: mf_x(mfx_lo(1):mfx_hi(1), &
            &                                   mfx_lo(2):mfx_hi(2), &
            &                                   mfx_lo(3):mfx_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: mfy_lo, mfy_hi
        real(amrex_real), intent(inout) :: mf_y(mfy_lo(1):mfy_hi(1), &
            &                                   mfy_lo(2):mfy_hi(2), &
            &                                   mfy_lo(3):mfy_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: mfz_lo, mfz_hi
        real(amrex_real), intent(inout) :: mf_z(mfz_lo(1):mfz_hi(1), &
            &                                   mfz_lo(2):mfz_hi(2), &
            &                                   mfz_lo(3):mfz_hi(3))

        ! ** OUT: total spread weights to staggered (face-centered) faces
        !         corresponding to the MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: wfx_lo, wfx_hi
        real(amrex_real), intent(inout) :: weight_x(wfx_lo(1):wfx_hi(1), &
            &                                       wfx_lo(2):wfx_hi(2), &
            &                                       wfx_lo(3):wfx_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: wfy_lo, wfy_hi
        real(amrex_real), intent(inout) :: weight_y(wfy_lo(1):wfy_hi(1), &
            &                                       wfy_lo(2):wfy_hi(2), &
            &                                       wfy_lo(3):wfy_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: wfz_lo, wfz_hi
        real(amrex_real), intent(inout) :: weight_z(wfz_lo(1):wfz_hi(1), &
            &                                       wfz_lo(2):wfz_hi(2), &
            &                                       wfz_lo(3):wfz_hi(3))

        ! ** IN:  coordinates of staggered (face-centered) faces corresponding to the
        !         MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: cx_lo, cx_hi
        real(amrex_real), intent(in   ) :: coords_x(cx_lo(1):cx_hi(1), &
            &                                       cx_lo(2):cx_hi(2), &
            &                                       cx_lo(3):cx_hi(3), AMREX_SPACEDIM)

        integer(c_int), dimension(3), intent(in   ) :: cy_lo, cy_hi
        real(amrex_real), intent(in   ) :: coords_y(cy_lo(1):cy_hi(1), &
            &                                       cy_lo(2):cy_hi(2), &
            &                                       cy_lo(3):cy_hi(3), AMREX_SPACEDIM)

        integer(c_int), dimension(3), intent(in   ) :: cz_lo, cz_hi
        real(amrex_real), intent(in   ) :: coords_z(cz_lo(1):cz_hi(1), &
            &                                       cz_lo(2):cz_hi(2), &
            &                                       cz_lo(3):cz_hi(3), AMREX_SPACEDIM)

        ! ** IN:  quantity to spread (v_spread), given kernel position (pos), and the
        !         fluid grid discretization (dx)
        real(amrex_real), dimension(AMREX_SPACEDIM), intent(in   ) :: pos, v_spread, dx


        !________________________________________________________________________
        ! i, j, k   => face-centered indices
        integer :: i, j, k, ilo, ihi, jlo, jhi, klo, khi, gs
        ! ll        => loop counter over AMREX_SPACEDIM
        integer :: ll
        ! invvol    => 1/dx^AMREX_SPACEDIM
        ! weight    => kernel weight function
        real(amrex_real) :: invvol, weight
        ! pos_grid  => (pos - (cell position))/dx
        ! invdx     => 1/dx
        real(amrex_real), dimension(AMREX_SPACEDIM) :: pos_grid, invdx


        !________________________________________________________________________
        ! Using function pointer to specify kernel type - some question as to
        ! optimal approach here. DRL.

        abstract interface
          function kernel_np (r_in)
             use amrex_fort_module,      only: amrex_real
             real(amrex_real) :: kernel_np
             real(amrex_real), intent (in) :: r_in
          end function kernel_np
        end interface

        procedure (kernel_np), pointer :: kernel_ptr => null()

        ! TODO: It makes sense to specify this function pointer globally (as the
        ! pkernel_fluid parameter is also a global variable). I.e. when the
        ! pkernel_fluid value is set => we save ourselves the followint IF
        ! branch point. JPB.

        if(pkernel_fluid_in(1) .eq. 3) then
          kernel_ptr => kernel_3p
          gs = 2
        else
          kernel_ptr => kernel_6p
          gs = 4
        endif

        !procedure (func), pointer :: f_ptr => null ()


        !________________________________________________________________________
        ! compute geometric quantities : 1/dx and 1/dx^AMREX_SPACEDIM
        invdx(:) = 1d0/dx(:)
        invvol = 1d0
        do ll = 1, AMREX_SPACEDIM
            invvol = invvol * invdx(ll)
        end do

        if(ghost .eq. 0) then
          ilo = max(lo(1), int(pos(1) * invdx(1) - gs))
          ihi = min(hi(1), int(pos(1) * invdx(1) + gs))
          jlo = max(lo(2), int(pos(2) * invdx(2) - gs))
          jhi = min(hi(2), int(pos(2) * invdx(2) + gs))
          klo = max(lo(3), int(pos(3) * invdx(3) - gs))
          khi = min(hi(3), int(pos(3) * invdx(3) + gs))
        else
          ilo = int(pos(1) * invdx(1) - gs)
          ihi = int(pos(1) * invdx(1) + gs)
          jlo = int(pos(2) * invdx(2) - gs)
          jhi = int(pos(2) * invdx(2) + gs)
          klo = int(pos(3) * invdx(3) - gs)
          khi = int(pos(3) * invdx(3) + gs)
        endif
        !________________________________________________________________________
        ! x-components
        ! do k = lo(3), hi(3)
        !     do j = lo(2), hi(2)
        !         do i = lo(1), hi(1) + 1
        do k = klo, khi
            do j = jlo, jhi
                do i = ilo, ihi + 1

                    pos_grid(:) = pos(:) - coords_x(i, j, k, :)
                    pos_grid(:) = pos_grid(:) * invdx(:)

                    weight = 1d0
                    do ll = 1, AMREX_SPACEDIM
                        weight = weight * kernel_ptr(pos_grid(ll));
                    end do

                    mf_x(i, j, k)     = mf_x(i, j, k) + v_spread(1) * weight * invvol

                    weight_x(i, j, k) = weight_x(i, j, k) + weight
                end do
            end do
        end do

        !________________________________________________________________________
        ! y-components
        ! do k = lo(3), hi(3)
        !     do j = lo(2), hi(2) + 1
        !         do i = lo(1), hi(1)
        do k = klo, khi
            do j = jlo, jhi + 1
                do i = ilo, ihi

                    pos_grid(:) = pos(:) - coords_y(i, j, k, :)
                    pos_grid(:) = pos_grid(:) * invdx(:)

                    weight = 1d0
                    do ll = 1, AMREX_SPACEDIM
                        weight = weight * kernel_ptr(pos_grid(ll));
                    end do

                    mf_y(i, j, k)     = mf_y(i, j, k) + v_spread(2) * weight * invvol

                    weight_y(i, j, k) = weight_y(i, j, k) + weight

                end do
            end do
        end do


        !________________________________________________________________________
        ! z-components
        ! do k = lo(3), hi(3) + 1
        !     do j = lo(2), hi(2)
        !         do i = lo(1), hi(1)
        do k = klo, khi + 1
            do j = jlo, jhi
                do i = ilo, ihi

                    pos_grid(:) = pos(:) - coords_z(i, j, k, :)
                    pos_grid(:) = pos_grid(:) * invdx(:)

                    weight = 1d0
                    do ll = 1, AMREX_SPACEDIM
                        weight = weight * kernel_ptr(pos_grid(ll))
                    end do

                    mf_z(i, j, k) = mf_z(i, j, k) + v_spread(3) * weight * invvol

                    weight_z(i, j, k) = weight_z(i, j, k) + weight

                end do
            end do
        end do

    end subroutine spread_kernel

    subroutine spread_kernel_cc(lo,       hi,               &
                &                mf,     mf_lo,   mf_hi, &
                &                weights, wf_lo,   wf_hi, &
                &                coords, c_lo,    c_hi,  &
                &                pos, v_spread, dx, ghost, pkernel_fluid_in)


        !________________________________________________________________________
        ! ** work region
        integer(c_int), dimension(3), intent(in   ) :: lo, hi
        integer(c_int), intent(in   ) :: ghost
        integer(c_int), intent(in   ) :: pkernel_fluid_in(1)

        ! ** OUT: vector quantity `v_spread` is spread to MultiFab `mf`
        integer(c_int), dimension(3), intent(in   ) :: mf_lo, mf_hi
        real(amrex_real), intent(inout) :: mf(mf_lo(1):mf_hi(1), &
            &                                   mf_lo(2):mf_hi(2), &
            &                                   mf_lo(3):mf_hi(3))

        ! ** OUT: total spread weights to cell centers
        !         corresponding to the multiFab `mf`
        integer(c_int), dimension(3), intent(in   ) :: wf_lo, wf_hi
        real(amrex_real), intent(inout) :: weights(wf_lo(1):wf_hi(1), &
            &                                       wf_lo(2):wf_hi(2), &
            &                                       wf_lo(3):wf_hi(3))

        ! ** IN:  coordinates of cell centered multifab `mf`
        integer(c_int), dimension(3), intent(in   ) :: c_lo, c_hi
        real(amrex_real), intent(in   ) :: coords(c_lo(1):c_hi(1), &
            &                                       c_lo(2):c_hi(2), &
            &                                       c_lo(3):c_hi(3), AMREX_SPACEDIM)

        ! ** IN:  quantity to spread (v_spread), given kernel position (pos), and the
        !         fluid grid discretization (dx)
        real(amrex_real), dimension(AMREX_SPACEDIM), intent(in   ) :: pos, v_spread, dx


        !________________________________________________________________________
        ! i, j, k   => cell-centered indices
        integer :: i, j, k, ilo, ihi, jlo, jhi, klo, khi, gs
        ! ll        => loop counter over AMREX_SPACEDIM
        integer :: ll
        ! invvol    => 1/dx^AMREX_SPACEDIM
        ! weight    => kernel weight function
        real(amrex_real) :: invvol, weight
        ! pos_grid  => (pos - (cell position))/dx
        ! invdx     => 1/dx
        real(amrex_real), dimension(AMREX_SPACEDIM) :: pos_grid, invdx


        !________________________________________________________________________
        ! Using function pointer to specify kernel type - some question as to
        ! optimal approach here. DRL.

        abstract interface
          function kernel_np (r_in)
             use amrex_fort_module,      only: amrex_real
             real(amrex_real) :: kernel_np
             real(amrex_real), intent (in) :: r_in
          end function kernel_np
        end interface

        procedure (kernel_np), pointer :: kernel_ptr => null()

        ! TODO: It makes sense to specify this function pointer globally (as the
        ! pkernel_fluid parameter is also a global variable). I.e. when the
        ! pkernel_fluid value is set => we save ourselves the followint IF
        ! branch point. JPB.

        if(pkernel_fluid_in(1) .eq. 3) then
          kernel_ptr => kernel_3p
          gs = 2
        else
          kernel_ptr => kernel_6p
          gs = 4
        endif

        !procedure (func), pointer :: f_ptr => null ()


        !________________________________________________________________________
        ! compute geometric quantities : 1/dx and 1/dx^AMREX_SPACEDIM
        invdx(:) = 1d0/dx(:)
        invvol = 1d0
        do ll = 1, AMREX_SPACEDIM
            invvol = invvol * invdx(ll)
        end do

        if(ghost .eq. 0) then
          ilo = max(lo(1), int(pos(1) * invdx(1) - gs))
          ihi = min(hi(1), int(pos(1) * invdx(1) + gs))
          jlo = max(lo(2), int(pos(2) * invdx(2) - gs))
          jhi = min(hi(2), int(pos(2) * invdx(2) + gs))
          klo = max(lo(3), int(pos(3) * invdx(3) - gs))
          khi = min(hi(3), int(pos(3) * invdx(3) + gs))
        else
          ilo = int(pos(1) * invdx(1) - gs)
          ihi = int(pos(1) * invdx(1) + gs)
          jlo = int(pos(2) * invdx(2) - gs)
          jhi = int(pos(2) * invdx(2) + gs)
          klo = int(pos(3) * invdx(3) - gs)
          khi = int(pos(3) * invdx(3) + gs)
        endif
        !________________________________________________________________________
        ! x-components
        ! do k = lo(3), hi(3)
        !     do j = lo(2), hi(2)
        !         do i = lo(1), hi(1) + 1
        do k = klo, khi
            do j = jlo, jhi
                do i = ilo, ihi

                    pos_grid(:) = pos(:) - coords(i, j, k, :)
                    pos_grid(:) = pos_grid(:) * invdx(:)

                    weight = 1d0
                    do ll = 1, AMREX_SPACEDIM
                        weight = weight * kernel_ptr(pos_grid(ll));
                    end do

                    mf(i, j, k)     = mf(i, j, k) + v_spread(1) * weight * invvol

                    weights(i, j, k) = weights(i, j, k) + weight
                end do
            end do
        end do



    end subroutine spread_kernel_cc



    subroutine spread_markers(lo,         hi,                  &
            &                 tile_lo,    tile_hi,             &
            &                 mf_x,       mfx_lo,   mfx_hi,    &
            &                 mf_y,       mfy_lo,   mfy_hi,    &
            &                 mf_z,       mfz_lo,   mfz_hi,    &
            &                 weight_x,   wfx_lo,   wfx_hi,    &
            &                 weight_y,   wfy_lo,   wfy_hi,    &
            &                 weight_z,   wfz_lo,   wfz_hi,    &
            &                 coords_x,   cx_lo,    cx_hi,     &
            &                 coords_y,   cy_lo,    cy_hi,     &
            &                 coords_z,   cz_lo,    cz_hi,     &
            &                 pos_marker, v_marker, n_marker,  &
            &                 dx, ghost, pkernel_fluid_in ) &
            bind(C, name="spread_markers")

        !________________________________________________________________________
        ! ** work region
        integer(c_int), dimension(3), intent(in   ) :: lo, hi, tile_lo, tile_hi

        ! ** OUT: vector quantity `v_marker` is spread to (staggered) MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: mfx_lo, mfx_hi
        real(amrex_real), intent(inout) :: mf_x(mfx_lo(1):mfx_hi(1), &
            &                                   mfx_lo(2):mfx_hi(2), &
            &                                   mfx_lo(3):mfx_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: mfy_lo, mfy_hi
        real(amrex_real), intent(inout) :: mf_y(mfy_lo(1):mfy_hi(1), &
            &                                   mfy_lo(2):mfy_hi(2), &
            &                                   mfy_lo(3):mfy_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: mfz_lo, mfz_hi
        real(amrex_real), intent(inout) :: mf_z(mfz_lo(1):mfz_hi(1), &
            &                                   mfz_lo(2):mfz_hi(2), &
            &                                   mfz_lo(3):mfz_hi(3))

        ! ** OUT: total spread weights to staggered (face-centered) faces
        !         corresponding to the MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: wfx_lo, wfx_hi
        real(amrex_real), intent(inout) :: weight_x(wfx_lo(1):wfx_hi(1), &
            &                                       wfx_lo(2):wfx_hi(2), &
            &                                       wfx_lo(3):wfx_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: wfy_lo, wfy_hi
        real(amrex_real), intent(inout) :: weight_y(wfy_lo(1):wfy_hi(1), &
            &                                       wfy_lo(2):wfy_hi(2), &
            &                                       wfy_lo(3):wfy_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: wfz_lo, wfz_hi
        real(amrex_real), intent(inout) :: weight_z(wfz_lo(1):wfz_hi(1), &
            &                                       wfz_lo(2):wfz_hi(2), &
            &                                       wfz_lo(3):wfz_hi(3))

        ! ** IN:  coordinates of staggered (face-centered) faces corresponding to the
        !         MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: cx_lo, cx_hi
        real(amrex_real), intent(in   ) :: coords_x(cx_lo(1):cx_hi(1), &
            &                                       cx_lo(2):cx_hi(2), &
            &                                       cx_lo(3):cx_hi(3), AMREX_SPACEDIM)

        integer(c_int), dimension(3), intent(in   ) :: cy_lo, cy_hi
        real(amrex_real), intent(in   ) :: coords_y(cy_lo(1):cy_hi(1), &
            &                                       cy_lo(2):cy_hi(2), &
            &                                       cy_lo(3):cy_hi(3), AMREX_SPACEDIM)

        integer(c_int), dimension(3), intent(in   ) :: cz_lo, cz_hi
        real(amrex_real), intent(in   ) :: coords_z(cz_lo(1):cz_hi(1), &
            &                                       cz_lo(2):cz_hi(2), &
            &                                       cz_lo(3):cz_hi(3), AMREX_SPACEDIM)


        ! ** IN:  `n_marker`-many marker positions (pos_marker). The `pos_marker` and
        !         `v_marker` arrays are `Vector<RealVecr>` in c-land.
        integer(c_int), intent(in   ) :: n_marker, ghost
        integer(c_int), intent(in   ) :: pkernel_fluid_in(1)
        real(amrex_real), intent(in   ) :: pos_marker(AMREX_SPACEDIM, n_marker);
        ! ** IN:  each marker spreads a vector qunantity `v_marker` to the stagged
        !         MultiFabs `mf_*`
        real(amrex_real), intent(in   ) :: v_marker(AMREX_SPACEDIM, n_marker);
        ! ** IN:  fluid grid discretization
        real(amrex_real), dimension(AMREX_SPACEDIM), intent(in   ) :: dx


        !________________________________________________________________________
        ! i        <= loop counter over markers
        integer :: i
        ! pos      <= position of current marker
        ! v_spread <= current markers vector quantity to spread
        real(amrex_real), dimension(AMREX_SPACEDIM) :: pos, v_spread

        do i =  1, n_marker

            pos = pos_marker(:, i)
            if(ghost .eq. 0) then
              ! skip marker if outside tile box (prevent double-counting)
              if (pos(1) .lt. tile_lo(1)*dx(1) ) cycle
              if (pos(1) .ge. (tile_hi(1)+1)*dx(1) ) cycle

              if (pos(2) .lt. tile_lo(2)*dx(2) ) cycle
              if (pos(2) .ge. (tile_hi(2)+1)*dx(2) ) cycle

              if (pos(3) .lt. tile_lo(3)*dx(3) ) cycle
              if (pos(3) .ge. (tile_hi(3)+1)*dx(3) ) cycle

            endif

            v_spread = v_marker(:, i)

            call spread_kernel(lo,       hi,               &
                &              mf_x,     mfx_lo,   mfx_hi, &
                &              mf_y,     mfy_lo,   mfy_hi, &
                &              mf_z,     mfz_lo,   mfz_hi, &
                &              weight_x, wfx_lo,   wfx_hi, &
                &              weight_y, wfy_lo,   wfy_hi, &
                &              weight_z, wfz_lo,   wfz_hi, &
                &              coords_x, cx_lo,    cx_hi,  &
                &              coords_y, cy_lo,    cy_hi,  &
                &              coords_z, cz_lo,    cz_hi,  &
                &              pos,  v_spread, dx, ghost,  &
                &              pkernel_fluid_in)

        end do

    end subroutine spread_markers



    subroutine interpolate_kernel(lo,       hi,               &
            &                     mf_x,     mfx_lo,   mfx_hi, &
            &                     mf_y,     mfy_lo,   mfy_hi, &
            &                     mf_z,     mfz_lo,   mfz_hi, &
            &                     weight_x, wfx_lo,   wfx_hi, &
            &                     weight_y, wfy_lo,   wfy_hi, &
            &                     weight_z, wfz_lo,   wfz_hi, &
            &                     coords_x, cx_lo,    cx_hi,  &
            &                     coords_y, cy_lo,    cy_hi,  &
            &                     coords_z, cz_lo,    cz_hi,  &
            &                     pos,      v_interp, dx, pkernel_fluid_in)


        !________________________________________________________________________
        ! ** work region
      integer(c_int), dimension(3), intent(in   ) :: lo, hi
      integer(c_int), intent(in   ) :: pkernel_fluid_in(1)

        ! ** IN:  interpolating v_interp from staggered MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: mfx_lo, mfx_hi
        real(amrex_real), intent(in   ) :: mf_x(mfx_lo(1):mfx_hi(1), &
            &                                   mfx_lo(2):mfx_hi(2), &
            &                                   mfx_lo(3):mfx_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: mfy_lo, mfy_hi
        real(amrex_real), intent(in   ) :: mf_y(mfy_lo(1):mfy_hi(1), &
            &                                   mfy_lo(2):mfy_hi(2), &
            &                                   mfy_lo(3):mfy_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: mfz_lo, mfz_hi
        real(amrex_real), intent(in   ) :: mf_z(mfz_lo(1):mfz_hi(1), &
            &                                   mfz_lo(2):mfz_hi(2), &
            &                                   mfz_lo(3):mfz_hi(3))

        ! ** IN:  total spread weights to staggered (face-centered) faces
        !         corresponding to the MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: wfx_lo, wfx_hi
        real(amrex_real), intent(in   ) :: weight_x(wfx_lo(1):wfx_hi(1), &
            &                                       wfx_lo(2):wfx_hi(2), &
            &                                       wfx_lo(3):wfx_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: wfy_lo, wfy_hi
        real(amrex_real), intent(in   ) :: weight_y(wfy_lo(1):wfy_hi(1), &
            &                                       wfy_lo(2):wfy_hi(2), &
            &                                       wfy_lo(3):wfy_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: wfz_lo, wfz_hi
        real(amrex_real), intent(in   ) :: weight_z(wfz_lo(1):wfz_hi(1), &
            &                                       wfz_lo(2):wfz_hi(2), &
            &                                       wfz_lo(3):wfz_hi(3))

        ! ** IN:  coordinates of staggered (face-centered) faces corresponding to the
        !         MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: cx_lo, cx_hi
        real(amrex_real), intent(in   ) :: coords_x(cx_lo(1):cx_hi(1), &
            &                                       cx_lo(2):cx_hi(2), &
            &                                       cx_lo(3):cx_hi(3), AMREX_SPACEDIM)

        integer(c_int), dimension(3), intent(in   ) :: cy_lo, cy_hi
        real(amrex_real), intent(in   ) :: coords_y(cy_lo(1):cy_hi(1), &
            &                                       cy_lo(2):cy_hi(2), &
            &                                       cy_lo(3):cy_hi(3), AMREX_SPACEDIM)

        integer(c_int), dimension(3), intent(in   ) :: cz_lo, cz_hi
        real(amrex_real), intent(in   ) :: coords_z(cz_lo(1):cz_hi(1), &
            &                                       cz_lo(2):cz_hi(2), &
            &                                       cz_lo(3):cz_hi(3), AMREX_SPACEDIM)

        ! ** OUT: vector to interpolate (v_interp) from staggered `mf_*`
        real(amrex_real), dimension(AMREX_SPACEDIM), intent(inout) :: v_interp
        ! ** IN:  kernel position (pos), and the fluid grid discretization (dx)
        real(amrex_real), dimension(AMREX_SPACEDIM), intent(in   ) :: pos, dx


        !________________________________________________________________________
        ! i,    j,    k    => face-centered indices
        integer :: i, j, k, ilo, ihi, jlo, jhi, klo, khi, gs
        ! ll               => loop counter over AMREX_SPACEDIM
        integer :: ll
        ! weight           => kernel weight function
        ! wfrac            => weight/total_weight (weight_{x,y,z})
        real(amrex_real) :: weight, wfrac
        ! pos_grid         => (pos - (cell position))/dx
        ! invdx            => 1/dx
        real(amrex_real), dimension(AMREX_SPACEDIM) :: pos_grid, invdx

        !________________________________________________________________________
        ! Using function pointer to specify kernel type - some question as to
        ! optimal approach here. DRL.

        abstract interface
          function kernel_np (r_in)
             use amrex_fort_module, only: amrex_real
             real(amrex_real) :: kernel_np
             real(amrex_real), intent (in) :: r_in
          end function kernel_np
        end interface

        procedure (kernel_np), pointer :: kernel_ptr => null()

        ! TODO: It makes sense to specify this function pointer globally (as the
        ! pkernel_fluid parameter is also a global variable). I.e. when the
        ! pkernel_fluid value is set => we save ourselves the followint IF
        ! branch point. JPB.

        if(pkernel_fluid_in(1) .eq. 3) then
          kernel_ptr => kernel_3p
          gs = 2
        else
          kernel_ptr => kernel_6p
          gs = 4
        endif

        !________________________________________________________________________
        ! compute geometric quantity 1/dx
        invdx(:) = 1d0/dx(:)

        ilo = max(lo(1), int(pos(1) * invdx(1) - gs))
        ihi = min(hi(1), int(pos(1) * invdx(1) + gs))
        jlo = max(lo(2), int(pos(2) * invdx(2) - gs))
        jhi = min(hi(2), int(pos(2) * invdx(2) + gs))
        klo = max(lo(3), int(pos(3) * invdx(3) - gs))
        khi = min(hi(3), int(pos(3) * invdx(3) + gs))



        !________________________________________________________________________
        ! x-components
        ! do k = lo(3), hi(3)
        !     do j = lo(2), hi(2)
        !         do i = lo(1), hi(1) + 1
        do k = klo, khi
            do j = jlo, jhi
                do i = ilo, ihi + 1

                    pos_grid(:) = pos(:) - coords_x(i, j, k, :)
                    pos_grid(:) = pos_grid(:) * invdx(:)

                    weight = 1d0
                    do ll = 1, AMREX_SPACEDIM
                        weight = weight * kernel_ptr(pos_grid(ll));
                    end do

                    if (weight_x(i, j, k) .gt. 0) then
                        wfrac = weight/weight_x(i, j, k)
                    else
                        wfrac = 1d0
                    end if

                    v_interp(1) = v_interp(1) + mf_x(i, j, k) * wfrac * weight

                end do
            end do
        end do


        !________________________________________________________________________
        ! y-components
        ! do k = lo(3), hi(3)
        !     do j = lo(2), hi(2) + 1
        !         do i = lo(1), hi(1)
        do k = klo, khi
            do j = jlo, jhi + 1
                do i = ilo, ihi

                    pos_grid(:) = pos(:) - coords_y(i, j, k, :)
                    pos_grid(:) = pos_grid(:) * invdx(:)

                    weight = 1d0
                    do ll = 1, AMREX_SPACEDIM
                        weight = weight * kernel_ptr(pos_grid(ll));
                    end do

                    if (weight_y(i, j, k) .gt. 0) then
                        wfrac = weight/weight_y(i, j, k)
                    else
                        wfrac = 1d0
                    end if

                    v_interp(2) = v_interp(2) + mf_y(i, j, k) * wfrac * weight

                end do
            end do
        end do


        !________________________________________________________________________
        ! z-components
        ! do k = lo(3), hi(3) + 1
        !     do j = lo(2), hi(2)
        !         do i = lo(1), hi(1)
        do k = klo, khi + 1
            do j = jlo, jhi
                do i = ilo, ihi

                    pos_grid(:) = pos(:) - coords_z(i, j, k, :)
                    pos_grid(:) = pos_grid(:) * invdx(:)

                    weight = 1d0
                    do ll = 1, AMREX_SPACEDIM
                        weight = weight * kernel_ptr(pos_grid(ll));
                    end do

                    if (weight_z(i, j, k) .gt. 0) then
                        wfrac = weight/weight_z(i, j, k)
                    else
                        wfrac = 1d0
                    end if

                    v_interp(3) = v_interp(3) + mf_z(i, j, k) * wfrac * weight

                end do
            end do
        end do

    end subroutine interpolate_kernel



    subroutine interpolate_markers(lo,         hi,                  &
            &                      tile_lo,    tile_hi,             &
            &                      mf_x,       mfx_lo,   mfx_hi,    &
            &                      mf_y,       mfy_lo,   mfy_hi,    &
            &                      mf_z,       mfz_lo,   mfz_hi,    &
            &                      weight_x,   wfx_lo,   wfx_hi,    &
            &                      weight_y,   wfy_lo,   wfy_hi,    &
            &                      weight_z,   wfz_lo,   wfz_hi,    &
            &                      coords_x,   cx_lo,    cx_hi,     &
            &                      coords_y,   cy_lo,    cy_hi,     &
            &                      coords_z,   cz_lo,    cz_hi,     &
            &                      pos_marker, v_marker, n_marker,  &
            &                      dx, pkernel_fluid_in           ) &
            bind(C, name="interpolate_markers")

        !________________________________________________________________________
        ! ** work region
      integer(c_int), dimension(3), intent(in   ) :: lo, hi, tile_lo, tile_hi
      integer(c_int), intent(in   ) :: pkernel_fluid_in(1)

        ! ** IN:  vector quantity `v_marker` is interpolated from (staggered) MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: mfx_lo, mfx_hi
        real(amrex_real), intent(in   ) :: mf_x(mfx_lo(1):mfx_hi(1), &
            &                                   mfx_lo(2):mfx_hi(2), &
            &                                   mfx_lo(3):mfx_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: mfy_lo, mfy_hi
        real(amrex_real), intent(in   ) :: mf_y(mfy_lo(1):mfy_hi(1), &
            &                                   mfy_lo(2):mfy_hi(2), &
            &                                   mfy_lo(3):mfy_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: mfz_lo, mfz_hi
        real(amrex_real), intent(in   ) :: mf_z(mfz_lo(1):mfz_hi(1), &
            &                                   mfz_lo(2):mfz_hi(2), &
            &                                   mfz_lo(3):mfz_hi(3))

        ! ** IN:  total spread weights to staggered (face-centered) faces
        !         corresponding to the MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: wfx_lo, wfx_hi
        real(amrex_real), intent(in   ) :: weight_x(wfx_lo(1):wfx_hi(1), &
            &                                       wfx_lo(2):wfx_hi(2), &
            &                                       wfx_lo(3):wfx_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: wfy_lo, wfy_hi
        real(amrex_real), intent(in   ) :: weight_y(wfy_lo(1):wfy_hi(1), &
            &                                       wfy_lo(2):wfy_hi(2), &
            &                                       wfy_lo(3):wfy_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: wfz_lo, wfz_hi
        real(amrex_real), intent(in   ) :: weight_z(wfz_lo(1):wfz_hi(1), &
            &                                       wfz_lo(2):wfz_hi(2), &
            &                                       wfz_lo(3):wfz_hi(3))

        ! ** IN:  coordinates of staggered (face-centered) faces corresponding to the
        !         MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: cx_lo, cx_hi
        real(amrex_real), intent(in   ) :: coords_x(cx_lo(1):cx_hi(1), &
            &                                       cx_lo(2):cx_hi(2), &
            &                                       cx_lo(3):cx_hi(3), AMREX_SPACEDIM)

        integer(c_int), dimension(3), intent(in   ) :: cy_lo, cy_hi
        real(amrex_real), intent(in   ) :: coords_y(cy_lo(1):cy_hi(1), &
            &                                       cy_lo(2):cy_hi(2), &
            &                                       cy_lo(3):cy_hi(3), AMREX_SPACEDIM)

        integer(c_int), dimension(3), intent(in   ) :: cz_lo, cz_hi
        real(amrex_real), intent(in   ) :: coords_z(cz_lo(1):cz_hi(1), &
            &                                       cz_lo(2):cz_hi(2), &
            &                                       cz_lo(3):cz_hi(3), AMREX_SPACEDIM)

        ! ** IN:  `n_marker`-many marker positions (pos_marker). The `pos_marker` and
        !         `v_marker` arrays are `Vector<RealVecr>` in c-land.
        integer(c_int), intent(in   ) :: n_marker
        real(amrex_real), intent(in   ) :: pos_marker(AMREX_SPACEDIM, n_marker);
        ! ** OUT: each marker vector `v_marker` is interpolated from the stagged
        !         MultiFabs `mf_*`
        real(amrex_real), intent(inout) :: v_marker(AMREX_SPACEDIM, n_marker);
        ! ** IN:  fluid grid discretization
        real(amrex_real), dimension(AMREX_SPACEDIM), intent(in   ) :: dx


        !________________________________________________________________________
        ! i        <= loop counter over markers
        integer :: i
        ! pos      <= position of current marker
        ! v_spread <= current markers vector quantity to spread
        real(amrex_real), dimension(AMREX_SPACEDIM) :: pos, v_spread


        do i =  1, n_marker
            pos = pos_marker(:, i)

            ! skip marker if outside tile box (prevent double-counting)
            if (pos(1) .lt. tile_lo(1)*dx(1) ) cycle
            if (pos(1) .ge. (tile_hi(1)+1)*dx(1) ) cycle

            if (pos(2) .lt. tile_lo(2)*dx(2) ) cycle
            if (pos(2) .ge. (tile_hi(2)+1)*dx(2) ) cycle

            if (pos(3) .lt. tile_lo(3)*dx(3) ) cycle
            if (pos(3) .ge. (tile_hi(3)+1)*dx(3) ) cycle


            v_spread = (/0d0, 0d0, 0d0/) ! v_marker(:, i)


            call interpolate_kernel(lo,       hi,               &
                &                   mf_x,     mfx_lo,   mfx_hi, &
                &                   mf_y,     mfy_lo,   mfy_hi, &
                &                   mf_z,     mfz_lo,   mfz_hi, &
                &                   weight_x, wfx_lo,   wfx_hi, &
                &                   weight_y, wfy_lo,   wfy_hi, &
                &                   weight_z, wfz_lo,   wfz_hi, &
                &                   coords_x, cx_lo,    cx_hi,  &
                &                   coords_y, cy_lo,    cy_hi,  &
                &                   coords_z, cz_lo,    cz_hi,  &
                &                   pos,      v_spread, dx, pkernel_fluid_in)

            v_marker(:, i) = v_spread(:)

        end do

    end subroutine interpolate_markers




    subroutine inv_interpolate_kernel(lo,       hi,               &
            &                         mf_x,     mfx_lo,   mfx_hi, &
            &                         mf_y,     mfy_lo,   mfy_hi, &
            &                         mf_z,     mfz_lo,   mfz_hi, &
            &                         coords_x, cx_lo,    cx_hi,  &
            &                         coords_y, cy_lo,    cy_hi,  &
            &                         coords_z, cz_lo,    cz_hi,  &
            &                         pos,      v_spread, dx, pkernel_fluid_in)


        !________________________________________________________________________
        ! ** work region
      integer(c_int), dimension(3), intent(in   ) :: lo, hi
      integer(c_int), intent(in   ) :: pkernel_fluid_in(1)

        ! ** OUT: vector quantity `v_spread` is spread to (staggered) MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: mfx_lo, mfx_hi
        real(amrex_real), intent(inout) :: mf_x(mfx_lo(1):mfx_hi(1), &
            &                                   mfx_lo(2):mfx_hi(2), &
            &                                   mfx_lo(3):mfx_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: mfy_lo, mfy_hi
        real(amrex_real), intent(inout) :: mf_y(mfy_lo(1):mfy_hi(1), &
            &                                   mfy_lo(2):mfy_hi(2), &
            &                                   mfy_lo(3):mfy_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: mfz_lo, mfz_hi
        real(amrex_real), intent(inout) :: mf_z(mfz_lo(1):mfz_hi(1), &
            &                                   mfz_lo(2):mfz_hi(2), &
            &                                   mfz_lo(3):mfz_hi(3))

        ! ** IN:  coordinates of staggered (face-centered) faces corresponding to the
        !         MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: cx_lo, cx_hi
        real(amrex_real), intent(in   ) :: coords_x(cx_lo(1):cx_hi(1), &
            &                                       cx_lo(2):cx_hi(2), &
            &                                       cx_lo(3):cx_hi(3), AMREX_SPACEDIM)

        integer(c_int), dimension(3), intent(in   ) :: cy_lo, cy_hi
        real(amrex_real), intent(in   ) :: coords_y(cy_lo(1):cy_hi(1), &
            &                                       cy_lo(2):cy_hi(2), &
            &                                       cy_lo(3):cy_hi(3), AMREX_SPACEDIM)

        integer(c_int), dimension(3), intent(in   ) :: cz_lo, cz_hi
        real(amrex_real), intent(in   ) :: coords_z(cz_lo(1):cz_hi(1), &
            &                                       cz_lo(2):cz_hi(2), &
            &                                       cz_lo(3):cz_hi(3), AMREX_SPACEDIM)

        ! ** IN:  quantity to spread (v_spread), given kernel position (pos), and the
        !         fluid grid discretization (dx)
        real(amrex_real), dimension(AMREX_SPACEDIM), intent(in   ) :: pos, v_spread, dx


        !________________________________________________________________________
        ! i, j, k   => face-centered indices
        integer :: i, j, k, ilo, ihi, jlo, jhi, klo, khi, gs
        ! ll        => loop counter over AMREX_SPACEDIM
        integer :: ll, ct_face
        ! invvol    => 1/dx^AMREX_SPACEDIM
        ! weight    => kernel weight function
        real(amrex_real) :: invvol, weight, w_threshold
        ! pos_grid  => (pos - (cell position))/dx
        ! invdx     => 1/dx
        real(amrex_real), dimension(AMREX_SPACEDIM) :: pos_grid, invdx, v_scaled

        !________________________________________________________________________
        ! Using function pointer to specify kernel type - some question as to
        ! optimal approach here. DRL.

        abstract interface
          function kernel_np (r_in)
             use amrex_fort_module, only: amrex_real
             real(amrex_real) :: kernel_np
             real(amrex_real), intent (in) :: r_in
          end function kernel_np
        end interface

        procedure (kernel_np), pointer :: kernel_ptr => null()

        ! TODO: It makes sense to specify this function pointer globally (as the
        ! pkernel_fluid parameter is also a global variable). I.e. when the
        ! pkernel_fluid value is set => we save ourselves the followint IF
        ! branch point. JPB.

        if(pkernel_fluid_in(1) .eq. 3) then
          kernel_ptr => kernel_3p
          gs = 2
        else
          kernel_ptr => kernel_6p
          gs = 4
        endif

        w_threshold = 1e-4

        !________________________________________________________________________
        ! compute geometric quantities : 1/dx and 1/dx^AMREX_SPACEDIM
        invdx(:) = 1d0/dx(:)
        invvol = 1d0
        do ll = 1, AMREX_SPACEDIM
            invvol = invvol * invdx(ll)
        end do


        ilo = max(lo(1), int(pos(1) * invdx(1) - gs))
        ihi = min(hi(1), int(pos(1) * invdx(1) + gs))
        jlo = max(lo(2), int(pos(2) * invdx(2) - gs))
        jhi = min(hi(2), int(pos(2) * invdx(2) + gs))
        klo = max(lo(3), int(pos(3) * invdx(3) - gs))
        khi = min(hi(3), int(pos(3) * invdx(3) + gs))


        !________________________________________________________________________
        ! x-components
        ! do k = lo(3), hi(3)
        !     do j = lo(2), hi(2)
        !         do i = lo(1), hi(1) + 1
        ct_face = 0
        do k = klo, khi
            do j = jlo, jhi
                do i = ilo, ihi + 1

                    pos_grid(:) = pos(:) - coords_x(i, j, k, :)
                    pos_grid(:) = pos_grid(:) * invdx(:)

                    weight = 1d0
                    do ll = 1, AMREX_SPACEDIM
                        weight = weight * kernel_ptr(pos_grid(ll));
                    end do

                    if (weight .gt. w_threshold) then
                        ct_face = ct_face + 1
                    end if

                end do
            end do
        end do

        v_scaled(1) = v_spread(1)/ct_face

        do k = klo, khi
            do j = jlo, jhi
                do i = ilo, ihi + 1

                    pos_grid(:) = pos(:) - coords_x(i, j, k, :)
                    pos_grid(:) = pos_grid(:) * invdx(:)

                    weight = 1d0
                    do ll = 1, AMREX_SPACEDIM
                        weight = weight * kernel_ptr(pos_grid(ll));
                    end do

                    if (weight .gt. w_threshold) then
                        mf_x(i, j, k) = mf_x(i, j, k) + v_scaled(1) / weight
                    end if

                end do
            end do
        end do


        !________________________________________________________________________
        ! y-components
        ! do k = lo(3), hi(3)
        !     do j = lo(2), hi(2) + 1
        !         do i = lo(1), hi(1)
        ct_face = 0
        do k = klo, khi
            do j = jlo, jhi + 1
                do i = ilo, ihi

                    pos_grid(:) = pos(:) - coords_y(i, j, k, :)
                    pos_grid(:) = pos_grid(:) * invdx(:)

                    weight = 1d0
                    do ll = 1, AMREX_SPACEDIM
                        weight = weight * kernel_ptr(pos_grid(ll));
                    end do

                    if (weight .gt. w_threshold) then
                        ct_face = ct_face + 1
                    end if

                end do
            end do
        end do

        v_scaled(2) = v_spread(2)/ct_face

        do k = klo, khi
            do j = jlo, jhi + 1
                do i = ilo, ihi

                    pos_grid(:) = pos(:) - coords_y(i, j, k, :)
                    pos_grid(:) = pos_grid(:) * invdx(:)

                    weight = 1d0
                    do ll = 1, AMREX_SPACEDIM
                        weight = weight * kernel_ptr(pos_grid(ll));
                    end do

                    if (weight .gt. w_threshold) then
                        mf_y(i, j, k) = mf_y(i, j, k) + v_scaled(2) / weight
                    end if

                end do
            end do
        end do



        !________________________________________________________________________
        ! z-components
        ! do k = lo(3), hi(3) + 1
        !     do j = lo(2), hi(2)
        !         do i = lo(1), hi(1)
        ct_face = 0
        do k = klo, khi + 1
            do j = jlo, jhi
                do i = ilo, ihi

                    pos_grid(:) = pos(:) - coords_z(i, j, k, :)
                    pos_grid(:) = pos_grid(:) * invdx(:)

                    weight = 1d0
                    do ll = 1, AMREX_SPACEDIM
                        weight = weight * kernel_ptr(pos_grid(ll));
                    end do

                    if (weight .gt. w_threshold) then
                        ct_face = ct_face + 1
                    end if

                end do
            end do
        end do

        v_scaled(3) = v_spread(3)/ct_face

        do k = klo, khi + 1
            do j = jlo, jhi
                do i = ilo, ihi

                    pos_grid(:) = pos(:) - coords_z(i, j, k, :)
                    pos_grid(:) = pos_grid(:) * invdx(:)

                    weight = 1d0
                    do ll = 1, AMREX_SPACEDIM
                        weight = weight * kernel_ptr(pos_grid(ll));
                    end do

                    if (weight .gt. w_threshold) then
                        mf_z(i, j, k) = mf_z(i, j, k) + v_scaled(3) / weight
                    end if

                end do
            end do
        end do


    end subroutine inv_interpolate_kernel



    subroutine inv_interpolate_markers(lo,         hi,                  &
            &                          tile_lo,    tile_hi,             &
            &                          mf_x,       mfx_lo,   mfx_hi,    &
            &                          mf_y,       mfy_lo,   mfy_hi,    &
            &                          mf_z,       mfz_lo,   mfz_hi,    &
            &                          coords_x,   cx_lo,    cx_hi,     &
            &                          coords_y,   cy_lo,    cy_hi,     &
            &                          coords_z,   cz_lo,    cz_hi,     &
            &                          pos_marker, v_marker, n_marker,  &
            &                          dx, pkernel_fluid_in           ) &
            bind(C, name="inv_interpolate_markers")

        !________________________________________________________________________
        ! ** work region
      integer(c_int), dimension(3), intent(in   ) :: lo, hi, tile_lo, tile_hi
      integer(c_int), intent(in   ) :: pkernel_fluid_in(1)

        ! ** OUT: vector quantity `v_marker` is spread to (staggered) MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: mfx_lo, mfx_hi
        real(amrex_real), intent(inout) :: mf_x(mfx_lo(1):mfx_hi(1), &
            &                                   mfx_lo(2):mfx_hi(2), &
            &                                   mfx_lo(3):mfx_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: mfy_lo, mfy_hi
        real(amrex_real), intent(inout) :: mf_y(mfy_lo(1):mfy_hi(1), &
            &                                   mfy_lo(2):mfy_hi(2), &
            &                                   mfy_lo(3):mfy_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: mfz_lo, mfz_hi
        real(amrex_real), intent(inout) :: mf_z(mfz_lo(1):mfz_hi(1), &
            &                                   mfz_lo(2):mfz_hi(2), &
            &                                   mfz_lo(3):mfz_hi(3))

        ! ** IN:  coordinates of staggered (face-centered) faces corresponding to the
        !         MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: cx_lo, cx_hi
        real(amrex_real), intent(in   ) :: coords_x(cx_lo(1):cx_hi(1), &
            &                                       cx_lo(2):cx_hi(2), &
            &                                       cx_lo(3):cx_hi(3), AMREX_SPACEDIM)

        integer(c_int), dimension(3), intent(in   ) :: cy_lo, cy_hi
        real(amrex_real), intent(in   ) :: coords_y(cy_lo(1):cy_hi(1), &
            &                                       cy_lo(2):cy_hi(2), &
            &                                       cy_lo(3):cy_hi(3), AMREX_SPACEDIM)

        integer(c_int), dimension(3), intent(in   ) :: cz_lo, cz_hi
        real(amrex_real), intent(in   ) :: coords_z(cz_lo(1):cz_hi(1), &
            &                                       cz_lo(2):cz_hi(2), &
            &                                       cz_lo(3):cz_hi(3), AMREX_SPACEDIM)

        ! ** IN:  `n_marker`-many marker positions (pos_marker). The `pos_marker` and
        !         `v_marker` arrays are `Vector<RealVecr>` in c-land.
        integer(c_int), intent(in   ) :: n_marker
        real(amrex_real), intent(in   ) :: pos_marker(AMREX_SPACEDIM, n_marker);
        ! ** IN:  each marker spreads a vector qunantity `v_marker` to the stagged
        !         MultiFabs `mf_*`
        real(amrex_real), intent(in   ) :: v_marker(AMREX_SPACEDIM, n_marker);
        ! ** IN:  fluid grid discretization
        real(amrex_real), dimension(AMREX_SPACEDIM), intent(in   ) :: dx


        !________________________________________________________________________
        ! i        <= loop counter over markers
        integer :: i
        ! pos      <= position of current marker
        ! v_spread <= current markers vector quantity to spread
        real(amrex_real), dimension(AMREX_SPACEDIM) :: pos, v_spread


        do i =  1, n_marker
            pos = pos_marker(:, i)
            v_spread = v_marker(:, i)

            call inv_interpolate_kernel(lo,       hi,               &
                &                       mf_x,     mfx_lo,   mfx_hi, &
                &                       mf_y,     mfy_lo,   mfy_hi, &
                &                       mf_z,     mfz_lo,   mfz_hi, &
                &                       coords_x, cx_lo,    cx_hi,  &
                &                       coords_y, cy_lo,    cy_hi,  &
                &                       coords_z, cz_lo,    cz_hi,  &
                &                       pos,      v_spread, dx, pkernel_fluid_in )

        end do

    end subroutine inv_interpolate_markers


    subroutine inv_spread_kernel(lo,       hi,               &
            &                    mf_x,     mfx_lo,   mfx_hi, &
            &                    mf_y,     mfy_lo,   mfy_hi, &
            &                    mf_z,     mfz_lo,   mfz_hi, &
            &                    coords_x, cx_lo,    cx_hi,  &
            &                    coords_y, cy_lo,    cy_hi,  &
            &                    coords_z, cz_lo,    cz_hi,  &
            &                    pos,      v_interp, dx,     &
            &                    pkernel_fluid_in)


        !________________________________________________________________________
        ! ** work region
        integer(c_int), dimension(3), intent(in   ) :: lo, hi
        integer(c_int), intent(in   ) :: pkernel_fluid_in(1)

        ! ** IN:  interpolating v_interp from staggered MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: mfx_lo, mfx_hi
        real(amrex_real), intent(in   ) :: mf_x(mfx_lo(1):mfx_hi(1), &
            &                                   mfx_lo(2):mfx_hi(2), &
            &                                   mfx_lo(3):mfx_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: mfy_lo, mfy_hi
        real(amrex_real), intent(in   ) :: mf_y(mfy_lo(1):mfy_hi(1), &
            &                                   mfy_lo(2):mfy_hi(2), &
            &                                   mfy_lo(3):mfy_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: mfz_lo, mfz_hi
        real(amrex_real), intent(in   ) :: mf_z(mfz_lo(1):mfz_hi(1), &
            &                                   mfz_lo(2):mfz_hi(2), &
            &                                   mfz_lo(3):mfz_hi(3))

        ! ** IN:  coordinates of staggered (face-centered) faces corresponding to the
        !         MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: cx_lo, cx_hi
        real(amrex_real), intent(in   ) :: coords_x(cx_lo(1):cx_hi(1), &
            &                                       cx_lo(2):cx_hi(2), &
            &                                       cx_lo(3):cx_hi(3), AMREX_SPACEDIM)

        integer(c_int), dimension(3), intent(in   ) :: cy_lo, cy_hi
        real(amrex_real), intent(in   ) :: coords_y(cy_lo(1):cy_hi(1), &
            &                                       cy_lo(2):cy_hi(2), &
            &                                       cy_lo(3):cy_hi(3), AMREX_SPACEDIM)

        integer(c_int), dimension(3), intent(in   ) :: cz_lo, cz_hi
        real(amrex_real), intent(in   ) :: coords_z(cz_lo(1):cz_hi(1), &
            &                                       cz_lo(2):cz_hi(2), &
            &                                       cz_lo(3):cz_hi(3), AMREX_SPACEDIM)

        ! ** OUT: vector to interpolate (v_interp) from staggered `mf_*`
        real(amrex_real), dimension(AMREX_SPACEDIM), intent(inout) :: v_interp
        ! ** IN:  kernel position (pos), and the fluid grid discretization (dx)
        real(amrex_real), dimension(AMREX_SPACEDIM), intent(in   ) :: pos, dx


        !________________________________________________________________________
        ! i,    j,    k    => face-centered indices
        integer :: i, j, k, ilo, ihi, jlo, jhi, klo, khi, gs
        ! ll               => loop counter over AMREX_SPACEDIM
        integer :: ll, ct_face
        ! weight           => kernel weight function
        ! wfrac            => weight/total_weight (weight_{x,y,z})
        real(amrex_real) :: weight, wfrac, invvol, w_threshold
        ! pos_grid         => (pos - (cell position))/dx
        ! invdx            => 1/dx
        real(amrex_real), dimension(AMREX_SPACEDIM) :: pos_grid, invdx, interp_scaled

        !________________________________________________________________________
        ! Using function pointer to specify kernel type - some question as to
        ! optimal approach here. DRL.

        abstract interface
          function kernel_np (r_in)
             use amrex_fort_module, only: amrex_real
             real(amrex_real) :: kernel_np
             real(amrex_real), intent (in) :: r_in
          end function kernel_np
        end interface

        procedure (kernel_np), pointer :: kernel_ptr => null()

        ! TODO: It makes sense to specify this function pointer globally (as the
        ! pkernel_fluid parameter is also a global variable). I.e. when the
        ! pkernel_fluid value is set => we save ourselves the followint IF
        ! branch point. JPB.

        if(pkernel_fluid_in(1) .eq. 3) then
          kernel_ptr => kernel_3p
          gs = 2
        else
          kernel_ptr => kernel_6p
          gs = 4
        endif


        w_threshold = 1e-4


        !________________________________________________________________________
        ! compute geometric quantity 1/dx
        invdx(:) = 1d0/dx(:)
        invvol = 1d0
        do ll = 1, AMREX_SPACEDIM
            invvol = invvol * invdx(ll)
        end do


        ! ilo = max(lo(1), int(pos(1) * invdx(1) - 3))
        ! ihi = min(hi(1), int(pos(1) * invdx(1) + 3))
        ! jlo = max(lo(2), int(pos(2) * invdx(2) - 3))
        ! jhi = min(hi(2), int(pos(2) * invdx(2) + 3))
        ! klo = max(lo(3), int(pos(3) * invdx(3) - 3))
        ! khi = min(hi(3), int(pos(3) * invdx(3) + 3))
        ilo = int(pos(1) * invdx(1) - gs)
        ihi = int(pos(1) * invdx(1) + gs)
        jlo = int(pos(2) * invdx(2) - gs)
        jhi = int(pos(2) * invdx(2) + gs)
        klo = int(pos(3) * invdx(3) - gs)
        khi = int(pos(3) * invdx(3) + gs)



        interp_scaled(:) = 0

        !________________________________________________________________________
        ! x-components
        ! do k = lo(3), hi(3)
        !     do j = lo(2), hi(2)
        !         do i = lo(1), hi(1) + 1
        ct_face = 0
        do k = klo, khi
            do j = jlo, jhi
                do i = ilo, ihi + 1

                    pos_grid(:) = pos(:) - coords_x(i, j, k, :)
                    pos_grid(:) = pos_grid(:) * invdx(:)

                    weight = 1d0
                    do ll = 1, AMREX_SPACEDIM
                        weight = weight * kernel_ptr(pos_grid(ll));
                    end do

                    if (weight .gt. 0) then
                        interp_scaled(1) = interp_scaled(1) + mf_x(i, j, k) / weight / invvol
                        ct_face = ct_face + 1
                    end if

                end do
            end do
        end do

        v_interp(1) = v_interp(1) + interp_scaled(1)/ct_face

        !________________________________________________________________________
        ! y-components
        ! do k = lo(3), hi(3)
        !     do j = lo(2), hi(2) + 1
        !         do i = lo(1), hi(1)
        ct_face = 0
        do k = klo, khi
            do j = jlo, jhi + 1
                do i = ilo, ihi

                    pos_grid(:) = pos(:) - coords_y(i, j, k, :)
                    pos_grid(:) = pos_grid(:) * invdx(:)

                    weight = 1d0
                    do ll = 1, AMREX_SPACEDIM
                        weight = weight * kernel_ptr(pos_grid(ll));
                    end do

                    if (weight .gt. 0) then
                        interp_scaled(2) = interp_scaled(2) + mf_y(i, j, k) / weight / invvol
                        ct_face = ct_face + 1
                    end if

                end do
            end do
        end do

        v_interp(2) = v_interp(2) + interp_scaled(2)/ct_face

        !________________________________________________________________________
        ! z-components
        ! do k = lo(3), hi(3) + 1
        !     do j = lo(2), hi(2)
        !         do i = lo(1), hi(1)
        ct_face = 0;
        do k = klo, khi + 1
            do j = jlo, jhi
                do i = ilo, ihi

                    pos_grid(:) = pos(:) - coords_z(i, j, k, :)
                    pos_grid(:) = pos_grid(:) * invdx(:)

                    weight = 1d0
                    do ll = 1, AMREX_SPACEDIM
                        weight = weight * kernel_ptr(pos_grid(ll));
                    end do

                    if (weight .gt. 0) then
                        interp_scaled(3) = interp_scaled(3) + mf_z(i, j, k) / weight / invvol
                        ct_face = ct_face + 1
                    end if

                end do
            end do
        end do

        v_interp(3) = v_interp(3) + interp_scaled(3)/ct_face


    end subroutine inv_spread_kernel



    subroutine inv_spread_markers(lo,         hi,                  &
            &                     tile_lo,    tile_hi,             &
            &                     mf_x,       mfx_lo,   mfx_hi,    &
            &                     mf_y,       mfy_lo,   mfy_hi,    &
            &                     mf_z,       mfz_lo,   mfz_hi,    &
            &                     coords_x,   cx_lo,    cx_hi,     &
            &                     coords_y,   cy_lo,    cy_hi,     &
            &                     coords_z,   cz_lo,    cz_hi,     &
            &                     pos_marker, v_marker, n_marker,  &
            &                     dx, pkernel_fluid_in           ) &
            bind(C, name="inv_spread_markers")

        !________________________________________________________________________
        ! ** work region
        integer(c_int), dimension(3), intent(in   ) :: lo, hi, tile_lo, tile_hi
        integer(c_int), intent(in   ) :: pkernel_fluid_in(1)

        ! ** IN:  vector quantity `v_marker` is interpolated from (staggered) MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: mfx_lo, mfx_hi
        real(amrex_real), intent(in   ) :: mf_x(mfx_lo(1):mfx_hi(1), &
            &                                   mfx_lo(2):mfx_hi(2), &
            &                                   mfx_lo(3):mfx_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: mfy_lo, mfy_hi
        real(amrex_real), intent(in   ) :: mf_y(mfy_lo(1):mfy_hi(1), &
            &                                   mfy_lo(2):mfy_hi(2), &
            &                                   mfy_lo(3):mfy_hi(3))

        integer(c_int), dimension(3), intent(in   ) :: mfz_lo, mfz_hi
        real(amrex_real), intent(in   ) :: mf_z(mfz_lo(1):mfz_hi(1), &
            &                                   mfz_lo(2):mfz_hi(2), &
            &                                   mfz_lo(3):mfz_hi(3))

        ! ** IN:  coordinates of staggered (face-centered) faces corresponding to the
        !         MultiFabs `mf_*`
        integer(c_int), dimension(3), intent(in   ) :: cx_lo, cx_hi
        real(amrex_real), intent(in   ) :: coords_x(cx_lo(1):cx_hi(1), &
            &                                       cx_lo(2):cx_hi(2), &
            &                                       cx_lo(3):cx_hi(3), AMREX_SPACEDIM)

        integer(c_int), dimension(3), intent(in   ) :: cy_lo, cy_hi
        real(amrex_real), intent(in   ) :: coords_y(cy_lo(1):cy_hi(1), &
            &                                       cy_lo(2):cy_hi(2), &
            &                                       cy_lo(3):cy_hi(3), AMREX_SPACEDIM)

        integer(c_int), dimension(3), intent(in   ) :: cz_lo, cz_hi
        real(amrex_real), intent(in   ) :: coords_z(cz_lo(1):cz_hi(1), &
            &                                       cz_lo(2):cz_hi(2), &
            &                                       cz_lo(3):cz_hi(3), AMREX_SPACEDIM)

        ! ** IN:  `n_marker`-many marker positions (pos_marker). The `pos_marker` and
        !         `v_marker` arrays are `Vector<RealVecr>` in c-land.
        integer(c_int),   intent(in   ) :: n_marker
        real(amrex_real), intent(in   ) :: pos_marker(AMREX_SPACEDIM, n_marker);
        ! ** OUT: each marker vector `v_marker` is interpolated from the stagged
        !         MultiFabs `mf_*`
        real(amrex_real), intent(inout) :: v_marker(AMREX_SPACEDIM, n_marker);
        ! ** IN:  fluid grid discretization
        real(amrex_real), dimension(AMREX_SPACEDIM), intent(in   ) :: dx


        !________________________________________________________________________
        ! i        <= loop counter over markers
        integer :: i
        ! pos      <= position of current marker
        ! v_spread <= current markers vector quantity to spread
        real(amrex_real), dimension(AMREX_SPACEDIM) :: pos, v_spread


        do i =  1, n_marker
            pos = pos_marker(:, i)
            v_spread = (/0d0, 0d0, 0d0/) ! v_marker(:, i)

            call inv_spread_kernel(lo,       hi,               &
                &                   mf_x,     mfx_lo,   mfx_hi, &
                &                   mf_y,     mfy_lo,   mfy_hi, &
                &                   mf_z,     mfz_lo,   mfz_hi, &
                &                   coords_x, cx_lo,    cx_hi,  &
                &                   coords_y, cy_lo,    cy_hi,  &
                &                   coords_z, cz_lo,    cz_hi,  &
                &                   pos,      v_spread, dx, pkernel_fluid_in )

            v_marker(:, i) = v_spread(:)

        end do

    end subroutine inv_spread_markers



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

        !________________________________________________________________________
        ! ** work region
        integer(c_int), dimension(3), intent(in   ) :: lo, hi

        ! ** OUT: (staggered) Immersed-boundary velocity components
        integer(c_int), dimension(3), intent(in   ) :: udlo, udhi, vdlo, vdhi, wdlo, wdhi
        real(amrex_real), intent(  out) :: u_d(udlo(1):udhi(1), &
            &                                  udlo(2):udhi(2), &
            &                                  udlo(3):udhi(3))
        real(amrex_real), intent(  out) :: v_d(vdlo(1):vdhi(1), &
            &                                  vdlo(2):vdhi(2), &
            &                                  vdlo(3):vdhi(3))
        real(amrex_real), intent(  out) :: w_d(wdlo(1):wdhi(1), &
            &                                  wdlo(2):wdhi(2), &
            &                                  wdlo(3):wdhi(3))

        ! ** IN:  (staggered) RHS terms w/o IBM force terms
        integer(c_int), dimension(3), intent(in   ) :: uslo, ushi, vslo, vshi, wslo, wshi
        real(amrex_real), intent(in   ) :: u_s(uslo(1):ushi(1), &
            &                                  uslo(2):ushi(2), &
            &                                  uslo(3):ushi(3))
        real(amrex_real), intent(in   ) :: v_s(vslo(1):vshi(1), &
            &                                  vslo(2):vshi(2), &
            &                                  vslo(3):vshi(3))
        real(amrex_real), intent(in   ) :: w_s(wslo(1):wshi(1), &
            &                                  wslo(2):wshi(2), &
            &                                  wslo(3):wshi(3))

        ! ** OUT: (cell-centered) effective tag
        integer(c_int), dimension(3), intent(in   ) :: etlo, ethi
        integer(c_int),   intent(  out) :: et(etlo(1):ethi(1), &
            &                                 etlo(2):ethi(2), &
            &                                 etlo(3):ethi(3))

        integer(c_int), dimension(3), intent(in   ) :: philo, phihi
        real(amrex_real), intent(in   ) :: phi(philo(1):phihi(1), &
            &                                  philo(2):phihi(2), &
            &                                  philo(3):phihi(3))

        ! ** IN:  (cell-centered) tags
        integer(c_int), dimension(3), intent(in   ) :: taglo, taghi
        integer(c_int),   intent(in   ) :: tag(taglo(1):taghi(1), &
            &                                  taglo(2):taghi(2), &
            &                                  taglo(3):taghi(3))

        ! ** IN:  (staggered) IBM (desiered) veloctity
        integer(c_int), dimension(3), intent(in   ) :: vello, velhi
        real(amrex_real), intent(in   ) :: vel(vello(1):velhi(1), &
            &                                  vello(2):velhi(2), &
            &                                  vello(3):velhi(3), 3)


        !________________________________________________________________________
        ! i,    j,    k    => face-centered indices
        ! i_cc, j_cc, k_cc => cell-centered indices
        integer :: i, j, k, i_cc, j_cc, k_cc
        integer :: eff_tag

        !________________________________________________________________________
        ! x-components
        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1) + 1

                    i_cc = min(i, hi(1))
                    et(i_cc, j, k) = 0

                    ! eff_tag = effective_tag(i, j, k, 0, tag, taglo, taghi)
                    ! Tag interface and internal cells
                    if (tag(i_cc, j, k) .ge. 1) then
                        eff_tag = 1
                    else
                        eff_tag = 0
                    end if

                    if ( eff_tag .ge. 1 ) then
                        u_d(i, j, k) = vel(i, j, k, 1)
                        et(i_cc, j, k)  = 1
                    else
                        u_d(i, j, k) = 0. ! u_s(i, j, k)
                    end if

                end do
            end do
        end do


        !________________________________________________________________________
        ! y-components
        do k = lo(3), hi(3)
            do j = lo(2), hi(2) + 1
                do i = lo(1), hi(1)

                    j_cc = min(j, hi(2))
                    et(i, j_cc, k) = 0

                    ! eff_tag = effective_tag(i, j, k, 0, tag, taglo, taghi)
                    ! Tag interface and internal cells
                    if (tag(i, j_cc, k) .ge. 1) then
                        eff_tag = 1
                    else
                        eff_tag = 0
                    end if

                    if ( eff_tag .ge. 1 ) then
                        v_d(i, j, k) = vel(i, j, k, 2)
                        et(i, j_cc, k)  = 1
                    else
                        v_d(i, j, k) = 0. ! v_s(i, j, k)
                    end if

                end do
            end do
        end do


        !________________________________________________________________________
        ! z-components
        do k = lo(3), hi(3) + 1
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    k_cc = min(k, hi(3))
                    et(i, j, k_cc) = 0

                    ! eff_tag = effective_tag(i, j, k, 0, tag, taglo, taghi)
                    ! Tag interface and internal cells
                    if (tag(i, j, k_cc) .ge. 1) then
                        eff_tag = 1
                    else
                        eff_tag = 0
                    end if

                    if ( eff_tag .ge. 1 ) then
                        w_d(i, j, k) = vel(i, j, k, 3)
                        et(i, j, k_cc)  = 1
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

        !________________________________________________________________________
        ! ** work region
        integer(c_int), dimension(3), intent(in   ) :: lo, hi


        ! ** OUT: (staggered) Immersed-boundary spreading/interpolation weights
        integer(c_int), dimension(3), intent(in   ) :: fulo, fuhi, fvlo, fvhi, fwlo, fwhi
        real(amrex_real), intent(  out) :: f_u(fulo(1):fuhi(1), &
            &                                  fulo(2):fuhi(2), &
            &                                  fulo(3):fuhi(3))
        real(amrex_real), intent(  out) :: f_v(fvlo(1):fvhi(1), &
            &                                  fvlo(2):fvhi(2), &
            &                                  fvlo(3):fvhi(3))
        real(amrex_real), intent(  out) :: f_w(fwlo(1):fwhi(1), &
            &                                  fwlo(2):fwhi(2), &
            &                                  fwlo(3):fwhi(3))

        ! ** IN:  (cell-centered) effective tag
        integer(c_int), dimension(3), intent(in   ) :: etlo, ethi
        integer(c_int), intent(in   ) :: et(etlo(1):ethi(1), &
            &                               etlo(2):ethi(2), &
            &                               etlo(3):ethi(3))


        !________________________________________________________________________
        ! i,    j,    k    => face-centered indices
        ! i_cc, j_cc, k_cc => cell-centered indices
        integer :: i, j, k, i_cc, j_cc, k_cc


        !________________________________________________________________________
        ! x-components
        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1) + 1

                    i_cc = min(i, hi(1))

                    if (  et(i_cc, j, k)  .ge. 1 ) then
                        f_u(i, j, k) = 1.
                    else
                        f_u(i, j, k) = 0.
                    end if
                end do
            end do
        end do


        !________________________________________________________________________
        ! y-components
        do k = lo(3), hi(3)
            do j = lo(2), hi(2) + 1
                do i = lo(1), hi(1)

                    j_cc = min(j, hi(2))

                    if (  et(i, j_cc, k)  .ge. 1 ) then
                        f_v(i, j, k) = 1.
                    else
                        f_v(i, j, k) = 0.
                    end if
                end do
            end do
        end do


        !________________________________________________________________________
        ! z-components
        do k = lo(3), hi(3) + 1
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                    k_cc = min(k, hi(3))

                    if (  et(i, j, k_cc)  .ge. 1 ) then
                        f_w(i, j, k) = 1.
                    else
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

        !________________________________________________________________________
        ! ** work region
        integer(c_int), dimension(3), intent(in   ) :: lo, hi

        ! ** OUT: (staggered) explicit IBM forces
        integer(c_int), dimension(3), intent(in   ) :: fulo, fuhi, fvlo, fvhi, fwlo, fwhi
        real(amrex_real), intent(  out) :: f_u(fulo(1):fuhi(1), &
            &                                  fulo(2):fuhi(2), &
            &                                  fulo(3):fuhi(3))
        real(amrex_real), intent(  out) :: f_v(fvlo(1):fvhi(1), &
            &                                  fvlo(2):fvhi(2), &
            &                                  fvlo(3):fvhi(3))
        real(amrex_real), intent(  out) :: f_w(fwlo(1):fwhi(1), &
            &                                  fwlo(2):fwhi(2), &
            &                                  fwlo(3):fwhi(3))

        ! ** IN:  (staggered) RHS terms w/o IBM forcer terms
        integer(c_int), dimension(3), intent(in   ) :: uglo, ughi, vglo, vghi, wglo, wghi
        real(amrex_real), intent(in   ) :: u_g(uglo(1):ughi(1), &
            &                                  uglo(2):ughi(2), &
            &                                  uglo(3):ughi(3))
        real(amrex_real), intent(in   ) :: v_g(vglo(1):vghi(1), &
            &                                  vglo(2):vghi(2), &
            &                                  vglo(3):vghi(3))
        real(amrex_real), intent(in   ) :: w_g(wglo(1):wghi(1), &
            &                                  wglo(2):wghi(2), &
            &                                  wglo(3):wghi(3))

        ! ** IN:  (staggered) IBM (desired) velocity
        integer(c_int), dimension(3), intent(in   ) :: udlo, udhi, vdlo, vdhi, wdlo, wdhi
        real(amrex_real), intent(in   ) :: u_d(udlo(1):udhi(1), &
            &                                  udlo(2):udhi(2), &
            &                                  udlo(3):udhi(3))
        real(amrex_real), intent(in   ) :: v_d(vdlo(1):vdhi(1), &
            &                                  vdlo(2):vdhi(2), &
            &                                  vdlo(3):vdhi(3))
        real(amrex_real), intent(in   ) :: w_d(wdlo(1):wdhi(1), &
            &                                  wdlo(2):wdhi(2), &
            &                                  wdlo(3):wdhi(3))

        ! ** IN:  (cell_centered) effective tags
        integer(c_int), dimension(3), intent(in   ) :: etlo, ethi
        integer(c_int),   intent(in   ) :: et(etlo(1):ethi(1), &
            &                                 etlo(2):ethi(2), &
            &                                 etlo(3):ethi(3))

        ! ** IN: time-step
        real(amrex_real), intent(in   ) :: dt


        !________________________________________________________________________
        ! i,    j,    k    => face-centered indices
        ! i_cc, j_cc, k_cc => cell-centered indices
        integer :: i, j, k, i_cc, j_cc, k_cc



        !________________________________________________________________________
        ! x-components
        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1) + 1

                    i_cc = min(i, hi(1))

                    if ( et(i_cc, j, k) .ge. 1 ) then
                        f_u(i, j, k) = (u_d(i, j, k) - u_g(i, j, k)) / dt
                    else
                        f_u(i, j, k) = 0.
                    end if
                end do
            end do
        end do


        !________________________________________________________________________
        ! y-components
        do k = lo(3), hi(3)
            do j = lo(2), hi(2) + 1
                do i = lo(1), hi(1)

                    j_cc = min(j, hi(2))

                    if ( et(i, j_cc, k) .ge. 1 ) then
                        f_v(i, j, k) = (v_d(i, j, k) - v_g(i, j, k)) / dt
                    else
                        f_v(i, j, k) = 0.
                    end if
                end do
            end do
        end do


        !________________________________________________________________________
        ! z-components
        do k = lo(3), hi(3) + 1
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                    k_cc = min(k, hi(3))

                    if ( et(i, j, k_cc) .ge. 1 ) then
                        f_w(i, j, k) = (w_d(i, j, k) - w_g(i, j, k)) / dt
                    else
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
