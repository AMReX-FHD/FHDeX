module particle_functions_module

  use amrex_fort_module, only: amrex_real
  use iso_c_binding, only: c_ptr, c_int, c_f_pointer
  use cell_sorted_particle_module, only : particle_t, remove_particle_from_cell
  use species_type_module, only: species_t
  use rng_functions_module
  use paramplane_type_module
  use common_namelist_module, only: k_B, T_init, permittivity, eepsilon, images, pkernel_es, &
                                    prob_lo, prob_hi, bc_es_lo, bc_es_hi, rmin, p_int_tog, &
                                    fixed_dt, graphene_tog, mass, particle_n0, particle_neff, &
                                    visc_type, variance_coef_mom, pkernel_fluid, dry_move_tog, &
                                    move_tog, nspecies, drag_tog, es_tog, rfd_tog, qval, &
                                    visc_coef, bc_vel_lo, bc_vel_hi, rfd_delta, sigma

  implicit none

  private
  
contains

  subroutine repulsive_force(part1,part2,dx, dr2) &
                             bind(c,name="repulsive_force")
        
    type(particle_t), intent(inout) :: part1 
    type(particle_t), intent(inout) :: part2
    real(amrex_real), intent(in) :: dx(3), dr2

    real(amrex_real) :: ff

    !This needs to be updated to handle multispecies.
    ff = part1%eepsilon*4*(1./(dr2*dr2*dr2*dr2))* &
         (-12*((part1%sigma+part2%sigma)/4d0)**12/(dr2*dr2*dr2)+ &
          6*((part1%sigma+part2%sigma)/4d0)**6)

    part1%force = part1%force - dx*ff

  end subroutine repulsive_force

  subroutine force_function2(part1,part2,domsize) &
                             bind(c,name="force_function2")

    type(particle_t), intent(inout) :: part1 !is this defined correctly?
    type(particle_t), intent(inout) :: part2
    real(amrex_real), intent(in) :: domsize(3)

    integer :: i,j,k, bound, ii, jj, kk, imagecounter, xswitch, partno, n, pairs, imag
    real(amrex_real) :: dx(3), dx0(3), dr, dr2, rtdr2, maxdist, ee

    ee = (1d0/(permittivity*4*3.142))

    dx0 = part1%pos-part2%pos

    ii=0
    jj=0
    kk=0

    pairs = 0

    if(images .eq. 0) then
       imag = 1
    else
       imag = images
    endif

    maxdist = (imag)*domsize(1)

    if(imag*domsize(2) .lt. maxdist) then
       maxdist = (imag)*domsize(2)
    endif

    if((imag*domsize(3) .lt. imag*domsize(2)) .or. (imag*domsize(3) .lt. imag*domsize(1))) then
       maxdist = (imag)*domsize(3)
    endif

    maxdist = 0.99*maxdist 


    do ii = -images, images
       do jj = -images, images 
          do kk = -images, images

             dx(1) = dx0(1) - ii*domsize(1)
             dx(2) = dx0(2) - jj*domsize(2)
             dx(3) = dx0(3) - kk*domsize(3)

             dr2 = dot_product(dx,dx)

             if(dr2 .ne. 0 ) then !might find a more careful way of doing this

                rtdr2 = sqrt(dr2)

                if(rtdr2 .lt. maxdist) then
                   part1%force = part1%force + ee*(dx/rtdr2)*part1%q*part2%q/dr2

                   pairs = pairs + 1

                endif
             endif

          end do
       end do
    end do

  end subroutine force_function2


  ! check if particle is near a boundary--assumes boundary is parallel to either xy,xz,yz plane
  ! inputs: particle p     : just need the particles position and p3m radius 
  !         int dir        : either +-1,+-2,+-3
  !         int near_wall  : output, 0=part's p3m_radius/2 does not overlap with wall, 1=yes    
  subroutine near_wall_check(part, dir, near_wall) &
                             bind(c, name="near_wall_check")

    type(particle_t), intent(in   ) :: part
    integer,          intent(in   ) :: dir
    integer,          intent(inout) :: near_wall

    integer :: abs_dir
    abs_dir = abs(dir)

    if ((abs_dir.eq.1) .or. (abs_dir.eq.2) .or. (abs_dir.eq.3)) then 
       if (dir .gt. 0) then
          if ((part%pos(abs_dir) + 0.5*part%p3m_radius) .gt. prob_hi(abs_dir)) then 
             near_wall = 1
          endif
       else
          if ((part%pos(abs_dir) - 0.5*part%p3m_radius) .lt. prob_lo(abs_dir)) then 
             near_wall = 1
          endif
       endif
    else 
       ! throw error!
       print*, "Incorrect function call!!!!" 
    endif

  end subroutine near_wall_check

  ! calculate location of an image charge reflected across some plane--assumes boundary
  ! is parallel to either xy,xz,yz plane
  ! inputs: particle p             : just need the particles position 
  !         int dir                : either +-1,+-2,+-3
  !         real(3) im_charge_pos  : output, 0=part's p3m_radius does not overlap with wall, 1=yes    
  !
  !         abs(dir) > 0 means reflect part pos across plane of x(dir) = prob_hi(dir)
  !                  < 0 means reflect part pos across plane of x(dir) = prob_lo(dir)
  subroutine calc_im_charge_loc(part, dir, im_charge_pos) &
                                bind(c, name="calc_im_charge_loc")

    type(particle_t), intent(in   ) :: part
    integer,          intent(in   ) :: dir
    real(amrex_real), intent(inout) :: im_charge_pos(3)

    integer :: abs_dir
    abs_dir = abs(dir)

    if (abs_dir.eq.1) then 
       if (dir .gt. 0) then 
          ! place image charge to the right
          im_charge_pos(2) = part%pos(2)
          im_charge_pos(3) = part%pos(3)
          im_charge_pos(1) = part%pos(1) + 2.d0*(prob_hi(1)-part%pos(1))
       else 
          ! place image charge to the left
          im_charge_pos(2) = part%pos(2)
          im_charge_pos(3) = part%pos(3)
          im_charge_pos(1) = part%pos(1) - 2.d0*(part%pos(1)-prob_lo(1))
       endif
    else if (abs_dir.eq.2) then 
       if (dir .gt. 0) then 
          ! place image charge above
          im_charge_pos(1) = part%pos(1)
          im_charge_pos(3) = part%pos(3)
          im_charge_pos(2) = part%pos(2) + 2.d0*(prob_hi(2)-part%pos(2))
       else 
          ! place image charge below
          im_charge_pos(1) = part%pos(1)
          im_charge_pos(3) = part%pos(3)
          im_charge_pos(2) = part%pos(2) - 2.d0*(part%pos(2)-prob_lo(2))
       endif
    else if (abs_dir.eq.3) then 
       if (dir .gt. 0) then 
          ! place image charge behind
          im_charge_pos(1) = part%pos(1)
          im_charge_pos(2) = part%pos(2)
          im_charge_pos(3) = part%pos(3) + 2.d0*(prob_hi(3)-part%pos(3))
       else 
          ! place image charge in front
          im_charge_pos(1) = part%pos(1)
          im_charge_pos(2) = part%pos(2)
          im_charge_pos(3) = part%pos(3) - 2.d0*(part%pos(3)-prob_lo(3))
       endif
    else 
       print*, "Incorrect function call!!!!" 
       ! throw error!
    endif

  end subroutine calc_im_charge_loc

  subroutine compute_p3m_force_mag(r, mag, dx) &
                                   bind(c, name="compute_p3m_force_mag")

    real(amrex_real), intent(in   ) :: r  
    real(amrex_real), intent(inout) :: mag
    real(amrex_real), intent(in   ) :: dx(3)

    real(amrex_real) :: vals6(70), points6(70), vals4(50), points4(50),r_cell_frac, m, r_norm
    integer          :: r_cell

    !!!!!!!!!!!!!!!!!!!!!!!!!!
    ! p6 tables:
    !!!!!!!!!!!!!!!!!!!!!!!!!!

    ! this is force per dx**2
    vals6 =(/0., 0.0130335, 0.0259347, 0.0384891, 0.0507132, 0.0624416, 0.0735093, &
         0.0837511, 0.0931669, 0.101922, 0.109521, 0.116293, 0.121745, &
         0.126535, 0.130335, 0.132978, 0.134134, 0.135621, 0.13529, 0.134629, &
         0.133308, 0.130995, 0.129343, 0.126205, 0.122075, 0.118276, 0.115137, &
         0.110677, 0.106547, 0.102748, 0.0979574, 0.0936624, 0.0898631, &
         0.0860637, 0.0816036, 0.0782998, 0.0745005, 0.0707011, 0.0680581, &
         0.0645891, 0.0616157, 0.0586423, 0.0561644, 0.0535214, 0.0513739, &
         0.0487309, 0.0469138, 0.0450967, 0.0431145, 0.0414626, 0.0399759, &
         0.038324, 0.0366721, 0.0355157, 0.034029, 0.0328727, 0.0317164, &
         0.0307252, 0.0295689, 0.0285778, 0.0275866, 0.0267607, 0.0259347, &
         0.0251088, 0.0242829, 0.0236221, 0.0229613, 0.0221354, 0.0214746, &
         0.0209791/)
    ! these are fractions of cell size dx
    points6 =(/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, &
         1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, &
         2.8, 2.9, 3., 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4., 4.1, &
         4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5., 5.1, 5.2, 5.3, 5.4, 5.5, &
         5.6, 5.7, 5.8, 5.9, 6., 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9/)

    vals4 =(/0._16, 0.0191993_16, 0.0384104_16, 0.0561438_16, 0.0745141_16, 0.0914394_16, 0.106869_16, &
            0.121719_16, 0.134211_16, 0.145262_16, 0.154958_16, 0.162921_16, 0.168841_16, 0.173264_16, &
            0.176417_16, 0.17854_16, 0.177695_16, 0.175969_16, 0.174664_16, 0.17083_16, 0.165755_16, &
            0.162663_16, 0.156037_16, 0.148626_16, 0.144648_16, 0.1381_16, 0.131473_16, 0.12615_16, &
            0.120057_16, 0.112581_16, 0.106808_16, 0.100863_16, 0.0947845_16, 0.0905315_16, 0.0856055_16, &
            0.0818209_16, 0.0768459_16, 0.0728331_16, 0.0694836_16, 0.0650267_16, 0.0627171_16, 0.0591556_16, &
            0.0561682_16, 0.0539003_16, 0.0511791_16, 0.0492405_16, 0.0470903_16, 0.0453286_16, 0.043644_16, 0.0417381_16/)

    ! these are fractions of cell size dx
    points4 =(/0._16, 0.1_16, 0.2_16, 0.3_16, 0.4_16, 0.5_16, 0.6_16, 0.7_16, 0.8_16, 0.9_16, 1.0_16, 1.1_16, 1.2_16, 1.3_16, &
              1.4_16, 1.5_16, 1.6_16, 1.7_16, 1.8_16, 1.9_16, 2.0_16, 2.1_16, 2.2_16, 2.3_16, 2.4_16, 2.5_16, 2.6_16, 2.7_16, &
              2.8_16, 2.9_16, 3.0_16, 3.1_16, 3.2_16, 3.3_16, 3.4_16, 3.5_16, 3.6_16, 3.7_16, 3.8_16, 3.9_16, 4.0_16, 4.1_16, &
              4.2_16, 4.3_16, 4.4_16, 4.5_16, 4.6_16, 4.7_16, 4.8_16, 4.9_16/)

    !test = (/0.1_16,0.2_16,0.3_16,0.4_16/)

    ! Divison/multiplication by 0.1 below is bc the spacing in the table above is 0.1
    r_norm = r/dx(1)                ! separation dist in units of dx=dy=dz
    r_cell = floor(r_norm/0.1)     ! scaling by 10 allows r_cell to index the points/val arrays above 
    r_cell_frac = r_norm/0.1d0-r_cell !
    r_cell = r_cell + 1             ! shift simply bc points/val array indices begin at 1, but first val=0.0
    r_cell_frac = r_cell_frac*0.1d0   ! for use in point-slope formula below

    !print *, "RCELLFRAC: ", r_cell_frac, r_norm
    !print *, "r: ", r_cell_frac
    !print *, "cr: ", r_cell, r
    !print *, "r_cell_frac: ", r_cell_frac
    if (pkernel_es .eq. 6) then 

       ! do linear interpolation of force between vals(i+1) and val(i) 
       m = (vals6(r_cell+1)-vals6(r_cell))/(points6(r_cell+1)-points6(r_cell))

       mag  = m*r_cell_frac + vals6(r_cell)

    elseif (pkernel_es .eq. 4) then 

       ! do linear interpolation of force between vals(i+1) and val(i) 
       m = (vals4(r_cell+1)-vals4(r_cell))/(points4(r_cell+1)-points4(r_cell))

       mag  = m*r_cell_frac + vals4(r_cell)

        !print *, "MAG: ", m, r_cell_frac, vals4(r_cell), r_cell

    else
       print*, "P3M implemented only for pkernel 4 and 6!"
       ! throw error
    endif

  end subroutine compute_p3m_force_mag

  subroutine calculate_force(particles, np, lo, hi, cell_part_ids, cell_part_cnt, clo, chi, &
                             plo, phi, partno) &
                             bind(c,name="calculate_force")

    !this is everything we pass in
    type(particle_t), intent(inout), target :: particles(np)
    integer(c_int), intent(in) :: np, partno 
    integer(c_int), intent(in) :: lo(3), hi(3)
    integer(c_int), intent(in) :: clo(3), chi(3) 
    type(c_ptr), intent(inout) :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int), intent(inout) :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    real(amrex_real), intent(in) :: plo(3), phi(3) 

    integer :: i, j, k, p, n
    type(particle_t), pointer :: part
    type(particle_t), pointer :: part2 !added this
    real(amrex_real) :: dx(3), domsize(3)

    domsize = phi - plo

    !calculate N^2 interaction

    part => particles(partno) !this defines one particle--we can access all the data by doing part%something

    do n = 1, np

       part2 => particles(n) !this defines one particle--we can access all the data by doing part%something

       if(n .ne. partno) then

          call force_function2(part,part2,domsize)

       endif

    end do

  end subroutine calculate_force

  subroutine amrex_compute_p3m_sr_correction_nl(rparticles, np, neighbors, &
                                                nn, nl, size, rcount, charge, chargelo, chargehi, &
                                                coords, coordslo, coordshi, lo, hi, dx) &
                                                bind(c,name='amrex_compute_p3m_sr_correction_nl')

    integer,          intent(in   ) :: np, nn, size, chargelo(3), chargehi(3), coordslo(3), coordshi(3)
    integer,          intent(in   ) :: lo(3), hi(3)
    real(amrex_real), intent(in   ) :: dx(3)
    real(amrex_real), intent(inout) :: rcount
    type(particle_t), intent(inout) :: rparticles(np)
    type(particle_t), intent(inout) :: neighbors(nn)
    integer,          intent(in   ) :: nl(size)

    real(amrex_real), intent(in   ) :: charge(chargelo(1):chargehi(1),chargelo(2):chargehi(2),chargelo(3):chargehi(3))
    real(amrex_real), intent(in   ) :: coords(coordslo(1):coordshi(1),coordslo(2):coordshi(2),coordslo(3):coordshi(3),1:AMREX_SPACEDIM)

    real(amrex_real) :: dr(3), r2, r, coef, mass, correction_force_mag, ee, vals(70), points(70)
    real(amrex_real) :: dx2_inv, r_cell_frac,m, r_norm
    real(amrex_real) :: p3m_radius, im_charge_pos(3)
    integer :: i, j, index, nneighbors, store, ks, lookup_idx, k,  r_cell
    integer :: near_wall_below, near_wall_above
    integer :: near_wall_below_NL_part, near_wall_above_NL_part

    type(particle_t)                    :: particles(np+nn)

    double precision, allocatable :: weights(:,:,:,:)
    integer, allocatable :: indicies(:,:,:,:,:)

    ! initialize to 0 in case pkernel doesn't equal 6
    correction_force_mag=0.d0

    ! so that we only reference the p3m radius once (assuming it's the same forall particles)
    p3m_radius = particles(1)%p3m_radius

    ee = 1.d0/(permittivity*4_16*3.141592653589793238_16) 
    dx2_inv = 1.d0/(dx(1)*dx(1)) ! assumes isotropic cells


    particles(    1:np) = rparticles
    particles(np+1:   ) = neighbors

!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! p6 tables:
!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! this is force per dx**2
    vals =(/0., 0.0130335, 0.0259347, 0.0384891, 0.0507132, 0.0624416, 0.0735093, &
         0.0837511, 0.0931669, 0.101922, 0.109521, 0.116293, 0.121745, &
         0.126535, 0.130335, 0.132978, 0.134134, 0.135621, 0.13529, 0.134629, &
         0.133308, 0.130995, 0.129343, 0.126205, 0.122075, 0.118276, 0.115137, &
         0.110677, 0.106547, 0.102748, 0.0979574, 0.0936624, 0.0898631, &
         0.0860637, 0.0816036, 0.0782998, 0.0745005, 0.0707011, 0.0680581, &
         0.0645891, 0.0616157, 0.0586423, 0.0561644, 0.0535214, 0.0513739, &
         0.0487309, 0.0469138, 0.0450967, 0.0431145, 0.0414626, 0.0399759, &
         0.038324, 0.0366721, 0.0355157, 0.034029, 0.0328727, 0.0317164, &
         0.0307252, 0.0295689, 0.0285778, 0.0275866, 0.0267607, 0.0259347, &
         0.0251088, 0.0242829, 0.0236221, 0.0229613, 0.0221354, 0.0214746, &
         0.0209791/)
    ! these are fractions of cell size dx
    points =(/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, &
         1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, &
         2.8, 2.9, 3., 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4., 4.1, &
         4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5., 5.1, 5.2, 5.3, 5.4, 5.5, &
         5.6, 5.7, 5.8, 5.9, 6., 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9/)

    index = 1
    do i = 1, np
       !print*, "particle w charge ", int(particles(i)%q/abs(particles(i)%q)), " pos: " , particles(i)%pos
       near_wall_below = 0 ! reset for each particle
       near_wall_above = 0 ! reset for each particle

!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Peform correction in the presence of walls. 
       ! 
       ! Current implementation only works if a particle can be in the nbhd of one wall,
       ! in other words NO corners are allowed
       !
       ! Furthermore, it is assumed that the walls are in the x-z plane. 
!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! check if wall below 
       if (bc_es_lo(2) .ne. (-1)) then 
          call near_wall_check(particles(i), -2, near_wall_below)
       endif
       ! check if wall above
       if (bc_es_hi(2) .ne. (-1)) then 
          call near_wall_check(particles(i), 2, near_wall_above)
       endif
       ! print*, 'near wall below: ',   near_wall_below

       ! image charge interactions for wall below
       if (near_wall_below.eq.1) then
          ! compute image charge location
          call calc_im_charge_loc(particles(i), -2, im_charge_pos)
          !print*, 'im charge loc = ', im_charge_pos

          ! compute sep vector btwn charge and its image
          dr(1) = particles(i)%pos(1) - im_charge_pos(1)
          dr(2) = particles(i)%pos(2) - im_charge_pos(2)
          dr(3) = particles(i)%pos(3) - im_charge_pos(3)

          r2 = dot_product(dr,dr) 
          r = sqrt(r2)                    ! separation dist

          ! perform coulomb and p3m interaction with image charge
          if ((bc_es_lo(2) .eq. 1) .and. (r .lt. (particles(i)%p3m_radius))) then                      ! hom. dirichlet--image charge opposite that of particle

             ! coulomb
             particles(i)%force = particles(i)%force + ee*(dr/r)*particles(i)%q*(-1.d0*particles(i)%q)/r2
             !print*, 'coulomb interaction w self image: ', ee*(dr/r)*particles(i)%q*(-1.d0*particles(i)%q)/r2
             ! p3m 
             call compute_p3m_force_mag(r, correction_force_mag, dx)

             particles(i)%force = particles(i)%force - ee*particles(i)%q*(-1.d0*particles(i)%q)*(dr/r)*correction_force_mag*dx2_inv

          else if ((bc_es_lo(2) .eq. 2) .and. (r .lt. (particles(i)%p3m_radius))) then                 ! hom. neumann  --image charge equal that of particle

             ! coulomb
             particles(i)%force = particles(i)%force + ee*(dr/r)*particles(i)%q*(1.d0*particles(i)%q)/r2
             !print*, 'coulomb interaction w self image: ', ee*(dr/r)*particles(i)%q*(1.d0*particles(i)%q)/r2
             ! p3m 
             call compute_p3m_force_mag(r, correction_force_mag, dx)
             particles(i)%force = particles(i)%force - ee*particles(i)%q*(1.d0*particles(i)%q)*(dr/r)*correction_force_mag*dx2_inv

          endif
       endif

       ! image charge interactions for wall above
       if (near_wall_above.eq.1) then
          ! compute image charge location
          call calc_im_charge_loc(particles(i), 2, im_charge_pos)
          !print*, 'im charge loc = ', im_charge_pos
          ! compute sep vector btwn charge and its image
          dr(1) = particles(i)%pos(1) - im_charge_pos(1)
          dr(2) = particles(i)%pos(2) - im_charge_pos(2)
          dr(3) = particles(i)%pos(3) - im_charge_pos(3)

          r2 = dot_product(dr,dr) 
          r = sqrt(r2)                    ! separation dist


          ! perform coulomb and p3m interaction with image charge
          if ((bc_es_hi(2) .eq. 1) .and. (r .lt. (particles(i)%p3m_radius))) then                      ! hom. dirichlet--image charge opposite that of particle
             ! coulomb
             particles(i)%force = particles(i)%force + ee*(dr/r)*particles(i)%q*(-1.d0*particles(i)%q)/r2
             !print*, 'coulomb interaction w self image: ', ee*(dr/r)*particles(i)%q*(-1.d0*particles(i)%q)/r2
             ! p3m 
             call compute_p3m_force_mag(r, correction_force_mag, dx)
             particles(i)%force = particles(i)%force - ee*particles(i)%q*(-1.d0*particles(i)%q)*(dr/r)*correction_force_mag*dx2_inv

          else if ((bc_es_hi(2) .eq. 2) .and. (r .lt. (particles(i)%p3m_radius))) then                 ! hom. neumann  --image charge equal that of particle
             ! coulomb
             particles(i)%force = particles(i)%force + ee*(dr/r)*particles(i)%q*(1.d0*particles(i)%q)/r2
             !print*, 'coulomb interaction w self image: ', ee*(dr/r)*particles(i)%q*(1.d0*particles(i)%q)/r2
             ! p3m
             call compute_p3m_force_mag(r, correction_force_mag, dx)
             particles(i)%force = particles(i)%force - ee*particles(i)%q*(1.d0*particles(i)%q)*(dr/r)*correction_force_mag*dx2_inv
          endif

       endif

       nneighbors = nl(index)
       index = index + 1

       !print *, "particle ", i, " has ", nneighbors, " neighbours. Force ", particles(i)%force


       ! loop through neighbor list
       do j = index, index + nneighbors - 1
          near_wall_above_NL_part = 0 
          near_wall_below_NL_part = 0 

          dr(1) = particles(i)%pos(1) - particles(nl(j))%pos(1)
          dr(2) = particles(i)%pos(2) - particles(nl(j))%pos(2)
          dr(3) = particles(i)%pos(3) - particles(nl(j))%pos(3)

          r2 = dot_product(dr,dr) 
          r = sqrt(r2)                    ! separation dist


          if (r .lt. (particles(i)%p3m_radius)) then 

!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! Do local (short range) coulomb interaction within coulombRadiusFactor
!!!!!!!!!!!!!!!!!!!!!!!!!!
             particles(i)%force = particles(i)%force + ee*(dr/r)*particles(i)%q*particles(nl(j))%q/r2

             ! print*, 'Coulomb interction with NL part: ', ee*(dr/r)*particles(i)%q*particles(nl(j))%q/r2

             !print *, "particle ", i, " force ", particles(i)%force, particles(i)%pos, particles(nl(j))%pos

!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! Compute correction for fact that the above, sr coulomb interactions accounted for in poisson solve
             ! and hence are double counted 
             !
             ! Currently only implemented for pkernel=6
!!!!!!!!!!!!!!!!!!!!!!!!!!
                !print *, "calling with ", r, (particles(i)%p3m_radius)

             call compute_p3m_force_mag(r, correction_force_mag, dx)

             !print *, correction_force_mag, ee, particles(i)%q, dx2_inv

             ! force correction is negative: F_tot_electrostatic = F_sr_coulomb + F_poisson - F_correction
             particles(i)%force = particles(i)%force - ee*particles(i)%q*particles(nl(j))%q*(dr/r)*correction_force_mag*dx2_inv

!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! If part near wall, check if NL part also near wall
!!!!!!!!!!!!!!!!!!!!!!!!!!
             if (near_wall_below.eq.1) then
                call near_wall_check(particles(nl(j)), -2, near_wall_below_NL_part)
                if (near_wall_below_NL_part.eq.1) then 
                   ! compute NL particle's image charge location
                   call calc_im_charge_loc(particles(nl(j)), -2, im_charge_pos)
                   !print*, 'im charge loc = ', im_charge_pos

                   ! compute sep vector btwn charge and NL's image
                   dr(1) = particles(i)%pos(1) - im_charge_pos(1)
                   dr(2) = particles(i)%pos(2) - im_charge_pos(2)
                   dr(3) = particles(i)%pos(3) - im_charge_pos(3)

                   r2 = dot_product(dr,dr) 
                   r = sqrt(r2)                    ! separation dist

                   ! perform coulomb and p3m interaction with NL particle's image charge
                   if ((bc_es_lo(2) .eq. 1) .and. (r .lt. (particles(i)%p3m_radius))) then                      ! hom. dirichlet--image charge opposite that of particle

                      ! coulomb
                      !print*, 'Pre force: ', particles(i)%force

                      particles(i)%force = particles(i)%force + ee*(dr/r)*particles(i)%q*(-1.d0*particles(nl(j))%q)/r2
                      !print*, 'Coulomb: ', ee*(dr/r)*particles(i)%q*(-1.d0*particles(nl(j))%q)/r2

                      ! p3m
                      call compute_p3m_force_mag(r, correction_force_mag, dx)
                      particles(i)%force = particles(i)%force - ee*particles(i)%q*(-1.d0*particles(nl(j))%q)*(dr/r)*correction_force_mag*dx2_inv

                      !print *, "Below: ", particles(i)%id, correction_force_mag, particles(i)%force

                   else if ((bc_es_lo(2) .eq. 2) .and. (r .lt. (particles(i)%p3m_radius))) then                 ! hom. neumann  --image charge equal that of particle

                      ! coulomb
                      particles(i)%force = particles(i)%force + ee*(dr/r)*particles(i)%q*(1.d0*particles(nl(j))%q)/r2
                      !print*, 'Coulomb interaction w NL im part: ', ee*(dr/r)*particles(i)%q*(1.d0*particles(nl(j))%q)/r2

                      ! p3m 
                      call compute_p3m_force_mag(r, correction_force_mag, dx)
                      particles(i)%force = particles(i)%force - ee*particles(i)%q*(1.d0*particles(nl(j))%q)*(dr/r)*correction_force_mag*dx2_inv

                   endif

                endif
             endif

             if (near_wall_above.eq.1) then
                call near_wall_check(particles(nl(j)), 2, near_wall_above_NL_part)
                if (near_wall_above_NL_part.eq.1) then 
                   ! compute NL particle's image charge location
                   call calc_im_charge_loc(particles(nl(j)), 2, im_charge_pos)
                   !print*, 'im charge loc = ', im_charge_pos

                   ! compute sep vector btwn charge and NL's image
                   dr(1) = particles(i)%pos(1) - im_charge_pos(1)
                   dr(2) = particles(i)%pos(2) - im_charge_pos(2)
                   dr(3) = particles(i)%pos(3) - im_charge_pos(3)

                   r2 = dot_product(dr,dr) 
                   r = sqrt(r2)                    ! separation dist

                   ! perform coulomb and p3m interaction with NL particle's image charge
                   if ((bc_es_hi(2) .eq. 1) .and. (r .lt. (particles(i)%p3m_radius))) then                      ! hom. dirichlet--image charge opposite that of particle

                      ! coulomb
                      particles(i)%force = particles(i)%force + ee*(dr/r)*particles(i)%q*(-1.d0*particles(nl(j))%q)/r2
                      !print*, 'Coulomb interaction w NL im part: ', ee*(dr/r)*particles(i)%q*(-1.d0*particles(nl(j))%q)/r2
                      ! p3m 
                      call compute_p3m_force_mag(r, correction_force_mag, dx)
                      particles(i)%force = particles(i)%force - ee*particles(i)%q*(-1.d0*particles(nl(j))%q)*(dr/r)*correction_force_mag*dx2_inv

                   else if ((bc_es_hi(2) .eq. 2) .and. (r .lt. (particles(i)%p3m_radius))) then                 ! hom. neumann  --image charge equal that of particle
                      !print *, "Pre ", particles(i)%force
                      ! coulomb
                      particles(i)%force = particles(i)%force + ee*(dr/r)*particles(i)%q*(1.d0*particles(nl(j))%q)/r2
                      !print*, 'Coulomb interaction w NL im part: ', ee*(dr/r)*particles(i)%q*(1.d0*particles(nl(j))%q)/r2
                      ! p3m 
                      call compute_p3m_force_mag(r, correction_force_mag, dx)
                      particles(i)%force = particles(i)%force - ee*particles(i)%q*(1.d0*particles(nl(j))%q)*(dr/r)*correction_force_mag*dx2_inv

                      !print *, "Post ", correction_force_mag, particles(i)%force, dr, r

                   endif
                endif
             endif

             !print *, "Corr:", ee*particles(i)%q*particles(nl(j))%q*correction_force_mag*dx2_inv
             !print *, "norm:", correction_force_mag

          end if


       end do

       index = index + nneighbors

        !print *, particles(i)%id, particles(i)%force


    end do

    rparticles(:) = particles(1:np)
    neighbors(:)  = particles(np+1:)

  end subroutine amrex_compute_p3m_sr_correction_nl

  subroutine amrex_compute_forces_nl(rparticles, np, neighbors, & 
                                     nn, nl, size, rcount) &
                                     bind(c,name='amrex_compute_forces_nl')

    integer,          intent(in   ) :: np, nn, size
    real(amrex_real), intent(inout) :: rcount
    type(particle_t), intent(inout) :: rparticles(np)
    type(particle_t), intent(inout) :: neighbors(nn)
    integer,          intent(in   ) :: nl(size)

    real(amrex_real) :: dx(3), r2, r, coef, mass
    integer :: i, j, index, nneighbors

    type(particle_t)                    :: particles(np+nn)

    particles(    1:np) = rparticles
    particles(np+1:   ) = neighbors


    index = 1
    do i = 1, np

       !We need to check this properly

       nneighbors = nl(index)
       index = index + 1

       !print *, "particle ", i, " has ", nneighbors, " neighbours."

       do j = index, index + nneighbors - 1

          if((p_int_tog(particles(i)%species) .ne. 0) .and. (p_int_tog(particles(nl(j))%species) .ne. 0)) then

             dx(1) = particles(i)%pos(1) - particles(nl(j))%pos(1)
             dx(2) = particles(i)%pos(2) - particles(nl(j))%pos(2)
             dx(3) = particles(i)%pos(3) - particles(nl(j))%pos(3)

             r2 = dx(1) * dx(1) + dx(2) * dx(2) + dx(3) * dx(3)
             r2 = max(r2, rmin*rmin) 
             r = sqrt(r2)

             !repulsive interaction

             !print *, r , (1.122*particles(i)%sigma/2.0)
             if (r .lt. (1.122*particles(i)%sigma/2.0)) then ! NOTE! Should be able to set neighbor cell list with cut_off distance in mind

                !print *, "Repulsing, ", i, r, dx
                rcount = rcount + 1

                call repulsive_force(particles(i),particles(nl(j)),dx,r2) 
             end if
          endif

       end do

       index = index + nneighbors

    end do

    rparticles(:) = particles(1:np)
    neighbors(:)  = particles(np+1:)

  end subroutine amrex_compute_forces_nl

  subroutine move_particles_dsmc(particles, np, lo, hi, &
                                 cell_part_ids, cell_part_cnt, clo, chi, plo, phi, dx, dt, &
                                 paramplanes, ns, time, flux) &
                                 bind(c,name="move_particles_dsmc")

    type(particle_t), intent(inout), target :: particles(np)
    type(particle_t) :: toppart
    type(paramplane_t), intent(in), target :: paramplanes(ns)
    integer(c_int), intent(in) :: np, ns
    integer(c_int), intent(in) :: lo(3), hi(3)
    integer(c_int), intent(in) :: clo(3), chi(3)
    type(c_ptr), intent(inout) :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int), intent(inout) :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    real(amrex_real), intent(in) :: plo(3), phi(3), dx(3)
    real(amrex_real), intent(in) :: dt, time 
    integer(c_int), intent(inout) :: flux(2)

    integer :: i, j, k, p, cell_np, new_np, intsurf, intside, push, intcount, ii, fluxL, fluxR, count, numcoll
    integer :: cell(3)
    integer(c_int), pointer :: cell_parts(:)
    type(particle_t), pointer :: part
    type(paramplane_t), pointer :: surf
    real(amrex_real) inv_dx(3), runtime, inttime, adjalt, adj, inv_dt, domsize(3), posalt(3), prex, postx, radius, radius1, interval, omega, bessj0, dbessj0, bJ1, prefact, pi, t


    adj = 0.9999999
    adjalt = 2d0*(1d0 - adj)
    fluxL = 0
    fluxR = 0

    inv_dx = 1.d0/dx
    inv_dt = 1.d0/dt

    domsize = phi - plo

    do p = 1, ns

       surf => paramplanes(p)  

       surf%fxleft = 0
       surf%fyleft = 0
       surf%fzleft = 0

       surf%fxright = 0
       surf%fyright = 0
       surf%fzright = 0 

    enddo

    intcount = 0

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             cell_np = cell_part_cnt(i,j,k)
             call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])

             new_np = cell_np
             p = 1

             do while (p <= new_np)

                !if(cell_parts(p) .eq. 1320) then

                !print *, "Accessing ", p, " which is ", cell_parts(p), " and has id ", i,j,k

                !  endif
                part => particles(cell_parts(p))

                ! print*, part%id, part%vel(3)

                runtime = dt


                do while (runtime .gt. 0)

                   call find_intersect(part,runtime, paramplanes, ns, intsurf, inttime, intside, phi, plo)

                   !print *, runtime, inttime

                   posalt(1) = inttime*part%vel(1)*adjalt
                   posalt(2) = inttime*part%vel(2)*adjalt
#if (BL_SPACEDIM == 3)
                   posalt(3) = inttime*part%vel(3)*adjalt
#endif

                   ! if(intsurf .eq.  5 .or. intsurf .eq.  6) then
                   !    write(*,*) "YAY"

                   ! endif

                   if(intsurf .eq. 7) then
                      if(part%vel(1) .lt. 0) then
                         fluxL = fluxL + 1
                      else
                         fluxR = fluxR + 1
                      endif
                   endif




                   ! move the particle in a straight line, adj factor prevents double detection of boundary intersection
                   part%pos(1) = part%pos(1) + inttime*part%vel(1)*adj
                   part%pos(2) = part%pos(2) + inttime*part%vel(2)*adj
#if (BL_SPACEDIM == 3)
                   part%pos(3) = part%pos(3) + inttime*part%vel(3)*adj
#endif
                   runtime = runtime - inttime



                   !  if(intsurf .eq.  5) then
                   !     write(*,*) "5", part%pos(1), part%pos(2), part%pos(3)
                   !     count5=count5+1
                   !  elseif(intsurf .eq. 6) then
                   !     write(*,*)  "6", part%pos(1), part%pos(2), part%pos(3)
                   !    count6= count6+1
                   ! endif
                   if(intsurf .gt. 0) then

                      surf => paramplanes(intsurf)


                      call apply_bc(surf, part, intside, domsize, push, time, inttime)                      

                      if(push .eq. 1) then

                         part%pos(1) = part%pos(1) + posalt(1)
                         part%pos(2) = part%pos(2) + posalt(2)
#if (BL_SPACEDIM == 3)
                         part%pos(3) = part%pos(3) + posalt(3)
#endif 
                      endif

                   endif

                end do



                ! if it has changed cells, remove from vector.
                ! otherwise continue

                cell(1) = floor((part%pos(1) - plo(1))*inv_dx(1))              
                cell(2) = floor((part%pos(2) - plo(2))*inv_dx(2))
#if (BL_SPACEDIM == 3)
                cell(3) = floor((part%pos(3) - plo(3))*inv_dx(3))
#else
                cell(3) = 0
#endif
                ! if(intsurf .eq.  5) then
                !      write(*,*) "5", part%pos(1), part%pos(2), part%pos(3)
                !      count5=count5+1
                !   elseif(intsurf .eq. 6) then
                !      write(*,*)  "6", part%pos(1), part%pos(2), part%pos(3)
                !     count6= count6+1
                !  endif
                if ((cell(1) /= i) .or. (cell(2) /= j) .or. (cell(3) /= k)) then
                   part%sorted = 0

                   call remove_particle_from_cell(cell_parts, cell_np, new_np, p)

                else
                   p = p + 1
                end if
             end do
             cell_part_cnt(i,j,k) = new_np           
          end do

       end do
    end do

    !print *, "intcount: ", intcount

    do p = 1, ns

       surf => paramplanes(p)  

       surf%fxleft = surf%fxleft*inv_dt
       surf%fyleft = surf%fyleft*inv_dt
       surf%fzleft = surf%fzleft*inv_dt

       surf%fxright = surf%fxright*inv_dt
       surf%fyright = surf%fyright*inv_dt
       surf%fzright = surf%fzright*inv_dt

    enddo

    flux(1) = fluxL
    flux(2) = fluxR

    if(graphene_tog .eq. 1) then

       surf=>paramplanes(6)

       pi=3.1415926535897932
       numcoll=floor(pi*(prob_hi(1)**2)*fixed_dt*(particle_n0(1)/particle_neff)*sqrt((k_b*t_init(1))/(2*pi*mass(1))))
       !print *, "topnum: ", numcoll, prob_hi(1), fixed_dt

       do count=1,  numcoll
          call topparticle(surf, time, inttime)
       end do

       call laser(surf, time)
       surf%a0graph=surf%agraph
       surf%b0graph=surf%bgraph

       interval=prob_hi(1)/100
       radius=0
       bJ1 = bessel_jn(1,2.4048)
       prefact = 9144**2/(prob_hi(1)*prob_hi(1)*3.14159*bJ1**2)
       omega=12.5*(10**6)*2*3.1415926535897932
       surf=>paramplanes(6)
       do ii=1, 1
          radius=interval*ii
          radius=radius*2.4048/prob_hi(1)
          bessj0 =-prefact*bessel_jn(0, radius)*(surf%a0graph*sin(omega*time)+surf%b0graph*cos(time*omega))/omega
          !surf%besslist(ii)=bessj0
          dbessj0=prefact*bessel_jn(0, radius)*(surf%a0graph*sin(omega*time)+surf%b0graph*cos(time*omega))
          surf%dbesslist(ii)=bessj0
       enddo

       ! print*,'position',part%pos
       ! print*,'vel',part%vel
       ! print*,'fortran move', surf%velz, part%id
       ! print*, 'a', surf%a0graph
       ! print*, sin((time*omega))
       ! print*, surf%velz
       ! print*, 'velocity',part%vel
       ! print*, "hack"
    endif

  end subroutine move_particles_dsmc

  ! subroutine redirect(part)

  !   type(particle_t), intent(inout) :: part

  !   double precision speed

  !   speed = sqrt(part%vel(1)**2 + part%vel(2)**2 + part%vel(3)**2)
    
  !   if(speed .ne. 0) then
  !      part%dir = part%vel/speed
  !   endif

  ! end subroutine redirect
  
  subroutine get_interpolation_weights(cc, rr, ixf, onemdxf)

    double precision, intent(in)  :: ixf(3), onemdxf(3)

    double precision, intent(inout)  :: cc(0:7)
    double precision, intent(inout)  :: rr(0:7)

#if (BL_SPACEDIM == 3)
    cc(0) = onemdxf(1)*onemdxf(2)*onemdxf(3)
    cc(4) = onemdxf(1)*onemdxf(2)*ixf(3)
    cc(2) = onemdxf(1)*onemdxf(3)*ixf(2)
    cc(6) = onemdxf(1)*ixf(2)*ixf(3)
    cc(5) = onemdxf(2)*onemdxf(3)*ixf(1)
    cc(1) = onemdxf(2)*ixf(1)*ixf(3)
    cc(3) = onemdxf(3)*ixf(1)*ixf(2)
    cc(7) = ixf(1)*ixf(2)*ixf(3)

    rr = cc/sum(cc)
#endif

#if (BL_SPACEDIM == 2)
    cc(0) = onemdxf(1)*onemdxf(2)
    cc(2) = onemdxf(1)*ixf(2)
    cc(1) = ixf(1)*onemdxf(2)
    cc(3) = ixf(1)*ixf(2)

    rr = cc/(sum(cc)*2)
#endif

  end subroutine get_interpolation_weights

  subroutine get_local_properties(cc, rr, fi, velx, velxlo, velxhi, vely, velylo, velyhi, &
#if (BL_SPACEDIM == 3)     
                                  velz, velzlo, velzhi, &
#endif
                                  beta, localvel, localbeta, betalo, betahi)
    
    integer,          intent(in   ) :: fi(3), velxlo(3), velxhi(3), velylo(3), velyhi(3), betalo(3), betahi(3)
#if (AMREX_SPACEDIM == 3)
    integer,          intent(in   ) :: velzlo(3), velzhi(3)
#endif
    double precision, intent(inout) :: rr(0:7), localvel(3), localbeta, cc(0:7)

    double precision, intent(in   ) :: velx(velxlo(1):velxhi(1),velxlo(2):velxhi(2),velxlo(3):velxhi(3))
    double precision, intent(in   ) :: vely(velylo(1):velyhi(1),velylo(2):velyhi(2),velylo(3):velyhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: velz(velzlo(1):velzhi(1),velzlo(2):velzhi(2),velzlo(3):velzhi(3))
#endif

    double precision, intent(in   ) :: beta(betalo(1):betahi(1),betalo(2):betahi(2),betalo(3):betahi(3))

#if (BL_SPACEDIM == 3)

    if (visc_type .gt. 0) then
       localbeta = beta(fi(1),fi(2),fi(3))
    else
       !3d visc
       localbeta = rr(0)*beta(fi(1),fi(2),fi(3)) + rr(4)*beta(fi(1),fi(2),fi(3)+1) + rr(2)*beta(fi(1),fi(2)+1,fi(3)) + rr(6)*beta(fi(1),fi(2)+1,fi(3)+1) + rr(1)*beta(fi(1)+1,fi(2),fi(3)) + rr(5)*beta(fi(1)+1,fi(2),fi(3)+1) + rr(3)*beta(fi(1)+1,fi(2)+1,fi(3)) + rr(7)*beta(fi(1)+1,fi(2)+1,fi(3)+1)
    endif

    localvel(1) = rr(0)*velx(fi(1),fi(2),fi(3)) + rr(4)*velx(fi(1),fi(2),fi(3)+1) + rr(2)*velx(fi(1),fi(2)+1,fi(3)) + rr(6)*velx(fi(1),fi(2)+1,fi(3)+1) + rr(1)*velx(fi(1)+1,fi(2),fi(3)) + rr(5)*velx(fi(1)+1,fi(2),fi(3)+1) + rr(3)*velx(fi(1)+1,fi(2)+1,fi(3)) + rr(7)*velx(fi(1)+1,fi(2)+1,fi(3)+1)
    localvel(2) = rr(0)*vely(fi(1),fi(2),fi(3)) + rr(4)*vely(fi(1),fi(2),fi(3)+1) + rr(2)*vely(fi(1),fi(2)+1,fi(3)) + rr(6)*vely(fi(1),fi(2)+1,fi(3)+1) + rr(1)*vely(fi(1)+1,fi(2),fi(3)) + rr(5)*vely(fi(1)+1,fi(2),fi(3)+1) + rr(3)*vely(fi(1)+1,fi(2)+1,fi(3)) + rr(7)*vely(fi(1)+1,fi(2)+1,fi(3)+1)
    localvel(3) = rr(0)*velz(fi(1),fi(2),fi(3)) + rr(4)*velz(fi(1),fi(2),fi(3)+1) + rr(2)*velz(fi(1),fi(2)+1,fi(3)) + rr(6)*velz(fi(1),fi(2)+1,fi(3)+1) + rr(1)*velz(fi(1)+1,fi(2),fi(3)) + rr(5)*velz(fi(1)+1,fi(2),fi(3)+1) + rr(3)*velz(fi(1)+1,fi(2)+1,fi(3)) + rr(7)*velz(fi(1)+1,fi(2)+1,fi(3)+1)
#endif

#if (BL_SPACEDIM == 2)

    !CHECK THIS!

    if (visc_type .gt. 0) then
       localbeta = beta(fi(1),fi(2),fi(3))
    else
       localbeta = beta(fi(1),fi(2),fi(3))*rr(0) + beta(fi(1),fi(2)+1,fi(3))*rr(2) + beta(fi(1)+1,fi(2),fi(3))*rr(1) + beta(fi(1)+1,fi(2)+1,fi(3))*rr(3)
    endif
    !2d xvel
    localvel(1) = velx(fi(1),fi(2),fi(3))*rr(0) + velx(fi(1),fi(2)+1,fi(3))*rr(2) + velx(fi(1)+1,fi(2),fi(3))*rr(1) + velx(fi(1)+1,fi(2)+1,fi(3))*rr(3)
    localvel(2) = vely(fi(1),fi(2),fi(3))*rr(0) + vely(fi(1),fi(2)+1,fi(3))*rr(2) + vely(fi(1)+1,fi(2),fi(3))*rr(1) + vely(fi(1)+1,fi(2)+1,fi(3))*rr(3)
#endif

  end subroutine get_local_properties

  subroutine distribute_momentum(deltap, rr, fi ,sourcex, sourcexlo, sourcexhi, sourcey, sourceylo, sourceyhi &
#if (BL_SPACEDIM == 3)
                                 , sourcez, sourcezlo, sourcezhi &
#endif
                                 )
    
    integer,          intent(in   ) :: fi(3), sourcexlo(3), sourcexhi(3), sourceylo(3), sourceyhi(3)
#if (AMREX_SPACEDIM == 3)
    integer,          intent(in   ) :: sourcezlo(3), sourcezhi(3)
#endif

    double precision, intent(in   ) :: rr(0:7)

    double precision, intent(inout) :: deltap(3)

    double precision, intent(inout) :: sourcex(sourcexlo(1):sourcexhi(1),sourcexlo(2):sourcexhi(2),sourcexlo(3):sourcexhi(3))
    double precision, intent(inout) :: sourcey(sourceylo(1):sourceyhi(1),sourceylo(2):sourceyhi(2),sourceylo(3):sourceyhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(inout) :: sourcez(sourcezlo(1):sourcezhi(1),sourcezlo(2):sourcezhi(2),sourcezlo(3):sourcezhi(3))
#endif

    double precision nodalp

#if (BL_SPACEDIM == 3)
    !distribute x momentum change 
    nodalp = rr(0)*deltap(1)
    sourcex(fi(1),fi(2),fi(3)) = sourcex(fi(1),fi(2),fi(3)) + nodalp
    sourcex(fi(1),fi(2)-1,fi(3)) = sourcex(fi(1),fi(2)-1,fi(3)) + nodalp
    sourcex(fi(1),fi(2),fi(3)-1) = sourcex(fi(1),fi(2),fi(3)-1) + nodalp
    sourcex(fi(1),fi(2)-1,fi(3)-1) = sourcex(fi(1),fi(2)-1,fi(3)-1) + nodalp

    nodalp = rr(4)*deltap(1)
    sourcex(fi(1),fi(2),fi(3)+1) = sourcex(fi(1),fi(2),fi(3)+1) + nodalp
    sourcex(fi(1),fi(2)-1,fi(3)+1) = sourcex(fi(1),fi(2)-1,fi(3)+1) + nodalp
    sourcex(fi(1),fi(2),fi(3)) = sourcex(fi(1),fi(2),fi(3)) + nodalp
    sourcex(fi(1),fi(2)-1,fi(3)) = sourcex(fi(1),fi(2)-1,fi(3)) + nodalp

    nodalp = rr(2)*deltap(1)
    sourcex(fi(1),fi(2)+1,fi(3)) = sourcex(fi(1),fi(2)+1,fi(3)) + nodalp
    sourcex(fi(1),fi(2),fi(3)) = sourcex(fi(1),fi(2),fi(3)) + nodalp
    sourcex(fi(1),fi(2)+1,fi(3)-1) = sourcex(fi(1),fi(2)+1,fi(3)-1) + nodalp
    sourcex(fi(1),fi(2),fi(3)-1) = sourcex(fi(1),fi(2),fi(3)-1) + nodalp

    nodalp = rr(6)*deltap(1)
    sourcex(fi(1),fi(2)+1,fi(3)+1) = sourcex(fi(1),fi(2)+1,fi(3)+1) + nodalp
    sourcex(fi(1),fi(2),fi(3)+1) = sourcex(fi(1),fi(2),fi(3)+1) + nodalp
    sourcex(fi(1),fi(2)+1,fi(3)) = sourcex(fi(1),fi(2)+1,fi(3)) + nodalp
    sourcex(fi(1),fi(2),fi(3)) = sourcex(fi(1),fi(2),fi(3)) + nodalp

    nodalp = rr(1)*deltap(1)
    sourcex(fi(1)+1,fi(2),fi(3)) = sourcex(fi(1)+1,fi(2),fi(3)) + nodalp
    sourcex(fi(1)+1,fi(2)-1,fi(3)) = sourcex(fi(1)+1,fi(2)-1,fi(3)) + nodalp
    sourcex(fi(1)+1,fi(2),fi(3)-1) = sourcex(fi(1)+1,fi(2),fi(3)-1) + nodalp
    sourcex(fi(1)+1,fi(2)-1,fi(3)-1) = sourcex(fi(1)+1,fi(2)-1,fi(3)-1) + nodalp

    nodalp = rr(5)*deltap(1)
    sourcex(fi(1)+1,fi(2),fi(3)+1) = sourcex(fi(1)+1,fi(2),fi(3)+1) + nodalp
    sourcex(fi(1)+1,fi(2)-1,fi(3)+1) = sourcex(fi(1)+1,fi(2)-1,fi(3)+1) + nodalp
    sourcex(fi(1)+1,fi(2),fi(3)) = sourcex(fi(1)+1,fi(2),fi(3)) + nodalp
    sourcex(fi(1)+1,fi(2)-1,fi(3)) = sourcex(fi(1)+1,fi(2)-1,fi(3)) + nodalp

    nodalp = rr(3)*deltap(1)
    sourcex(fi(1)+1,fi(2)+1,fi(3)) = sourcex(fi(1)+1,fi(2)+1,fi(3)) + nodalp
    sourcex(fi(1)+1,fi(2),fi(3)) = sourcex(fi(1)+1,fi(2),fi(3)) + nodalp
    sourcex(fi(1)+1,fi(2)+1,fi(3)-1) = sourcex(fi(1)+1,fi(2)+1,fi(3)-1) + nodalp
    sourcex(fi(1)+1,fi(2),fi(3)-1) = sourcex(fi(1)+1,fi(2),fi(3)-1) + nodalp

    nodalp = rr(7)*deltap(1)
    sourcex(fi(1)+1,fi(2)+1,fi(3)+1) = sourcex(fi(1)+1,fi(2)+1,fi(3)+1) + nodalp
    sourcex(fi(1)+1,fi(2),fi(3)+1) = sourcex(fi(1)+1,fi(2),fi(3)+1) + nodalp
    sourcex(fi(1)+1,fi(2)+1,fi(3)) = sourcex(fi(1)+1,fi(2)+1,fi(3)) + nodalp
    sourcex(fi(1)+1,fi(2),fi(3)) = sourcex(fi(1)+1,fi(2),fi(3)) + nodalp

    !distribute y momentum change 
    nodalp = rr(0)*deltap(2)
    sourcey(fi(1),fi(2),fi(3)) = sourcey(fi(1),fi(2),fi(3)) + nodalp
    sourcey(fi(1)-1,fi(2),fi(3)) = sourcey(fi(1)-1,fi(2),fi(3)) + nodalp
    sourcey(fi(1),fi(2),fi(3)-1) = sourcey(fi(1),fi(2),fi(3)-1) + nodalp
    sourcey(fi(1)-1,fi(2),fi(3)-1) = sourcey(fi(1)-1,fi(2),fi(3)-1) + nodalp

    nodalp = rr(4)*deltap(2)
    sourcey(fi(1),fi(2),fi(3)+1) = sourcey(fi(1),fi(2),fi(3)+1) + nodalp
    sourcey(fi(1)-1,fi(2),fi(3)+1) = sourcey(fi(1)-1,fi(2),fi(3)+1) + nodalp
    sourcey(fi(1),fi(2),fi(3)) = sourcey(fi(1),fi(2),fi(3)) + nodalp
    sourcey(fi(1)-1,fi(2),fi(3)) = sourcey(fi(1)-1,fi(2),fi(3)) + nodalp

    nodalp = rr(2)*deltap(2)
    sourcey(fi(1),fi(2)+1,fi(3)) = sourcey(fi(1),fi(2)+1,fi(3)) + nodalp
    sourcey(fi(1)-1,fi(2)+1,fi(3)) = sourcey(fi(1)-1,fi(2)+1,fi(3)) + nodalp
    sourcey(fi(1),fi(2)+1,fi(3)-1) = sourcey(fi(1),fi(2)+1,fi(3)-1) + nodalp
    sourcey(fi(1)-1,fi(2)+1,fi(3)-1) = sourcey(fi(1)-1,fi(2)+1,fi(3)-1) + nodalp

    nodalp = rr(6)*deltap(2)
    sourcey(fi(1),fi(2)+1,fi(3)+1) = sourcey(fi(1),fi(2)+1,fi(3)+1) + nodalp
    sourcey(fi(1)-1,fi(2)+1,fi(3)+1) = sourcey(fi(1)-1,fi(2)+1,fi(3)+1) + nodalp
    sourcey(fi(1),fi(2)+1,fi(3)) = sourcey(fi(1),fi(2)+1,fi(3)) + nodalp
    sourcey(fi(1)-1,fi(2)+1,fi(3)) = sourcey(fi(1)-1,fi(2)+1,fi(3)) + nodalp

    nodalp = rr(1)*deltap(2)
    sourcey(fi(1)+1,fi(2),fi(3)) = sourcey(fi(1)+1,fi(2),fi(3)) + nodalp
    sourcey(fi(1),fi(2),fi(3)) = sourcey(fi(1),fi(2),fi(3)) + nodalp
    sourcey(fi(1)+1,fi(2),fi(3)-1) = sourcey(fi(1)+1,fi(2),fi(3)-1) + nodalp
    sourcey(fi(1),fi(2),fi(3)-1) = sourcey(fi(1),fi(2),fi(3)-1) + nodalp

    nodalp = rr(5)*deltap(2)
    sourcey(fi(1)+1,fi(2),fi(3)+1) = sourcey(fi(1)+1,fi(2),fi(3)+1) + nodalp
    sourcey(fi(1),fi(2),fi(3)+1) = sourcey(fi(1),fi(2),fi(3)+1) + nodalp
    sourcey(fi(1)+1,fi(2),fi(3)) = sourcey(fi(1)+1,fi(2),fi(3)) + nodalp
    sourcey(fi(1),fi(2),fi(3)) = sourcey(fi(1),fi(2),fi(3)) + nodalp

    nodalp = rr(3)*deltap(2)
    sourcey(fi(1)+1,fi(2)+1,fi(3)) = sourcey(fi(1)+1,fi(2)+1,fi(3)) + nodalp
    sourcey(fi(1),fi(2)+1,fi(3)) = sourcey(fi(1),fi(2)+1,fi(3)) + nodalp
    sourcey(fi(1)+1,fi(2)+1,fi(3)-1) = sourcey(fi(1)+1,fi(2)+1,fi(3)-1) + nodalp
    sourcey(fi(1),fi(2)+1,fi(3)-1) = sourcey(fi(1),fi(2)+1,fi(3)-1) + nodalp

    nodalp = rr(7)*deltap(2)
    sourcey(fi(1)+1,fi(2)+1,fi(3)+1) = sourcey(fi(1)+1,fi(2)+1,fi(3)+1) + nodalp
    sourcey(fi(1),fi(2)+1,fi(3)+1) = sourcey(fi(1),fi(2)+1,fi(3)+1) + nodalp
    sourcey(fi(1)+1,fi(2)+1,fi(3)) = sourcey(fi(1)+1,fi(2)+1,fi(3)) + nodalp
    sourcey(fi(1),fi(2)+1,fi(3)) = sourcey(fi(1),fi(2)+1,fi(3)) + nodalp

    !distribute z momentum change 
    nodalp = rr(0)*deltap(3)
    sourcez(fi(1),fi(2),fi(3)) = sourcez(fi(1),fi(2),fi(3)) + nodalp
    sourcez(fi(1)-1,fi(2),fi(3)) = sourcez(fi(1)-1,fi(2),fi(3)) + nodalp
    sourcez(fi(1),fi(2)-1,fi(3)) = sourcez(fi(1),fi(2)-1,fi(3)) + nodalp
    sourcez(fi(1)-1,fi(2)-1,fi(3)) = sourcez(fi(1)-1,fi(2)-1,fi(3)) + nodalp

    nodalp = rr(4)*deltap(3)
    sourcez(fi(1),fi(2),fi(3)+1) = sourcez(fi(1),fi(2),fi(3)+1) + nodalp
    sourcez(fi(1)-1,fi(2),fi(3)+1) = sourcez(fi(1)-1,fi(2),fi(3)+1) + nodalp
    sourcez(fi(1),fi(2)-1,fi(3)+1) = sourcez(fi(1),fi(2)-1,fi(3)+1) + nodalp
    sourcez(fi(1)-1,fi(2)-1,fi(3)+1) = sourcez(fi(1)-1,fi(2)-1,fi(3)+1) + nodalp

    nodalp = rr(2)*deltap(3)
    sourcez(fi(1),fi(2)+1,fi(3)) = sourcez(fi(1),fi(2)+1,fi(3)) + nodalp
    sourcez(fi(1)-1,fi(2)+1,fi(3)) = sourcez(fi(1)-1,fi(2)+1,fi(3)) + nodalp
    sourcez(fi(1),fi(2),fi(3)) = sourcez(fi(1),fi(2),fi(3)) + nodalp
    sourcez(fi(1)-1,fi(2),fi(3)) = sourceZ(fi(1)-1,fi(2),fi(3)) + nodalp

    nodalp = rr(6)*deltap(3)
    sourcez(fi(1),fi(2)+1,fi(3)+1) = sourcez(fi(1),fi(2)+1,fi(3)+1) + nodalp
    sourcez(fi(1)-1,fi(2)+1,fi(3)+1) = sourcez(fi(1)-1,fi(2)+1,fi(3)+1) + nodalp
    sourcez(fi(1),fi(2),fi(3)+1) = sourcez(fi(1),fi(2),fi(3)+1) + nodalp
    sourcez(fi(1)-1,fi(2),fi(3)+1) = sourcez(fi(1)-1,fi(2),fi(3)+1) + nodalp

    nodalp = rr(1)*deltap(3)
    sourcez(fi(1)+1,fi(2),fi(3)) = sourcez(fi(1)+1,fi(2),fi(3)) + nodalp
    sourcez(fi(1),fi(2),fi(3)) = sourcez(fi(1),fi(2),fi(3)) + nodalp
    sourcez(fi(1)+1,fi(2)-1,fi(3)) = sourcez(fi(1)+1,fi(2)-1,fi(3)) + nodalp
    sourcez(fi(1),fi(2)-1,fi(3)) = sourcez(fi(1),fi(2)-1,fi(3)) + nodalp

    nodalp = rr(5)*deltap(3)
    sourcez(fi(1)+1,fi(2),fi(3)+1) = sourcez(fi(1)+1,fi(2),fi(3)+1) + nodalp
    sourcez(fi(1),fi(2),fi(3)+1) = sourcez(fi(1),fi(2),fi(3)+1) + nodalp
    sourcez(fi(1)+1,fi(2)-1,fi(3)+1) = sourcez(fi(1)+1,fi(2)-1,fi(3)+1) + nodalp
    sourcez(fi(1),fi(2)-1,fi(3)+1) = sourcez(fi(1),fi(2)-1,fi(3)+1) + nodalp

    nodalp = rr(3)*deltap(3)
    sourcez(fi(1)+1,fi(2)+1,fi(3)) = sourcez(fi(1)+1,fi(2)+1,fi(3)) + nodalp
    sourcez(fi(1),fi(2)+1,fi(3)) = sourcez(fi(1),fi(2)+1,fi(3)) + nodalp
    sourcez(fi(1)+1,fi(2),fi(3)) = sourcez(fi(1)+1,fi(2),fi(3)) + nodalp
    sourcez(fi(1),fi(2),fi(3)) = sourcez(fi(1),fi(2),fi(3)) + nodalp

    nodalp = rr(7)*deltap(3)
    sourcez(fi(1)+1,fi(2)+1,fi(3)+1) = sourcez(fi(1)+1,fi(2)+1,fi(3)+1) + nodalp
    sourcez(fi(1),fi(2)+1,fi(3)+1) = sourcez(fi(1),fi(2)+1,fi(3)+1) + nodalp
    sourcez(fi(1)+1,fi(2),fi(3)+1) = sourcez(fi(1)+1,fi(2),fi(3)+1) + nodalp
    sourcez(fi(1),fi(2),fi(3)+1) = sourcez(fi(1),fi(2),fi(3)+1) + nodalp
#endif

#if (BL_SPACEDIM == 2)
    !distribute x momentum change 
    nodalp = rr(0)*deltap(1)
    sourcex(fi(1),fi(2),fi(3)) = sourcex(fi(1),fi(2),fi(3)) + nodalp
    sourcex(fi(1),fi(2)-1,fi(3)) = sourcex(fi(1),fi(2)-1,fi(3)) + nodalp

    nodalp = rr(2)*deltap(1)
    sourcex(fi(1),fi(2)+1,fi(3)) = sourcex(fi(1),fi(2)+1,fi(3)) + nodalp
    sourcex(fi(1),fi(2),fi(3)) = sourcex(fi(1),fi(2),fi(3)) + nodalp

    nodalp = rr(1)*deltap(1)
    sourcex(fi(1)+1,fi(2),fi(3)) = sourcex(fi(1)+1,fi(2),fi(3)) + nodalp
    sourcex(fi(1)+1,fi(2)-1,fi(3)) = sourcex(fi(1)+1,fi(2)-1,fi(3)) + nodalp

    nodalp = rr(3)*deltap(1)
    sourcex(fi(1)+1,fi(2)+1,fi(3)) = sourcex(fi(1)+1,fi(2)+1,fi(3)) + nodalp
    sourcex(fi(1)+1,fi(2),fi(3)) = sourcex(fi(1)+1,fi(2),fi(3)) + nodalp

    !distribute y momentum change
    nodalp = rr(0)*deltap(2)
    sourcey(fi(1),fi(2),fi(3)) = sourcey(fi(1),fi(2),fi(3)) + nodalp
    sourcey(fi(1)-1,fi(2),fi(3)) = sourcey(fi(1)-1,fi(2),fi(3)) + nodalp

    nodalp = rr(2)*deltap(2)
    sourcey(fi(1),fi(2)+1,fi(3)) = sourcey(fi(1),fi(2)+1,fi(3)) + nodalp
    sourcey(fi(1)-1,fi(2)+1,fi(3)) = sourcey(fi(1)-1,fi(2)+1,fi(3)) + nodalp

    nodalp = rr(1)*deltap(2)
    sourcey(fi(1)+1,fi(2),fi(3)) = sourcey(fi(1)+1,fi(2),fi(3)) + nodalp
    sourcey(fi(1),fi(2),fi(3)) = sourcey(fi(1),fi(2),fi(3)) + nodalp

    nodalp = rr(3)*deltap(2)
    sourcey(fi(1)+1,fi(2)+1,fi(3)) = sourcey(fi(1)+1,fi(2)+1,fi(3)) + nodalp
    sourcey(fi(1),fi(2)+1,fi(3)) = sourcey(fi(1),fi(2)+1,fi(3)) + nodalp
#endif

  end subroutine distribute_momentum


 ! extra diffusion term when 
  subroutine dry(dt,part,dry_terms, mb) &
                           bind(c,name="dry")

    type(particle_t), intent(inout) :: part 
    double precision, intent(inout) :: dry_terms(3)
    double precision, intent(in   ) :: dt, mb(3)
    real(amrex_real) runtime, normalrand(3),std(3),bfac(3)

    ! Brownian forcing
    call get_particle_normal(normalrand(1))
    call get_particle_normal(normalrand(2))
    call get_particle_normal(normalrand(3))

!    normalrand(1) = 0.5d0
!    normalrand(2) = 0.5d0
!    normalrand(3) = 0.5d0

    !std = sqrt(part%dry_diff*k_B*2d0*t_init(1))
    std(1) = sqrt(2.0*mb(1)*part%dry_diff)
    std(2) = sqrt(2.0*mb(2)*part%dry_diff)
    std(3) = sqrt(2.0*mb(3)*part%dry_diff)

    !DRL: dry diffusion coef: part%dry_diff, temperature: t_init(1)

    bfac(1) = variance_coef_mom*std(1)*normalrand(1)/sqrt(dt)
    bfac(2) = variance_coef_mom*std(2)*normalrand(2)/sqrt(dt)
    bfac(3) = variance_coef_mom*std(3)*normalrand(3)/sqrt(dt)

    !KK does this have all the forces in it already? need to check
    dry_terms(1) = mb(1)*part%dry_diff*part%force(1)/(k_B*t_init(1))+bfac(1)
    dry_terms(2) = mb(2)*part%dry_diff*part%force(2)/(k_B*t_init(1))+bfac(2)
    dry_terms(3) = mb(3)*part%dry_diff*part%force(3)/(k_B*t_init(1))+bfac(3)

  end subroutine dry
  
  subroutine peskin_3pt(r,w)

    double precision, intent(in   ) :: r
    double precision, intent(inout) :: w

    double precision r2, r1

    r1 = abs(r)
    r2 = r1*r1

    if(r1 .le. 0.5) then

       w = (1 + sqrt(1-3*r2))/3

    elseif(r1 .le. 1.5) then

       w = (5 - 3*r1 - sqrt(-3*(1-r1)*(1-r1)+1)) / 6.

    else

       w = 0

    endif


  end subroutine peskin_3pt


  subroutine peskin_4pt(r,w)

    double precision, intent(in   ) :: r
    double precision, intent(inout) :: w

    double precision rr 
    rr = r*r

    if(r .le. -2) then

       w = 0

    elseif(r .le. -1) then

       w = 0.125*(5 + 2*r - sqrt(-7 - 12*r - 4*rr))

    elseif(r .le. 0) then

       w = 0.125*(3 + 2*r + sqrt(1 - 4*r - 4*rr))

    elseif(r .le. 1) then

       w = 0.125*(3 - 2*r + sqrt(1 + 4*r - 4*rr))

    elseif(r .le. 2) then

       w = 0.125*(5 - 2*r - sqrt(-7 + 12*r - 4*rr))

    else

       w = 0

    endif


  end subroutine peskin_4pt

  subroutine bspline_6pt(r,w)

    !Doesn't work?

    double precision, intent(in   ) :: r
    double precision, intent(inout) :: w

    double precision r1, r2, r3, r4, r5

    r1 = abs(r)
    r2 = r1*r1
    r3 = r2*r1
    r4 = r3*r1
    r5 = r4*r1

    if(r1 .le. 1) then

       w = 0.55 - 0.5*r2 + 0.25*r4 - (1d0/12d0)*r5

    elseif(r1 .le. 2) then

       w = 0.425 + 0.625*r - 1.75*r2 + 1.25*r3 - 0.375*r4 + (1d0/24d0)*r5

    elseif(r1 .le. 3) then

       w = 2.025 - 3.375*r + 2.25*r2 - 0.75*r3 + 0.125*r4 - (1d0/120d0)*r5;

    else

       w = 0

    endif


  end subroutine bspline_6pt

  subroutine peskin_6pt(r,w)

    double precision, intent(in   ) :: r
    double precision, intent(inout) :: w

    double precision alpha, beta, gamm, K, R1, R2, R3, discr
    integer sgn

    K = 59d0/60 - sqrt(29d0)/20

    R1 = r - ceiling(r) + 1
    R2 = R1*R1
    R3 = R2*R1
    alpha = 28
    beta  = 9d0/4 - 1.5 * (K + R2) + (22./3-7*K)*R1 - 7./3*R3;
    gamm = 0.25 * ( 0.5*(161./36 - 59./6*K + 5*K*K)*R2 + 1./3*(-109./24 + 5*K)*R2*R2 + 5./18*R3*R3  );
    discr = beta*beta - 4 * alpha * gamm;

    if((3/2 - K) .gt. 0) then

       sgn = 1

    else

       sgn = -1

    endif

    if(r .le. -3) then

       w = 0

    elseif(r .le. -2) then

       w = 1./(2*alpha) * ( -beta + sgn * sqrt(discr) );

    elseif(r .le. -1) then

       w = -3./(2*alpha) * ( -beta + sgn * sqrt(discr) ) - 1./16 + 1./8*( K+(r+2)*(r+2) ) + 1./12*(3*K-1)*(r+2) + 1./12*(r+2)*(r+2)*(r+2); 

    elseif(r .le. 0) then

       w = 2./(2*alpha) * ( -beta + sgn * sqrt(discr) ) + 1./4 + 1./6*(4-3*K)*(r+1) - 1./6*(r+1)*(r+1)*(r+1);

    elseif(r .le. 1) then

       w = 2./(2*alpha) * ( -beta + sgn * sqrt(discr) ) + 5./8 - 1./4 * ( K+r*r );

    elseif(r .le. 2) then

       w = -3./(2*alpha) * ( -beta + sgn * sqrt(discr) ) + 1./4 - 1./6*(4-3*K)*(r-1) + 1./6*(r-1)*(r-1)*(r-1);

    elseif(r .le. 3) then

       w = 1./(2*alpha) * ( -beta + sgn * sqrt(discr) ) - 1./16 + 1./8*(K+(r-2)*(r-2)) - 1./12*(3*K-1)*(r-2) - 1./12*(r-2)*(r-2)*(r-2); 	

    else

       w = 0

    endif

  end subroutine peskin_6pt

  subroutine get_weights(dxf, dxfinv, weights, indicies, &
                         coordsu, coordsulo, coordsuhi, &
                         coordsv, coordsvlo, coordsvhi, &
#if (BL_SPACEDIM == 3)
                         coordsw, coordswlo, coordswhi, &
#endif
                         part, ks, plof, rejected)
    
    double precision, intent(in   ) :: dxf(3), dxfinv(3), plof(3)
    integer,          intent(in   ) :: ks, coordsulo(3), coordsvlo(3), coordswlo(3), coordsuhi(3), coordsvhi(3), coordswhi(3)
    double precision, intent(inout) :: rejected
    type(particle_t), intent(in   ) :: part
    double precision, intent(inout) :: weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3)
    integer         , intent(inout) :: indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3)

    double precision, intent(in   ) :: coordsu(coordsulo(1):coordsuhi(1),coordsulo(2):coordsuhi(2),coordsulo(3):coordsuhi(3),1:AMREX_SPACEDIM)
    double precision, intent(in   ) :: coordsv(coordsvlo(1):coordsvhi(1),coordsvlo(2):coordsvhi(2),coordsvlo(3):coordsvhi(3),1:AMREX_SPACEDIM)
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: coordsw(coordswlo(1):coordswhi(1),coordswlo(2):coordswhi(2),coordswlo(3):coordswhi(3),1:AMREX_SPACEDIM)
#endif

    integer :: fi(3), fn(3),i, j, k
    double precision :: xx,yy,zz, w1, w2, w3, fr(3), fd(3), wcheck(3)

    !find fluid cell

    fr(1) = (part%pos(1) - plof(1))*dxfinv(1)
    fr(2) = (part%pos(2) - plof(2))*dxfinv(2)
    fr(3) = (part%pos(3) - plof(3))*dxfinv(3)

    !print*, "Real pos: ", fr, part%q

    fi(1) = floor(fr(1))
    fi(2) = floor(fr(2))
    fi(3) = floor(fr(3))

    fd(1) = fr(1) - fi(1)
    fd(2) = fr(2) - fi(2)
    fd(3) = fr(3) - fi(3)

    if(fd(1) .lt. 0.5) then
       fn(1) = -1
    else
       fn(1) = 0
    endif

    if(fd(2) .lt. 0.5) then
       fn(2) = -1
    else
       fn(2) = 0
    endif

    if(fd(3) .lt. 0.5) then
       fn(3) = -1
    else
       fn(3) = 0
    endif

    wcheck = 0

    rejected = 0

    if(  (((fi(1)-(ks-1)) .ge. coordsulo(1)) .and. ((fi(1)+ks) .le. coordsuhi(1)))  .and. (((fi(2)-(ks-1)+fn(2)) .ge. coordsulo(2)) .and. ((fi(2)+ks+fn(2)) .le. coordsuhi(2))) .and. (((fi(3)-(ks-1)+fn(3)) .ge. coordsulo(3)) .and. ((fi(3)+ks+fn(3)) .le. coordsuhi(3))) ) then

        if( ((((fi(1)-(ks-1)+fn(1)) .ge. coordsvlo(1)) .and. (fi(1)+ks+fn(1)) .le. coordsvhi(1)))  .and. (((fi(2)-(ks-1)) .ge. coordsvlo(2)) .and. ((fi(2)+ks) .le. coordsvhi(2))) .and. (((fi(3)-(ks-1)+fn(3)) .ge. coordsvlo(3)) .and. ((fi(3)+ks+fn(3)) .le. coordsvhi(3))) ) then

          if( ((((fi(1)-(ks-1)+fn(1)) .ge. coordswlo(1)) .and. (fi(1)+ks+fn(1)) .le. coordswhi(1)))  .and. (((fi(2)-(ks-1)+fn(2)) .ge. coordswlo(2)) .and. ((fi(2)+ks+fn(2)) .le. coordswhi(2))) .and. (((fi(3)-(ks-1)) .ge. coordswlo(3)) .and. ((fi(3)+ks) .le. coordswhi(3))) ) then

      do k = -(ks-1), ks
         do j = -(ks-1), ks
            do i = -(ks-1), ks

           !print *, "cord: ", fi(1)+i,fi(2)+j+fn(2),fi(3)+k+fn(3)      

               xx = part%pos(1) - coordsu(fi(1)+i,fi(2)+j+fn(2),fi(3)+k+fn(3),1)
               yy = part%pos(2) - coordsu(fi(1)+i,fi(2)+j+fn(2),fi(3)+k+fn(3),2)
               zz = part%pos(3) - coordsu(fi(1)+i,fi(2)+j+fn(2),fi(3)+k+fn(3),3)

               if(pkernel_fluid .eq. 3) then
                  call peskin_3pt(xx*dxfinv(1),w1)
                  call peskin_3pt(yy*dxfinv(2),w2)
                  call peskin_3pt(zz*dxfinv(3),w3)
               elseif(pkernel_fluid .eq. 4) then
                  call peskin_4pt(xx*dxfinv(1),w1)
                  call peskin_4pt(yy*dxfinv(2),w2)
                  call peskin_4pt(zz*dxfinv(3),w3)
               elseif(pkernel_fluid .eq. 6) then
                  call peskin_6pt(xx*dxfinv(1),w1)
                  call peskin_6pt(yy*dxfinv(2),w2)
                  call peskin_6pt(zz*dxfinv(3),w3)
               endif

               weights(i,j,k,1) = w1*w2*w3

               indicies(i,j,k,1,1) = fi(1)+i
               indicies(i,j,k,1,2) = fi(2)+j+fn(2)
               indicies(i,j,k,1,3) = fi(3)+k+fn(3)

               wcheck(1) = wcheck(1) + weights(i,j,k,1)

               !print *, fi(1)+i,fi(2)+j+fn(2),fi(3)+k+fn(3), coordsu(fi(1)+i,fi(2)+j+fn(2),fi(3)+k+fn(3),1), coordsu(fi(1)+i,fi(2)+j+fn(2),fi(3)+k+fn(3),2), coordsu(fi(1)+i,fi(2)+j+fn(2),fi(3)+k+fn(3),3)

               !print*, "xw: ", w1, "I: ", indicies(i,j,k,1,1), indicies(i,j,k,1,2), "D: ", xx*dxfinv(1), yy*dxfinv(2)

               !print *, i,j,k,w1*w2*w3

               !print *, "Accessing: ", fi(1)+i+fn(1),fi(2)+j,fi(3)+k+fn(3)

               xx = part%pos(1) - coordsv(fi(1)+i+fn(1),fi(2)+j,fi(3)+k+fn(3),1)
               yy = part%pos(2) - coordsv(fi(1)+i+fn(1),fi(2)+j,fi(3)+k+fn(3),2)
               zz = part%pos(3) - coordsv(fi(1)+i+fn(1),fi(2)+j,fi(3)+k+fn(3),3)

               if(pkernel_fluid .eq. 3) then
                  call peskin_3pt(xx*dxfinv(1),w1)
                  call peskin_3pt(yy*dxfinv(2),w2)
                  call peskin_3pt(zz*dxfinv(3),w3)
               elseif(pkernel_fluid .eq. 4) then
                  call peskin_4pt(xx*dxfinv(1),w1)
                  call peskin_4pt(yy*dxfinv(2),w2)
                  call peskin_4pt(zz*dxfinv(3),w3)
               elseif(pkernel_fluid .eq. 6) then
                  call peskin_6pt(xx*dxfinv(1),w1)
                  call peskin_6pt(yy*dxfinv(2),w2)
                  call peskin_6pt(zz*dxfinv(3),w3)
               endif

               weights(i,j,k,2) = w1*w2*w3

               indicies(i,j,k,2,1) = fi(1)+i+fn(1)
               indicies(i,j,k,2,2) = fi(2)+j
               indicies(i,j,k,2,3) = fi(3)+k+fn(3)

               wcheck(2) = wcheck(2) + weights(i,j,k,2)


               xx = part%pos(1) - coordsw(fi(1)+i+fn(1),fi(2)+j+fn(2),fi(3)+k,1)
               yy = part%pos(2) - coordsw(fi(1)+i+fn(1),fi(2)+j+fn(2),fi(3)+k,2)
               zz = part%pos(3) - coordsw(fi(1)+i+fn(1),fi(2)+j+fn(2),fi(3)+k,3)

               if(pkernel_fluid .eq. 3) then
                  call peskin_3pt(xx*dxfinv(1),w1)
                  call peskin_3pt(yy*dxfinv(2),w2)
                  call peskin_3pt(zz*dxfinv(3),w3)
               elseif(pkernel_fluid .eq. 4) then
                  call peskin_4pt(xx*dxfinv(1),w1)
                  call peskin_4pt(yy*dxfinv(2),w2)
                  call peskin_4pt(zz*dxfinv(3),w3)
               elseif(pkernel_fluid .eq. 6) then
                  call peskin_6pt(xx*dxfinv(1),w1)
                  call peskin_6pt(yy*dxfinv(2),w2)
                  call peskin_6pt(zz*dxfinv(3),w3)
               endif

               weights(i,j,k,3) = w1*w2*w3

               indicies(i,j,k,3,1) = fi(1)+i+fn(1)
               indicies(i,j,k,3,2) = fi(2)+j+fn(2)
               indicies(i,j,k,3,3) = fi(3)+k

               wcheck(3) = wcheck(3) + weights(i,j,k,3)

               !print *, i, j, k,  weights(i,j,k,:)

            enddo
         enddo
      enddo

      rejected =  0      

    endif
    endif
    endif


    !print*, "Total: ", wcheck

  end subroutine get_weights

  subroutine get_weights_scalar_cc(dx, dxinv, weights, indicies, &
                                   coords, coordslo, coordshi, &
                                   part, ks, lo, hi, plof, store)

    double precision, intent(in   ) :: dx(3), dxinv(3), plof(3)
    integer,          intent(in   ) :: ks, coordslo(3), coordshi(3), lo(3), hi(3), store
    type(particle_t), intent(in   ) :: part
    double precision, intent(inout) :: weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3)
    integer         , intent(inout) :: indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3)

    double precision, intent(in   ) :: coords(coordslo(1):coordshi(1),coordslo(2):coordshi(2),coordslo(3):coordshi(3),1:AMREX_SPACEDIM)

    integer :: fi(3), fn(3),i, j, k, wcount
    double precision :: xx,yy,zz, w1, w2, w3, fr(3), fd(3), wcheck

    !find scalar cell

    fr(1) = (part%pos(1) - plof(1))*dxinv(1)
    fr(2) = (part%pos(2) - plof(2))*dxinv(2)
    fr(3) = (part%pos(3) - plof(3))*dxinv(3)

    fi(1) = floor(fr(1))
    fi(2) = floor(fr(2))
    fi(3) = floor(fr(3))

    fd = fr - fi

    if(fd(1) .lt. 0.5) then
       fn(1) = -1
    else
       fn(1) = 0
    endif

    if(fd(2) .lt. 0.5) then
       fn(2) = -1
    else
       fn(2) = 0
    endif

    if(fd(3) .lt. 0.5) then
       fn(3) = -1
    else
       fn(3) = 0
    endif

    wcheck = 0
    wcount = 0

    do k = -(ks-1), ks
       do j = -(ks-1), ks
          do i = -(ks-1), ks

             xx = part%pos(1) - coords(fi(1)+i+fn(1),fi(2)+j+fn(2),fi(3)+k+fn(3),1)
             yy = part%pos(2) - coords(fi(1)+i+fn(1),fi(2)+j+fn(2),fi(3)+k+fn(3),2)
             zz = part%pos(3) - coords(fi(1)+i+fn(1),fi(2)+j+fn(2),fi(3)+k+fn(3),3)

             !        print *, "Relpos: ", xx*dxinv(1), yy*dxinv(2), zz*dxinv(3), xx, yy, zz
             !        print *, "Pos: ", part%pos
             !        print *, "Cell: ", fi
             !        print *, "coords: ", fi(1)+i+fn(1), fi(2)+j+fn(2), fi(3)+k+fn(3)
             !        print *, "Realcoords: ", coords(fi(1)+i+fn(1),fi(2)+j+fn(2),fi(3)+k+fn(3),1), coords(fi(1)+i+fn(1),fi(2)+j+fn(2),fi(3)+k+fn(3),2), coords(fi(1)+i+fn(1),fi(2)+j+fn(2),fi(3)+k+fn(3),3)

             if(pkernel_es .eq. 3) then
                call peskin_3pt(xx*dxinv(1),w1)
                call peskin_3pt(yy*dxinv(2),w2)
                call peskin_3pt(zz*dxinv(3),w3)
             elseif(pkernel_es .eq. 4) then
                call peskin_4pt(xx*dxinv(1),w1)
                call peskin_4pt(yy*dxinv(2),w2)
                call peskin_4pt(zz*dxinv(3),w3)
             elseif(pkernel_es .eq. 6) then
                call peskin_6pt(xx*dxinv(1),w1)
                call peskin_6pt(yy*dxinv(2),w2)
                call peskin_6pt(zz*dxinv(3),w3)
             endif

             weights(i,j,k,store) = w1*w2*w3

             if(weights(i,j,k,store) .ne. 0) then
                wcount = wcount +1
             endif
             indicies(i,j,k,store,1) = fi(1)+i+fn(1)
             indicies(i,j,k,store,2) = fi(2)+j+fn(2)
             indicies(i,j,k,store,3) = fi(3)+k+fn(3)

             wcheck = wcheck + weights(i,j,k,store)

             !print *, "Indicies: ", indicies(i,j,k,store,:)

          enddo
       enddo
    enddo

    !xx=0.5001
    !call peskin_3pt(xx,w1)

    !print *, "w: ", w1
    !print*, "Total: ", wcheck, "count: ", wcount, "kernel: ", pkernel_es

  end subroutine get_weights_scalar_cc

  subroutine spread_op(weights, indicies, &
                       sourceu, sourceulo, sourceuhi, &
                       sourcev, sourcevlo, sourcevhi, &
#if (BL_SPACEDIM == 3)
                       sourcew, sourcewlo, sourcewhi, &
#endif
                       part, ks, dxf)

    integer,          intent(in   ) :: ks, sourceulo(3), sourcevlo(3), sourcewlo(3), sourceuhi(3), sourcevhi(3), sourcewhi(3)
    double precision, intent(in   ) :: dxf(3)
    type(particle_t), intent(in   ) :: part
    double precision, intent(inout) :: weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3)
    integer         , intent(inout) :: indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3)

    double precision, intent(inout) :: sourceu(sourceulo(1):sourceuhi(1),sourceulo(2):sourceuhi(2),sourceulo(3):sourceuhi(3))
    double precision, intent(inout) :: sourcev(sourcevlo(1):sourcevhi(1),sourcevlo(2):sourcevhi(2),sourcevlo(3):sourcevhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(inout) :: sourcew(sourcewlo(1):sourcewhi(1),sourcewlo(2):sourcewhi(2),sourcewlo(3):sourcewhi(3))
#endif

    integer :: i, j, k, ii1, jj1, kk1, ii2, jj2, kk2, ii3, jj3, kk3
    double precision :: volinv, spreadcheck(3)

    volinv = 1/(dxf(1)*dxf(2)*dxf(3))
    spreadcheck = 0

    !print *, "Spreading ", part%force, ", volinv: ", volinv

    do k = -(ks-1), ks
       do j = -(ks-1), ks
          do i = -(ks-1), ks

             ii1 = indicies(i,j,k,1,1)
             jj1 = indicies(i,j,k,1,2)
             kk1 = indicies(i,j,k,1,3)

             

             !print *, "Touching pre: ", ii1,jj1,kk1,sourceu(ii1,jj1,kk1), part%force(1),weights(i,j,k,1)
             sourceu(ii1,jj1,kk1) = sourceu(ii1,jj1,kk1) + part%force(1)*weights(i,j,k,1)*volinv

             !print *, "Touching: ", ii1,jj1,kk1,weights(i,j,k,1), sourceu(ii1,jj1,kk1), part%force(1)

             !print *, "Touching post: ", ii1,jj1,kk1,sourceu(ii1,jj1,kk1), part%force(1),weights(i,j,k,1)

             spreadcheck(1) = spreadcheck(1) + sourceu(ii1,jj1,kk1)
             !print*, "S: ", sourceu(ii1,jj1,kk1)
             ii2 = indicies(i,j,k,2,1)
             jj2 = indicies(i,j,k,2,2)
             kk2 = indicies(i,j,k,2,3)

             sourcev(ii2,jj2,kk2) = sourcev(ii2,jj2,kk2) + part%force(2)*weights(i,j,k,2)*volinv

             spreadcheck(2) = spreadcheck(2) + sourcev(ii2,jj2,kk2)

             ii3 = indicies(i,j,k,3,1)
             jj3 = indicies(i,j,k,3,2)
             kk3 = indicies(i,j,k,3,3)

             sourcew(ii3,jj3,kk3) = sourcew(ii3,jj3,kk3) + part%force(3)*weights(i,j,k,3)*volinv

             spreadcheck(3) = spreadcheck(3) + sourcew(ii3,jj3,kk3)

          enddo
       enddo
    enddo

    !print *, "Spread1 ", spreadcheck/volinv

    spreadcheck = 0;

    do k = sourcewlo(3), sourcewhi(3)
       do j = sourcewlo(2), sourcewhi(2)
          do i = sourcewlo(1), sourcewhi(1)

            spreadcheck(3) = spreadcheck(3) + sourcew(i,j,k)

          enddo
       enddo
     enddo
    !  part => particles(1)
    ! part2 => particles(2)

    !print *, "Spread2 ", spreadcheck

  end subroutine spread_op

  subroutine spread_op_scalar_cc(weights, indicies, &
                                 source, sourcelo, sourcehi, &
                                 part, ks, dx, mq, store)

    integer,          intent(in   ) :: ks, sourcelo(3), sourcehi(3), store
    double precision, intent(in   ) :: dx(3), mq
    type(particle_t), intent(in   ) :: part
    double precision, intent(inout) :: weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3)
    integer         , intent(inout) :: indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3)

    double precision, intent(inout) :: source(sourcelo(1):sourcehi(1),sourcelo(2):sourcehi(2),sourcelo(3):sourcehi(3))

    integer :: i, j, k, ii, jj, kk
    double precision :: volinv, qm, pvol

    volinv = 1/(dx(1)*dx(2)*dx(3))

    if(pkernel_es .eq. 3) then 
       !pvol = 6.28319
       pvol = 1
    elseif(pkernel_es .eq. 4) then  
       pvol = 1
    elseif(pkernel_es .eq. 6) then  
       pvol = 1
    endif



    if(mq .eq. 0) then
       qm = pvol*part%q/permittivity

    else
       qm = part%mass
    endif

    do k = -(ks-1), ks
       do j = -(ks-1), ks
          do i = -(ks-1), ks

             ii = indicies(i,j,k,1,1)
             jj = indicies(i,j,k,1,2)
             kk = indicies(i,j,k,1,3)


             source(ii,jj,kk) = source(ii,jj,kk) + qm*weights(i,j,k,store)*volinv

          enddo
       enddo
    enddo

  end subroutine spread_op_scalar_cc

  subroutine inter_op(weights, indicies, &
                      velu, velulo, veluhi, &
                      velv, velvlo, velvhi, &
#if (BL_SPACEDIM == 3)
                      velw, velwlo, velwhi, &
#endif
                      part, ks, dxf, boundflag, midpoint, rejected)
    
    integer,          intent(in   ) :: ks, velulo(3), velvlo(3), velwlo(3), veluhi(3), velvhi(3), velwhi(3), midpoint
    integer,          intent(inout) :: boundflag
    double precision, intent(inout) :: rejected
    double precision, intent(in   ) :: dxf(3)
    type(particle_t), intent(inout) :: part
    double precision, intent(inout) :: weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3)
    integer         , intent(inout) :: indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3)

    double precision, intent(in   ) :: velu(velulo(1):veluhi(1),velulo(2):veluhi(2),velulo(3):veluhi(3))
    double precision, intent(in   ) :: velv(velvlo(1):velvhi(1),velvlo(2):velvhi(2),velvlo(3):velvhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: velw(velwlo(1):velwhi(1),velwlo(2):velwhi(2),velwlo(3):velwhi(3))
#endif

    integer :: i, j, k, ii, jj, kk
    double precision :: oldvel(3)


    boundflag = 0

    if(midpoint .eq. 0) then

       part%vel = 0

!print *, "Dims:", -(ks-1), ks

       do k = -(ks-1), ks
          do j = -(ks-1), ks
             do i = -(ks-1), ks

                ii = indicies(i,j,k,1,1)
                jj = indicies(i,j,k,1,2)
                kk = indicies(i,j,k,1,3)

                part%vel(1) = part%vel(1) + weights(i,j,k,1)*velu(ii,jj,kk)     

                ii = indicies(i,j,k,2,1)
                jj = indicies(i,j,k,2,2)
                kk = indicies(i,j,k,2,3)

                part%vel(2) = part%vel(2) + weights(i,j,k,2)*velv(ii,jj,kk)

                ii = indicies(i,j,k,3,1)
                jj = indicies(i,j,k,3,2)
                kk = indicies(i,j,k,3,3)

                part%vel(3) = part%vel(3) + weights(i,j,k,3)*velw(ii,jj,kk)

      ! print *, ii, jj, kk, velu(ii,jj,kk), weights(i,j,k,1)

             enddo
          enddo
       enddo

    else

       oldvel = part%vel
       part%vel = 0

       do k = -(ks-1), ks
          do j = -(ks-1), ks
             do i = -(ks-1), ks

                ii = indicies(i,j,k,1,1)
                jj = indicies(i,j,k,1,2)
                kk = indicies(i,j,k,1,3)

                if((ii .gt. veluhi(1)) .or. (ii .lt. velulo(1)) .or. (jj .gt. veluhi(2)) .or. (jj .lt. velulo(2)) .or. (kk .gt. veluhi(3)) .or. (kk .lt. velulo(3))) then
                   boundflag = 1
                else
                   part%vel(1) = part%vel(1) + weights(i,j,k,1)*velu(ii,jj,kk)
                endif

                !print*, "V: ", velu(ii,jj,kk), "I: ", i, j, k, "W: ", weights(i,j,k,1)

                ii = indicies(i,j,k,2,1)
                jj = indicies(i,j,k,2,2)
                kk = indicies(i,j,k,2,3)

                if((ii .gt. velvhi(1)) .or. (ii .lt. velvlo(1)) .or. (jj .gt. velvhi(2)) .or. (jj .lt. velvlo(2)) .or. (kk .gt. velvhi(3)) .or. (kk .lt. velvlo(3))) then
                   boundflag = 1
                else
                   part%vel(2) = part%vel(2) + weights(i,j,k,2)*velv(ii,jj,kk)
                endif

                ii = indicies(i,j,k,3,1)
                jj = indicies(i,j,k,3,2)
                kk = indicies(i,j,k,3,3)

                if((ii .gt. velwhi(1)) .or. (ii .lt. velwlo(1)) .or. (jj .gt. velwhi(2)) .or. (jj .lt. velwlo(2)) .or. (kk .gt. velwhi(3)) .or. (kk .lt. velwlo(3))) then
                   boundflag = 1
                else
                   part%vel(3) = part%vel(3) + weights(i,j,k,3)*velw(ii,jj,kk)
                endif
                !print *, ii, jj, kk, ": ", velu(ii,jj,kk), velv(ii,jj,kk), velw(ii,jj,kk)
                !print *, "weights: ", weights(i,j,k,:)

             enddo
          enddo
       enddo

       if(boundflag .eq. 1) then
          part%vel = oldvel
          rejected = rejected + 1
          print *, "Midpoint interpolation rejected."
       endif

    endif


    !print*, "Intervel: ", part%vel
    !print*, "a_rel: ", (1.0/(6*3.142*part%vel(1)*visc_coef))/dxf(1)

    part%multi = part%vel(1)

  end subroutine inter_op

  subroutine rfd(weights, indicies, &
                 sourceu, sourceulo, sourceuhi, &
                 sourcev, sourcevlo, sourcevhi, &
#if (BL_SPACEDIM == 3)
                 sourcew, sourcewlo, sourcewhi, &
#endif
                 coordsx, coordsxlo, coordsxhi, &
                 coordsy, coordsylo, coordsyhi, &
#if (BL_SPACEDIM == 3)
                 coordsz, coordszlo, coordszhi, &
#endif
                 part, ks, dxf, plof)

    integer,          intent(in   ) :: ks, sourceulo(3), sourcevlo(3), sourcewlo(3), sourceuhi(3), sourcevhi(3), sourcewhi(3)
    integer,          intent(in   ) :: coordsxlo(3), coordsylo(3), coordszlo(3), coordsxhi(3), coordsyhi(3), coordszhi(3)
    double precision, intent(in   ) :: dxf(3), plof(3)
    type(particle_t), intent(inout) :: part
    double precision, intent(inout) :: weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3)
    integer         , intent(inout) :: indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3)

    double precision, intent(inout) :: sourceu(sourceulo(1):sourceuhi(1),sourceulo(2):sourceuhi(2),sourceulo(3):sourceuhi(3))
    double precision, intent(inout) :: sourcev(sourcevlo(1):sourcevhi(1),sourcevlo(2):sourcevhi(2),sourcevlo(3):sourcevhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(inout) :: sourcew(sourcewlo(1):sourcewhi(1),sourcewlo(2):sourcewhi(2),sourcewlo(3):sourcewhi(3))
#endif
    double precision, intent(in   ) :: coordsx(coordsxlo(1):coordsxhi(1),coordsxlo(2):coordsxhi(2),coordsxlo(3):coordsxhi(3))
    double precision, intent(in   ) :: coordsy(coordsylo(1):coordsyhi(1),coordsylo(2):coordsyhi(2),coordsylo(3):coordsyhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: coordsz(coordszlo(1):coordszhi(1),coordszlo(2):coordszhi(2),coordszlo(3):coordszhi(3))
#endif

    integer :: i, j, k, ii, jj, kk
    double precision :: uloc, vloc, wloc, volinv, normalrand(3), delta, norm, dxfinv(3), rejected

    volinv = 1/(dxf(1)*dxf(2)*dxf(3))
    dxfinv = 1/dxf

    !delta = 1d-4*dxf(1)
    delta = rfd_delta*dxf(1)

    !print*, "Fluid vel: ", uloc, wloc, vloc

    call get_particle_normal(normalrand(1))
    call get_particle_normal(normalrand(2))
    call get_particle_normal(normalrand(3))

    part%pos = part%pos + delta*normalrand/2

    part%force(1) = variance_coef_mom*k_B*T_init(1)*normalrand(1)/(delta)
    part%force(2) = variance_coef_mom*k_B*T_init(1)*normalrand(2)/(delta)
    part%force(3) = variance_coef_mom*k_B*T_init(1)*normalrand(3)/(delta)

    !print *, "F: ", part%force

    call get_weights(dxf, dxfinv, weights, indicies, &
                     coordsx, coordsxlo, coordsxhi, &
                     coordsy, coordsylo, coordsyhi, &
#if (BL_SPACEDIM == 3)
                     coordsz, coordszlo, coordszhi, &
#endif
                     part, ks, plof, rejected)

    call spread_op(weights, indicies, &
                   sourceu, sourceulo, sourceuhi, &
                   sourcev, sourcevlo, sourcevhi, &
#if (BL_SPACEDIM == 3)
                   sourcew, sourcewlo, sourcewhi, &
#endif
                   part, ks, dxf)

    part%pos = part%pos - delta*normalrand

    call get_weights(dxf, dxfinv, weights, indicies, &
                     coordsx, coordsxlo, coordsxhi, &
                     coordsy, coordsylo, coordsyhi, &
#if (BL_SPACEDIM == 3)
                     coordsz, coordszlo, coordszhi, &
#endif
                     part, ks, plof, rejected)

    part%force(1) = -k_B*T_init(1)*normalrand(1)/(delta)
    part%force(2) = -k_B*T_init(1)*normalrand(2)/(delta)
    part%force(3) = -k_B*T_init(1)*normalrand(3)/(delta)

    call spread_op(weights, indicies, &
                   sourceu, sourceulo, sourceuhi, &
                   sourcev, sourcevlo, sourcevhi, &
#if (BL_SPACEDIM == 3)
                   sourcew, sourcewlo, sourcewhi, &
#endif
                   part, ks, dxf)

    part%pos = part%pos + delta*normalrand/2

  end subroutine rfd

  subroutine drag(weights, indicies, &
                  sourceu, sourceulo, sourceuhi, &
                  sourcev, sourcevlo, sourcevhi, &
#if (BL_SPACEDIM == 3)
                  sourcew, sourcewlo, sourcewhi, &
#endif
                  velu, velulo, veluhi, &
                  velv, velvlo, velvhi, &
#if (BL_SPACEDIM == 3)
                  velw, velwlo, velwhi, &
#endif
                  part, ks, dxf)

    integer,          intent(in   ) :: ks, sourceulo(3), sourcevlo(3), sourcewlo(3), sourceuhi(3), sourcevhi(3), sourcewhi(3), velulo(3), velvlo(3), velwlo(3), veluhi(3), velvhi(3), velwhi(3)
    double precision, intent(in   ) :: dxf(3)
    type(particle_t), intent(inout) :: part
    double precision, intent(inout) :: weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3)
    integer         , intent(inout) :: indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3)

    double precision, intent(inout) :: sourceu(sourceulo(1):sourceuhi(1),sourceulo(2):sourceuhi(2),sourceulo(3):sourceuhi(3))
    double precision, intent(inout) :: sourcev(sourcevlo(1):sourcevhi(1),sourcevlo(2):sourcevhi(2),sourcevlo(3):sourcevhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(inout) :: sourcew(sourcewlo(1):sourcewhi(1),sourcewlo(2):sourcewhi(2),sourcewlo(3):sourcewhi(3))
#endif
    double precision, intent(in   ) :: velu(velulo(1):veluhi(1),velulo(2):veluhi(2),velulo(3):veluhi(3))
    double precision, intent(in   ) :: velv(velvlo(1):velvhi(1),velvlo(2):velvhi(2),velvlo(3):velvhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: velw(velwlo(1):velwhi(1),velwlo(2):velwhi(2),velwlo(3):velwhi(3))
#endif

    integer :: i, j, k, ii, jj, kk, boundflag, midpoint
    double precision :: uloc, vloc, wloc, volinv, normalrand(3), delta, norm, rejected


    uloc = part%vel(1)
    vloc = part%vel(2)
    wloc = part%vel(3)

    midpoint = 0

    call inter_op(weights, indicies, &
         velu, velulo, veluhi, &
         velv, velvlo, velvhi, &
#if (BL_SPACEDIM == 3)
         velw, velwlo, velwhi, &
#endif
         part, ks, dxf, boundflag, midpoint, rejected)

    part%force(1) = part%force(1) + (uloc-part%vel(1))*k_B*T_init(1)/part%wet_diff
    part%force(2) = part%force(2) + (vloc-part%vel(2))*k_B*T_init(1)/part%wet_diff
    part%force(3) = part%force(3) + (wloc-part%vel(3))*k_B*T_init(1)/part%wet_diff

    part%vel(1) = uloc
    part%vel(2) = vloc
    part%vel(3) = wloc

  end subroutine drag

  subroutine emf(weights, indicies, &
                 sourceu, sourceulo, sourceuhi, &
                 sourcev, sourcevlo, sourcevhi, &
#if (BL_SPACEDIM == 3)
                 sourcew, sourcewlo, sourcewhi, &
#endif
                 fieldu, fieldulo, fielduhi, &
                 fieldv, fieldvlo, fieldvhi, &
#if (BL_SPACEDIM == 3)
                 fieldw, fieldwlo, fieldwhi, &
#endif
                 part, ks, dxp)

    integer,          intent(in   ) :: ks, sourceulo(3), sourcevlo(3), sourcewlo(3), sourceuhi(3), sourcevhi(3), sourcewhi(3), fieldulo(3), fieldvlo(3), fieldwlo(3), fielduhi(3), fieldvhi(3), fieldwhi(3)
    double precision, intent(in   ) :: dxp(3)
    type(particle_t), intent(inout) :: part
    double precision, intent(inout) :: weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3)
    integer         , intent(inout) :: indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3)

    double precision, intent(inout) :: sourceu(sourceulo(1):sourceuhi(1),sourceulo(2):sourceuhi(2),sourceulo(3):sourceuhi(3))
    double precision, intent(inout) :: sourcev(sourcevlo(1):sourcevhi(1),sourcevlo(2):sourcevhi(2),sourcevlo(3):sourcevhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(inout) :: sourcew(sourcewlo(1):sourcewhi(1),sourcewlo(2):sourcewhi(2),sourcewlo(3):sourcewhi(3))
#endif
    double precision, intent(in   ) :: fieldu(fieldulo(1):fielduhi(1),fieldulo(2):fielduhi(2),fieldulo(3):fielduhi(3))
    double precision, intent(in   ) :: fieldv(fieldvlo(1):fieldvhi(1),fieldvlo(2):fieldvhi(2),fieldvlo(3):fieldvhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: fieldw(fieldwlo(1):fieldwhi(1),fieldwlo(2):fieldwhi(2),fieldwlo(3):fieldwhi(3))
#endif

    integer :: i, j, k, ii, jj, kk, boundflag, midpoint
    double precision :: uloc(3), volinv, normalrand(3), delta, norm, rejected


    uloc = part%vel

    midpoint = 0

    !reusing the velocity interpolator here, so force is being temporarily stored in particle velocity - clean this up later.

    call inter_op(weights, indicies, &
         fieldu, fieldulo, fielduhi, &
         fieldv, fieldvlo, fieldvhi, &
#if (BL_SPACEDIM == 3)
         fieldw, fieldwlo, fieldwhi, &
#endif
         part, ks, dxp, boundflag, midpoint, rejected)

    !KK look here!
    !print *, "Poisson force: ", part%vel*part%q
    part%force = part%force + part%vel*part%q


    !print *, "FORCE: ", part%force
    part%vel = uloc

  end subroutine emf

  subroutine set_pos(part1, part2, dx, sep)

    type(particle_t), intent(inout) :: part1, part2
    double precision, intent(in   ) :: dx(3), sep

    double precision cosTheta, sinTheta, cosPhi, sinPhi, fx, fy, fz

    call get_angles(costheta, sintheta, cosphi, sinphi)
    call get_uniform(fx)
    call get_uniform(fy)
    call get_uniform(fz)

    part1%pos(1) = prob_hi(1)/2.0 + fx*dx(1) - (sep/2.0)*sintheta*cosphi*dx(1)
    part1%pos(2) = prob_hi(2)/2.0 + fy*dx(2) - (sep/2.0)*sintheta*sinphi*dx(2)
    part1%pos(3) = prob_hi(3)/2.0 + fz*dx(3) - (sep/2.0)*costheta*dx(3)

    part2%pos(1) = prob_hi(1)/2.0 + fx*dx(1) + (sep/2.0)*sintheta*cosphi*dx(1)
    part2%pos(2) = prob_hi(2)/2.0 + fy*dx(2) + (sep/2.0)*sintheta*sinphi*dx(2)
    part2%pos(3) = prob_hi(3)/2.0 + fz*dx(3) + (sep/2.0)*costheta*dx(3)

  end subroutine set_pos

  subroutine move_ions_fhd(particles, np, rejected, moves, maxspeed, maxdist, diffinst, &
                           lo, hi, &
                           cell_part_ids, cell_part_cnt, clo, chi, plo, phi, &
                           dx, dt, plof, dxf, dxe, &
                           velx, velxlo, velxhi, &
                           vely, velylo, velyhi, &
#if (BL_SPACEDIM == 3)
                           velz, velzlo, velzhi, &
#endif
                           efx, efxlo, efxhi, &
                           efy, efylo, efyhi, &
#if (BL_SPACEDIM == 3)
                           efz, efzlo, efzhi, &
#endif
                           coordsx, coordsxlo, coordsxhi, &
                           coordsy, coordsylo, coordsyhi, &
#if (BL_SPACEDIM == 3)
                           coordsz, coordszlo, coordszhi, &
#endif
                           
                           sourcex, sourcexlo, sourcexhi, &
                           sourcey, sourceylo, sourceyhi, &
#if (BL_SPACEDIM == 3)
                           sourcez, sourcezlo, sourcezhi, &
#endif
                           mobility, mlo, mhi, &
                           paramplanes, ns, kinetic, sw) &
                           bind(c,name="move_ions_fhd")

    integer,          intent(in   ) :: np, ns
    integer,          intent(in   ) :: lo(3), hi(3), clo(3), chi(3), velxlo(3), velxhi(3), velylo(3), velyhi(3)
    integer,          intent(in   ) :: efylo(3), efyhi(3), efxlo(3), efxhi(3), sw, mlo(3), mhi(3)
    integer,          intent(in   ) :: sourcexlo(3), sourcexhi(3), sourceylo(3), sourceyhi(3)
    integer,          intent(in   ) :: coordsxlo(3), coordsxhi(3), coordsylo(3), coordsyhi(3)
#if (AMREX_SPACEDIM == 3)
    integer,          intent(in   ) :: velzlo(3), velzhi(3), efzlo(3), efzhi(3), sourcezlo(3), sourcezhi(3), coordszlo(3), coordszhi(3)
#endif
    type(particle_t), intent(inout), target :: particles(np)
    type(paramplane_t),  intent(in),    target :: paramplanes(ns)
    double precision, intent(in   ) :: dx(3), dxf(3), dxe(3), dt, plo(3), phi(3), plof(3)
    double precision, intent(inout) :: kinetic
    double precision, intent(inout) :: rejected, moves, maxspeed, maxdist, diffinst
    double precision, intent(in   ) :: velx(velxlo(1):velxhi(1),velxlo(2):velxhi(2),velxlo(3):velxhi(3))
    double precision, intent(in   ) :: vely(velylo(1):velyhi(1),velylo(2):velyhi(2),velylo(3):velyhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: velz(velzlo(1):velzhi(1),velzlo(2):velzhi(2),velzlo(3):velzhi(3))
#endif

    double precision, intent(in   ) :: efx(efxlo(1):efxhi(1),efxlo(2):efxhi(2),efxlo(3):efxhi(3))
    double precision, intent(in   ) :: efy(efylo(1):efyhi(1),efylo(2):efyhi(2),efylo(3):efyhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: efz(efzlo(1):efzhi(1),efzlo(2):efzhi(2),efzlo(3):efzhi(3))
#endif

    double precision, intent(in   ) :: coordsx(coordsxlo(1):coordsxhi(1),coordsxlo(2):coordsxhi(2),coordsxlo(3):coordsxhi(3),1:AMREX_SPACEDIM)
    double precision, intent(in   ) :: coordsy(coordsylo(1):coordsyhi(1),coordsylo(2):coordsyhi(2),coordsylo(3):coordsyhi(3),1:AMREX_SPACEDIM)
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: coordsz(coordszlo(1):coordszhi(1),coordszlo(2):coordszhi(2),coordszlo(3):coordszhi(3),1:AMREX_SPACEDIM)
#endif

    double precision, intent(inout) :: sourcex(sourcexlo(1):sourcexhi(1),sourcexlo(2):sourcexhi(2),sourcexlo(3):sourcexhi(3))
    double precision, intent(inout) :: sourcey(sourceylo(1):sourceyhi(1),sourceylo(2):sourceyhi(2),sourceylo(3):sourceyhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(inout) :: sourcez(sourcezlo(1):sourcezhi(1),sourcezlo(2):sourcezhi(2),sourcezlo(3):sourcezhi(3))
#endif

    double precision, intent(inout) :: mobility(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3),1:(nspecies*AMREX_SPACEDIM))

    type(c_ptr),      intent(inout) :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int),   intent(inout) :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))

    integer :: i, j, k, p, cell_np, new_np, intside, intsurf, push, loopcount, pointcount, ks, boundflag, midpoint
    integer :: ni(3), fi(3)
    integer(c_int), pointer :: cell_parts(:)
    type(particle_t), pointer :: part
    type(paramplane_t), pointer :: surf
    real(amrex_real) dxinv(3), dxfinv(3), dxeinv(3), onemdxf(3), ixf(3), localvel(3), deltap(3), std
    real(amrex_real) normalrand(3), tempvel(3), intold, inttime, runerr, runtime, adj, adjalt
    real(amrex_real) domsize(3), posalt(3), propvec(3), norm(3), dry_terms(3)
    real(amrex_real) diffest, diffav, distav, veltest, posold(3), velold(3)
    real(amrex_real) speed, mb(3), dist, wcheck

    double precision, allocatable :: weights(:,:,:,:)
    integer, allocatable :: indicies(:,:,:,:,:)

    if(pkernel_fluid .eq. 3) then
       ks = 2
    elseif(pkernel_fluid .eq. 4) then
       ks = 3
    elseif(pkernel_fluid .eq. 6) then
       ks = 4
    endif

    allocate(weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3))
    allocate(indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3))

    domsize = phi - plo

    adj = 0.99999d0
    adjalt = 2d0*(1d0 - adj)

    dxinv = 1.d0/dx

    dxeinv = 1.d0/dxe

    dxfinv = 1.d0/dxf
    onemdxf = 1.d0 - dxf

    diffav = 0
    distav = 0
    diffinst = 0
    veltest = 0
    moves = 0
    rejected = 0

    maxdist = 0
    dist = 0 
    maxspeed = 0
    speed = 0

    kinetic = 0

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             cell_np = cell_part_cnt(i,j,k)
             call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])

             new_np = cell_np
             p = 1

             do while (p <= new_np)

                part => particles(cell_parts(p))

                !print *, "Moving", part%id
 
                !Get peskin kernel weights. Weights are stored in 'weights', indicies contains the indicies to which the weights are applied.

                call get_weights(dxf, dxfinv, weights, indicies, &
                     coordsx, coordsxlo, coordsxhi, &
                     coordsy, coordsylo, coordsyhi, &
#if (BL_SPACEDIM == 3)
                     coordsz, coordszlo, coordszhi, &
#endif
                     part, ks, plof, rejected)

                !use weights and indicies to interpolate velocity fields onto particle

                midpoint = 0
                !print *, "Pos1: ", part%pos

                call inter_op(weights, indicies, &
                     velx, velxlo, velxhi, &
                     vely, velylo, velyhi, &
#if (BL_SPACEDIM == 3)
                     velz, velzlo, velzhi, &
#endif
                     part, ks, dxf, boundflag, midpoint, rejected)

                   !print *, part%id, "Vel", part%vel(1), "Pos", part%pos(1)

                if(move_tog .eq. 2) then !mid point time stepping - First step 1/2 a time step then interpolate velocity field

                   posold = part%pos
                   velold = part%vel

                   runtime = dt*0.5

                   do while (runtime .gt. 0)

                      !check 

                      call find_intersect(part,runtime, paramplanes, ns, intsurf, inttime, intside, phi, plo)
                      !intsurf = 0
                      !inttime = runtime


                      posalt(1) = inttime*part%vel(1)*adjalt
                      posalt(2) = inttime*part%vel(2)*adjalt
#if (BL_SPACEDIM == 3)
                      posalt(3) = inttime*part%vel(3)*adjalt
#endif

                      ! move the particle in a straight line, adj factor prevents double detection of boundary intersection
                      part%pos(1) = part%pos(1) + inttime*part%vel(1)*adj
                      part%pos(2) = part%pos(2) + inttime*part%vel(2)*adj
#if (BL_SPACEDIM == 3)
                      part%pos(3) = part%pos(3) + inttime*part%vel(3)*adj
#endif

                      if(intsurf .gt. 0) then

                         surf => paramplanes(intsurf)

                         if(surf%periodicity .eq. 0) then


                           call apply_bc(surf, part, intside, domsize, push, 1, 1)


                           runtime = runtime - inttime

                         else

                          runtime = runtime - inttime

                          part%pos(1) = part%pos(1) + runtime*part%vel(1)
                          part%pos(2) = part%pos(2) + runtime*part%vel(2)
#if (BL_SPACEDIM == 3)
                          part%pos(3) = part%pos(3) + runtime*part%vel(3)
#endif
                          runtime = 0

                        endif

                      else
                       runtime = 0

                      endif


                   end do



                   midpoint = 1
                   moves = moves + 1

                   call get_weights(dxf, dxfinv, weights, indicies, &
                        coordsx, coordsxlo, coordsxhi, &
                        coordsy, coordsylo, coordsyhi, &
#if (BL_SPACEDIM == 3)
                        coordsz, coordszlo, coordszhi, &
#endif
                        part, ks, plof, wcheck)


                   !if(rejected .ne. 0) then

                        !print *, "Interpolating"
                     call inter_op(weights, indicies, &
                          velx, velxlo, velxhi, &
                          vely, velylo, velyhi, &
#if (BL_SPACEDIM == 3)
                          velz, velzlo, velzhi, &
#endif
                          part, ks, dxf, boundflag, midpoint, rejected)

                    !endif


                    !print *, "Moved", part%pos(1)-posold(1)

                    part%pos = posold

        

                endif
                
                runtime = dt

                if (dry_move_tog .eq. 1) then

!                   !Get fluid cell - possibly replace this with peskin interp
!                   fi(1) = floor((part%pos(1) - plof(1))*dxfinv(1))
!                   fi(2) = floor((part%pos(2) - plof(2))*dxfinv(2))
!#if (BL_SPACEDIM == 3)
!                   fi(3) = floor((part%pos(3) - plof(3))*dxfinv(3))
!#else
!                   fi(3) = 0
!#endif
!                   mb(1) = mobility(fi(1),fi(2),fi(3),(part%species-1)*AMREX_SPACEDIM + 1)
!                   mb(2) = mobility(fi(1),fi(2),fi(3),(part%species-1)*AMREX_SPACEDIM + 2)
!#if (BL_SPACEDIM == 3)
!                   mb(3) = mobility(fi(1),fi(2),fi(3),(part%species-1)*AMREX_SPACEDIM + 3)
!#endif
                   call get_explicit_mobility(mb, part, plo, phi)
                   call dry(dt,part,dry_terms, mb)

                   part%vel = part%vel + dry_terms
                endif

                !print *, "adding", dry_terms(1), mb

                speed = part%vel(1)**2 + part%vel(2)**2 + part%vel(3)**2              

                if(speed .gt. maxspeed) then
                   maxspeed = speed
                endif

                do while (runtime .gt. 0)

                   call find_intersect(part,runtime, paramplanes, ns, intsurf, inttime, intside, phi, plo)

                   posalt(1) = inttime*part%vel(1)*adjalt
                   posalt(2) = inttime*part%vel(2)*adjalt
#if (BL_SPACEDIM == 3)
                   posalt(3) = inttime*part%vel(3)*adjalt
#endif

                   ! move the particle in a straight line, adj factor prevents double detection of boundary intersection

                   part%pos(1) = part%pos(1) + inttime*(part%vel(1))*adj
                   part%pos(2) = part%pos(2) + inttime*(part%vel(2))*adj
#if (BL_SPACEDIM == 3)
                   part%pos(3) = part%pos(3) + inttime*(part%vel(3))*adj
#endif
                   runtime = runtime - inttime


                   if(intsurf .gt. 0) then

                      surf => paramplanes(intsurf)

                      call apply_bc(surf, part, intside, domsize, push, 1, 1)


                      if(push .eq. 1) then

                         part%pos(1) = part%pos(1) + posalt(1)
                         part%pos(2) = part%pos(2) + posalt(2)
#if (BL_SPACEDIM == 3)
                         part%pos(3) = part%pos(3) + posalt(3)
#endif
                      endif

                   endif

                end do

                !print *, part%id, "pre Vel", part%vel(1), "Pos", part%pos(1)

!!!!!!!!!! Mean square displacement measurer.

                !print *, "Before: ", part%abspos

                part%abspos = part%abspos + dt*part%vel

                !print *, "After: ", part%abspos

                dist = sqrt(dot_product(dt*part%vel,dt*part%vel))/part%radius

                if(dist .gt. maxdist) then

                   maxdist = dist

                endif

                !              distav = distav + dt*sqrt(part%vel(1)**2+part%vel(2)**2+part%vel(3)**2)

                part%travel_time = part%travel_time + dt

                !print *, "tt: ", part%travel_time

                norm = part%abspos

                diffest = (norm(1)**2 + norm(2)**2 + norm(3)**2)/(6*part%travel_time)

                diffinst = diffinst + diffest
 

                !print *, "Diffest: ", diffest

                !              part%step_count = part%step_count + 1

                ! if it has changed cells, remove from vector.
                ! otherwise continue
!                print *, "HERE1!"
                ni(1) = floor((part%pos(1) - plo(1))*dxinv(1))
                ni(2) = floor((part%pos(2) - plo(2))*dxinv(2))
#if (BL_SPACEDIM == 3)
                ni(3) = floor((part%pos(3) - plo(3))*dxinv(3))
#else
                ni(3) = 0
#endif

                if ((ni(1) /= i) .or. (ni(2) /= j) .or. (ni(3) /= k)) then
                   part%sorted = 0
                   call remove_particle_from_cell(cell_parts, cell_np, new_np, p)  
                else
                   p = p + 1
                end if
             end do



             cell_part_cnt(i,j,k) = new_np

          end do
       end do
    end do

    ! write out grid/tile diagnostics
    if (.false.) then
       print *, "On grid/tile ", lo, hi
       print *, "I see ", np, " particles." 
       if(move_tog .eq. 2) then
          print *, "Fraction of midpoint moves rejected: ", rejected/moves
       endif
       print *, "Maximum observed speed: ", sqrt(maxspeed)
       print *, "Maximum observed displacement (fraction of radius): ", maxdist
       print *, "Average diffusion coefficient: ", diffinst/np
    end if

    deallocate(weights)
    deallocate(indicies)

  end subroutine move_ions_fhd

  subroutine spread_ions_fhd(particles, np, lo, hi, &
                             cell_part_ids, cell_part_cnt, clo, chi, plo, phi, &
                             dx, dt, plof, dxf, dxe, &
                             velx, velxlo, velxhi, &
                             vely, velylo, velyhi, &
#if (BL_SPACEDIM == 3)
                             velz, velzlo, velzhi, &
#endif
                             efx, efxlo, efxhi, &
                             efy, efylo, efyhi, &
#if (BL_SPACEDIM == 3)
                             efz, efzlo, efzhi, &
#endif
                             charge, chargelo, chargehi, &
                             coordsx, coordsxlo, coordsxhi, &
                             coordsy, coordsylo, coordsyhi, &
#if (BL_SPACEDIM == 3)
                             coordsz, coordszlo, coordszhi, &
#endif
                             cellcenters, cellcenterslo, cellcentershi, &
                             sourcex, sourcexlo, sourcexhi, &
                             sourcey, sourceylo, sourceyhi, &
#if (BL_SPACEDIM == 3)
                             sourcez, sourcezlo, sourcezhi, &
#endif
                             paramplanes, ns, potential, sw) &
                             bind(c,name="spread_ions_fhd")

    integer,          intent(in   )         :: np, ns, lo(3), hi(3), clo(3), chi(3), velxlo(3), velxhi(3), velylo(3), velyhi(3), efylo(3), efyhi(3), efxlo(3), efxhi(3), sw
    integer,          intent(in   )         :: sourcexlo(3), sourcexhi(3), sourceylo(3), sourceyhi(3)
    integer,          intent(in   )         :: coordsxlo(3), coordsxhi(3), coordsylo(3), coordsyhi(3)
#if (AMREX_SPACEDIM == 3)
    integer,          intent(in   )         :: velzlo(3), velzhi(3), efzlo(3), efzhi(3), sourcezlo(3), sourcezhi(3), coordszlo(3), coordszhi(3), chargelo(3), chargehi(3)
#endif
    type(particle_t), intent(inout), target :: particles(np)
    type(paramplane_t),  intent(in),    target :: paramplanes(ns)

    integer,          intent(in   )         :: cellcenterslo(3), cellcentershi(3)

    double precision, intent(in   )         :: dx(3), dxf(3), dxe(3), dt, plo(3), phi(3), plof(3)

    double precision, intent(in   ) :: velx(velxlo(1):velxhi(1),velxlo(2):velxhi(2),velxlo(3):velxhi(3))
    double precision, intent(in   ) :: vely(velylo(1):velyhi(1),velylo(2):velyhi(2),velylo(3):velyhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: velz(velzlo(1):velzhi(1),velzlo(2):velzhi(2),velzlo(3):velzhi(3))
#endif

    double precision, intent(in   ) :: efx(efxlo(1):efxhi(1),efxlo(2):efxhi(2),efxlo(3):efxhi(3))
    double precision, intent(in   ) :: efy(efylo(1):efyhi(1),efylo(2):efyhi(2),efylo(3):efyhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: efz(efzlo(1):efzhi(1),efzlo(2):efzhi(2),efzlo(3):efzhi(3))
#endif

    double precision, intent(in   ) :: charge(chargelo(1):chargehi(1),chargelo(2):chargehi(2),chargelo(3):chargehi(3))

    double precision, intent(in   ) :: coordsx(coordsxlo(1):coordsxhi(1),coordsxlo(2):coordsxhi(2),coordsxlo(3):coordsxhi(3),1:AMREX_SPACEDIM)
    double precision, intent(in   ) :: coordsy(coordsylo(1):coordsyhi(1),coordsylo(2):coordsyhi(2),coordsylo(3):coordsyhi(3),1:AMREX_SPACEDIM)
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: coordsz(coordszlo(1):coordszhi(1),coordszlo(2):coordszhi(2),coordszlo(3):coordszhi(3),1:AMREX_SPACEDIM)
#endif

    double precision, intent(in   ) :: cellcenters(cellcenterslo(1):cellcentershi(1),cellcenterslo(2):cellcentershi(2),cellcenterslo(3):cellcentershi(3),1:AMREX_SPACEDIM)

    double precision, intent(inout) :: sourcex(sourcexlo(1):sourcexhi(1),sourcexlo(2):sourcexhi(2),sourcexlo(3):sourcexhi(3))
    double precision, intent(inout) :: sourcey(sourceylo(1):sourceyhi(1),sourceylo(2):sourceyhi(2),sourceylo(3):sourceyhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(inout) :: sourcez(sourcezlo(1):sourcezhi(1),sourcezlo(2):sourcezhi(2),sourcezlo(3):sourcezhi(3))
#endif

    type(c_ptr),      intent(inout) :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int),   intent(inout) :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))

    double precision, intent(inout) :: potential

    integer :: i, j, k, p, cell_np, new_np, intside, intsurf, push, loopcount, pointcount, ks, boundflag, midpoint, store
    integer :: ni(3), fi(3)
    integer(c_int), pointer :: cell_parts(:)
    type(particle_t), pointer :: part, part2
    type(paramplane_t), pointer :: surf
    real(amrex_real) dxinv(3), dxfinv(3), dxeinv(3), onemdxf(3), ixf(3), localvel(3), deltap(3), std, normalrand(3), tempvel(3), intold, inttime, runerr, runtime, adj, adjalt, domsize(3), posalt(3), propvec(3), norm(3), &
         diffest, diffav, distav, diffinst, veltest, posold(3), delta, volinv, sep, rejected

    double precision, allocatable :: weights(:,:,:,:)
    integer, allocatable :: indicies(:,:,:,:,:)
    real(amrex_real) poisson_force(3)

    if(pkernel_fluid .gt. pkernel_es) then
       if(pkernel_fluid .eq. 3) then
          ks = 2
       elseif(pkernel_fluid .eq. 4) then
          ks = 3
       elseif(pkernel_fluid .eq. 6) then
          ks = 4
       endif
    else
       if(pkernel_es .eq. 3) then
          ks = 2
       elseif(pkernel_es .eq. 4) then
          ks = 3
       elseif(pkernel_es .eq. 6) then
          ks = 4
       endif
    endif

    allocate(weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3))
    allocate(indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3))

    domsize = phi - plo

    adj = 0.999999
    adjalt = 2d0*(1d0 - adj)

    dxinv = 1.d0/dx

    dxeinv = 1.d0/dxe

    dxfinv = 1.d0/dxf
    onemdxf = 1.d0 - dxf

    diffav = 0
    distav = 0
    diffinst = 0
    veltest = 0
    potential = 0

    p = 1

    do while (p <= np)

       part => particles(p)

       if(es_tog .eq. 2) then

          call calculate_force(particles, np, lo, hi, cell_part_ids, cell_part_cnt, clo, chi, plo, phi, p) !pairwise coulomb calc

          !print *, "Coulomb force: ", part%force
          potential = potential + part%potential
       endif
       !print *, "start weights", p
       !Get peskin kernel weights. Weights are stored in 'weights', indicies contains the indicies to which the weights are applied.
       call get_weights(dxf, dxfinv, weights, indicies, &
            coordsx, coordsxlo, coordsxhi, &
            coordsy, coordsylo, coordsyhi, &
#if (BL_SPACEDIM == 3)
            coordsz, coordszlo, coordszhi, &
#endif
            part, ks, plof, rejected)
       !print *, "finish weights", p
       if(drag_tog .eq. 1) then

          call drag(weights, indicies, &
               sourcex, sourcexlo, sourcexhi, &
               sourcey, sourceylo, sourceyhi, &
#if (BL_SPACEDIM == 3)
               sourcez, sourcezlo, sourcezhi, &
#endif
               velx, velxlo, velxhi, &
               vely, velylo, velyhi, &
#if (BL_SPACEDIM == 3)
               velz, velzlo, velzhi, &
#endif
               part, ks, dxf)
       endif


       !-------------es part, look at more efficient way of doing this

       store = 1
       call get_weights_scalar_cc(dxe, dxeinv, weights, indicies, &
            cellcenters, cellcenterslo, cellcentershi, &
            part, ks, lo, hi, plof, store)

       store = 2
       call get_weights_scalar_cc(dxe, dxeinv, weights, indicies, &
            cellcenters, cellcenterslo, cellcentershi, &
            part, ks, lo, hi, plof, store)

       store = 3
       call get_weights_scalar_cc(dxe, dxeinv, weights, indicies, &
            cellcenters, cellcenterslo, cellcentershi, &
            part, ks, lo, hi, plof, store)
       !print*, "particle ", p, " force before emf:   ", part%force
       poisson_force = -1.d0*part%force
       call emf(weights, indicies, &
            sourcex, sourcexlo, sourcexhi, &
            sourcey, sourceylo, sourceyhi, &
#if (BL_SPACEDIM == 3)
            sourcez, sourcezlo, sourcezhi, &
#endif
            efx, efxlo, efxhi, &
            efy, efylo, efyhi, &
#if (BL_SPACEDIM == 3)
            efz, efzlo, efzhi, &
#endif
            part, ks, dxe)

       poisson_force = poisson_force + part%force
      !print*, "Poisson force on particle ", p, "is: ", poisson_force
       !print*, "particle ", p, " force after emf:    ", part%force
       !------------------

       call get_weights(dxf, dxfinv, weights, indicies, &
            coordsx, coordsxlo, coordsxhi, &
            coordsy, coordsylo, coordsyhi, &
#if (BL_SPACEDIM == 3)
            coordsz, coordszlo, coordszhi, &
#endif
            part, ks, plof, rejected)

       !  print*, "SPREAD"

       !      if(mod(int(part%step_count),100) .eq. 0) then

       !        part%step_count = 0
       !        part%diff_av = 0
       !        part%travel_time = 0

       !      endif


       !      part%step_count = part%step_count + 1
       !      part%diff_av = part%diff_av + norm2(part%force)
       !      part%travel_time = part%travel_time + (norm2(part%force) - part%diff_av/(part%step_count))**2


       !      if(mod(int(part%step_count),10) .eq. 0) then
       !        
       !        print *, "Sep: ", part%drag_factor
       !        print *, "Average: ", part%diff_av/part%step_count
       !        !print *, "% SD: ", 100+100*(sqrt(part%travel_time)/part%step_count-part%diff_av/part%step_count)/(part%diff_av/part%step_count)
       !        part%drag_factor = part%drag_factor + 0.1

       !      endif
       ! print *, "Part force: ", norm2(part%force)

    tempvel(3) = 0;

    do k = sourcezlo(3), sourcezhi(3)
       do j = sourcezlo(2), sourcezhi(2)
          do i = sourcezlo(1), sourcezhi(1)

            tempvel(3) = tempvel(3) + sourcez(i,j,k)

          enddo
       enddo
     enddo
    !  part => particles(1)
    ! part2 => particles(2)

    !print *, "SOURCE1:", tempvel(3)*(prob_hi(1)-prob_lo(1))*(prob_hi(2)-prob_lo(2))*(prob_hi(3)-prob_lo(3))

       call spread_op(weights, indicies, &
            sourcex, sourcexlo, sourcexhi, &
            sourcey, sourceylo, sourceyhi, &
#if (BL_SPACEDIM == 3)
            sourcez, sourcezlo, sourcezhi, &
#endif
            part, ks, dxf)
       !print*, 'Part force at end of spread_ions: ', part%force, part%q
       p = p + 1

    end do

    tempvel(3) = 0;

    do k = sourcezlo(3), sourcezhi(3)
       do j = sourcezlo(2), sourcezhi(2)
          do i = sourcezlo(1), sourcezhi(1)

            tempvel(3) = tempvel(3) + sourcez(i,j,k)

          enddo
       enddo
     enddo
    !  part => particles(1)
    ! part2 => particles(2)

!    print *, "SOURCE2:", tempvel(3)*(prob_hi(1)-prob_lo(1))*(prob_hi(2)-prob_lo(2))*(prob_hi(3)-prob_lo(3))
    !print *, "SOURCE2:", tempvel(3)

    !call set_pos(part, part2,dxe,part%drag_factor)

    deallocate(weights)
    deallocate(indicies)

  end subroutine spread_ions_fhd

  subroutine do_rfd(particles, np, lo, hi, &
                    cell_part_ids, cell_part_cnt, clo, chi, plo, phi, &
                    dx, dt, plof, dxf, dxe, &
                    velx, velxlo, velxhi, &
                    vely, velylo, velyhi, &
#if (BL_SPACEDIM == 3)
                    velz, velzlo, velzhi, &
#endif
                    efx, efxlo, efxhi, &
                    efy, efylo, efyhi, &
#if (BL_SPACEDIM == 3)
                    efz, efzlo, efzhi, &
#endif
                    coordsx, coordsxlo, coordsxhi, &
                    coordsy, coordsylo, coordsyhi, &
#if (BL_SPACEDIM == 3)
                    coordsz, coordszlo, coordszhi, &
#endif
                    cellcenters, cellcenterslo, cellcentershi, &
                    sourcex, sourcexlo, sourcexhi, &
                    sourcey, sourceylo, sourceyhi, &
#if (BL_SPACEDIM == 3)
                    sourcez, sourcezlo, sourcezhi, &
#endif
                    paramplanes, ns, sw) &
                    bind(c,name="do_rfd")

    integer,          intent(in   )         :: np, ns, lo(3), hi(3), clo(3), chi(3), velxlo(3), velxhi(3), velylo(3), velyhi(3), efylo(3), efyhi(3), efxlo(3), efxhi(3), sw
    integer,          intent(in   )         :: sourcexlo(3), sourcexhi(3), sourceylo(3), sourceyhi(3)
    integer,          intent(in   )         :: coordsxlo(3), coordsxhi(3), coordsylo(3), coordsyhi(3)
#if (AMREX_SPACEDIM == 3)
    integer,          intent(in   )         :: velzlo(3), velzhi(3), efzlo(3), efzhi(3), sourcezlo(3), sourcezhi(3), coordszlo(3), coordszhi(3)
#endif
    type(particle_t), intent(inout), target :: particles(np)
    type(paramplane_t),  intent(in),    target :: paramplanes(ns)

    integer,          intent(in   )         :: cellcenterslo(3), cellcentershi(3)

    double precision, intent(in   )         :: dx(3), dxf(3), dxe(3), dt, plo(3), phi(3), plof(3)

    double precision, intent(in   ) :: velx(velxlo(1):velxhi(1),velxlo(2):velxhi(2),velxlo(3):velxhi(3))
    double precision, intent(in   ) :: vely(velylo(1):velyhi(1),velylo(2):velyhi(2),velylo(3):velyhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: velz(velzlo(1):velzhi(1),velzlo(2):velzhi(2),velzlo(3):velzhi(3))
#endif

    double precision, intent(in   ) :: efx(efxlo(1):efxhi(1),efxlo(2):efxhi(2),efxlo(3):efxhi(3))
    double precision, intent(in   ) :: efy(efylo(1):efyhi(1),efylo(2):efyhi(2),efylo(3):efyhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: efz(efzlo(1):efzhi(1),efzlo(2):efzhi(2),efzlo(3):efzhi(3))
#endif

    double precision, intent(in   ) :: coordsx(coordsxlo(1):coordsxhi(1),coordsxlo(2):coordsxhi(2),coordsxlo(3):coordsxhi(3),1:AMREX_SPACEDIM)
    double precision, intent(in   ) :: coordsy(coordsylo(1):coordsyhi(1),coordsylo(2):coordsyhi(2),coordsylo(3):coordsyhi(3),1:AMREX_SPACEDIM)
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: coordsz(coordszlo(1):coordszhi(1),coordszlo(2):coordszhi(2),coordszlo(3):coordszhi(3),1:AMREX_SPACEDIM)
#endif

    double precision, intent(in   ) :: cellcenters(cellcenterslo(1):cellcentershi(1),cellcenterslo(2):cellcentershi(2),cellcenterslo(3):cellcentershi(3),1:AMREX_SPACEDIM)

    double precision, intent(inout) :: sourcex(sourcexlo(1):sourcexhi(1),sourcexlo(2):sourcexhi(2),sourcexlo(3):sourcexhi(3))
    double precision, intent(inout) :: sourcey(sourceylo(1):sourceyhi(1),sourceylo(2):sourceyhi(2),sourceylo(3):sourceyhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(inout) :: sourcez(sourcezlo(1):sourcezhi(1),sourcezlo(2):sourcezhi(2),sourcezlo(3):sourcezhi(3))
#endif

    type(c_ptr),      intent(inout) :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int),   intent(inout) :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))

    integer :: i, j, k, p, cell_np, new_np, intside, intsurf, push, loopcount, pointcount, ks, boundflag, midpoint, store
    integer :: ni(3), fi(3)
    integer(c_int), pointer :: cell_parts(:)
    type(particle_t), pointer :: part
    type(paramplane_t), pointer :: surf
    real(amrex_real) dxinv(3), dxfinv(3), dxeinv(3), onemdxf(3), ixf(3), localvel(3), deltap(3), std, normalrand(3), tempvel(3), intold, inttime, runerr, runtime, domsize(3), posalt(3), propvec(3), norm(3), &
         diffest, diffav, distav, diffinst, veltest, posold(3), delta, volinv

    double precision, allocatable :: weights(:,:,:,:)
    integer, allocatable :: indicies(:,:,:,:,:)

    if(pkernel_fluid .gt. pkernel_es) then
       if(pkernel_fluid .eq. 3) then
          ks = 2
       elseif(pkernel_fluid .eq. 4) then
          ks = 3
       elseif(pkernel_fluid .eq. 6) then
          ks = 3
       endif
    else
       if(pkernel_es .eq. 3) then
          ks = 2
       elseif(pkernel_es .eq. 4) then
          ks = 3
       elseif(pkernel_es .eq. 6) then
          ks = 3
       endif
    endif

    allocate(weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3))
    allocate(indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3))

    domsize = phi - plo


    dxinv = 1.d0/dx

    dxeinv = 1.d0/dxe

    dxfinv = 1.d0/dxf
    onemdxf = 1.d0 - dxf

    diffav = 0
    distav = 0
    diffinst = 0
    veltest = 0


    p = 1

    do while (p <= np)

       part => particles(p)

       if((rfd_tog .eq. 1) .and. (variance_coef_mom .ne. 0)) then
          part%force = 0

          call rfd(weights, indicies, &
               sourcex, sourcexlo, sourcexhi, &
               sourcey, sourceylo, sourceyhi, &
#if (BL_SPACEDIM == 3)
               sourcez, sourcezlo, sourcezhi, &
#endif
               coordsx, coordsxlo, coordsxhi, &
               coordsy, coordsylo, coordsyhi, &
#if (BL_SPACEDIM == 3)
               coordsz, coordszlo, coordszhi, &
#endif
               part, ks, dxf, plof)

       endif

       part%force = 0

       p = p + 1


    end do

    deallocate(weights)
    deallocate(indicies)

  end subroutine do_rfd

  
  subroutine collect_charge(particles, np, lo, hi, &
                            cell_part_ids, cell_part_cnt, clo, chi, plo, phi, dx, dt, ploes, dxes, &
                            cellcenters, cellcenterslo, cellcentershi, &
                            charge, chargelo, chargehi) &
                            bind(c,name="collect_charge")

    integer,          intent(in   )         :: np, lo(3), hi(3), clo(3), chi(3)
    integer,          intent(in   )         :: chargelo(3), chargehi(3)
    integer,          intent(in   )         :: cellcenterslo(3), cellcentershi(3)

    type(particle_t), intent(inout), target :: particles(np)

    double precision, intent(in   )         :: dx(3), dxes(3), dt, plo(3), phi(3), ploes(3)

    double precision, intent(in   ) :: cellcenters(cellcenterslo(1):cellcentershi(1),cellcenterslo(2):cellcentershi(2),cellcenterslo(3):cellcentershi(3),1:AMREX_SPACEDIM)
    double precision, intent(inout) :: charge(chargelo(1):chargehi(1),chargelo(2):chargehi(2),chargelo(3):chargehi(3))

    type(c_ptr),      intent(inout) :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int),   intent(inout) :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))

    integer :: i, j, k, p, cell_np, new_np, intside, intsurf, push, loopcount, pointcount, ks, store
    integer :: ni(3), fi(3)
    integer(c_int), pointer :: cell_parts(:)
    type(particle_t), pointer :: part
    type(paramplane_t), pointer :: surf
    real(amrex_real) dxinv(3), dxesinv(3), onemdxf(3), ixf(3), diffav, distav, domsize(3), qm, diffinst, volinv

    double precision, allocatable :: weights(:,:,:,:)
    integer, allocatable :: indicies(:,:,:,:,:)

    if(pkernel_es .eq. 3) then
       ks = 2
    elseif(pkernel_es .eq. 4) then
       ks = 3
    elseif(pkernel_es .eq. 6) then
       ks = 4
    endif

    allocate(weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3))
    allocate(indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3))

    domsize = phi - plo

    dxinv = 1.d0/dx

    dxesinv = 1.d0/dxes

    volinv = 1d0/(dxes(1)*dxes(2)*dxes(3))


    diffav = 0
    distav = 0
    diffinst = 0

    store = 1
    qm = 0

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             cell_np = cell_part_cnt(i,j,k)

             call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])

             new_np = cell_np
             p = 1
             do while (p <= new_np)

                !              runtime = dt
                part => particles(cell_parts(p))


                call get_weights_scalar_cc(dxes, dxesinv, weights, indicies, &
                     cellcenters, cellcenterslo, cellcentershi, &
                     part, ks, lo, hi, ploes, store)


                call spread_op_scalar_cc(weights, indicies, &
                     charge, chargelo, chargehi, &
                     part, ks, dxes, qm, store)

                p = p + 1

             enddo

          end do
       end do
    end do

    deallocate(weights)
    deallocate(indicies)

  end subroutine collect_charge

  subroutine get_mobility(nmob, tmob, spec, z)

    type(species_t), intent(in   ) :: spec
    real(amrex_real),intent(in   ) :: z
    real(amrex_real),intent(inout) :: nmob, tmob

    real(amrex_real) a

    a = k_b*t_init(1)/(spec%dry_diff*visc_coef*3.142*6)

    !print *, "z, a, nmob: ", z, a, nmob

    nmob = max(1 - 9*a/(8*z) + (a**3)/(2*z**3) - (a**5)/(8*(z**5)),0d0)
    tmob = max(1 - 9*a/(16*z) + 2*(a**3)/(16*(z**3)) - (a**5)/(16*(z**5)),0d0)

  end subroutine get_mobility


  subroutine mob_interp(h, tmob, nmob)

    real(amrex_real),intent(in   ) :: h
    real(amrex_real),intent(inout) :: nmob, tmob

    !0.4706e-6, diff 1.1708e-05 1.326e-05
!    tmob = max(0.935656798937413 + 10.48743922905448/(1.456211069295241 + h)**5 - 7.45512854038917/(1.456211069295241 + h)**3 - 0.1666926153450568/(1.456211069295241 + h),0d0)
!   nmob = max(0.8784836695101248 + 16.178193275967942/(1.8846739602964013 + h)**3 - 13.450951287145715/(1.8846739602964013 + h)**2 + 0.9286418269991527/(1.8846739602964013 + h),0d0)

    !Single plane approximation
    nmob = max(1 - 9/(8*h) + 1/(2*h**3) - 1/(8*(h**5)),0d0)
    tmob = max(1 - 9/(16*h) + 2/(16*(h**3)) -1/(16*(h**5)),0d0)

    !1e-6, diff 1.1708e-05 1.326e-05
!    tmob = max(0.975240720746649 + 0.6618096094545308/(0.6399549837528644 + h)**3 - 1.261848923090074/(0.6399549837528644 + h)**2 -  0.2682891286963005/(0.6399549837528644 + h),0d0)
!    nmob = max(0.9707837307169467 + 6.987793458468398/(1.5507333020341336 + h)**3 - 6.037411430364329/(1.5507333020341336 + h)**2 -  0.5177318798737919/(1.5507333020341336 + h),0d0)

!     print *, h, nmob, tmob

  end subroutine mob_interp

  subroutine get_mobility_diff(nmob, tmob, part, z)

    real(amrex_real),intent(in   ) :: z
    real(amrex_real),intent(inout) :: nmob, tmob
    type(particle_t),intent(in   ), target :: part

    real(amrex_real) awet, atotal, hwet, htotal, tmobwet, nmobwet, tmobtotal, nmobtotal, h

    awet = k_b*t_init(1)/(part%wet_diff*visc_coef*3.142*6)
    atotal = k_b*t_init(1)/(part%total_diff*visc_coef*3.142*6)

    hwet = z/awet
    htotal = z/atotal

    call mob_interp(hwet, tmobwet, nmobwet)
    call mob_interp(htotal, tmobtotal, nmobtotal)

    tmob = max((tmobtotal*part%total_diff - tmobwet*part%wet_diff)/part%dry_diff,0.0)
    nmob = max((nmobtotal*part%total_diff - nmobwet*part%wet_diff)/part%dry_diff,0.0)

  end subroutine get_mobility_diff

  subroutine get_explicit_mobility(mob, part, plo, phi) &
                           bind(c,name="get_explicit_mobility")

    real(amrex_real),intent(in   ) :: plo(3), phi(3)
    real(amrex_real),intent(inout) :: mob(3)
    type(particle_t),intent(in   ), target :: part

    real(amrex_real) nmob, tmob, z

    mob(1:3) = 1

    if((bc_vel_lo(1) .eq. 2) .and. (bc_vel_hi(1) .eq. 2)) then

       z = part%pos(1)

       if(z .gt. (phi(1)-plo(1))/2.0) then
          z = phi(1) - z
       endif

       call get_mobility_diff(nmob, tmob, part, z)

       mob(1) = nmob
       mob(2) = tmob
#if (BL_SPACEDIM == 3)               
       mob(3) = tmob
#endif
    endif

    if((bc_vel_lo(2) .eq. 2) .and. (bc_vel_hi(2) .eq. 2)) then

       z = part%pos(2)

       if(z .gt. (phi(2)-plo(2))/2.0) then
          z = phi(2) - z
       endif

       call get_mobility_diff(nmob,tmob, part, z)

       mob(1) = tmob
       mob(2) = nmob
#if (BL_SPACEDIM == 3)               
       mob(3) = tmob
#endif
    endif


#if (BL_SPACEDIM == 3)               
    if((bc_vel_lo(3) .eq. 2) .and. (bc_vel_hi(3) .eq. 2)) then

       z = part%pos(3)

       if(z .gt. (phi(3)-plo(3))/2.0) then
          z = phi(3) - z
       endif

       call get_mobility_diff(nmob,tmob, part, z)

       mob(1) = tmob
       mob(2) = tmob
#if (BL_SPACEDIM == 3)               
       mob(3) = nmob
#endif

    endif
#endif

   ! print *, "z, nmob, pos: ", z, nmob, part%pos(3)

   !call sleep(100)
  end subroutine get_explicit_mobility

  subroutine compute_dry_mobility(lo, hi, mobility, mlo, mhi, dx, plo, phi, ngc, species) &
                                  bind(c,name="compute_dry_mobility")

    integer,          intent(in   )         :: lo(3), hi(3), mlo(3), mhi(3), ngc
    double precision, intent(in   )         :: plo(3), phi(3), dx(3)

    double precision, intent(inout)         :: mobility(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3),1:(nspecies*AMREX_SPACEDIM))
    type(species_t),  intent(in  ),  target :: species(nspecies)

    type(species_t), pointer :: spec

    integer :: i, j, k, l
    double precision :: z, nmob, tmob

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             do l = 1, nspecies

                mobility(i,j,k,(l-1)*AMREX_SPACEDIM + 1) = 1
                mobility(i,j,k,(l-1)*AMREX_SPACEDIM + 2) = 1
#if (BL_SPACEDIM == 3)               
                mobility(i,j,k,(l-1)*AMREX_SPACEDIM + 3) = 1
#endif
                spec => species(l)

                if((bc_vel_lo(1) .eq. 2) .and. (bc_vel_hi(1) .eq. 2)) then

                   z = dx(1)*(i+0.5)

                   if(z .gt. (phi(1)-plo(1))/2.0) then
                      z = phi(1) - z
                   endif

                   call get_mobility(nmob, tmob, spec, z)

                   mobility(i,j,k,(l-1)*AMREX_SPACEDIM + 1) = nmob
                   mobility(i,j,k,(l-1)*AMREX_SPACEDIM + 2) = tmob
#if (BL_SPACEDIM == 3)               
                   mobility(i,j,k,(l-1)*AMREX_SPACEDIM + 3) = tmob
#endif
                endif

                if((bc_vel_lo(2) .eq. 2) .and. (bc_vel_hi(2) .eq. 2)) then

                   z = dx(2)*(j+0.5)

                   if(z .gt. (phi(2)-plo(2))/2.0) then
                      z = phi(2) - z
                   endif

                   call get_mobility(nmob, tmob, spec, z)

                   mobility(i,j,k,(l-1)*AMREX_SPACEDIM + 1) = tmob
                   mobility(i,j,k,(l-1)*AMREX_SPACEDIM + 2) = nmob
#if (BL_SPACEDIM == 3)               
                   mobility(i,j,k,(l-1)*AMREX_SPACEDIM + 3) = tmob
#endif
                endif


#if (BL_SPACEDIM == 3)               
                if((bc_vel_lo(3) .eq. 2) .and. (bc_vel_hi(3) .eq. 2)) then

                   z = dx(3)*(k+0.5)

                   if(z .gt. (phi(3)-plo(3))/2.0) then
                      z = phi(3) - z
                   endif

                   call get_mobility(nmob, tmob, spec, z)

                   mobility(i,j,k,(l-1)*AMREX_SPACEDIM + 1) = tmob
                   mobility(i,j,k,(l-1)*AMREX_SPACEDIM + 2) = tmob

                   mobility(i,j,k,(l-1)*AMREX_SPACEDIM + 3) = nmob

                endif
#endif

             enddo

          end do
       end do
    end do

  end subroutine compute_dry_mobility

  subroutine sync_particles(spec3xPos, spec3yPos, spec3zPos, spec3xForce, spec3yForce, spec3zForce, &
                            particles, np, length) &
                            bind(c,name="sync_particles")

    integer,          intent(in   ) :: length, np
    double precision, intent(inout) :: spec3xPos(length), spec3yPos(length), spec3zPos(length), spec3xForce(length), spec3yForce(length), spec3zForce(length)

    type(particle_t), intent(inout), target :: particles(np)

    integer i, id
    type(particle_t), pointer :: part

    do i = 1, length

       spec3xPos(i) = 0
       spec3yPos(i) = 0
       spec3zPos(i) = 0
    enddo

    do i = 1, np

       part => particles(i)

       !if(part%species .eq. 1) then

       id = part%id

       spec3xPos(id) = part%pos(1)
       spec3yPos(id) = part%pos(2)
       spec3zPos(id) = part%pos(3)

       !print *, "id: ", id, "zpos: ", spec3zPos(id), "real: ", part%pos(3)

       !endif

    enddo

  end subroutine sync_particles

  subroutine force_particles(spec3xPos, spec3yPos, spec3zPos, spec3xForce, spec3yForce, spec3zForce, &
                             particles, np, length) &
                             bind(c,name="force_particles")

    integer,          intent(in   ) :: length, np
    double precision, intent(inout) :: spec3xPos(length), spec3yPos(length), spec3zPos(length), spec3xForce(length), spec3yForce(length), spec3zForce(length)

    type(particle_t), intent(inout), target :: particles(np)

    integer i, id
    type(particle_t), pointer :: part

    do i = 1, np

       part => particles(i)

       !if(part%species .eq. 1) then

       id = part%id

       part%force(1) = part%force(1) + spec3xForce(id)
       part%force(2) = part%force(2) + spec3yForce(id)
       part%force(3) = part%force(3) + spec3zForce(id)

       !endif

    enddo

  end subroutine force_particles

end module particle_functions_module
