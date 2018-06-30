  module collision_functions_module

  use amrex_error_module
  use amrex_fort_module
  use amrex_constants_module
  use common_namelist_module
  use rng_functions_module
  use cell_sorted_particle_module

contains

  subroutine init_cells(particles, lo, hi, cell_part_ids, cell_part_cnt, clo, chi, cellpairs, cplo, cphi, cellfactor, cflo, cfhi, cellvols, cvlo, cvhi, np, neff, cp, d, delt) bind(c,name='init_cells')

    use iso_c_binding, only: c_ptr, c_int, c_f_pointer

    implicit none

    integer,          intent(in      ) :: np, lo(3), hi(3), clo(3), chi(3), cplo(3), cphi(3), cflo(3), cfhi(3), cvlo(3), cvhi(3)
    double precision, intent(in      ) :: neff, delt, cp, d

    double precision, intent(inout   ) :: cellpairs(cplo(1):cphi(1),cplo(2):cphi(2),cplo(3):cphi(3))
    double precision, intent(inout   ) :: cellfactor(cplo(1):cphi(1),cplo(2):cphi(2),cplo(3):cphi(3))
    double precision, intent(inout   ) :: cellvols(cvlo(1):cvhi(1),cvlo(2):cvhi(2),cvlo(3):cvhi(3))


    type(c_ptr), intent(inout)      :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int), intent(inout)   :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))

    integer(c_int), pointer :: cell_parts(:)

    type(particle_t), intent(in), target :: particles(np)

    type(particle_t), pointer :: p1
    type(particle_t), pointer :: p2


    double precision fac
    integer i,j,k,l,m, cell_np

    !fac1 = delt*neff/4d0

    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)

          cell_np = cell_part_cnt(i,j,k)
          call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])

          if (cell_np .lt. 2) then

            !if there less than two particles in the cell, use the modal thermal speed

            cellfactor(i, j, k) = 3.14159265359*d*d*cp;

          else

            !fac2 = fac1*(cell_np**2)/cellvols(i,j,k)
            
            do l = 1, cell_np
              do m = 1, cell_np
              
                p1 => particles(cell_parts(l))
                p2 => particles(cell_parts(m))

                fac = sqrt((p1%vel(1)-p2%vel(1))**2 + (p1%vel(2)-p2%vel(2))**2 + (p1%vel(3)-p2%vel(3))**2)*3.14159265359*(p1%radius + p2%radius)**2

                if(fac .gt. cellfactor(i,j,k)) then

                  cellfactor(i,j,k) = fac

                endif         

              enddo
            enddo
          
          endif

        enddo
      enddo
    enddo

  end subroutine init_cells

  subroutine collide_cells(particles, lo, hi, cell_part_ids, cell_part_cnt, clo, chi, cellpairs, cplo, cphi, cellfactor, cflo, cfhi, cellvols, cvlo, cvhi, np, neff, cp, d, delt) bind(c,name='collide_cells')

    use iso_c_binding, only: c_ptr, c_int, c_f_pointer

    implicit none

    integer,          intent(in      ) :: np, lo(3), hi(3), clo(3), chi(3), cplo(3), cphi(3), cflo(3), cfhi(3), cvlo(3), cvhi(3)
    double precision, intent(in      ) :: neff, delt, cp, d

    double precision, intent(inout   ) :: cellpairs(cplo(1):cphi(1),cplo(2):cphi(2),cplo(3):cphi(3))
    double precision, intent(inout   ) :: cellfactor(cplo(1):cphi(1),cplo(2):cphi(2),cplo(3):cphi(3))
    double precision, intent(inout   ) :: cellvols(cvlo(1):cvhi(1),cvlo(2):cvhi(2),cvlo(3):cvhi(3))

    type(c_ptr),    intent(in)         :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int), intent(in)         :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))

    integer(c_int), pointer :: cell_parts(:)

    type(particle_t), intent(in), target :: particles(np)

    type(particle_t), pointer :: p1
    type(particle_t), pointer :: p2

    !Currently assuming single species, although mass is explicitly included for future extension
    !Go through this and optimise later


    !double precision fac1, fac2, fac3, test, pairfrac
    integer i,j,k,l,m,n,pairs, cell_np
    double precision fac1, fac2(3), fac3, sintheta, costheta, sinphi, cosphi, cmvel(3), totalmass

    fac3 = delt*neff/2d0

    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)

          cell_np = cell_part_cnt(i,j,k)
          call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])

          if (cell_np .gt. 1) then

            cellpairs(i,j,k) = fac3*cellfactor(i,j,k)*cell_np*cell_np/cellvols(i,j,k)

            pairs = floor(cellpairs(i,j,k))

            if( get_uniform_func() .lt. (cellpairs(i,j,k) - pairs) ) then

              pairs = pairs + 1;

            endif

            !print *, "Attempting ", pairs, " pairs. Cell factor: ", cellfactor(i,j,k)

            do n = 1, pairs

              call get_selector(l,cell_np)
              m = l;

              do while(m .eq. l)

                call get_selector(m,cell_np)

              enddo

              p1 => particles(cell_parts(l))
              p2 => particles(cell_parts(m))

              fac1 = sqrt((p1%vel(1)-p2%vel(1))**2 + (p1%vel(2)-p2%vel(2))**2 + (p1%vel(3)-p2%vel(3))**2)*(3.14159265359*(p1%radius + p2%radius)**2)

              if(fac1 .gt. cellfactor(i,j,k)) then

                cellfactor(i,j,k) = 2d0*fac1

                !print *, "Maxrel updated in cell ", i, j, k

              endif

              if (get_uniform_func()*cellfactor(i,j,k) .lt. fac1) then
                !print *, "Pre energy sqr: ", p1%vel(1)**2 + p1%vel(2)**2 + p1%vel(3)**2 + p2%vel(1)**2 + p2%vel(2)**2 + p2%vel(3)**2

                totalmass = p1%mass + p2%mass
                cmvel = (p1%mass*p1%vel + p2%mass*p2%vel)/totalmass

                call get_angles(costheta, sintheta, cosphi, sinphi)

                fac2 = (fac1/(3.14159265359*(p1%radius + p2%radius)**2))*(/sintheta*cosphi,sintheta*sinphi,costheta/)

                p1%vel = cmvel + p2%mass*fac2/totalmass
                p2%vel = cmvel - p1%mass*fac2/totalmass

                !print *, "Post energy sqr: ", p1%vel(1)**2 + p1%vel(2)**2 + p1%vel(3)**2 + p2%vel(1)**2 + p2%vel(2)**2 + p2%vel(3)**2

              endif
            
            enddo

          endif

        enddo
      enddo
    enddo

  end subroutine collide_cells

  subroutine evaluate_fields(particles, lo, hi, cell_part_ids, cell_part_cnt, clo, chi, members, mlo, mhi, density, dlo, dhi, velx, velxlo, velxhi, vely, velylo, velyhi, velz, velzlo, velzhi, &
                             temp, templo, temphi, speed, speedlo, speedhi, cellvols, cvlo, cvhi, neff, np) bind(c,name='evaluate_fields')


    use amrex_fort_module, only: amrex_real
    use iso_c_binding, only: c_ptr, c_int, c_f_pointer
    use cell_sorted_particle_module, only: particle_t

    implicit none

    integer,          intent(in      ) :: clo(3), chi(3), cvlo(3), cvhi(3), mlo(3), mhi(3), lo(3), hi(3), np
    integer,          intent(in      ) :: dlo(3), dhi(3), velxlo(3), velxhi(3), velylo(3), velyhi(3), velzlo(3), velzhi(3), templo(3), temphi(3), speedlo(3), speedhi(3)
    double precision, intent(in      ) :: neff

    double precision, intent(inout   ) :: cellvols(cvlo(1):cvhi(1),cvlo(2):cvhi(2),cvlo(3):cvhi(3))
    double precision, intent(inout   ) :: density(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
    double precision, intent(inout   ) :: velx(velxlo(1):velxhi(1),velxlo(2):velxhi(2),velxlo(3):velxhi(3))
    double precision, intent(inout   ) :: vely(velylo(1):velyhi(1),velylo(2):velyhi(2),velylo(3):velyhi(3))
    double precision, intent(inout   ) :: velz(velzlo(1):velzhi(1),velzlo(2):velzhi(2),velzlo(3):velzhi(3))
    double precision, intent(inout   ) :: temp(templo(1):temphi(1),templo(2):temphi(2),templo(3):temphi(3))
    double precision, intent(inout   ) :: speed(speedlo(1):speedhi(1),speedlo(2):speedhi(2),speedlo(3):speedhi(3))
    double precision, intent(inout   ) :: members(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    type(c_ptr), intent(inout)      :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int), intent(inout)   :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))

    type(particle_t), intent(inout), target :: particles(np)

    type(particle_t), pointer :: part
    integer(c_int), pointer :: cell_parts(:)


    !Go through this and optimise later

    !double precision fac1, fac2, fac3, test, pairfrac
    integer i,j,k,p,cell_np
    double precision membersinv

    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)

          members(i,j,k) = 0;
          density(i,j,k) = 0;
          velx(i,j,k) = 0;
          vely(i,j,k) = 0;
          velz(i,j,k) = 0;
          temp(i,j,k) = 0;

          cell_np = cell_part_cnt(i,j,k)
          call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])

          if(cell_np .ne. 0) then 
            membersinv = 1d0/cell_np
          else
            membersinv = 0d0
          endif
     

          members(i,j,k) = cell_np

          do p = 1, cell_np

            part => particles(cell_parts(p))

            density(i,j,k) = density(i,j,k) + part%mass

            velx(i,j,k) = velx(i,j,k) + part%vel(1)
            vely(i,j,k) = vely(i,j,k) + part%vel(2)
            velz(i,j,k) = velz(i,j,k) + part%vel(3)

          enddo

          density(i,j,k) = density(i,j,k)*neff/cellvols(i,j,k)
        
          velx(i,j,k) = velx(i,j,k)*membersinv
          vely(i,j,k) = vely(i,j,k)*membersinv
          velz(i,j,k) = velz(i,j,k)*membersinv

          do p = 1, cell_np

            part => particles(cell_parts(p))

            temp(i,j,k) = temp(i,j,k) + (1d0/part%R)*( (velx(i,j,k)-part%vel(1))**2 + (vely(i,j,k)-part%vel(2))**2 + (velz(i,j,k)-part%vel(3))**2 )

          enddo

          temp(i,j,k) = temp(i,j,k)*membersinv

        enddo
      enddo
    enddo

  end subroutine evaluate_fields





  subroutine evaluate_stats(particles, lo, hi, cell_part_ids, cell_part_cnt, clo, chi, &

                             members, mlo, mhi, &
                             density, dlo, dhi, & 
                             velx, velxlo, velxhi, & 
                             vely, velylo, velyhi, &
                             velz, velzlo, velzhi, &
                             temp, templo, temphi, & 
                          
                             membersmean, mmlo, mmhi, &
                             densitymean, dmlo, dmhi, &
                             velxmean, velxmlo, velxmhi, &
                             velymean, velymlo, velymhi, &
                             velzmean, velzmlo, velzmhi, &
                             tempmean, tempmlo, tempmhi, &
                             speedmean, speedmlo, speedmhi, &

                             membersvar, mvlo, mvhi, &
                             densityvar, dvlo, dvhi, &
                             velxvar, velxvlo, velxvhi, &
                             velyvar, velyvlo, velyvhi, &
                             velzvar, velzvlo, velzvhi, &
                             tempvar, tempvlo, tempvhi, &
                             speedvar, speedvlo, speedvhi, &

                             cellvols, cvlo, cvhi, np, neff, delt, steps) bind(c,name='evaluate_stats')

    use iso_c_binding, only: c_ptr, c_int, c_f_pointer
    use cell_sorted_particle_module, only: particle_t

    implicit none

    integer,          intent(in      ) :: np, steps, lo(3), hi(3), clo(3), chi(3), cvlo(3), cvhi(3)
    integer,          intent(in      ) :: mlo(3), mhi(3), dlo(3), dhi(3), velxlo(3), velxhi(3), velylo(3), velyhi(3), velzlo(3), velzhi(3), templo(3), temphi(3)
    integer,          intent(in      ) :: dmlo(3), dmhi(3), mmlo(3), mmhi(3), velxmlo(3), velxmhi(3), velymlo(3), velymhi(3), velzmlo(3), velzmhi(3), tempmlo(3), tempmhi(3), speedmlo(3), speedmhi(3)
    integer,          intent(in      ) :: dvlo(3), dvhi(3), mvlo(3), mvhi(3), velxvlo(3), velxvhi(3), velyvlo(3), velyvhi(3), velzvlo(3), velzvhi(3), tempvlo(3), tempvhi(3), speedvlo(3), speedvhi(3)
    double precision, intent(in      ) :: neff, delt

    double precision, intent(inout   ) :: cellvols(cvlo(1):cvhi(1),cvlo(2):cvhi(2),cvlo(3):cvhi(3))
    double precision, intent(inout   ) :: members(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
    double precision, intent(inout   ) :: density(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
    double precision, intent(inout   ) :: velx(velxlo(1):velxhi(1),velxlo(2):velxhi(2),velxlo(3):velxhi(3))
    double precision, intent(inout   ) :: vely(velylo(1):velyhi(1),velylo(2):velyhi(2),velylo(3):velyhi(3))
    double precision, intent(inout   ) :: velz(velzlo(1):velzhi(1),velzlo(2):velzhi(2),velzlo(3):velzhi(3))
    double precision, intent(inout   ) :: temp(templo(1):temphi(1),templo(2):temphi(2),templo(3):temphi(3))

    double precision, intent(inout   ) :: densitymean(dmlo(1):dmhi(1),dmlo(2):dmhi(2),dmlo(3):dmhi(3))
    double precision, intent(inout   ) :: membersmean(mmlo(1):mmhi(1),mmlo(2):mmhi(2),mmlo(3):mmhi(3))
    double precision, intent(inout   ) :: velxmean(velxmlo(1):velxmhi(1),velxmlo(2):velxmhi(2),velxmlo(3):velxmhi(3))
    double precision, intent(inout   ) :: velymean(velymlo(1):velymhi(1),velymlo(2):velymhi(2),velymlo(3):velymhi(3))
    double precision, intent(inout   ) :: velzmean(velzmlo(1):velzmhi(1),velzmlo(2):velzmhi(2),velzmlo(3):velzmhi(3))
    double precision, intent(inout   ) :: tempmean(tempmlo(1):tempmhi(1),tempmlo(2):tempmhi(2),tempmlo(3):tempmhi(3))
    double precision, intent(inout   ) :: speedmean(speedmlo(1):speedmhi(1),speedmlo(2):speedmhi(2),speedmlo(3):speedmhi(3))

    double precision, intent(inout   ) :: densityvar(dvlo(1):dvhi(1),dvlo(2):dvhi(2),dvlo(3):dvhi(3))
    double precision, intent(inout   ) :: membersvar(mvlo(1):mvhi(1),mvlo(2):mvhi(2),mvlo(3):mvhi(3))
    double precision, intent(inout   ) :: velxvar(velxvlo(1):velxvhi(1),velxvlo(2):velxvhi(2),velxvlo(3):velxvhi(3))
    double precision, intent(inout   ) :: velyvar(velyvlo(1):velyvhi(1),velyvlo(2):velyvhi(2),velyvlo(3):velyvhi(3))
    double precision, intent(inout   ) :: velzvar(velzvlo(1):velzvhi(1),velzvlo(2):velzvhi(2),velzvlo(3):velzvhi(3))
    double precision, intent(inout   ) :: tempvar(tempvlo(1):tempvhi(1),tempvlo(2):tempvhi(2),tempvlo(3):tempvhi(3))
    double precision, intent(inout   ) :: speedvar(speedvlo(1):speedvhi(1),speedvlo(2):speedvhi(2),speedvlo(3):speedvhi(3))

    type(c_ptr), intent(inout)      :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int), intent(inout)   :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))

    type(particle_t), intent(inout), target :: particles(np)

    !type(particle_t), pointer :: part
    !integer(c_int), pointer :: cell_parts(:)


    !Go through this and optimise later

    !double precision fac1, fac2, fac3, test, pairfrac
    integer i,j,k
    double precision stepsminusone, stepsinv

    stepsminusone = steps - 1
    stepsinv = 1d0/steps

    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)

          velxmean(i,j,k) = (velxmean(i,j,k)*stepsminusone + velx(i,j,k))*stepsinv
          velymean(i,j,k) = (velymean(i,j,k)*stepsminusone + vely(i,j,k))*stepsinv
          velzmean(i,j,k) = (velzmean(i,j,k)*stepsminusone + velz(i,j,k))*stepsinv

          densitymean(i,j,k) = (densitymean(i,j,k)*stepsminusone + density(i,j,k))*stepsinv
          tempmean(i,j,k) = (tempmean(i,j,k)*stepsminusone + temp(i,j,k))*stepsinv

          membersmean(i,j,k) = (membersmean(i,j,k)*stepsminusone + members(i,j,k))*stepsinv

        enddo
      enddo
    enddo

    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)

          speedvar(i,j,k) = (speedvar(i,j,k)*stepsminusone + (velx(i,j,k) - velxmean(i,j,k))**2 + (vely(i,j,k) - velymean(i,j,k))**2 + (velz(i,j,k) - velzmean(i,j,k))**2)/steps

        enddo
      enddo
    enddo

    if((lo(1) .eq. 0) .and. (lo(2) .eq. 0) .and. (lo(3) .eq. 0)) then

      !print *, speedvar(cmlo(1),cmlo(2),cmlo(3)), particles(1)%R*tempmean(cmlo(1),cmlo(2),cmlo(3))/(membersmean(cmlo(1),cmlo(2),cmlo(3)))

    endif

  end subroutine evaluate_stats
  
end module collision_functions_module

