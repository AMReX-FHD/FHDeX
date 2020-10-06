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
    !type(particle_t), pointer :: p2

    double precision fac
    integer i,j,k,l,cell_np
    !integer m

    !fac1 = delt*neff/4d0

    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)

          cell_np = cell_part_cnt(i,j,k)
          call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])	

          if (cell_np .lt. 2) then

            !if there less than two particles in the cell, use the modal thermal speed

            cellfactor(i, j, k) = 3.14159265359*d*d*cp;

            !print *, "update cell ", i,j,k

          else

            !fac2 = fac1*(cell_np**2)/cellvols(i,j,k)
            
            !do l = 1, cell_np
            !  do m = 1, cell_np
            !  
            !    p1 => particles(cell_parts(l))
            !    p2 => particles(cell_parts(m))
            !
            !    fac = sqrt((p1%vel(1)-p2%vel(1))**2 + (p1%vel(2)-p2%vel(2))**2 + (p1%vel(3)-p2%vel(3))**2)*3.14159265359*(p1%radius + p2%radius)**2

            !    if(fac .gt. cellfactor(i,j,k)) then

             !     cellfactor(i,j,k) = 2d0*fac

               ! endif         

           !   enddo
         !   enddo


            do l = 1, cell_np
             
              
                p1 => particles(cell_parts(l))

                fac = 2d0*sqrt((p1%vel(1))**2 + (p1%vel(2))**2 + (p1%vel(3))**2)*3.14159265359*(2*p1%radius)**2

                if(fac .gt. cellfactor(i,j,k)) then

                  cellfactor(i,j,k) = 1.5*fac

                endif                      
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

            

          !print *, "parts: ", cell_np
          call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])

          if (cell_np .gt. 1) then

            cellpairs(i,j,k) = fac3*cellfactor(i,j,k)*cell_np*cell_np/cellvols(i,j,k)

            pairs = floor(cellpairs(i,j,k))

            if( get_uniform_func() .lt. (cellpairs(i,j,k) - pairs) ) then

              pairs = pairs + 1;

            endif

            !pairs = floor(pairs*0.1)

            !print *, "Attempting ", pairs, " pairs. Cell factor: ", cellfactor(i,j,k), cellvols(i,j,k)

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

                cellfactor(i,j,k) = 1.5d0*fac1

                print *, "Maxrel updated in cell ", i, j, k

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

  subroutine evaluate_fields(particles, lo, hi, cell_part_ids, cell_part_cnt, clo, chi, instant, ilo, ihi, cellvols, cvlo, cvhi, neff, np) bind(c,name='evaluate_fields')


    use amrex_fort_module, only: amrex_real
    use iso_c_binding, only: c_ptr, c_int, c_f_pointer
    use cell_sorted_particle_module, only: particle_t

    implicit none

    integer,          intent(in      ) :: clo(3), chi(3), cvlo(3), cvhi(3), ilo(3), ihi(3), lo(3), hi(3), np
    double precision, intent(in      ) :: neff

    double precision, intent(inout   ) :: cellvols(cvlo(1):cvhi(1),cvlo(2):cvhi(2),cvlo(3):cvhi(3))
    double precision, intent(inout   ) :: instant(ilo(1):ihi(1),ilo(2):ihi(2),ilo(3):ihi(3),14)

    type(c_ptr), intent(inout)      :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int), intent(inout)   :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))

    type(particle_t), intent(inout), target :: particles(np)

    type(particle_t), pointer :: part
    integer(c_int), pointer :: cell_parts(:)


    !Go through this and optimise later

    !double precision fac1, fac2, fac3, test, pairfrac
    integer i,j,k,p,cell_np, ti, tj, tk, totalparticles
    double precision membersinv, nrg, totalpx, totalpy, totalpz, rmean

    totalpx = 0
    totalpy = 0
    totalpz = 0
    totalparticles = 0

    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)

          cell_np = cell_part_cnt(i,j,k)
          call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])

          if(cell_np .ne. 0) then 
            membersinv = 1d0/cell_np
          else
            membersinv = 0d0
          endif
    
          instant(i,j,k,1) = cell_np

 

          do p = 1, cell_np

            part => particles(cell_parts(p))
            
            instant(i,j,k,2) = instant(i,j,k,2) + part%mass

            instant(i,j,k,3) = instant(i,j,k,3) + part%vel(1)
            instant(i,j,k,4) = instant(i,j,k,4) + part%vel(2)
            instant(i,j,k,5) = instant(i,j,k,5) + part%vel(3)

            instant(i,j,k,7) = instant(i,j,k,7) + part%vel(1)*part%mass
            instant(i,j,k,8) = instant(i,j,k,8) + part%vel(2)*part%mass
            instant(i,j,k,9) = instant(i,j,k,9) + part%vel(3)*part%mass

            !if(part%q .gt. 0) then

              instant(i,j,k,12) = instant(i,j,k,12) + part%vel(1)*part%q
              instant(i,j,k,13) = instant(i,j,k,13) + part%vel(2)*part%q
              instant(i,j,k,14) = instant(i,j,k,14) + part%vel(3)*part%q

            !endif

            nrg = 0.5*part%mass*(part%vel(1)*part%vel(1) + part%vel(2)*part%vel(2) + part%vel(3)*part%vel(3))

            instant(i,j,k,10) = instant(i,j,k,10) + nrg

            totalparticles = totalparticles + 1;

            totalpx = totalpx +  part%vel(1)
            totalpy = totalpy +  part%vel(2)
            totalpz = totalpz +  part%vel(3)

          enddo

          instant(i,j,k,2) = instant(i,j,k,2)*neff/cellvols(i,j,k)
   
          instant(i,j,k,3) = instant(i,j,k,3)*membersinv
          instant(i,j,k,4) = instant(i,j,k,4)*membersinv
          instant(i,j,k,5) = instant(i,j,k,5)*membersinv

          instant(i,j,k,12) = instant(i,j,k,12)*neff/cellvols(i,j,k)
          instant(i,j,k,13) = instant(i,j,k,13)*neff/cellvols(i,j,k)
          instant(i,j,k,14) = instant(i,j,k,14)*neff/cellvols(i,j,k)

          instant(i,j,k,7) = instant(i,j,k,7)*neff/cellvols(i,j,k)
          instant(i,j,k,8) = instant(i,j,k,8)*neff/cellvols(i,j,k)
          instant(i,j,k,9) = instant(i,j,k,9)*neff/cellvols(i,j,k)

          instant(i,j,k,10) = instant(i,j,k,10)*neff/cellvols(i,j,k)

          rmean = 0

          do p = 1, cell_np

            part => particles(cell_parts(p))

            rmean = rmean + part%R

            instant(i,j,k,6) = instant(i,j,k,6) + &
                 (1d0/part%R)*( (instant(i,j,k,3)-part%vel(1))**2 &
                               +(instant(i,j,k,4)-part%vel(2))**2 &
                               +(instant(i,j,k,5)-part%vel(3))**2 )

          enddo

          instant(i,j,k,6) = instant(i,j,k,6)*membersinv*0.33333333333333333

          instant(i,j,k,11) = instant(i,j,k,2)*rmean*instant(i,j,k,5)

        enddo
      enddo
    enddo

  end subroutine evaluate_fields

  subroutine evaluate_fields_pp(particles, lo, hi, cell_part_ids, cell_part_cnt, clo, chi, instant, ilo, ihi, cellvols, cvlo, cvhi, neff, np, dx) bind(c,name='evaluate_fields_pp')


    use amrex_fort_module, only: amrex_real
    use iso_c_binding, only: c_ptr, c_int, c_f_pointer
    use cell_sorted_particle_module, only: particle_t

    implicit none

    integer,          intent(in      ) :: clo(3), chi(3), cvlo(3), cvhi(3), ilo(3), ihi(3), lo(3), hi(3), np
    double precision, intent(in      ) :: neff, dx(3)

    double precision, intent(inout   ) :: cellvols(cvlo(1):cvhi(1),cvlo(2):cvhi(2),cvlo(3):cvhi(3))
    double precision, intent(inout   ) :: instant(ilo(1):ihi(1),ilo(2):ihi(2),ilo(3):ihi(3),14)

    type(c_ptr), intent(inout)      :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int), intent(inout)   :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))

    type(particle_t), intent(inout), target :: particles(np)

    type(particle_t), pointer :: part
    integer(c_int), pointer :: cell_parts(:)


    !Go through this and optimise later

    !double precision fac1, fac2, fac3, test, pairfrac
    integer i,j,k,p,cell_np, ti, tj, tk, totalparticles
    double precision membersinv, nrg, totalpx, totalpy, totalpz, rmean, dxinv(3)

    totalpx = 0
    totalpy = 0
    totalpz = 0
    totalparticles = 0

    dxinv=1d0/dx

    do p = 1, np

      part => particles(p)

      i = floor(part%pos(1)*dxinv(1))
      j = floor(part%pos(2)*dxinv(2))
      k = floor(part%pos(3)*dxinv(3))

      instant(i,j,k,1) = instant(i,j,k,1) + 1

      instant(i,j,k,2) = instant(i,j,k,2) + part%mass

      instant(i,j,k,3) = instant(i,j,k,3) + part%vel(1)
      instant(i,j,k,4) = instant(i,j,k,4) + part%vel(2)
      instant(i,j,k,5) = instant(i,j,k,5) + part%vel(3)

      instant(i,j,k,7) = instant(i,j,k,7) + part%vel(1)*part%mass
      instant(i,j,k,8) = instant(i,j,k,8) + part%vel(2)*part%mass
      instant(i,j,k,9) = instant(i,j,k,9) + part%vel(3)*part%mass

      !if(part%q .gt. 0) then

        instant(i,j,k,12) = instant(i,j,k,12) + part%vel(1)*part%q
        instant(i,j,k,13) = instant(i,j,k,13) + part%vel(2)*part%q
        instant(i,j,k,14) = instant(i,j,k,14) + part%vel(3)*part%q

      !endif

      nrg = 0.5*part%mass*(part%vel(1)*part%vel(1) + part%vel(2)*part%vel(2) + part%vel(3)*part%vel(3))

      instant(i,j,k,10) = instant(i,j,k,10) + nrg

      instant(i,j,k,11) = instant(i,j,k,11) + part%R

      totalparticles = totalparticles + 1;

      totalpx = totalpx +  part%vel(1)
      totalpy = totalpy +  part%vel(2)
      totalpz = totalpz +  part%vel(3)
        
    enddo

    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)

          membersinv = 1d0/instant(i,j,k,1)

          instant(i,j,k,2) = instant(i,j,k,2)*neff/cellvols(i,j,k)
   
          instant(i,j,k,3) = instant(i,j,k,3)*membersinv
          instant(i,j,k,4) = instant(i,j,k,4)*membersinv
          instant(i,j,k,5) = instant(i,j,k,5)*membersinv

          instant(i,j,k,12) = instant(i,j,k,12)*neff/cellvols(i,j,k)
          instant(i,j,k,13) = instant(i,j,k,13)*neff/cellvols(i,j,k)
          instant(i,j,k,14) = instant(i,j,k,14)*neff/cellvols(i,j,k)

          instant(i,j,k,7) = instant(i,j,k,7)*neff/cellvols(i,j,k)
          instant(i,j,k,8) = instant(i,j,k,8)*neff/cellvols(i,j,k)
          instant(i,j,k,9) = instant(i,j,k,9)*neff/cellvols(i,j,k)

          instant(i,j,k,10) = instant(i,j,k,10)*neff/cellvols(i,j,k)

          instant(i,j,k,6) = instant(i,j,k,6)*membersinv*0.33333333333333333

          instant(i,j,k,11) = instant(i,j,k,11)*instant(i,j,k,2)*instant(i,j,k,5)*membersinv

        enddo
      enddo
    enddo

    do p = 1, np

      part => particles(p)

      i = floor(part%pos(1)*dxinv(1))
      j = floor(part%pos(2)*dxinv(2))
      k = floor(part%pos(3)*dxinv(3))

      instant(i,j,k,6) = instant(i,j,k,6) + &
           (1d0/part%R)*( (instant(i,j,k,3)-part%vel(1))**2 &
                         +(instant(i,j,k,4)-part%vel(2))**2 &
                         +(instant(i,j,k,5)-part%vel(3))**2 )


    enddo

  end subroutine evaluate_fields_pp

  subroutine evaluate_means(lo, hi, &
                            instant, ilo, ihi, &
                            means, mlo, mhi, & 
                            vars, vlo, vhi, & 
                            cellvols, cvlo, cvhi, &
                            n0, T0, delt, steps, &
                            avcurrent) bind(c,name='evaluate_means')

    use iso_c_binding, only: c_ptr, c_int, c_f_pointer
    use cell_sorted_particle_module, only: particle_t

    implicit none

    integer,          intent(in      ) :: steps, lo(3), hi(3), cvlo(3), cvhi(3), ilo(3), ihi(3), mlo(3), mhi(3), vlo(3), vhi(3)
    double precision, intent(in      ) :: delt, n0, T0

    double precision, intent(inout   ) :: cellvols(cvlo(1):cvhi(1),cvlo(2):cvhi(2),cvlo(3):cvhi(3))
    double precision, intent(inout   ) :: instant(ilo(1):ihi(1),ilo(2):ihi(2),ilo(3):ihi(3),14)
    double precision, intent(inout   ) :: means(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3),14)
    double precision, intent(inout   ) :: vars(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),18)

    double precision, intent(inout) :: avcurrent(1:3)

    !double precision fac1, fac2, fac3, test, pairfrac
    integer i,j,k, ti, tj, tk, kc, jc
    double precision stepsminusone, stepsinv, cv, cvinv, delg, qmean, delpx, delpy, delpz, delrho, delvelx, delvely, delvelz, delenergy, densitymeaninv

    stepsminusone = steps - 1
    stepsinv = 1d0/steps

    cvinv = 2.0/(3.0)
    cv = 1.0/cvinv

    do k = mlo(3), mhi(3)
      do j = mlo(2), mhi(2)
        do i = mlo(1), mhi(1)

          !Means      

          means(i,j,k,1) = (means(i,j,k,1)*stepsminusone + instant(i,j,k,1))*stepsinv !member density
          means(i,j,k,2) = (means(i,j,k,2)*stepsminusone + instant(i,j,k,2))*stepsinv !mass density

!          if(means(i,j,k,1) .ne. 0) then
!            print *, i, j, k, means(i,j,k,1)
!          endif

          !densitymean(i,j,k) = 1.79233
          if(means(i,j,k,2) .gt. 0) then
            densitymeaninv = 1.0/means(i,j,k,2)
          else
            densitymeaninv = 0
          endif

          means(i,j,k,7) = (means(i,j,k,7)*stepsminusone + instant(i,j,k,7))*stepsinv !momentum density
          means(i,j,k,8) = (means(i,j,k,8)*stepsminusone + instant(i,j,k,8))*stepsinv
          means(i,j,k,9) = (means(i,j,k,9)*stepsminusone + instant(i,j,k,9))*stepsinv

          means(i,j,k,12) = (means(i,j,k,12)*stepsminusone + instant(i,j,k,12))*stepsinv !total current
          means(i,j,k,13) = (means(i,j,k,13)*stepsminusone + instant(i,j,k,13))*stepsinv
          means(i,j,k,14) = (means(i,j,k,14)*stepsminusone + instant(i,j,k,14))*stepsinv

          

          !pxmean(i,j,k) = 0
          !pymean(i,j,k) = 0
          !pzmean(i,j,k) = 0

          means(i,j,k,10) = (means(i,j,k,10)*stepsminusone + instant(i,j,k,10))*stepsinv
          !energymean(i,j,k) = 152835.9985

          means(i,j,k,3) = means(i,j,k,7)*densitymeaninv   !velocity
          means(i,j,k,4) = means(i,j,k,8)*densitymeaninv
          means(i,j,k,5) = means(i,j,k,9)*densitymeaninv

          means(i,j,k,6) = cvinv*densitymeaninv*(means(i,j,k,10) - 0.5*densitymeaninv*(means(i,j,k,7)*means(i,j,k,7) + means(i,j,k,8)*means(i,j,k,8) + means(i,j,k,9)*means(i,j,k,9)) )  !temperature - wrong for multispec

          means(i,j,k,1) = (means(i,j,k,1)*stepsminusone + means(i,j,k,1))*stepsinv !members

          means(i,j,k,11) = cvinv*(means(i,j,k,10) -0.5*densitymeaninv*(means(i,j,k,7)*means(i,j,k,7) + means(i,j,k,8)*means(i,j,k,8) + means(i,j,k,9)*means(i,j,k,9))  ) !pressure - wrong for multispec

        enddo
      enddo
    enddo

    avcurrent = 0

    do k = mlo(3), mhi(3)
      do j = mlo(2), mhi(2)
        do i = mlo(1), mhi(1)
           avcurrent(1) = avcurrent(1) + means(i,j,k,12)
           avcurrent(2) = avcurrent(2) + means(i,j,k,13)
           avcurrent(3) = avcurrent(3) + means(i,j,k,14)
        enddo
      enddo
    enddo

  end subroutine evaluate_means

  subroutine evaluate_corrs(particles, lo, hi, cell_part_ids, cell_part_cnt, clo, chi, &
                            instant, ilo, ihi, &
                            means, mlo, mhi, & 
                            vars, vlo, vhi, & 
                            cellvols, cvlo, cvhi, np, neff, n0, T0,delt, steps, &
                            varcurrent) bind(c,name='evaluate_corrs')

    use iso_c_binding, only: c_ptr, c_int, c_f_pointer
    use cell_sorted_particle_module, only: particle_t

    implicit none

    integer,          intent(in      ) :: np, steps, lo(3), hi(3), clo(3), chi(3), cvlo(3), cvhi(3), ilo(3), ihi(3), mlo(3), mhi(3), vlo(3), vhi(3)
    double precision, intent(in      ) :: neff, delt, n0, T0

    double precision, intent(inout   ) :: cellvols(cvlo(1):cvhi(1),cvlo(2):cvhi(2),cvlo(3):cvhi(3))
    double precision, intent(inout   ) :: instant(ilo(1):ihi(1),ilo(2):ihi(2),ilo(3):ihi(3),14)
    double precision, intent(inout   ) :: means(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3),14)
    double precision, intent(inout   ) :: vars(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),18)

    type(c_ptr), intent(inout)      :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int), intent(inout)   :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))

    type(particle_t), intent(inout), target :: particles(np)

    double precision, intent(inout) :: varcurrent(1:3)

    !Go through this and optimise later

    !double precision fac1, fac2, fac3, test, pairfrac
    integer i,j,k, it, jt, kt
    double precision stepsminusone, stepsinv, lhs, rhs, tempcm, cv, ncm, velcm, momentumcm, cvinv
    double precision delg, qmean, delpx, delpy, delpz, delrho, delvelx, delvely, delvelz, delenergy
    double precision densitymeaninv, deltemp, delix, deliy, deliz

    stepsminusone = steps - 1
    stepsinv = 1d0/steps

    cvinv = 2.0/(3.0*particles(1)%R)
    cv = 1.0/cvinv


    do k = mlo(3), mhi(3)
      do j = mlo(2), mhi(2)
        do i = mlo(1), mhi(1)

          !Vars
          qmean = cv*means(i,j,k,6)-0.5*(means(i,j,k,3)**2 + means(i,j,k,4)**2 + means(i,j,k,5)**2)

          if(means(i,j,k,2) .gt. 0) then
            densitymeaninv = 1.0/means(i,j,k,2)
          else
            densitymeaninv = 0
          endif

          delrho = instant(i,j,k,2) - means(i,j,k,2)

          delpx = instant(i,j,k,7) - means(i,j,k,7)
          delpy = instant(i,j,k,8) - means(i,j,k,8)
          delpz = instant(i,j,k,9) - means(i,j,k,9)

          delix = instant(i,j,k,12) - means(i,j,k,12)
          deliy = instant(i,j,k,13) - means(i,j,k,13)
          deliz = instant(i,j,k,14) - means(i,j,k,14)

          delenergy = instant(i,j,k,10) - means(i,j,k,10)

          delvelx = (delpx - means(i,j,k,3)*delrho)*densitymeaninv
          delvely = (delpy - means(i,j,k,4)*delrho)*densitymeaninv
          delvelz = (delpz - means(i,j,k,5)*delrho)*densitymeaninv

          vars(i,j,k,2) = (vars(i,j,k,2)*stepsminusone + delrho**2)*stepsinv

          vars(i,j,k,7) = (vars(i,j,k,7)*stepsminusone + delpx**2)*stepsinv
          vars(i,j,k,8) = (vars(i,j,k,8)*stepsminusone + delpy**2)*stepsinv
          vars(i,j,k,9) = (vars(i,j,k,9)*stepsminusone + delpz**2)*stepsinv

          vars(i,j,k,10) = (vars(i,j,k,10)*stepsminusone + delenergy**2)*stepsinv      

          vars(i,j,k,3) = (vars(i,j,k,3)*stepsminusone + delvelx**2)*stepsinv
          vars(i,j,k,4) = (vars(i,j,k,4)*stepsminusone + delvely**2)*stepsinv
          vars(i,j,k,5) = (vars(i,j,k,5)*stepsminusone + delvelz**2)*stepsinv

          delg = means(i,j,k,3)*delpx + means(i,j,k,4)*delpy + means(i,j,k,5)*delpz
          !delg = 0

          vars(i,j,k,12) = (vars(i,j,k,12)*stepsminusone + delg**2)*stepsinv

          vars(i,j,k,13) = (vars(i,j,k,13)*stepsminusone + delg*delenergy)*stepsinv
          vars(i,j,k,14) = (vars(i,j,k,14)*stepsminusone + delrho*delenergy)*stepsinv
          vars(i,j,k,15) = (vars(i,j,k,15)*stepsminusone + delrho*delg)*stepsinv

          vars(i,j,k,6) = (vars(i,j,k,6)*stepsminusone + cvinv*cvinv*densitymeaninv*densitymeaninv*(vars(i,j,k,10) + vars(i,j,k,12) - 2*vars(i,j,k,13) + qmean*(qmean*vars(i,j,k,2) - 2*vars(i,j,k,14) + 2*vars(i,j,k,15))))*stepsinv

          deltemp = (delenergy - delg - qmean*delrho)*cvinv*densitymeaninv

          vars(i,j,k,16) = (vars(i,j,k,16)*stepsminusone + delix**2)*stepsinv
          vars(i,j,k,17) = (vars(i,j,k,17)*stepsminusone + deliy**2)*stepsinv
          vars(i,j,k,18) = (vars(i,j,k,18)*stepsminusone + deliz**2)*stepsinv



          !print *, delrho
        enddo
      enddo
    enddo

    varcurrent = 0

    do k = vlo(3), vhi(3)
      do j = vlo(2), vhi(2)
        do i = vlo(1), vhi(1)
            varcurrent(1) = varcurrent(1) + vars(i,j,k,16)
            varcurrent(2) = varcurrent(2) + vars(i,j,k,17)
            varcurrent(3) = varcurrent(3) + vars(i,j,k,18)
         enddo
      enddo
    enddo

  end subroutine evaluate_corrs


  subroutine initialize_fields(particles, lo, hi, cell_part_ids, cell_part_cnt, clo, chi, &
                               instant, plo, phi, cellvols, cvlo, cvhi, neff, np, r, t) bind(c,name='initialize_fields')


    use amrex_fort_module, only: amrex_real
    use iso_c_binding, only: c_ptr, c_int, c_f_pointer
    use cell_sorted_particle_module, only: particle_t

    implicit none

    integer,          intent(in      ) :: clo(3), chi(3), cvlo(3), cvhi(3), plo(3), phi(3), lo(3), hi(3), np
    double precision, intent(in      ) :: neff, r, t

    double precision, intent(inout   ) :: cellvols(cvlo(1):cvhi(1),cvlo(2):cvhi(2),cvlo(3):cvhi(3))
    double precision, intent(inout   ) :: instant(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),11)

    type(c_ptr), intent(inout)      :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int), intent(inout)   :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))

    type(particle_t), intent(inout), target :: particles(np)

    type(particle_t), pointer :: part
    integer(c_int), pointer :: cell_parts(:)


    !double precision fac1, fac2, fac3, test, pairfrac
    integer i,j,k,p,cell_np
    double precision membersinv, varx, vary, varz, ratiox, ratioy, ratioz, totalenergy, totalvel(3)

    totalenergy = 0
    totalvel(1) = 0
    totalvel(2) = 0
    totalvel(3) = 0

    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)

          instant(i,j,k,1) = 0; !Members
          instant(i,j,k,2) = 0; !Density
          instant(i,j,k,3) = 0; !velx
          instant(i,j,k,4) = 0; !vely
          instant(i,j,k,5) = 0; !velz
          instant(i,j,k,6) = 0; !Temp

          cell_np = cell_part_cnt(i,j,k)
          call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])

          if(cell_np .ne. 0) then 
            membersinv = 1d0/cell_np
          else
            membersinv = 0d0
          endif

          instant(i,j,k,1) = cell_np

          do p = 1, cell_np

            part => particles(cell_parts(p))
            
            instant(i,j,k,3) = instant(i,j,k,3) + part%vel(1)
            instant(i,j,k,4) = instant(i,j,k,4) + part%vel(2)
            instant(i,j,k,5) = instant(i,j,k,5) + part%vel(3)

          enddo

          instant(i,j,k,3) = instant(i,j,k,3)*membersinv
          instant(i,j,k,4) = instant(i,j,k,4)*membersinv
          instant(i,j,k,5) = instant(i,j,k,5)*membersinv

          do p = 1, cell_np

            part => particles(cell_parts(p))

            part%vel(1) = part%vel(1) - instant(i,j,k,3)
            part%vel(2) = part%vel(2) - instant(i,j,k,4)
            part%vel(3) = part%vel(3) - instant(i,j,k,5)

          enddo

          instant(i,j,k,3) = 0;
          instant(i,j,k,4) = 0;
          instant(i,j,k,5) = 0;

          do p = 1, cell_np

            part => particles(cell_parts(p))

            instant(i,j,k,3) = instant(i,j,k,3) + part%vel(1)
            instant(i,j,k,4) = instant(i,j,k,4) + part%vel(2)
            instant(i,j,k,5) = instant(i,j,k,5) + part%vel(3)

          enddo

          instant(i,j,k,3) = instant(i,j,k,3)*membersinv
          instant(i,j,k,4) = instant(i,j,k,4)*membersinv
          instant(i,j,k,5) = instant(i,j,k,5)*membersinv

          varx = 0
          vary = 0
          varz = 0

          do p = 1, cell_np

            part => particles(cell_parts(p))
            varx = varx + (part%vel(1))**2
            vary = vary + (part%vel(2))**2
            varz = varz + (part%vel(3))**2

          enddo

          
          varx = varx*membersinv
          vary = vary*membersinv
          varz = varz*membersinv

          if(cell_np .gt. 1) then 
              ratiox = sqrt(r*t/varx)
              ratioy = sqrt(r*t/vary)
              ratioz = sqrt(r*t/varz)
          else
              ratiox = 1
              ratioy = 1
              ratioz = 1
          endif 

          do p = 1, cell_np

            part => particles(cell_parts(p))

            part%vel(1) = part%vel(1)*ratiox
            part%vel(2) = part%vel(2)*ratioy
            part%vel(3) = part%vel(3)*ratioz

            totalenergy = totalenergy + part%vel(1)**2 + part%vel(2)**2 + part%vel(3)**2

            totalvel = totalvel + part%vel

          enddo

          instant(i,j,k,3) = 0;
          instant(i,j,k,4) = 0;
          instant(i,j,k,5) = 0;

          do p = 1, cell_np

            part => particles(cell_parts(p))

            !density(i,j,k) = density(i,j,k) + part%mass

            instant(i,j,k,3) = instant(i,j,k,3) + part%vel(1)
            instant(i,j,k,4) = instant(i,j,k,4) + part%vel(2)
            instant(i,j,k,5) = instant(i,j,k,5) + part%vel(3)

          enddo

          instant(i,j,k,3) = instant(i,j,k,3)*membersinv
          instant(i,j,k,4) = instant(i,j,k,4)*membersinv
          instant(i,j,k,5) = instant(i,j,k,5)*membersinv


        enddo
      enddo
    enddo

    print *, "Corrected energy: ", totalenergy
    print *, "Corrected vel: ", totalvel

  end subroutine initialize_fields
  
end module collision_functions_module

