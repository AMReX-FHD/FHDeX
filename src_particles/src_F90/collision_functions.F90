  module collision_functions_module

  use amrex_error_module
  use amrex_fort_module
  use amrex_constants_module
  use common_namelist_module
  use rng_functions_module
  use particle_functions_module

contains

  subroutine init_cells(particles, cellmembers, cmlo, cmhi, celllists, cllo, clhi, cellpairs, cplo, cphi, cellfactor, cflo, cfhi, cellvols, cvlo, cvhi, ppc, np, neff, cp, d, delt) bind(c,name='init_cells')

    implicit none

    integer,          intent(in      ) :: np, ppc, cmlo(3), cmhi(3), cllo(3), clhi(3), cplo(3), cphi(3), cflo(3), cfhi(3), cvlo(3), cvhi(3)
    double precision, intent(in      ) :: neff, delt, cp, d

    integer         , intent(in      ) :: cellmembers(cmlo(1):cmhi(1),cmlo(2):cmhi(2),cmlo(3):cmhi(3))
    integer         , intent(in      ) :: celllists(cmlo(1):cmhi(1),cmlo(2):cmhi(2),cmlo(3):cmhi(3),1:ppc)
    double precision, intent(inout   ) :: cellpairs(cplo(1):cphi(1),cplo(2):cphi(2),cplo(3):cphi(3))
    double precision, intent(inout   ) :: cellfactor(cplo(1):cphi(1),cplo(2):cphi(2),cplo(3):cphi(3))
    double precision, intent(inout   ) :: cellvols(cvlo(1):cvhi(1),cvlo(2):cvhi(2),cvlo(3):cvhi(3))

    type(f_particle), intent(in), target :: particles(np)

    type(f_particle), pointer :: p1
    type(f_particle), pointer :: p2


    double precision fac1, fac2, fac3
    integer i,j,k,l,m

    fac1 = delt*neff/4d0

    do k = cplo(3), cphi(3)
      do j = cplo(2), cphi(2)
        do i = cplo(1), cphi(1)

          if (cellmembers(i,j,k) .lt. 2) then

            !if there less than two particles in the cell, use a fraction of modal thermal speed

            cellfactor(i, j, k) = 3.14159265359*d*d*cp/4d0;

          else

            fac2 = fac1*(cellmembers(i,j,k)**2)/cellvols(i,j,k)
            
            do l = 1, cellmembers(i,j,k)
              do m = 1, cellmembers(i,j,k)

                p1 => particles(l)
                p2 => particles(m)

                fac3 = sqrt((p1%vel(1)-p2%vel(1))**2 + (p1%vel(2)-p2%vel(2))**2 + (p1%vel(3)-p2%vel(3))**2)*3.14159265359*(p1%radius + p2%radius)**2

                if(fac3 .gt. cellfactor(i,j,k)) then

                  cellfactor(i,j,k) = fac3

                endif
                

              enddo
            enddo
          
          endif

        enddo
      enddo
    enddo

  end subroutine init_cells

  subroutine collide_cells(particles, cellmembers, cmlo, cmhi, celllists, cllo, clhi, cellpairs, cplo, cphi, cellfactor, cflo, cfhi, cellvols, cvlo, cvhi, ppc, np, neff, delt) bind(c,name='collide_cells')

    implicit none

    integer,          intent(in      ) :: np, ppc, cmlo(3), cmhi(3), cllo(3), clhi(3), cplo(3), cphi(3), cflo(3), cfhi(3), cvlo(3), cvhi(3)
    double precision, intent(in      ) :: neff, delt

    integer         , intent(in      ) :: cellmembers(cmlo(1):cmhi(1),cmlo(2):cmhi(2),cmlo(3):cmhi(3))
    integer         , intent(in      ) :: celllists(cmlo(1):cmhi(1),cmlo(2):cmhi(2),cmlo(3):cmhi(3),1:ppc)
    double precision, intent(inout   ) :: cellpairs(cplo(1):cphi(1),cplo(2):cphi(2),cplo(3):cphi(3))
    double precision, intent(inout   ) :: cellfactor(cplo(1):cphi(1),cplo(2):cphi(2),cplo(3):cphi(3))
    double precision, intent(inout   ) :: cellvols(cvlo(1):cvhi(1),cvlo(2):cvhi(2),cvlo(3):cvhi(3))

    type(f_particle), intent(in), target :: particles(np)

    type(f_particle), pointer :: p1
    type(f_particle), pointer :: p2


    !Currently assuming single species, although mass is explicitly included for future extension
    !Go through this and optimise later


    !double precision fac1, fac2, fac3, test, pairfrac
    integer i,j,k,l,m,n,pairs
    double precision fac1, fac2(3), fac3, sintheta, costheta, sinphi, cosphi, cmvel(3), totalmass

    fac3 = delt*neff/2d0

    do k = cplo(3), cphi(3)
      do j = cplo(2), cphi(2)
        do i = cplo(1), cphi(1)

          if (cellmembers(i,j,k) .gt. 1) then

            cellpairs(i,j,k) = fac3*cellfactor(i,j,k)*cellmembers(i,j,k)*cellmembers(i,j,k)/cellvols(i,j,k)

            pairs = floor(cellpairs(i,j,k))

            if( get_uniform_func() .lt. (cellpairs(i,j,k) - pairs) ) then

              pairs = pairs + 1;

            endif

            !print *, "Attempting ", pairs, " pairs. Cell factor: ", cellfactor(i,j,k)

            do n = 1, pairs

              call get_selector(l,cellmembers(i,j,k))
              m = l;

              do while(m .eq. l)

                call get_selector(m,cellmembers(i,j,k))

              enddo

              p1 => particles(l)
              p2 => particles(m)

              fac1 = sqrt((p1%vel(1)-p2%vel(1))**2 + (p1%vel(2)-p2%vel(2))**2 + (p1%vel(3)-p2%vel(3))**2)*(3.14159265359*(p1%radius + p2%radius)**2)

              if(fac1 .gt. cellfactor(i,j,k)) then

                cellfactor(i,j,k) = fac1

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
  
end module collision_functions_module





