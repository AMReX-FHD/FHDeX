module thermostat_module
  use amrex_fort_module, only: amrex_real, amrex_particle_real
  use iso_c_binding ,    only: c_int
  
  implicit none

  public thermostat, getVelocity, getTemp

contains        

  subroutine thermostat(particles, lo, hi, cell_part_ids, cell_part_cnt, clo, chi, &
    cellvols, cvlo, cvhi, neff, np, surfaces, ns, meanL, meanR, lC, rC)bind(c,name="thermostat")

    use iso_c_binding, only: c_ptr, c_int, c_f_pointer
    use cell_sorted_particle_module, only: particle_t
    use surfaces_module, only: surface_t

    integer,          intent(in      ) :: clo(3), chi(3), cvlo(3), cvhi(3), lo(3), hi(3), np
    double precision, intent(in      ) :: neff

    double precision, intent(inout   ) :: cellvols(cvlo(1):cvhi(1),cvlo(2):cvhi(2),cvlo(3):cvhi(3))

    type(c_ptr), intent(inout)      :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int), intent(inout)   :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))

    type(particle_t), intent(inout), target :: particles(np)

    type(particle_t), pointer :: part
    integer(c_int), pointer :: cell_parts(:)

    type(surface_t), intent(in), target :: surfaces(ns)
    integer(c_int), intent(in) :: ns
    double precision, intent(in) :: lC, rC, meanL, meanR


    integer i,j,k,p,cell_np
    double precision mid

    !mid = (hi(1)-lo(1))/2 !does not set correctly - fix later
    mid = 1.0
    !print *, "mid = ", mid

    !apply the scaling correction to the velocity
    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)

          cell_np = cell_part_cnt(i,j,k)
          call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])

          do p = 1, cell_np

            part => particles(cell_parts(p))

            !correct the velocities
            if(part%pos(1) < mid) then
              !part%vel(1) = part%vel(1)*sqrt(1d0/part%R)*lC
              part%vel(1) = (part%vel(1)-meanL)*lC
            else
              !part%vel(1) = part%vel(1)*sqrt(1d0/part%R)*rC
              part%vel(1) = (part%vel(1)-meanR)*rC
            endif

          enddo
        enddo
      enddo
    enddo


  end subroutine thermostat

  subroutine getVelocity(particles, lo, hi, cell_part_ids, cell_part_cnt, clo, chi, &
    cellvols, cvlo, cvhi, neff, np, surfaces, ns, pL, pR, vL, vR)bind(c,name="getVelocity")

    use iso_c_binding, only: c_ptr, c_int, c_f_pointer
    use cell_sorted_particle_module, only: particle_t
    use surfaces_module, only: surface_t

    integer,          intent(in      ) :: clo(3), chi(3), cvlo(3), cvhi(3), lo(3), hi(3), np
    double precision, intent(in      ) :: neff

    double precision, intent(inout   ) :: cellvols(cvlo(1):cvhi(1),cvlo(2):cvhi(2),cvlo(3):cvhi(3))

    type(c_ptr), intent(inout)      :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int), intent(inout)   :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))

    type(particle_t), intent(inout), target :: particles(np)

    type(particle_t), pointer :: part
    integer(c_int), pointer :: cell_parts(:)

    type(surface_t), intent(in), target :: surfaces(ns)
    integer(c_int), intent(in) :: ns
    integer(c_int), intent(inout) :: pL, pR
    double precision, intent(inout) :: vL, vR

    integer i,j,k,p,cell_np
    double precision mid

    !mid = (hi(1)-lo(1))/2 !does not set correctly - fix later
    mid = 1.0
    !print *, "mid = ", mid

    !get average velocities on left and right
    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)

          cell_np = cell_part_cnt(i,j,k)
          call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])

          do p = 1, cell_np

            part => particles(cell_parts(p))

            !print *, "Position, velocity : ", part%pos(1), part%vel(1)

            !increment the total velocity and member count
            if(part%pos(1) < mid) then
              vL = vL + part%vel(1)
              pL = pL + 1
            else
              vR = vR + part%vel(1)
              pR = pR + 1
            endif


          enddo
        enddo
      enddo
    enddo

  !print *, "pops, vels:", pL, pR, vL, vR

  end subroutine getVelocity

  subroutine getTemp(particles, lo, hi, cell_part_ids, cell_part_cnt, clo, chi, &
    cellvols, cvlo, cvhi, neff, np, surfaces, ns, meanL, meanR, varL, varR)bind(c,name="getTemp")

    use iso_c_binding, only: c_ptr, c_int, c_f_pointer
    use cell_sorted_particle_module, only: particle_t
    use surfaces_module, only: surface_t

    integer,          intent(in      ) :: clo(3), chi(3), cvlo(3), cvhi(3), lo(3), hi(3), np
    double precision, intent(in      ) :: neff

    double precision, intent(inout   ) :: cellvols(cvlo(1):cvhi(1),cvlo(2):cvhi(2),cvlo(3):cvhi(3))

    type(c_ptr), intent(inout)      :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int), intent(inout)   :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))

    type(particle_t), intent(inout), target :: particles(np)

    type(particle_t), pointer :: part
    integer(c_int), pointer :: cell_parts(:)

    type(surface_t), intent(in), target :: surfaces(ns)
    integer(c_int), intent(in) :: ns
    double precision, intent(in) :: meanL, meanR
    double precision, intent(inout) :: varL, varR


    integer i,j,k,p,cell_np
    double precision mid

    !mid = (hi(1)-lo(1))/2 !does not set correctly - fix later
    mid = 1.0
    !print *, "mid = ", mid
   
    !get the temperature on left and right through the velocity sample variance
    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)

          cell_np = cell_part_cnt(i,j,k)
          call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])

          do p = 1, cell_np

            part => particles(cell_parts(p))

            !measure temp through velocity sample variance
            if(part%pos(1) < mid) then
              varL = varL + (1d0/part%R) * (meanL-part%vel(1))**2
            else
              varR = varR + (1d0/part%R) * (meanR-part%vel(1))**2
            endif

          enddo
        enddo
      enddo
    enddo


  end subroutine getTemp



end module thermostat_module
