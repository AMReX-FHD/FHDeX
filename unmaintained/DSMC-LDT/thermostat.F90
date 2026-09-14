module thermostat_module
  use amrex_fort_module, only: amrex_real, amrex_particle_real
  use iso_c_binding ,    only: c_int

  implicit none

  public thermostat, getVelocity, getTemp

contains

  subroutine thermostat(particles, lo, hi, cell_part_ids, cell_part_cnt, clo, chi, &
    cellvols, cvlo, cvhi, neff, np, paramplanes, ns, meanLx, meanRx, meanLy, meanRy, &
    meanLz, meanRz, lC, rC)bind(c,name="thermostat")

    use iso_c_binding, only: c_ptr, c_int, c_f_pointer
    use cell_sorted_particle_module, only: particle_t
    use paramplane_module, only: paramplane_t

    integer,          intent(in      ) :: clo(3), chi(3), cvlo(3), cvhi(3), lo(3), hi(3), np
    double precision, intent(in      ) :: neff

    double precision, intent(inout   ) :: cellvols(cvlo(1):cvhi(1),cvlo(2):cvhi(2),cvlo(3):cvhi(3))

    type(c_ptr), intent(inout)      :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int), intent(inout)   :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))

    type(particle_t), intent(inout), target :: particles(np)

    type(particle_t), pointer :: part
    integer(c_int), pointer :: cell_parts(:)

    type(paramplane_t), intent(in), target :: paramplanes(ns)
    integer(c_int), intent(in) :: ns
    double precision, intent(in) :: lC, rC, meanLx, meanRx, meanLy, meanRy, meanLz, meanRz


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
              part%vel(1) = (part%vel(1)-meanLx)*lC
              part%vel(2) = (part%vel(2)-meanLy)*lC
              part%vel(3) = (part%vel(3)-meanLz)*lC
            else
              !part%vel(1) = part%vel(1)*sqrt(1d0/part%R)*rC
              part%vel(1) = (part%vel(1)-meanRx)*rC
              part%vel(2) = (part%vel(2)-meanRy)*rC
              part%vel(3) = (part%vel(3)-meanRz)*rC
            endif

          enddo
        enddo
      enddo
    enddo


  end subroutine thermostat

  subroutine getVelocity(particles, lo, hi, cell_part_ids, cell_part_cnt, clo, chi, &
    cellvols, cvlo, cvhi, neff, np, paramplanes, ns, pL, pR, vLx, vRx, vLy, vRy, &
    vLz, vRz)bind(c,name="getVelocity")

    use iso_c_binding, only: c_ptr, c_int, c_f_pointer
    use cell_sorted_particle_module, only: particle_t
    use paramplane_module, only: paramplane_t

    integer,          intent(in      ) :: clo(3), chi(3), cvlo(3), cvhi(3), lo(3), hi(3), np
    double precision, intent(in      ) :: neff

    double precision, intent(inout   ) :: cellvols(cvlo(1):cvhi(1),cvlo(2):cvhi(2),cvlo(3):cvhi(3))

    type(c_ptr), intent(inout)      :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int), intent(inout)   :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))

    type(particle_t), intent(inout), target :: particles(np)

    type(particle_t), pointer :: part
    integer(c_int), pointer :: cell_parts(:)

    type(paramplane_t), intent(in), target :: paramplanes(ns)
    integer(c_int), intent(in) :: ns
    integer(c_int), intent(inout) :: pL, pR
    double precision, intent(inout) :: vLx, vRx, vLy, vRy, vLz, vRz

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

            !print *, "Position, velocity : ",p, part%pos(1), part%vel(1)
            !print *,p,part%pos(1), part%vel(1), part%vel(2), part%vel(3)

            !increment the total velocity and member count
            if(part%pos(1) < mid) then
              vLx = vLx + part%vel(1)
              vLy = vLy + part%vel(2)
              vLz = vLz + part%vel(3)
              pL = pL + 1
            else
              vRx = vRx + part%vel(1)
              vRy = vRy + part%vel(2)
              vRz = vRz + part%vel(3)
              pR = pR + 1
            endif


          enddo
        enddo
      enddo
    enddo

    ! if (pR .lt. 10) then
    !   print *, pR, pL
    ! endif

  !print *, "pops, vels:", pL, pR, vL, vR

  end subroutine getVelocity

  subroutine getTemp(particles, lo, hi, cell_part_ids, cell_part_cnt, clo, chi, &
    cellvols, cvlo, cvhi, neff, np, paramplanes, ns, meanLx, meanRx, meanLy, meanRy, &
    meanLz, meanRz, varL, varR)bind(c,name="getTemp")

    use iso_c_binding, only: c_ptr, c_int, c_f_pointer
    use cell_sorted_particle_module, only: particle_t
    use paramplane_module, only: paramplane_t

    integer,          intent(in      ) :: clo(3), chi(3), cvlo(3), cvhi(3), lo(3), hi(3), np
    double precision, intent(in      ) :: neff

    double precision, intent(inout   ) :: cellvols(cvlo(1):cvhi(1),cvlo(2):cvhi(2),cvlo(3):cvhi(3))

    type(c_ptr), intent(inout)      :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
    integer(c_int), intent(inout)   :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))

    type(particle_t), intent(inout), target :: particles(np)

    type(particle_t), pointer :: part
    integer(c_int), pointer :: cell_parts(:)

    type(paramplane_t), intent(in), target :: paramplanes(ns)
    integer(c_int), intent(in) :: ns
    double precision, intent(in) :: meanLx, meanRx, meanLy, meanRy, meanLz, meanRz
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
              varL = varL + (1d0/part%R) * ((meanLx-part%vel(1))**2 + (meanLy-part%vel(2))**2 + (meanLz-part%vel(3))**2)
            else
              varR = varR + (1d0/part%R) * ((meanRx-part%vel(1))**2 + (meanRy-part%vel(2))**2 + (meanRz-part%vel(3))**2)
            endif

          enddo
        enddo
      enddo
    enddo


  end subroutine getTemp



end module thermostat_module
