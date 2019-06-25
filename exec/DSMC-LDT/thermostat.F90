module thermostat_module
  use amrex_fort_module, only: amrex_real, amrex_particle_real
  use iso_c_binding ,    only: c_int
  
  implicit none

  public thermostat

contains        

  subroutine thermostat(particles, lo, hi, cell_part_ids, cell_part_cnt, clo, chi, &
    cellvols, cvlo, cvhi, neff, np, surfaces, ns, tL, tR)bind(c,name="thermostat")

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
    double precision, intent(in) :: tL, tR


    integer i,j,k,p,cell_np, totalparticles
    double precision Lmembersinv, Rmembersinv, vLx, vRx, mid, mLt, mRt, lC, rC

    totalparticles = 0
    !mid = (hi(1)-lo(1))/2 !does not set correctly - fix later
    mid = 1.0
    Lmembersinv = 0
    Rmembersinv = 0
    vLx = 0
    vRx = 0
    mLt = 0
    mRt = 0
    print *, "mid = ", mid

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
              vLx = vLx + part%vel(1)
              Lmembersinv = Lmembersinv + 1
            else
              vRx = vRx + part%vel(1)
              Rmembersinv = Rmembersinv + 1
            endif

            totalparticles = totalparticles + 1;

          enddo
        enddo
      enddo
    enddo

    !print *, Lmembersinv, Rmembersinv

    Lmembersinv = (1d0)/Lmembersinv
    Rmembersinv = (1d0)/Rmembersinv

    print *, Lmembersinv, Rmembersinv


    !average velocities
    vLx = vLx * Lmembersinv
    vRx = vRx * Rmembersinv

    print *, vLx, vRx
   
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
              mLt = mLt + (1d0/part%R) * (vLx-part%vel(1))**2
            else
              mRt = mRt + (1d0/part%R) * (vRx-part%vel(1))**2
            endif

          enddo
        enddo
      enddo
    enddo

    !temperature left and right
    mRt = mRt * Rmembersinv
    mLt = mLt * Lmembersinv

    print *, mRt, mLt

    !apply a velocity scaling to correct the temperature
    lC = sqrt(tL / mLt)
    rC = sqrt(tR / mRt)

    !print *, lC, rC

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
              part%vel(1) = (part%vel(1)-vLx)*lC
            else
              !part%vel(1) = part%vel(1)*sqrt(1d0/part%R)*rC
              part%vel(1) = (part%vel(1)-vRx)*rC
            endif

          enddo
        enddo
      enddo
    enddo


  end subroutine thermostat



end module thermostat_module
