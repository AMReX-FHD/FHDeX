module precheck_module
  use amrex_fort_module, only: amrex_real, amrex_particle_real
  use iso_c_binding ,    only: c_int

  implicit none

  public precheck

contains

  subroutine precheck(part, paramplanes, ns, delt, flag, phi, plo)

    use iso_c_binding, only: c_int
    use cell_sorted_particle_module, only: particle_t
    use paramplane_module, only: paramplane_t

    type(particle_t), intent(inout) :: part
    type(paramplane_t), intent(in), target :: paramplanes(ns)
    real(amrex_real), intent(in) :: delt, phi(3), plo(3)
    integer(c_int), intent(inout) :: flag
    integer(c_int), intent(in) :: ns
#if (BL_SPACEDIM == 3)
    double precision proj(3), box1lo(3), box1hi(3), box2lo(3), box2hi(3)
#endif
#if (BL_SPACEDIM == 2)
    double precision proj(2), box1lo(2), box1hi(2), box2lo(2), box2hi(2)
#endif

#if (BL_SPACEDIM == 3)
    proj = part%pos + part%vel*delt

    box1lo = plo
    box1hi = (/ paramplanes(7)%x0 , phi(2) , phi(3) /)

    box2lo = (/ paramplanes(7)%x0 , plo(2) , plo(3) /)
    box2hi = phi


!    box1lo = plo
!    box1hi = (/ phi(1)*0.5 , phi(2) , phi(3) /)

!    box2lo = (/ phi(1)*0.5 , plo(2) , plo(3) /)
!    box2hi = phi
#endif

#if (BL_SPACEDIM == 2)
    proj(1) = part%pos(1) + part%vel(1)*delt
    proj(2) = part%pos(2) + part%vel(2)*delt

    box1lo(1) = plo(1)
    box1lo(2) = plo(2)

    box1hi = (/ paramplanes(5)%x0 , phi(2) /)

    box2lo = (/ paramplanes(5)%x0 , plo(2) /)

    box2hi(1) = phi(1)
    box2hi(2) = phi(2)
#endif

#if (BL_SPACEDIM == 3)
    if(  (part%pos(1) < box1hi(1)) .and. (part%pos(2) < box1hi(2)) .and. (part%pos(3) < box1hi(3)) .and. (part%pos(1) > box1lo(1)) .and. (part%pos(2) > box1lo(2)) .and. (part%pos(3) > box1lo(3))  ) then  !started in box1

      if(  (proj(1) < box1hi(1)) .and. (proj(2) < box1hi(2)) .and. (proj(3) < box1hi(3)) .and. (proj(1) > box1lo(1)) .and. (proj(2) > box1lo(2)) .and. (proj(3) > box1lo(3)) ) then  !ended in box 1

        flag = 1

      endif

    else  !started in box 2

      if(  (proj(1) < box2hi(1)) .and. (proj(2) < box2hi(2)) .and. (proj(3) < box2hi(3)) .and. (proj(1) > box2lo(1)) .and. (proj(2) > box2lo(2)) .and. (proj(3) > box2lo(3)) ) then  !ended in box 2

        flag = 1

      endif
    endif
#endif
#if (BL_SPACEDIM == 2)
    if(  (part%pos(1) < box1hi(1)) .and. (part%pos(2) < box1hi(2)) .and. (part%pos(1) > box1lo(1)) .and. (part%pos(2) > box1lo(2)) ) then  !started in box1

      if(  (proj(1) < box1hi(1)) .and. (proj(2) < box1hi(2)) .and. (proj(1) > box1lo(1)) .and. (proj(2) > box1lo(2)) ) then  !ended in box 1

        flag = 1

      endif

    else  !started in box 2

      if(  (proj(1) < box2hi(1)) .and. (proj(2) < box2hi(2)) .and. (proj(1) > box2lo(1)) .and. (proj(2) > box2lo(2))) then  !ended in box 2

        flag = 1

      endif
    endif
#endif
  end subroutine precheck

!  subroutine precheck(part, paramplanes, ns, delt, flag, phi, plo)

!    use iso_c_binding, only: c_int
!    use cell_sorted_particle_module, only: particle_t
!    use paramplane_module, only: paramplane_t

!    type(particle_t), intent(inout) :: part
!    type(paramplane_t), intent(in), target :: paramplanes(ns)
!    real(amrex_real), intent(in) :: delt, phi(3), plo(3)
!    integer(c_int), intent(inout) :: flag
!    integer(c_int), intent(in) :: ns

!    double precision proj(3), projrad, rad, height

!    proj = part%pos + part%vel*delt

!    projrad = proj(1)**2 + proj(2)**2 + proj(3)**2

!    flag = 1
!    if((projrad .ge. 0.95*(2.5e-5)**2) .or. (proj(3) .le. 0) .or. (proj(3) .ge. 0.5e-5)) then
!      flag = 0
!    endif
!
!  end subroutine precheck

end module precheck_module
