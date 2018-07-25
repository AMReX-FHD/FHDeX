module surfaces_module
  use amrex_fort_module, only: amrex_real, amrex_particle_real
  use iso_c_binding ,    only: c_int
  use rng_functions_module
  
  implicit none
  private
  
  public surface_t, find_intersect, apply_bc
  
  type, bind(C) :: surface_t
#if (BL_SPACEDIM == 3)
     real(amrex_particle_real) :: x0
     real(amrex_particle_real) :: y0
     real(amrex_particle_real) :: z0

     real(amrex_particle_real) :: ux
     real(amrex_particle_real) :: uy
     real(amrex_particle_real) :: uz

     real(amrex_particle_real) :: vx
     real(amrex_particle_real) :: vy
     real(amrex_particle_real) :: Vz

     real(amrex_particle_real) :: lnx
     real(amrex_particle_real) :: lny
     real(amrex_particle_real) :: lnz

     real(amrex_particle_real) :: rnx
     real(amrex_particle_real) :: rny
     real(amrex_particle_real) :: rnz

     real(amrex_particle_real) :: utop
     real(amrex_particle_real) :: vtop

     real(amrex_particle_real) :: costhetaleft
     real(amrex_particle_real) :: sinthetaleft
     real(amrex_particle_real) :: cosphileft
     real(amrex_particle_real) :: sinphileft

     real(amrex_particle_real) :: costhetaright
     real(amrex_particle_real) :: sinthetaright
     real(amrex_particle_real) :: cosphiright
     real(amrex_particle_real) :: sinphiright

#endif

#if (BL_SPACEDIM == 2)
     real(amrex_particle_real) :: x0
     real(amrex_particle_real) :: y0

     real(amrex_particle_real) :: ux
     real(amrex_particle_real) :: uy

     real(amrex_particle_real) :: lnx
     real(amrex_particle_real) :: lny

     real(amrex_particle_real) :: rnx
     real(amrex_particle_real) :: rny

     real(amrex_particle_real) :: utop

     real(amrex_particle_real) :: costhetaleft
     real(amrex_particle_real) :: sinthetaleft
     real(amrex_particle_real) :: costhetaright
     real(amrex_particle_real) :: sinthetaright

#endif

     real(amrex_particle_real) :: porosityleft
     real(amrex_particle_real) :: specularityleft
     real(amrex_particle_real) :: temperatureleft

     real(amrex_particle_real) :: porosityright
     real(amrex_particle_real) :: specularityright
     real(amrex_particle_real) :: temperatureright

     real(amrex_particle_real) :: periodicity

     real(amrex_particle_real) :: fxleft
     real(amrex_particle_real) :: fyleft
     real(amrex_particle_real) :: fzleft

     real(amrex_particle_real) :: fxright
     real(amrex_particle_real) :: fyright
     real(amrex_particle_real) :: fzright

     real(amrex_particle_real) :: fxleftav
     real(amrex_particle_real) :: fyleftav
     real(amrex_particle_real) :: fzleftav

     real(amrex_particle_real) :: fxrightav
     real(amrex_particle_real) :: fyrightav
     real(amrex_particle_real) :: fzrightav

     integer :: boundary

  end type surface_t

contains

#if (BL_SPACEDIM == 3)

  subroutine rotation(costheta, sintheta, cosphi, sinphi, cx, cy, cz)

    implicit none

    real(amrex_real), intent(inout) :: cx, cy, cz
    real(amrex_real), intent(in) :: costheta, sintheta, cosphi, sinphi
    
    real(amrex_real) cxnew, cynew, cznew

    cxnew = cx*costheta*cosphi + cz*cosphi*sintheta - cy*sinphi
    cynew = cy*cosphi + cx*costheta*sinphi + cz*sintheta*sinphi
    cznew = cz*costheta - cx*sintheta

    cx = cxnew
    cy = cynew
    cz = cznew

  end subroutine rotation

#endif

#if (BL_SPACEDIM == 2)

  subroutine rotation(costheta, sintheta, cx, cy)

    implicit none

    real(amrex_real), intent(inout) :: cx, cy
    real(amrex_real), intent(in) :: costheta, sintheta
    
    real(amrex_real) cxnew, cynew

    cxnew = cx*costheta + cy*sintheta
    cynew = -cx*sintheta + cy*costheta

    cx = cxnew
    cy = cynew

  end subroutine rotation

#endif

  subroutine precheck(part, surfaces, ns, delt, flag, phi, plo)

    use iso_c_binding, only: c_int
    use cell_sorted_particle_module, only: particle_t

    type(particle_t), intent(inout) :: part
    type(surface_t), intent(in), target :: surfaces(ns)
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
    box1hi = (/ surfaces(7)%x0 , phi(2) , phi(3) /)

    box2lo = (/ surfaces(7)%x0 , plo(2) , plo(3) /)
    box2hi = phi
#endif

#if (BL_SPACEDIM == 2)
    proj(1) = part%pos(1) + part%vel(1)*delt
    proj(2) = part%pos(2) + part%vel(2)*delt

    box1lo(1) = plo(1)
    box1lo(2) = plo(2)

    box1hi = (/ surfaces(7)%x0 , phi(2) /)

    box2lo = (/ surfaces(7)%x0 , plo(2) /)

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
  
  subroutine find_intersect(part, delt, surfaces, ns, intsurf, inttime, intside, phi, plo)
    
    use iso_c_binding, only: c_int
    use cell_sorted_particle_module, only: particle_t
    
    implicit none

    type(particle_t), intent(inout) :: part
    type(surface_t), intent(in), target :: surfaces(ns)
    real(amrex_real), intent(in) :: delt, phi(3), plo(3)
    real(amrex_real), intent(inout) :: inttime
    integer(c_int), intent(inout) :: intsurf, intside
    integer(c_int), intent(in) :: ns

    integer :: s, flag
    real(amrex_real) denominv, uval, vval, tval, dotprod
    type(surface_t), pointer :: surf

    inttime = delt
    intsurf = -1

    flag = 0    
    call precheck(part, surfaces, ns, delt, flag, phi, plo)
    
#if (BL_SPACEDIM == 3)

    if(flag .eq. 0) then

      do s = 1, ns

        surf => surfaces(s)

        denominv = 1d0/(part%vel(3)*surf%uy*surf%vx - part%vel(2)*surf%uz*surf%vx - part%vel(3)*surf%ux*surf%vy + part%vel(1)*surf%uz*surf%vy + part%vel(2)*surf%ux*surf%vz - part%vel(1)*surf%uy*surf%vz)

        uval = (part%vel(3)*part%pos(2)*surf%vx - part%vel(2)*part%pos(3)*surf%vx - part%vel(3)*surf%y0*surf%vx + part%vel(2)*surf%z0*surf%vx - part%vel(3)*part%pos(1)*surf%vy + part%vel(1)*part%pos(3)*surf%vy + part%vel(3)*surf%x0*surf%vy - part%vel(1)*surf%z0*surf%vy + part%vel(2)*part%pos(1)*surf%vz - part%vel(1)*part%pos(2)*surf%vz -  part%vel(2)*surf%x0*surf%vz + part%vel(1)*surf%y0*surf%vz)*denominv

        vval = (-part%vel(3)*part%pos(2)*surf%ux + part%vel(2)*part%pos(3)*surf%ux + part%vel(3)*surf%y0*surf%ux - part%vel(2)*surf%z0*surf%ux + part%vel(3)*part%pos(1)*surf%uy - part%vel(1)*part%pos(3)*surf%uy - part%vel(3)*surf%x0*surf%uy + part%vel(1)*surf%z0*surf%uy - part%vel(2)*part%pos(1)*surf%uz + part%vel(1)*part%pos(2)*surf%uz + part%vel(2)*surf%x0*surf%uz - part%vel(1)*surf%y0*surf%uz)*denominv

        tval = (-part%pos(3)*surf%uy*surf%vx + surf%z0*surf%uy*surf%vx + part%pos(2)*surf%uz*surf%vx - surf%y0*surf%uz*surf%vx + part%pos(3)*surf%ux*surf%vy - surf%z0*surf%ux*surf%vy - part%pos(1)*surf%uz*surf%vy + surf%x0*surf%uz*surf%vy - part%pos(2)*surf%ux*surf%vz + surf%y0*surf%ux*surf%vz + part%pos(1)*surf%uy*surf%vz - surf%x0*surf%uy*surf%vz)*denominv

        if(  ((uval .gt. 0) .and. (uval .lt. surf%utop))    .and.   ((vval .gt. 0) .and. (vval .lt. surf%vtop))    .and.     ((tval .gt. 0) .and. (tval .lt. inttime))   ) then

          inttime = tval
          intsurf = s

          dotprod = part%vel(1)*surf%lnx + part%vel(2)*surf%lny + part%vel(3)*surf%lnz

          if (dotprod .gt. 0) then
            intside = 1 !1 for rhs
          else
            intside = 0 !0 for lhs
          endif

          !print *, "Intersection! Intsurf: ", intsurf, ". Intside: ", intside, ". tval: ", tval

        endif
      enddo
    endif

#endif

#if (BL_SPACEDIM == 2)

    do s = 1, ns

      surf => surfaces(s)

      denominv = 1d0/(part%vel(2)*surf%ux - part%vel(1)*surf%uy)

      tval = (-part%pos(2)*surf%ux + surf%y0*surf%ux + part%pos(1)*surf%uy - surf%x0*surf%uy)*denominv

      uval = (part%vel(2)*part%pos(1) - part%vel(1)*part%pos(2) - part%vel(2)*surf%x0 + part%vel(1)*surf%y0)*denominv

      if(  ((uval .gt. 0) .and. (uval .lt. surf%utop)) ) then

        inttime = tval
        intsurf = s

        dotprod = part%vel(1)*surf%lnx + part%vel(2)*surf%lny

        if (dotprod .gt. 0) then
          intside = 1 !1 for rhs
        else
          intside = 0 !0 for lhs
        endif

      endif

    enddo

#endif
        
  end subroutine find_intersect

  subroutine apply_bc(surf, part, intside, domsize)
    
    use iso_c_binding, only: c_int
    use cell_sorted_particle_module, only: particle_t
    
    implicit none

    integer(c_int),   intent(in) :: intside
    double precision, intent(in) :: domsize(3)
    type(surface_t),  intent(inout) :: surf
    type(particle_t), intent(inout) :: part

    real(amrex_real) dotprod, srt

    if(intside .eq. 1) then
   
      if(get_uniform_func() > surf%porosityright) then

        surf%fxright = surf%fxright + part%mass*part%vel(1)
        surf%fyright = surf%fyright + part%mass*part%vel(1)
        surf%fzright = surf%fzright + part%mass*part%vel(1)

        if(get_uniform_func() < surf%specularityright) then

#if (BL_SPACEDIM == 2)
          dotprod = part%vel(1)*surf%rnx + part%vel(2)*surf%rny
#endif
#if (BL_SPACEDIM == 3)
          dotprod = part%vel(1)*surf%rnx + part%vel(2)*surf%rny + part%vel(3)*surf%rnz
#endif
          part%vel(1) = -2d0*dotprod*surf%rnx + part%vel(1)
          part%vel(2) = -2d0*dotprod*surf%rny + part%vel(2)
#if (BL_SPACEDIM == 3)
          part%vel(3) = -2d0*dotprod*surf%rnz + part%vel(3)
#endif
        else
          
          srt = sqrt(part%r*surf%temperatureright)

#if (BL_SPACEDIM == 3)
          part%vel(1) = srt*get_particle_normal_func()
          part%vel(2) = srt*get_particle_normal_func()
          part%vel(3) = 1.414213562*srt*sqrt(-log(get_uniform_func()))
        
          call rotation(surf%costhetaright, surf%sinthetaright, surf%cosphiright, surf%sinphiright, part%vel(1), part%vel(2), part%vel(3))
#endif

#if (BL_SPACEDIM == 2)
          part%vel(1) = 1.414213562*srt*sqrt(-log(get_uniform_func()))
          part%vel(2) = srt*get_particle_normal_func()
          part%vel(3) = srt*get_particle_normal_func()

          call rotation(surf%costhetaright, surf%sinthetaright, part%vel(1), part%vel(2))
#endif


        endif
      elseif(get_uniform_func() < surf%periodicity) then

        if(surf%boundary .eq. 1) then

          part%pos(1) = part%pos(1) + 0.9999*domsize(1)

        elseif(surf%boundary .eq. 2) then

          part%pos(1) = part%pos(1) - 0.9999*domsize(1)

        elseif(surf%boundary .eq. 3) then

          part%pos(2) = part%pos(2) + 0.9999*domsize(2)

        elseif(surf%boundary .eq. 4) then

          part%pos(2) = part%pos(2) - 0.9999*domsize(2)

#if (BL_SPACEDIM == 3)

        elseif(surf%boundary .eq. 5) then

          part%pos(3) = part%pos(3) + 0.9999*domsize(3)

        elseif(surf%boundary .eq. 6) then

          part%pos(3) = part%pos(3) - 0.9999*domsize(3)
#endif

        endif
      endif

    else

      if(get_uniform_func() > surf%porosityleft) then

        surf%fxleft = surf%fxleft + part%mass*part%vel(1)
        surf%fyleft = surf%fyleft + part%mass*part%vel(1)
        surf%fzleft = surf%fzleft + part%mass*part%vel(1)

        if(get_uniform_func() < surf%specularityleft) then

#if (BL_SPACEDIM == 2)
          dotprod = part%vel(1)*surf%lnx + part%vel(2)*surf%lny
#endif
#if (BL_SPACEDIM == 3)
          dotprod = part%vel(1)*surf%lnx + part%vel(2)*surf%lny + part%vel(3)*surf%lnz
#endif

          part%vel(1) = -2d0*dotprod*surf%lnx + part%vel(1)
          part%vel(2) = -2d0*dotprod*surf%lny + part%vel(2)
#if (BL_SPACEDIM == 3)
          part%vel(3) = -2d0*dotprod*surf%lnz + part%vel(3)
#endif
        else
          
          srt = sqrt(part%r*surf%temperatureleft)

#if (BL_SPACEDIM == 3)
          part%vel(1) = srt*get_particle_normal_func()
          part%vel(2) = srt*get_particle_normal_func()
          part%vel(3) = 1.414213562*srt*sqrt(-log(get_uniform_func()))

          call rotation(surf%costhetaleft, surf%sinthetaleft, surf%cosphileft, surf%sinphileft, part%vel(1), part%vel(2), part%vel(3))
#endif
#if (BL_SPACEDIM == 2)
          part%vel(1) = 1.414213562*srt*sqrt(-log(get_uniform_func()))
          part%vel(2) = srt*get_particle_normal_func()
          part%vel(3) = srt*get_particle_normal_func()

          call rotation(surf%costhetaleft, surf%sinthetaleft, part%vel(1), part%vel(2))
#endif
       
        endif
      elseif(get_uniform_func() < surf%periodicity) then

        if(surf%boundary .eq. 1) then

          part%pos(1) = part%pos(1) + 0.9999*domsize(1)

        elseif(surf%boundary .eq. 2) then

          part%pos(1) = part%pos(1) - 0.9999*domsize(1)

        elseif(surf%boundary .eq. 3) then

          part%pos(2) = part%pos(2) + 0.9999*domsize(2)

        elseif(surf%boundary .eq. 4) then

          part%pos(2) = part%pos(2) - 0.9999*domsize(2)

#if (BL_SPACEDIM == 3)
        elseif(surf%boundary .eq. 5) then

          part%pos(3) = part%pos(3) + 0.9999*domsize(3)

        elseif(surf%boundary .eq. 6) then

          part%pos(3) = part%pos(3) - 0.9999*domsize(3)
#endif
        endif

      endif
 
    endif
        
  end subroutine apply_bc
  
end module surfaces_module


































