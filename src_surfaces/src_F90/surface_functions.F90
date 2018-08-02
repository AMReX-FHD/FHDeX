


#if (BL_SPACEDIM == 3)

  subroutine rotation(costheta, sintheta, cosphi, sinphi, cx, cy, cz)

    use amrex_fort_module, only: amrex_real, amrex_particle_real
    use iso_c_binding ,    only: c_int
    implicit none

    real(amrex_real), intent(inout) :: cx, cy, cz
    real(amrex_real), intent(in) :: costheta, sintheta, cosphi, sinphi
    real(amrex_real) cxnew

    cxnew = cx*costheta*cosphi + cz*cosphi*sintheta - cy*sinphi
    cy = cy*cosphi + cx*costheta*sinphi + cz*sintheta*sinphi
    cz = cz*costheta - cx*sintheta

    cx = cxnew

  end subroutine rotation

  subroutine randomhemisphere(costheta, sintheta, cosphi, sinphi, cx, cy, cz)

    use rng_functions_module
    use amrex_fort_module, only: amrex_real, amrex_particle_real
    use iso_c_binding ,    only: c_int
    implicit none

    real(amrex_real), intent(inout) :: cx, cy, cz
    real(amrex_real), intent(in) :: costheta, sintheta, cosphi, sinphi
    
    real(amrex_real) mag, costhetanew, sinthetanew, cosphinew, sinphinew

    mag = sqrt(cx**2 + cy**2 + cz**2)

    call get_half_angles(costhetanew, sinthetanew, cosphinew, sinphinew)

    cx = mag*sinthetanew*cosphinew
    cy = mag*sinthetanew*sinphinew
    cz = mag*costhetanew

    call rotation(costheta, sintheta, cosphi, sinphi, cx, cy, cz)

  end subroutine randomhemisphere

#endif

#if (BL_SPACEDIM == 2)

  subroutine rotation(costheta, sintheta, cx, cy)

    use rng_functions_module
    use amrex_fort_module, only: amrex_real, amrex_particle_real
    use iso_c_binding ,    only: c_int
    implicit none

    real(amrex_real), intent(inout) :: cx, cy
    real(amrex_real), intent(in) :: costheta, sintheta
    
    real(amrex_real) cxnew

    cxnew = cx*costheta + cy*sintheta
    cy = -cx*sintheta + cy*costheta

    cx = cxnew

  end subroutine rotation

  subroutine randomhemisphere(costheta, sintheta, cx, cy, cz)

    use rng_functions_module
    use amrex_fort_module, only: amrex_real, amrex_particle_real
    use iso_c_binding ,    only: c_int
    implicit none

    real(amrex_real), intent(inout) :: cx, cy, cz
    real(amrex_real), intent(in) :: costheta, sintheta
    real(amrex_real) cxnew
    
    real(amrex_real) mag, costhetanew, sinthetanew, cosphinew, sinphinew

    mag = sqrt(cx**2 + cy**2 + cz**2)

    call get_half_angles(costhetanew, sinthetanew, cosphinew, sinphinew)

    cx = mag*costhetanew
    cy = mag*sinthetanew*cosphinew
    cz = mag*sinthetanew*sinphinew

    call rotation(costheta, sintheta, cx, cy)

  end subroutine randomhemisphere

#endif
  
  subroutine find_intersect(part, delt, surfaces, ns, intsurf, inttime, intside, phi, plo)
    
    use iso_c_binding, only: c_int
    use cell_sorted_particle_module, only: particle_t
    use surfaces_module
    use precheck_module
    
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

      if(  ((uval .gt. 0) .and. (uval .lt. surf%utop))    .and.     ((tval .gt. 0) .and. (tval .lt. inttime))  ) then

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

  subroutine apply_bc(surf, part, intside, domsize, push)
    
    use iso_c_binding, only: c_int
    use amrex_fort_module, only: amrex_real, amrex_particle_real
    use cell_sorted_particle_module, only: particle_t
    use surfaces_module
    use rng_functions_module
    
    implicit none

    integer(c_int),   intent(in) :: intside
    integer(c_int),   intent(inout) :: push
    double precision, intent(in) :: domsize(3)
    type(surface_t),  intent(inout) :: surf
    type(particle_t), intent(inout) :: part

    real(amrex_real) dotprod, srt

    if(intside .eq. 1) then
   
      if(get_uniform_func() > surf%porosityright) then

        push = 0

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

        push = 0

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
      else

        push = 1

        if(get_uniform_func() > surf%momentumright) then

#if (BL_SPACEDIM == 3)
          call randomhemisphere(surf%costhetaleft, surf%sinthetaleft, surf%cosphileft, surf%sinphileft, part%vel(1), part%vel(2), part%vel(3))
          !print *, "Scattering!"
#endif
#if (BL_SPACEDIM == 2)
          call randomhemisphere(surf%costhetaleft, surf%sinthetaleft, part%vel(1), part%vel(2), part%vel(3))
#endif

        endif

      endif

    else

      if(get_uniform_func() > surf%porosityleft) then

        push = 0

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
          !print *, "dot left: ", dotprod
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

        push = 0

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

      else
        push = 1      

        if(get_uniform_func() > surf%momentumleft) then
#if (BL_SPACEDIM == 3)
          call randomhemisphere(surf%costhetaright, surf%sinthetaright, surf%cosphiright, surf%sinphiright, part%vel(1), part%vel(2), part%vel(3))
#endif
#if (BL_SPACEDIM == 2)
          call randomhemisphere(surf%costhetaright, surf%sinthetaright, part%vel(1), part%vel(2), part%vel(3))
#endif

        endif

      endif
 
    endif
        
  end subroutine apply_bc
  


































