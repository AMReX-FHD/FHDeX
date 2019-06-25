


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

        !print *, "s: ", s, " vtop: ", surf%vtop, " utop ", surf%utop

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
    surf%coltime=tval
        
  end subroutine find_intersect

  subroutine apply_bc(surf, part, intside, domsize, push, time, inttime)
    
    use iso_c_binding, only: c_int
    use amrex_fort_module, only: amrex_real, amrex_particle_real
    use cell_sorted_particle_module, only: particle_t
    use surfaces_module
    use rng_functions_module
    use common_namelist_module, only: prob_hi, fixed_dt, graphene_tog
    
    implicit none

    integer(c_int),   intent(in) :: intside
    integer(c_int),   intent(inout) :: push
    double precision, intent(in) :: domsize(3)
    type(surface_t),  intent(inout) :: surf
    type(particle_t), intent(inout) :: part

    real(amrex_real) dotprod, srt, time, inttime
    real(amrex_real) :: normvel(3), j(3), oldvel(3)

    if(surf%boundary .eq. 6)then
       oldvel=part%vel
       ! write(*,*) "old", oldvel(3), part%id
      ! print*, part%id, part%vel(3)
    endif
    
    if(intside .eq. 1) then
   
      if(get_uniform_func() > surf%porosityright) then

        push = 0

        surf%fxright = surf%fxright + part%mass*part%vel(1)
        surf%fyright = surf%fyright + part%mass*part%vel(2)
        surf%fzright = surf%fzright + part%mass*part%vel(3)

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

          part%pos(1) = part%pos(1) + 0.9999999*domsize(1)

        elseif(surf%boundary .eq. 2) then

          part%pos(1) = part%pos(1) - 0.9999999*domsize(1)

        elseif(surf%boundary .eq. 3) then

          part%pos(2) = part%pos(2) + 0.9999999*domsize(2)

        elseif(surf%boundary .eq. 4) then

          part%pos(2) = part%pos(2) - 0.9999999*domsize(2)

#if (BL_SPACEDIM == 3)

        elseif(surf%boundary .eq. 5) then

          part%pos(3) = part%pos(3) + 0.9999999*domsize(3)

        elseif(surf%boundary .eq. 6) then

          part%pos(3) = part%pos(3) - 0.9999999*domsize(3)
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

          part%pos(1) = part%pos(1) + 0.9999999*domsize(1)

        elseif(surf%boundary .eq. 2) then

          part%pos(1) = part%pos(1) - 0.9999999*domsize(1)

        elseif(surf%boundary .eq. 3) then

          part%pos(2) = part%pos(2) + 0.9999999*domsize(2)

        elseif(surf%boundary .eq. 4) then

          part%pos(2) = part%pos(2) - 0.9999999*domsize(2)

#if (BL_SPACEDIM == 3)
        elseif(surf%boundary .eq. 5) then

          part%pos(3) = part%pos(3) + 0.9999999*domsize(3)

        elseif(surf%boundary .eq. 6) then

          part%pos(3) = part%pos(3) - 0.9999999*domsize(3)
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
   if(graphene_tog .eq. 1) then
   if(surf%boundary .eq. 6) then
      ! call test(part, surf, intside)
      call surf_velocity(surf, part, time, oldvel, inttime)
   endif
   endif
        
  end subroutine apply_bc


 subroutine test(part, surf, intside, dt)
    
    use iso_c_binding, only: c_int
    use amrex_fort_module, only: amrex_real, amrex_particle_real
    use cell_sorted_particle_module, only: particle_t
    use surfaces_module
    use rng_functions_module
    
    implicit none

    type(particle_t), intent(inout) :: part
    type(surface_t) :: surf
    integer(c_int) :: count5, count6, intside
    real(amrex_real) :: magnormvel, dt
    real(amrex_real), dimension(3):: rnorm, lnorm, j, normvel, surfvel

    rnorm=(/ surf%rnx, surf%rny, surf%rnz /)
    lnorm=(/ surf%lnx, surf%lny, surf%lnz /)
    surfvel=(/ surf%velx, surf%vely, surf%velz /)
    if(intside .eq. 1) then
       normvel=dot_product(part%vel, rnorm)*rnorm
    else
       normvel=dot_product(part%vel, lnorm)*lnorm
    endif
    j=normvel*part%mass
    magnormvel=norm2(normvel)
    normvel=normvel/magnormvel
   ! part%vel(3)=part%vel(3)+surf%velz
   ! part%vel=part%vel-surfvel
   ! part%vel=part%vel+surf%velz*normvel
       !write(*,*) surf%velz*normvel

             ! intsurf=surf%boundary
  
             !     if(intsurf .eq.  5) then
             !        write(*,*) "5", part%pos(1), part%pos(2), part%pos(3)
             !        count5=count5+1
             !     elseif(intsurf .eq. 6) then
             !        write(*,*)  "6", part%pos(1), part%pos(2), part%pos(3)
             !       count6= count6+1
    !    endif

    
end subroutine test
  
subroutine surf_velocity(surf, part, time, oldvel, inttime)
  
 use iso_c_binding, only: c_int
 use amrex_fort_module, only: amrex_real, amrex_particle_real
 use cell_sorted_particle_module, only: particle_t
 use surfaces_module
 use rng_functions_module
 use common_namelist_module, only: prob_hi, fixed_dt
 

 implicit none

 type(particle_t), intent(inout) :: part
 type(surface_t) :: surf 
 integer(c_int)  i, count, step, ii
 real(amrex_real) surfvel, r, f_x, a, r2, r3, time, bessj0, dbessj0, k, rho, tau, omega, dt, c, alpha, pi, graphi, grac, xvec, yvec, interval, radius, t, inttime
  real(amrex_real), dimension(3)::oldvel
 character (len=90) :: filename

    r=sqrt(part%pos(1)**2+part%pos(2)**2)
    c=91.4468
    alpha=0
    pi=3.1415926535897932
    a=(part%vel(3)+oldvel(3))*part%mass
    !parabola
    ! f_x=-a*r*r+a*d*r+100000
    bessj0=0
    t=time+inttime
        do i=1,1
           if(i .eq. 1)then
              k=2.4048
           elseif(i .eq. 2)then
              k=5.5201
           elseif(i .eq. 3)then
              k=8.6537
           elseif(i .eq. 4)then
              k=11.7915
           else
              k=14.9309
           endif
           r2=r*k/prob_hi(1)
           omega=sqrt(((c*(k**2))/prob_hi(1)**2)+alpha)/(3.141592653589793**2)
           bessj0 =a*surf%grac*bessel_jn(0,r2)*sin((t*omega)+surf%graphi)
           dbessj0=a*surf%grac*bessel_jn(0, r2)*omega
           graphi=omega*t
           grac=(c**2/(pi*prob_hi(1)**2))*bessel_jn(0,r2)/((bessel_jn(1, k)**2)*omega)
           xvec=surf%grac*cos(surf%graphi)+grac*cos(graphi)
           yvec=surf%grac*sin(surf%graphi)+grac*sin(graphi)
           surf%grac=sqrt(xvec**2+yvec**2)
           surf%graphi=atan2(yvec, xvec)
        enddo
       ! print*, surf%velz
      surf%velz=dbessj0*cos((t*omega)+surf%graphi)
      part%vel(3)=part%vel(3)+surf%velz
      step=time/fixed_dt
   
   !  if(step .eq. 300)then
    ! write(*,*) a
    ! write(*,*) "old", oldvel(3), part%id
    ! write(*,*) "new", part%vel(3)
   !  endif
  end subroutine surf_velocity

