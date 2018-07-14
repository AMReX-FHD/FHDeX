module surfaces_module
  use amrex_fort_module, only: amrex_real, amrex_particle_real
  use iso_c_binding ,    only: c_int
  use cell_sorted_particle_module
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
  
  subroutine find_intersect(cx,cy,cz,lx0,ly0,lz0,delt, surfaces, ns, intsurf, inttime, intside)
    
    use iso_c_binding, only: c_int
    
    implicit none

    type(surface_t), intent(in), target :: surfaces(ns)
    real(amrex_real), intent(in) :: delt, cx, cy, cz, lx0, ly0, lz0
    real(amrex_real), intent(inout) :: inttime
    integer(c_int), intent(inout) :: intsurf, intside
    integer(c_int), intent(in) :: ns


    integer :: s
    real(amrex_real) ux, uy, uz, vx, vy, vz, sx0, sy0, sz0, denominv, uval, vval, tval, umax, vmax, dotprod

    type(surface_t), pointer :: surf

    inttime = delt
    intsurf = -1

    do s = 1, ns

      surf => surfaces(s)
      ux = surf%ux
      uy = surf%uy
      uz = surf%uz

      vx = surf%vx
      vy = surf%vy
      vz = surf%vz

      sx0 = surf%x0
      sy0 = surf%y0
      sz0 = surf%z0

      umax = surf%utop
      vmax = surf%vtop

      denominv = 1d0/(cz*uy*vx - cy*uz*vx - cz*ux*vy + cx*uz*vy + cy*ux*vz - cx*uy*vz)

      uval = (cz*ly0*vx - cy*lz0*vx - cz*sy0*vx + cy*sz0*vx - cz*lx0*vy + cx*lz0*vy + cz*sx0*vy - cx*sz0*vy + cy*lx0*vz - cx*ly0*vz -  cy*sx0*vz + cx*sy0*vz)*denominv

      vval = (-cz*ly0*ux + cy*lz0*ux + cz*sy0*ux - cy*sz0*ux + cz*lx0*uy - cx*lz0*uy - cz*sx0*uy + cx*sz0*uy - cy*lx0*uz + cx*ly0*uz + cy*sx0*uz - cx*sy0*uz)*denominv

      tval = (-lz0*uy*vx + sz0*uy*vx + ly0*uz*vx - sy0*uz*vx + lz0*ux*vy - sz0*ux*vy - lx0*uz*vy + sx0*uz*vy - ly0*ux*vz + sy0*ux*vz + lx0*uy*vz - sx0*uy*vz)*denominv

      if(  ((uval .gt. 0) .and. (uval .lt. surf%utop))    .and.   ((vval .gt. 0) .and. (vval .lt. surf%vtop))    .and.     ((tval .gt. 0) .and. (tval .lt. inttime))   ) then

        inttime = tval
        intsurf = s

        dotprod = cx*surf%lnx + cy*surf%lny + cz*surf%lnz

        if (dotprod .gt. 0) then
          intside = 1 !1 for rhs
        else
          intside = 0 !0 for lhs
        endif

        !print *, "Intersection! Intsurf: ", intsurf, ". Intside: ", intside, ". tval: ", tval

      endif

    enddo
        
  end subroutine find_intersect


  subroutine apply_bc(cx, cy, cz, nx,ny,nz, costheta, sintheta, cosphi, sinphi, porosity, specularity, temperature, r)
    
    use iso_c_binding, only: c_int
    
    implicit none

    real(amrex_real), intent(in) :: nx, ny, nz, costheta, sintheta, cosphi, sinphi, porosity, specularity, temperature, r
    real(amrex_real), intent(inout) :: cx, cy, cz

    real(amrex_real) dotprod, srt

    !print *, "Norm: ", nx, ny, nz

    if(get_uniform_func() > porosity) then
      if(get_uniform_func() < specularity) then

        dotprod = cx*nx + cy*ny + cz*nz

        cx = -2d0*dotprod*nx + cx
        cy = -2d0*dotprod*ny + cy
        cz = -2d0*dotprod*nz + cz

      else
        
        srt = sqrt(r*temperature)

        cx = srt*get_particle_normal_func()
        cy = srt*get_particle_normal_func()
        cz = 1.414213562*srt*sqrt(-log(get_uniform_func()))

        !print *, "Pre rotation: ", cx, cy, cz
      
        call rotation(costheta, sintheta, cosphi, sinphi, cx, cy, cz)

      endif
    endif
        
  end subroutine apply_bc

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
  
  subroutine find_intersect(cx,cy,lx0,ly0,delt, surfaces, ns, intsurf, inttime, intside)
    
    use iso_c_binding, only: c_int
    
    implicit none

    type(surface_t), intent(in), target :: surfaces(ns)
    real(amrex_real), intent(in) :: delt, cx, cy, lx0, ly0
    real(amrex_real), intent(inout) :: inttime
    integer(c_int), intent(inout) :: intsurf, intside
    integer(c_int), intent(in) :: ns


    integer :: s
    real(amrex_real) ux, uy, sx0, sy0, denominv, uval, tval, umax, dotprod

    type(surface_t), pointer :: surf

    inttime = delt
    intsurf = -1

    do s = 1, ns

      surf => surfaces(s)
      ux = surf%ux
      uy = surf%uy

      sx0 = surf%x0
      sy0 = surf%y0

      umax = surf%utop

      denominv = 1d0/(cy*ux - cx*uy)

      tval = (-ly0*ux + sy0*ux + lx0*uy - sx0*uy)*denominv

      uval = (cy*lx0 - cx*ly0 - cy*sx0 + cx*sy0)*denominv

      if(  ((uval .gt. 0) .and. (uval .lt. surf%utop))  .and.     ((tval .gt. 0) .and. (tval .lt. inttime))   ) then

        inttime = tval
        intsurf = s

        dotprod = cx*surf%lnx + cy*surf%lny

        if (dotprod .gt. 0) then
          intside = 1 !1 for rhs
        else
          intside = 0 !0 for lhs
        endif

      endif

    enddo
        
  end subroutine find_intersect

  subroutine apply_bc(cx, cy,nx,ny, costheta, sintheta, porosity, specularity, temperature, r)
    
    use iso_c_binding, only: c_int
    
    implicit none

    real(amrex_real), intent(in) :: nx, ny, costheta, sintheta, porosity, specularity, temperature, r
    real(amrex_real), intent(inout) :: cx, cy

    real(amrex_real) dotprod

    if(get_uniform_func() > porosity) then
      if(get_uniform_func() < specularity) then

        dotprod = cx*nx + cy*ny 

        cx = -2d0*dotprod*nx + cx
        cy = -2d0*dotprod*ny + cy

      else

        cx = sqrt(r*temperature)*log(-get_uniform_func())
        cy = get_particle_normal_func()
      
        call rotation(costheta, sintheta, cx, cy)

      endif
    endif
        
  end subroutine apply_bc
#endif
  
end module surfaces_module







































