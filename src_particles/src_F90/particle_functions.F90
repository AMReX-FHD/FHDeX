module particle_functions_module

  use amrex_error_module
  use amrex_fort_module
  use amrex_constants_module
  use common_namelist_module

  implicit none
  private

  public  f_particle
  
  type, bind(C)  :: f_particle
#if (BL_SPACEDIM == 2)
     real(amrex_particle_real)   :: pos(2)
#endif
#if (BL_SPACEDIM == 3)
     real(amrex_particle_real)   :: pos(3)
#endif
     real(amrex_particle_real)   :: mass
     real(amrex_particle_real)   :: fluid_density

     real(amrex_particle_real)   :: temperature
     real(amrex_particle_real)   :: fluid_temperature

     real(amrex_particle_real)   :: fluid_viscosity
     real(amrex_particle_real)   :: radius
     real(amrex_particle_real)   :: drag_factor


     real(amrex_particle_real)   :: vel(3)
     real(amrex_particle_real)   :: angular_vel(2)
     real(amrex_particle_real)   :: fluid_vel(3)

     integer(c_int)              :: id
     integer(c_int)              :: cpu

     integer(c_int)              :: cell_reverse_index
     integer(c_int)              :: species
#if (BL_SPACEDIM == 2) 
     integer(c_int)              :: cell_index(2)
#endif
#if (BL_SPACEDIM == 3) 
     integer(c_int)              :: cell_index(3)
#endif


  end type f_particle

contains

  subroutine update_particles(particles, np, dt, dx, index_lo, index_hi, real_lo, real_hi, &
                                       velx, velxlo, velxhi, &
                                       vely, velylo, velyhi, &
#if (BL_SPACEDIM == 3)
                                       velz, velzlo, velzhi, &
#endif
                                       coordsx, coordsxlo, coordsxhi, &
                                       coordsy, coordsylo, coordsyhi, &
#if (BL_SPACEDIM == 3)
                                       coordsz, coordszlo, coordszhi, &
#endif
#if (AMREX_SPACEDIM == 2)
                                       beta, betalo, betahi, &
#endif
#if (AMREX_SPACEDIM == 3)
                                       betaxy, betaxylo, betaxyhi, &
                                       betaxz, betaxzlo, betaxzhi, &
                                       betayz, betayzlo, betayzhi, &
#endif
                                       rho, rholo, rhohi, &

                                       sourcex, sourcexlo, sourcexhi, &
                                       sourcey, sourceylo, sourceyhi &
#if (BL_SPACEDIM == 3)
                                       , sourcez, sourcezlo, sourcezhi &
#endif
                                       ) bind(c,name='update_particles')


    implicit none

    integer,          intent(in   )         :: np, index_lo(3), index_hi(3), velxlo(3), velxhi(3), velylo(3), velyhi(3)
    integer,          intent(in   )         :: sourcexlo(3), sourcexhi(3), sourceylo(3), sourceyhi(3), rholo(3), rhohi(3)
    integer,          intent(in   )         :: coordsxlo(3), coordsxhi(3), coordsylo(3), coordsyhi(3)
#if (AMREX_SPACEDIM == 2)
    integer,          intent(in   )         :: betalo(3), betahi(3)
#endif
#if (AMREX_SPACEDIM == 3)
    integer,          intent(in   )         :: velzlo(3), velzhi(3), sourcezlo(3), sourcezhi(3), coordszlo(3), coordszhi(3)
    integer,          intent(in   )         :: betaxylo(3), betaxyhi(3), betaxzlo(3), betaxzhi(3), betayzlo(3), betayzhi(3)
#endif
    type(f_particle), intent(inout), target :: particles(np)
    double precision, intent(in   )         :: dt
    double precision, intent(in   )         :: dx(3)

    double precision, intent(in   )         :: real_lo(3)
    double precision, intent(in   )         :: real_hi(3)

    double precision, intent(in   ) :: velx(velxlo(1):velxhi(1),velxlo(2):velxhi(2),velxlo(3):velxhi(3))
    double precision, intent(in   ) :: vely(velylo(1):velyhi(1),velylo(2):velyhi(2),velylo(3):velyhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: velz(velzlo(1):velzhi(1),velzlo(2):velzhi(2),velzlo(3):velzhi(3))
#endif

    double precision, intent(in   ) :: coordsx(coordsxlo(1):coordsxhi(1),coordsxlo(2):coordsxhi(2),coordsxlo(3):coordsxhi(3),1:AMREX_SPACEDIM)
    double precision, intent(in   ) :: coordsy(coordsylo(1):coordsyhi(1),coordsylo(2):coordsyhi(2),coordsylo(3):coordsyhi(3),1:AMREX_SPACEDIM)
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: coordsz(coordszlo(1):coordszhi(1),coordszlo(2):coordszhi(2),coordszlo(3):coordszhi(3),1:AMREX_SPACEDIM)
#endif

#if (AMREX_SPACEDIM == 2)
    double precision, intent(in   ) :: beta(betalo(1):betahi(1),betalo(2):betahi(2),betalo(3):betahi(3)) !nodal
#endif
#if (AMREX_SPACEDIM == 3)
    double precision, intent(in   ) :: betaxy(betaxylo(1):betaxyhi(1),betaxylo(2):betaxyhi(2),betaxylo(3):betaxyhi(3))
    double precision, intent(in   ) :: betaxz(betaxzlo(1):betaxzhi(1),betaxzlo(2):betaxzhi(2),betaxzlo(3):betaxzhi(3))
    double precision, intent(in   ) :: betayz(betayzlo(1):betayzhi(1),betayzlo(2):betayzhi(2),betayzlo(3):betayzhi(3))
#endif

    double precision, intent(in   ) :: rho(rholo(1):rhohi(1),rholo(2):rhohi(2),rholo(3):rhohi(3))

    double precision, intent(inout) :: sourcex(sourcexlo(1):sourcexhi(1),sourcexlo(2):sourcexhi(2),sourcexlo(3):sourcexhi(3))
    double precision, intent(inout) :: sourcey(sourceylo(1):sourceyhi(1),sourceylo(2):sourceyhi(2),sourceylo(3):sourceyhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(inout) :: sourcez(sourcezlo(1):sourcezhi(1),sourcezlo(2):sourcezhi(2),sourcezlo(3):sourcezhi(3))
#endif

    integer l,i,j,k, pindex(3)
    type(f_particle), pointer :: p
    double precision drag(3)
    double precision fluidvel(3)
    double precision dxinv(3)

    double precision cx(3)

    dxinv = 1.0d0/dx


    
    do l = 1, np
       
      p => particles(l)

#if (BL_SPACEDIM == 3)
      !Find cell index, and position within cell, of particle
      !This can probably be further optimised

      cx(1) = (p%pos(1) - real_lo(1))*dxinv(1)
      cx(2) = (p%pos(2) - real_lo(2))*dxinv(2)
#if (BL_SPACEDIM == 3)
      cx(3) = (p%pos(3) - real_lo(3))*dxinv(3)
#endif

      i = ceiling(cx(1))
      j = ceiling(cx(2))
#if (BL_SPACEDIM == 3)
     k = ceiling(cx(3))
#endif

      cx(1) = p%pos(1) - (i-1)*dx(1)
      cx(2) = p%pos(2) - (j-1)*dx(2)
#if (BL_SPACEDIM == 3)
      cx(3) = p%pos(3) - (k-1)*dx(3)
#endif
        
      !Abort if particle is out of bounds
      if (i .lt. index_lo(1) .or. i .gt. index_hi(1) .or. &
          j .lt. index_lo(2) .or. j .gt. index_hi(1) .or. &
          k .lt. index_lo(3) .or. k .gt. index_hi(1)) then
            print *,'PARTICLE ID ', p%id,'OUT OF BOUNDS: ',i,j,k
            print *,'Array bounds: ', index_lo(:), index_hi(:)
            print *,'Position: ', p%pos(1), p%pos(2), p%pos(3)
            call bl_error('Aborting in update_particles')
      end if

      !First order one directional interpolation
      fluidvel(1) = velx(i,j,k)
#endif
!       drag(1) = p%drag_factor*p%fluid_viscosity*(p%vel(1)-p%fluid_vel(1))
!       drag(2) = p%drag_factor*p%fluid_viscosity*(p%vel(2)-p%fluid_vel(2))
!       drag(3) = p%drag_factor*p%fluid_viscosity*(p%vel(3)-p%fluid_vel(3))

!       p%vel(1) = p%vel(1) + drag(1)*dt 
!       p%vel(2) = p%vel(2) + drag(2)*dt 
!       p%vel(3) = p%vel(3) + drag(3)*dt 

!       p%pos(1) = p%pos(1) + p%vel(1)*dt 
!       p%pos(2) = p%pos(2) + p%vel(2)*dt 
#if (BL_SPACEDIM == 3) 
!       p%pos(3) = p%pos(3) + p%vel(3)*dt
#endif

       !Remove division from this
       !Cell indicies use fortran convention

!       p%cell_index(1) = floor((p%pos(1) - real_lo(1))/dx(1)) + 1
!       p%cell_index(2) = floor((p%pos(2) - real_lo(2))/dx(2)) + 1
#if (BL_SPACEDIM == 3) 
!       p%cell_index(3) = floor((p%pos(3) - real_lo(3))/dx(3)) + 1
#endif

    end do

  end subroutine update_particles
  
end module particle_functions_module

