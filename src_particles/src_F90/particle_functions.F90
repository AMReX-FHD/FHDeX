module particle_functions_module

  use amrex_error_module
  use amrex_fort_module
  use amrex_constants_module
  use common_namelist_module
  !use bl_random_module

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
     real(amrex_particle_real)   :: accel_factor


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
                                       beta, betalo, betahi, &

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
    integer,          intent(in   )         :: betalo(3), betahi(3)
#if (AMREX_SPACEDIM == 3)
    integer,          intent(in   )         :: velzlo(3), velzhi(3), sourcezlo(3), sourcezhi(3), coordszlo(3), coordszhi(3)
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

    double precision, intent(in   ) :: beta(betalo(1):betahi(1),betalo(2):betahi(2),betalo(3):betahi(3))

    double precision, intent(in   ) :: rho(rholo(1):rhohi(1),rholo(2):rhohi(2),rholo(3):rhohi(3))

    double precision, intent(inout) :: sourcex(sourcexlo(1):sourcexhi(1),sourcexlo(2):sourcexhi(2),sourcexlo(3):sourcexhi(3))
    double precision, intent(inout) :: sourcey(sourceylo(1):sourceyhi(1),sourceylo(2):sourceyhi(2),sourceylo(3):sourceyhi(3))
#if (AMREX_SPACEDIM == 3)
    double precision, intent(inout) :: sourcez(sourcezlo(1):sourcezhi(1),sourcezlo(2):sourcezhi(2),sourcezlo(3):sourcezhi(3))
#endif

    integer l,i,j,k
    type(f_particle), pointer :: p
    double precision localvel(3)
    double precision localbeta, deltap(3), nodalp
    double precision dxinv(3)
    double precision onemxd(3)

#if (BL_SPACEDIM == 3)
    double precision c000,c001,c010,c011,c100,c101,c110,c111, ctotal
    double precision r000,r001,r010,r011,r100,r101,r110,r111
#endif

#if (BL_SPACEDIM == 2)
    double precision c00,c01,c10,c11, ctotal
    double precision r00,r01,r10,r11
#endif

    double precision xd(3)

    dxinv = 1.0d0/dx
    onemxd = 1.0d0-dx

    !Zero source terms
    do k = sourcexlo(3), sourcexhi(3)
      do j = sourcexlo(2), sourcexhi(2)
        do i = sourcexlo(1), sourcexhi(1)

          sourcex(i,j,k) = 0;

        enddo
      enddo
    enddo

    do k = sourceylo(3), sourceyhi(3)
      do j = sourceylo(2), sourceyhi(2)
        do i = sourceylo(1), sourceyhi(1)

          sourcey(i,j,k) = 0;

        enddo
      enddo
    enddo

#if (BL_SPACEDIM == 3)
    do k = sourcezlo(3), sourcezhi(3)
      do j = sourcezlo(2), sourcezhi(2)
        do i = sourcezlo(1), sourcezhi(1)

          sourcez(i,j,k) = 0;

        enddo
      enddo
    enddo
#endif

    
    do l = 1, np
       
      p => particles(l)

      !Find cell index, and position within cell, of particle
      !This can probably be further optimised

      i = floor((p%pos(1) - real_lo(1))*dxinv(1))
      j = floor((p%pos(2) - real_lo(2))*dxinv(2))
#if (BL_SPACEDIM == 3)
      k = floor((p%pos(3) - real_lo(3))*dxinv(3))
#else
      k = 0
#endif

      p%cell_index(1) = i;
      p%cell_index(2) = j;
#if (BL_SPACEDIM == 3)
      p%cell_index(3) = k;
#endif
        
#if (BL_SPACEDIM == 3)
      !Abort if particle is out of bounds
      !if (i .lt. index_lo(1) .or. i .gt. index_hi(1) .or. &
      !    j .lt. index_lo(2) .or. j .gt. index_hi(1) .or. &
      !    k .lt. index_lo(3) .or. k .gt. index_hi(1)) then
      !      print *,'PARTICLE ID ', p%id,'OUT OF BOUNDS: ',i,j,k
      !      print *,'Array bounds: ', index_lo(:), index_hi(:)
      !      print *,'Position: ', p%pos(1), p%pos(2), p%pos(3)
      !      call bl_error('Aborting in update_particles')
      !end if
#endif

      !Interpolate fields

      xd(1) = (p%pos(1) - coordsx(i,j,k,1))*dxInv(1)
      xd(2) = (p%pos(2) - coordsy(i,j,k,2))*dxInv(2)
#if (BL_SPACEDIM == 3)
      xd(3) = (p%pos(3) - coordsx(i,j,k,3))*dxInv(3)
#endif

#if (BL_SPACEDIM == 3)

      c000 = onemxd(1)*onemxd(2)*onemxd(3)
      c001 = onemxd(1)*onemxd(2)*xd(3)
      c010 = onemxd(1)*onemxd(3)*xd(2)
      c011 = onemxd(1)*xd(2)*xd(3)
      c101 = onemxd(2)*onemxd(3)*xd(1)
      c100 = onemxd(2)*xd(1)*xd(3)
      c110 = onemxd(3)*xd(1)*xd(2)
      c111 = xd(1)*xd(2)*xd(3)

      ctotal = (c000 + c001 + c010 + c011 + c100 + c101 + c110 + c111)*4.0

      r000 = c000/ctotal
      r001 = c001/ctotal
      r010 = c010/ctotal
      r011 = c011/ctotal
      r100 = c100/ctotal
      r101 = c101/ctotal
      r110 = c110/ctotal
      r111 = c111/ctotal

      if (visc_type .gt. 0) then
        localbeta = beta(i,j,k)
      else
        !3d visc
        localbeta = c000*beta(i,j,k) + c001*beta(i,j,k+1) + c010*beta(i,j+1,k) + c011*beta(i,j+1,k+1) + c100*beta(i+1,j,k) + c101*beta(i+1,j,k+1) + c110*beta(i+1,j+1,k) + c111*beta(i+1,j+1,k+1)
      endif

      localvel(1) = c000*velx(i,j,k) + c001*velx(i,j,k+1) + c010*velx(i,j+1,k) + c011*velx(i,j+1,k+1) + c100*velx(i+1,j,k) + c101*velx(i+1,j,k+1) + c110*velx(i+1,j+1,k) + c111*velx(i+1,j+1,k+1)
      localvel(2) = c000*vely(i,j,k) + c001*vely(i,j,k+1) + c010*vely(i,j+1,k) + c011*vely(i,j+1,k+1) + c100*vely(i+1,j,k) + c101*vely(i+1,j,k+1) + c110*vely(i+1,j+1,k) + c111*vely(i+1,j+1,k+1)
      localvel(3) = c000*velz(i,j,k) + c001*velz(i,j,k+1) + c010*velz(i,j+1,k) + c011*velz(i,j+1,k+1) + c100*velz(i+1,j,k) + c101*velz(i+1,j,k+1) + c110*velz(i+1,j+1,k) + c111*velz(i+1,j+1,k+1)

#endif

#if (BL_SPACEDIM == 2)

      c00 = onemxd(1)*onemxd(2)
      c01 = onemxd(1)*xd(2)
      c10 = xd(1)*onemxd(2)
      c11 = xd(1)*xd(2)

      ctotal = (c00 + c01 + c10 + c11)*2

      r00 = c00/ctotal
      r01 = c01/ctotal
      r10 = c10/ctotal
      r11 = c11/ctotal

      if (visc_type .gt. 0) then
        localbeta = beta(i,j,k)
      else
        localbeta = beta(i,j,k)*c00 + beta(i,j+1,k)*c01 + beta(i+1,j,k)*c10 + beta(i+1,j+1,k)*c11
      endif
      !2d xvel
      localvel(1) = velx(i,j,k)*c00 + velx(i,j+1,k)*c01 + velx(i+1,j,k)*c10 + velx(i+1,j+1,k)*c11
      localvel(2) = vely(i,j,k)*c00 + vely(i,j+1,k)*c01 + vely(i+1,j,k)*c10 + vely(i+1,j+1,k)*c11

#endif

      !Semi-implicit Euler velocity and position update

      deltap(1) = p%vel(1)
      deltap(2) = p%vel(2)
#if (BL_SPACEDIM == 3)
      deltap(3) = p%vel(3)
#endif

      !print *, velx(i,j,k)
      !print *, p%pos
      p%vel(1) = -p%accel_factor*localbeta*(p%vel(1)-localvel(1))*dt + p%vel(1)
      p%vel(2) = -p%accel_factor*localbeta*(p%vel(2)-localvel(2))*dt + p%vel(2)
#if (BL_SPACEDIM == 3)
      p%vel(3) = -p%accel_factor*localbeta*(p%vel(3)-localvel(3))*dt + p%vel(3)
#endif
      deltap(1) = p%mass*(p%vel(1) - deltap(1))
      deltap(2) = p%mass*(p%vel(2) - deltap(2))
#if (BL_SPACEDIM == 3)
      deltap(3) = p%mass*(p%vel(3) - deltap(3))
#endif

      p%pos(1) = p%pos(1) + p%vel(1)*dt 
      p%pos(2) = p%pos(2) + p%vel(2)*dt 
#if (BL_SPACEDIM == 3) 
      p%pos(3) = p%pos(3) + p%vel(3)*dt
#endif

#if (BL_SPACEDIM == 3)
      !distribute x momentum change 
       nodalp = r000*deltap(1)
       sourcex(i,j,k) = sourcex(i,j,k) + nodalp
       sourcex(i,j-1,k) = sourcex(i,j-1,k) + nodalp
       sourcex(i,j,k-1) = sourcex(i,j,k-1) + nodalp
       sourcex(i,j-1,k-1) = sourcex(i,j-1,k-1) + nodalp

       nodalp = r001*deltap(1)
       sourcex(i,j,k+1) = sourcex(i,j,k+1) + nodalp
       sourcex(i,j-1,k+1) = sourcex(i,j-1,k+1) + nodalp
       sourcex(i,j,k) = sourcex(i,j,k) + nodalp
       sourcex(i,j-1,k) = sourcex(i,j-1,k) + nodalp

       nodalp = r010*deltap(1)
       sourcex(i,j+1,k) = sourcex(i,j+1,k) + nodalp
       sourcex(i,j,k) = sourcex(i,j,k) + nodalp
       sourcex(i,j+1,k-1) = sourcex(i,j+1,k-1) + nodalp
       sourcex(i,j,k-1) = sourcex(i,j,k-1) + nodalp

       nodalp = r011*deltap(1)
       sourcex(i,j+1,k+1) = sourcex(i,j+1,k+1) + nodalp
       sourcex(i,j,k+1) = sourcex(i,j,k+1) + nodalp
       sourcex(i,j+1,k) = sourcex(i,j+1,k) + nodalp
       sourcex(i,j,k) = sourcex(i,j,k) + nodalp

       nodalp = r100*deltap(1)
       sourcex(i+1,j,k) = sourcex(i+1,j,k) + nodalp
       sourcex(i+1,j-1,k) = sourcex(i+1,j-1,k) + nodalp
       sourcex(i+1,j,k-1) = sourcex(i+1,j,k-1) + nodalp
       sourcex(i+1,j-1,k-1) = sourcex(i+1,j-1,k-1) + nodalp

       nodalp = r101*deltap(1)
       sourcex(i+1,j,k+1) = sourcex(i+1,j,k+1) + nodalp
       sourcex(i+1,j-1,k+1) = sourcex(i+1,j-1,k+1) + nodalp
       sourcex(i+1,j,k) = sourcex(i+1,j,k) + nodalp
       sourcex(i+1,j-1,k) = sourcex(i+1,j-1,k) + nodalp

       nodalp = r110*deltap(1)
       sourcex(i+1,j+1,k) = sourcex(i+1,j+1,k) + nodalp
       sourcex(i+1,j,k) = sourcex(i+1,j,k) + nodalp
       sourcex(i+1,j+1,k-1) = sourcex(i+1,j+1,k-1) + nodalp
       sourcex(i+1,j,k-1) = sourcex(i+1,j,k-1) + nodalp

       nodalp = r111*deltap(1)
       sourcex(i+1,j+1,k+1) = sourcex(i+1,j+1,k+1) + nodalp
       sourcex(i+1,j,k+1) = sourcex(i+1,j,k+1) + nodalp
       sourcex(i+1,j+1,k) = sourcex(i+1,j+1,k) + nodalp
       sourcex(i+1,j,k) = sourcex(i+1,j,k) + nodalp

      !distribute y momentum change 
       nodalp = r000*deltap(2)
       sourcey(i,j,k) = sourcey(i,j,k) + nodalp
       sourcey(i-1,j,k) = sourcey(i-1,j,k) + nodalp
       sourcey(i,j,k-1) = sourcey(i,j,k-1) + nodalp
       sourcey(i-1,j,k-1) = sourcey(i-1,j,k-1) + nodalp

       nodalp = r001*deltap(2)
       sourcey(i,j,k+1) = sourcey(i,j,k+1) + nodalp
       sourcey(i-1,j,k+1) = sourcey(i-1,j,k+1) + nodalp
       sourcey(i,j,k) = sourcey(i,j,k) + nodalp
       sourcey(i-1,j,k) = sourcey(i-1,j,k) + nodalp

       nodalp = r010*deltap(2)
       sourcey(i,j+1,k) = sourcey(i,j+1,k) + nodalp
       sourcey(i-1,j+1,k) = sourcey(i-1,j+1,k) + nodalp
       sourcey(i,j+1,k-1) = sourcey(i,j+1,k-1) + nodalp
       sourcey(i-1,j+1,k-1) = sourcey(i-1,j+1,k-1) + nodalp

       nodalp = r011*deltap(2)
       sourcey(i,j+1,k+1) = sourcey(i,j+1,k+1) + nodalp
       sourcey(i-1,j+1,k+1) = sourcey(i-1,j+1,k+1) + nodalp
       sourcey(i,j+1,k) = sourcey(i,j+1,k) + nodalp
       sourcey(i-1,j+1,k) = sourcey(i-1,j+1,k) + nodalp

       nodalp = r100*deltap(2)
       sourcey(i+1,j,k) = sourcey(i+1,j,k) + nodalp
       sourcey(i,j,k) = sourcey(i,j,k) + nodalp
       sourcey(i+1,j,k-1) = sourcey(i+1,j,k-1) + nodalp
       sourcey(i,j,k-1) = sourcey(i,j,k-1) + nodalp

       nodalp = r101*deltap(2)
       sourcey(i+1,j,k+1) = sourcey(i+1,j,k+1) + nodalp
       sourcey(i,j,k+1) = sourcey(i,j,k+1) + nodalp
       sourcey(i+1,j,k) = sourcey(i+1,j,k) + nodalp
       sourcey(i,j,k) = sourcey(i,j,k) + nodalp

       nodalp = r110*deltap(2)
       sourcey(i+1,j+1,k) = sourcey(i+1,j+1,k) + nodalp
       sourcey(i,j+1,k) = sourcey(i,j+1,k) + nodalp
       sourcey(i+1,j+1,k-1) = sourcey(i+1,j+1,k-1) + nodalp
       sourcey(i,j+1,k-1) = sourcey(i,j+1,k-1) + nodalp

       nodalp = r111*deltap(2)
       sourcey(i+1,j+1,k+1) = sourcey(i+1,j+1,k+1) + nodalp
       sourcey(i,j+1,k+1) = sourcey(i,j+1,k+1) + nodalp
       sourcey(i+1,j+1,k) = sourcey(i+1,j+1,k) + nodalp
       sourcey(i,j+1,k) = sourcey(i,j+1,k) + nodalp

      !distribute z momentum change 
       nodalp = r000*deltap(3)
       sourcez(i,j,k) = sourcez(i,j,k) + nodalp
       sourcez(i-1,j,k) = sourcez(i-1,j,k) + nodalp
       sourcez(i,j-1,k) = sourcez(i,j-1,k) + nodalp
       sourcez(i-1,j-1,k) = sourcez(i-1,j-1,k) + nodalp

       nodalp = r001*deltap(3)
       sourcez(i,j,k+1) = sourcez(i,j,k+1) + nodalp
       sourcez(i-1,j,k+1) = sourcez(i-1,j,k+1) + nodalp
       sourcez(i,j-1,k+1) = sourcez(i,j-1,k+1) + nodalp
       sourcez(i-1,j-1,k+1) = sourcez(i-1,j-1,k+1) + nodalp

       nodalp = r010*deltap(3)
       sourcez(i,j+1,k) = sourcez(i,j+1,k) + nodalp
       sourcez(i-1,j+1,k) = sourcez(i-1,j+1,k) + nodalp
       sourcez(i,j,k) = sourcez(i,j,k) + nodalp
       sourcez(i-1,j,k) = sourceZ(i-1,j,k) + nodalp

       nodalp = r011*deltap(3)
       sourcez(i,j+1,k+1) = sourcez(i,j+1,k+1) + nodalp
       sourcez(i-1,j+1,k+1) = sourcez(i-1,j+1,k+1) + nodalp
       sourcez(i,j,k+1) = sourcez(i,j,k+1) + nodalp
       sourcez(i-1,j,k+1) = sourcez(i-1,j,k+1) + nodalp

       nodalp = r100*deltap(3)
       sourcez(i+1,j,k) = sourcez(i+1,j,k) + nodalp
       sourcez(i,j,k) = sourcez(i,j,k) + nodalp
       sourcez(i+1,j-1,k) = sourcez(i+1,j-1,k) + nodalp
       sourcez(i,j-1,k) = sourcez(i,j-1,k) + nodalp

       nodalp = r101*deltap(3)
       sourcez(i+1,j,k+1) = sourcez(i+1,j,k+1) + nodalp
       sourcez(i,j,k+1) = sourcez(i,j,k+1) + nodalp
       sourcez(i+1,j-1,k+1) = sourcez(i+1,j-1,k+1) + nodalp
       sourcez(i,j-1,k+1) = sourcez(i,j-1,k+1) + nodalp

       nodalp = r110*deltap(3)
       sourcez(i+1,j+1,k) = sourcez(i+1,j+1,k) + nodalp
       sourcez(i,j+1,k) = sourcez(i,j+1,k) + nodalp
       sourcez(i+1,j,k) = sourcez(i+1,j,k) + nodalp
       sourcez(i,j,k) = sourcez(i,j,k) + nodalp

       nodalp = r111*deltap(3)
       sourcez(i+1,j+1,k+1) = sourcez(i+1,j+1,k+1) + nodalp
       sourcez(i,j+1,k+1) = sourcez(i,j+1,k+1) + nodalp
       sourcez(i+1,j,k+1) = sourcez(i+1,j,k+1) + nodalp
       sourcez(i,j,k+1) = sourcez(i,j,k+1) + nodalp
#endif

#if (BL_SPACEDIM == 2)
      !distribute x momentum change 
       nodalp = r00*deltap(1)
       sourcex(i,j,k) = sourcex(i,j,k) + nodalp
       sourcex(i,j-1,k) = sourcex(i,j-1,k) + nodalp

       nodalp = r01*deltap(1)
       sourcex(i,j+1,k) = sourcex(i,j+1,k) + nodalp
       sourcex(i,j,k) = sourcex(i,j,k) + nodalp

       nodalp = r10*deltap(1)
       sourcex(i+1,j,k) = sourcex(i+1,j,k) + nodalp
       sourcex(i+1,j-1,k) = sourcex(i+1,j-1,k) + nodalp

       nodalp = r11*deltap(1)
       sourcex(i+1,j+1,k) = sourcex(i+1,j+1,k) + nodalp
       sourcex(i+1,j,k) = sourcex(i+1,j,k) + nodalp

      !distribute y momentum change
       nodalp = r00*deltap(2)
       sourcey(i,j,k) = sourcey(i,j,k) + nodalp
       sourcey(i-1,j,k) = sourcey(i-1,j,k) + nodalp

       nodalp = r01*deltap(2)
       sourcey(i,j+1,k) = sourcey(i,j+1,k) + nodalp
       sourcey(i-1,j+1,k) = sourcey(i-1,j+1,k) + nodalp

       nodalp = r10*deltap(2)
       sourcey(i+1,j,k) = sourcey(i+1,j,k) + nodalp
       sourcey(i,j,k) = sourcey(i,j,k) + nodalp

       nodalp = r11*deltap(2)
       sourcey(i+1,j+1,k) = sourcey(i+1,j+1,k) + nodalp
       sourcey(i,j+1,k) = sourcey(i,j+1,k) + nodalp

#endif

    end do

  end subroutine update_particles
  
end module particle_functions_module































