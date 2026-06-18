      program GinzLand

      USE Random_Numbers
      implicit none

      integer ndim
      parameter(ndim = 256)

      double precision phi(0:ndim+1,0:ndim+1)
      double precision phin(0:ndim+1,0:ndim+1)
      double precision func(0:ndim+1,0:ndim+1)
      double precision rnums(0:ndim+1,0:ndim+1)

      integer i,j,n,imax,jmax,nstep,ntherm,nfreq

      double precision dx,dy,dt, value, cfl, time, factor

      double precision alpha,gamma, k, a,b,c,d
      double precision phi_init, phi_0, integral
      double precision mean,std, ave, pi, x , y, rad
      double precision maxphi,minphi


      write(6,*)" imax,jmax,cfl,nstep,ntherm,nfreq"
      read(5,*) imax,jmax,cfl,nstep,ntherm,nfreq
      write(6,*) imax,jmax,cfl,nstep,ntherm,nfreq
      write(6,*)"gamma, alpha, k"
      read(5,*)gamma, alpha, k
      write(6,*)gamma, alpha, k
      write(6,*)"phi_init, phi_0"
      read(5,*)phi_init, phi_0
      write(6,*)phi_init, phi_0
      write(6,*)"a,b,c,d"
      read(5,*)a,b,c,d
      write(6,*)a,b,c,d
      read(5,*)rad

      call RandomSeeds(1905849284)

      dx = 1.d0/dfloat(imax)
      dy = 1.d0/dfloat(jmax)

      dt = 0.25d0*cfl*dx**2 / gamma

      factor = sqrt(2.d0*alpha/(dt*dx*dy))

      phi = phi_init

      pi= 4.d0*atan2(1.d0,1.d0)

      do j=1,jmax
      do i=1,jmax

         x = dx*dfloat(i-1)
         y = dy*dfloat(j-1)

      !  phi(i,j) =  (sin(2.d0*pi*x)*sin(2.d0*pi*y)) **2
         if( (x-.5d0)**2 + (y-.5d0)**2 .lt. rad**2)then
            phi(i,j) = 1.d0
         else
            phi(i,j) = 0.d0
         endif

      enddo
      enddo


      mean = 0.d0
      std = 0.d0

      do n=1,ntherm + nstep

         time = n*dt

         call periodic(phi,imax,jmax,ndim)
         call integrate(phi, phi_0, integral,dx,dy,imax,jmax,ndim)

         ave = 0.d0
         maxphi = phi(1,1)
         minphi = phi(1,1)

         do j=1,jmax
            do i=1,imax

             func(i,j) = a+2.d0*b*phi(i,j)+3.d0*c*phi(i,j)**2+4.d0*d*phi(i,j)**3
             call RandomNormal(value)
             rnums(i,j) = value

             mean = mean + value
             std = std + value**2
             ave= ave+phi(i,j)
             maxphi = max(maxphi, phi(i,j))
             minphi = min(minphi, phi(i,j))

           enddo
         enddo

         mean = mean / dfloat(imax*jmax)
         std = std / dfloat(imax*jmax)
         ave = ave / dfloat(imax*jmax)
         std = sqrt(std - mean**2)

         write(6,*)" average, max, min", ave,maxphi,minphi

         do j=1,jmax
           do i= 1,imax

             phin(i,j) = phi(i,j) + dt*gamma*(phi(i+1,j)-2.d0*phi(i,j)+phi(i-1,j))/dx**2  &
              + dt*gamma*(phi(i,j+1)-2.d0*phi(i,j)+phi(i,j-1))/dy**2  &
              -dt*func(i,j)-dt*k*integral + dt*factor * rnums(i,j)

           enddo
         enddo

       do j=1,jmax
       do i=1,imax

          phi(i,j) = phin(i,j)

      enddo
      enddo

      call output(1,imax,1,jmax,phi,time,ndim)

      if(n.gt.ntherm .and. mod(n,nfreq).eq.0)then

         write(8,*)"time=",time
         do j=1,jmax
         write(8,*) (phi(i,j),i=1,imax)
         enddo
      endif

      enddo

      end

      subroutine periodic(phi,imax,jmax,ndim)

          integer imax,jmax,ndim
          double precision phi(0:ndim+1,0:ndim+1)

        do j=1,jmax

           phi(0,j) = phi(imax,j)
           phi(imax+1,j) = phi(1,j)

        enddo

        do i=0,imax+1

           phi(i,0) = phi(i,jmax)
           phi(i,jmax+1) = phi(i,1)

        enddo

       return

       end subroutine periodic


      subroutine integrate(phi, phi_0, integral,dx,dy,imax,jmax,ndim)


          integer imax,jmax,ndim
          double precision phi(0:ndim+1,0:ndim+1)

          double precision phi_0, integral, dx, dy

          integral = 0.d0

          do j=1,jmax
          do i=1,imax

              integral = integral + phi(i,j) - phi_0

          enddo
          enddo

          integral = integral * dx * dy

      return
      end subroutine integrate
