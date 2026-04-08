      subroutine compflux(flux,ranflux,u,npts,weight,
     1       dx,dt,dorand,ndim,iper,ires,ncoef)

      implicit none
      integer nloc,ndim,iper,npts,ires,ncoef

      parameter (nloc = 2000)

      double precision u(0:ndim+1,2), flux(1:ndim+1,2),
     1                 ranflux(1:ndim+1,2)
      double precision dx,dt,dorand,weight,varflux
      double precision uave,uplus,umin, coeff ,factor
      double precision sqrtcoef(0:nloc+2)

      double precision delreg, delnreg, sgnphi,phiabs


      integer i,j,k
      integer jend,jstart
      integer is_edge
      integer usereg

      double precision  uleft, uright, pi, uinit

      common /params/ uleft, uright, pi, uinit

      delreg = 4.d0*dx**2
      delnreg = uinit * delreg
c      delnreg = 10.d0
      usereg = 0


      do k=1,ncoef

        if(iper.eq.1)then

           u(0,k) = u(npts,k)
           u(npts+1,k) = u(1,k)

        elseif (ires.eq.1)then

             u(0,k) = uleft
             u(npts+1,k) = uright
        else
           u(0,k)      = u(1,k)
           u(npts+1,k) = u(npts,k)



         endif

         do j=1,npts+1

            flux(j,k) = 0.5d0*(u(j,k)-u(j-1,k))/dx

         enddo

         if(iper.eq.1)then

               flux(npts+1,k) = flux(1,k)

          elseif(ires.eq.0)then

              flux(1,k) = 0.d0
              flux(npts+1,k) = 0.d0

          else

              flux(1,k) = 2.d0*flux(1,k)
              flux(npts+1,k) = 2.d0*flux(npts+1,k)

          endif

      enddo

      if(usereg .eq. 0)then

        do j=0,npts +1

           sqrtcoef(j) = sqrt(max(u(j,ncoef),0.d0))

        enddo

      else

        do j=0,npts +1

           sgnphi = sign(1.d0,u(j,ncoef))
           phiabs = abs(u(j,ncoef))
           if(phiabs .le. delnreg/2.d0) then
               sqrtcoef(j) = u(j,ncoef)/sqrt(delnreg)
           elseif (phiabs <= delnreg) then
               sqrtcoef(j) = sgnphi*sqrt(delnreg)/2.d0
     1           -1.5*u(j,ncoef)/sqrt(delnreg)
     1           + sgnphi*4.d0*u(j,ncoef)**2/sqrt(delnreg**3)
     1         -2.d0*u(j,ncoef)**3/sqrt(delnreg**5)
           else
              sqrtcoef(j) = sgnphi* sqrt(phiabs)
           endif

c          if(u(j,ncoef) .lt.0.d0)then
c             sqrtcoef(j) = 0.d0
c          endif

        enddo

      endif


      do j=1,npts +1


c        umin =  max(u(j-1,ncoef),0.d0)
c        uplus = max(u(j,ncoef),0.d0)
c        uave = 2.d0*(uplus*umin)/(uplus+umin+1.d-16)

c        umin = sqrt( max(u(j-1,ncoef),0.d0))
c        uplus = sqrt(max(u(j,ncoef),0.d0))
         umin = sqrtcoef(j-1)
         uplus = sqrtcoef(j)
         uave = 0.5*(uplus+umin)

         factor = 1.d0

         if(j.eq.1.or.j.eq.npts+1)then
         if(iper.eq.0)then
             if(ires.eq.1)then
              factor = sqrt(2.d0)
             else
              factor = 0.d0
             endif
         endif
         endif

c        varflux = factor*sqrt(uave)
         varflux = factor*uave


c        flux(j) = flux(j)+varflux*(ranflux(j,1)+weight*ranflux(j,2))
         flux(j,1) = flux(j,1)+varflux*ranflux(j,1)
     1             *dorand/sqrt(dx*dt)

      enddo

        if(iper.eq.1)then

           do k=1,ncoef

             flux(npts+1,k) = flux(1,k)

           enddo

        endif


      return
      end
