      subroutine compflux(flux,ranflux,u,npts,weight,
     1       dx,dt,dorand,ndim,iper,ires)

      implicit none
      integer nloc,ndim,iper,npts,ires

      parameter (nloc = 2000)

      double precision u(0:ndim+1), flux(1:ndim+1),ranflux(1:ndim+1,2)
      double precision dx,dt,dorand,weight,varflux
      double precision uave,uplus,umin, coeff ,factor


      integer i,j
      integer jend,jstart
      integer is_edge

      double precision  uleft, uright, pi

      common /params/ uleft, uright, pi



        if(iper.eq.1)then

           u(0) = u(npts)
           u(npts+1) = u(1)

        elseif (ires.eq.1)then

             u(0) = uleft
             u(npts+1) = uright
        else
           u(0)      = u(1)
           u(npts+1) = u(npts)



         endif

         do j=1,npts+1

            flux(j) = 0.5d0*(u(j)-u(j-1))/dx

         enddo

         if(iper.eq.1)then

               flux(npts+1) = flux(1)

          elseif(ires.eq.0)then

              flux(1) = 0.d0
              flux(npts+1) = 0.d0

          endif



      do j=1,npts +1


         umin =  max(u(j-1),0.d0)
         uplus = max(u(j),0.d0)
c        uave = 2.d0*(uplus*umin)/(uplus+umin+1.d-16)
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

         varflux = factor*sqrt(uave)


c        flux(j) = flux(j)+varflux*(ranflux(j,1)+weight*ranflux(j,2))
         flux(j) = flux(j)+varflux*ranflux(j,1)
     1             *dorand/sqrt(dx*dt)

      enddo

        if(iper.eq.1)then

             flux(npts+1) = flux(1)

        endif


      return
      end
