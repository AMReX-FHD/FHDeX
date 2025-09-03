      subroutine init_part(u,part,xl,jpartl,jpartr,num_part,dx,
     1         ndim,pdim)

      implicit none
      integer nloc,ndim,pdim,jpartl,jpartr, num_part

      double precision u(0:ndim+1), part(pdim),xl(ndim+1)
      double precision dx, rn_unif

      integer j,npart,n,k

      parameter (nloc = 2000)


      num_part = 0

      do j=jpartl,jpartr

         call random_number(rn_unif)
c        npart = u(j)*dx+0.5d0
         npart = u(j)*dx+rn_unif

c        write(6,*)j,u(j),u(j)*dx,rn_unif,npart

         if(npart.gt.0)then
            do n=1,npart

               num_part = num_part + 1
               call random_number(rn_unif)
               part(num_part) = xl(j)+rn_unif*dx
c              write(6,*)n,num_part,xl(j),rn_unif,part(num_part)

            enddo
         endif

       enddo

      return
      end
      subroutine avg_down(u,part,xl,jpartl,jpartr,num_part,dx,
     1         ndim,pdim)

      implicit none
      integer nloc,ndim,pdim,jpartl,jpartr, num_part

      double precision u(0:ndim+1), part(pdim),xl(ndim+1)
      double precision dx

      integer j,npart,n,bin

      parameter (nloc = 2000)


      do j=jpartl,jpartr

         u(j) = 0.d0

      enddo

c     write(6,*) num_part,dx

      do n=1,num_part

          bin = part(n)/dx + 1

          if(bin.ge.jpartl .and. bin .le.jpartr)then

              u(bin) = u(bin) + 1.d0

          else

c             write(6,*)" particle ",n," at ", part(n)," bin ",bin,
c    1               jpartl,jpartr
          endif

       enddo

       do j=jpartl,jpartr

           u(j) = u(j)/dx


       enddo

      return
      end
      subroutine fill_part(u,part,xl,jpartl,jpartr,num_part,dx,
     1         npts,ndim,pdim,npghost)

      implicit none
      integer nloc,ndim,pdim,jpartl,jpartr, num_part
      integer npts,npghost

      double precision u(0:ndim+1), part(pdim),xl(ndim+1)
      double precision dx, rn_unif

      integer j,npart,n,k

      parameter (nloc = 2000)

c doesn't fill bc for reservoir.  do that separately

c     write(6,*)" interior particles ", num_part
c     k = num_part

      if(jpartl.gt.npghost .and. jpartr .le.npts-npghost)then

      do j=jpartl-npghost,jpartl-1

         call random_number(rn_unif)
         npart = u(j)*dx+rn_unif
c        npart = u(j)*dx+0.5d0

         if(npart.gt.0)then
            do n=1,npart

               num_part = num_part + 1
               call random_number(rn_unif)
               part(num_part) = xl(j)+rn_unif*dx

            enddo
         endif

       enddo

      do j=jpartr+1,jpartr+npghost

         call random_number(rn_unif)
         npart = u(j)*dx+rn_unif
c        npart = u(j)*dx+0.5d0

         if(npart.gt.0)then
            do n=1,npart

               num_part = num_part + 1
               call random_number(rn_unif)
               part(num_part) = xl(j)+rn_unif*dx

            enddo
         endif

       enddo

      endif

c     write(6,*)" number total, ghost ", num_part, num_part-k

      return
      end
      subroutine adv_part(part,part_new,xl,jpartl,jpartr,crossl,crossr,
     1    num_part,dx,dt,npts,ndim,pdim)

      implicit none
      integer nloc,ndim,pdim,jpartl,jpartr, num_part
      integer npts,npghost

      double precision part(pdim),part_new(pdim),xl(ndim+1)
      double precision dx,dt,crossl,crossr,inc,normal

      integer j,npart,n,ibad

      ibad = 0

      crossl = 0.d0
      crossr = 0.d0

      do n=1,num_part

          call blinvnormdist(normal)

          inc = normal*sqrt(dt)

          if(inc .gt.dx)then
              ibad = ibad + 1
c             write(6,*)"bad point ",ibad,inc,dx,inc/dx
              inc = dx
          endif
          part_new(n) = part(n) + inc

       enddo

       if(ibad.gt.0)then
           write(6,*)ibad," of ", num_part," took step gt dx"
       endif


      if(jpartl.gt.1 .and. jpartr.le.npts)then

          do n=1,num_part

              if(part_new(n).ge.xl(jpartl)
     1            .and. part(n).lt.xl(jpartl))then

                 crossl = crossl-1.d0

              endif
              if(part_new(n).lt.xl(jpartl)
     1            .and. part(n).ge.xl(jpartl))then

                 crossl = crossl+1.d0

              endif
              if(part_new(n).ge.xl(jpartr+1)
     1            .and. part(n).lt.xl(jpartr+1))then

                 crossr = crossr-1.d0

              endif
              if(part_new(n).lt.xl(jpartr+1)
     1            .and. part(n).ge.xl(jpartr+1))then

                 crossr = crossr+1.d0

              endif

          enddo

      endif

c     write(6,*)" number of particles crossing ", crossl, crossr


      return
      end
      subroutine resort(part,part_new,xl,jpartl,jpartr,
     1    num_part,num_part_new,dx,dt,npts,iper,ndim,pdim)

      implicit none
      integer nloc,ndim,pdim,jpartl,jpartr, num_part, num_part_new
      integer npts,npghost, iper

      double precision part(pdim),part_new(pdim),xl(ndim+1)
      double precision dx,dt,crossl,crossr

      integer j,npart,n,ibad

      num_part_new = 0

      do n=1,num_part

         if(part(n).ge.xl(jpartl) .and. part(n).lt.xl(jpartr+1))then

             num_part_new = num_part_new + 1
             part_new(num_part_new) = part(n)

         elseif(iper.eq.1 .and. jpartl .eq. 1 .and.
     1           part(n).lt.xl(1))then

             num_part_new = num_part_new + 1
             part_new(num_part_new) = part(n)+xl(npts+1)

         elseif(iper.eq.1 .and. jpartr .eq. npts .and.
     1           part(n).ge.xl(npts+1))then

             num_part_new = num_part_new + 1
             part_new(num_part_new) = part(n)-xl(npts+1)

         endif

       enddo

c      write(6,*)"before and after ", num_part,num_part_new,
c    1        num_part_new-num_part

       return
       end
       subroutine reflux(u,fluxregl,fluxregr,crossl,crossr,
     1          jpartl,jpartr,npts,dx, ndim)

      implicit none
      integer nloc,ndim,jpartl,jpartr,npts

      double precision u(0:ndim+1)
      double precision dx,crossl,crossr,fluxregl,fluxregr

      integer j,npart,n,bin

      if(jpartl .gt. 1)then

         u(jpartl-1) = u(jpartl-1)+ (crossl-fluxregl)/dx

      endif
      if(jpartr .lt. npts)then

         u(jpartr+1) = u(jpartr+1)- (crossr-fluxregr)/dx

      endif

      return
      end
