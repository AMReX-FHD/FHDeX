       program dk
c
c************************************************************
c
c     this program implements the muscl scheme for the buckley
c
c     by setting the flag igud to 1 the program sets
c     all slope to zero and implements godunov's
c     method
      implicit none

      integer ndim,pdim
      parameter (ndim = 2000)
      parameter (pdim = 2000000)

      double precision part(pdim),part_new(pdim)
c
       double precision flux(1:ndim+1)
       double precision ranflux(1:ndim+1,2)

       double precision u(0:ndim+1),unew(0:ndim+1)
      double precision x(1:ndim),xl(1:ndim+1)


      character*17 uvarfile
      character*10 label


      integer seed,ntherm,nstep,i,j,k,l,iper,iprint,istat,jstart
      integer icor,n,nout,npts,nstat,nsamples,ires
      integer jmidl,jmidr,jpartl,jpartr,is_hybrid
      integer num_part, npart, npghost, num_part_new, is_pure_part

      double precision xlen,dt,dx,factor,time,totmass,dorand,weight,cfl
      double precision minmass,midfact

      double precision crossl,crossr,fluxregl, fluxregr

      double precision  uleft, uright, pi, uinit

      common /params/  uleft, uright, pi

      namelist /input_param/ npts,xlen,iper,dorand,seed,ntherm,nstep,
     1   dt,nout,ires,nstat,icor,ires,uleft,uright, uinit, cfl,
     1   jmidl,jmidr,midfact,jpartl,jpartr,is_hybrid


       xlen = 1.d0
       iper = 1
       ntherm = 0
       dora nd = 1.d0
       seed = 43064
       cfl = -1.d0
       midfact = 1.d0
       jmidl = -5
       jmidr = -5
       jpartl = -5
       jpartr = -5
       is_hybrid = 0
       npghost = 1
       is_pure_part = 0


       read(5,input_param)
       write(6,input_param)

       if(is_hybrid.eq.1)then
          if(jpartl.eq.1 .and. jpartr.eq.npts)then
           is_pure_part = 1
            write(6,*)" pure particle simulation"
          elseif(jpartl.le.npghost .or. jpartr.gt.npts-npghost)then
            write(6,*)" patch to close to boundary"
            stop
          endif
       endif

       pi = 4.d0*atan2(1.d0,1.d0)

       time =0.0
        dx = xlen/dfloat(npts)

        do j=1,npts
          x(j) = (dfloat(j)-0.5d0)*dx
          xl(j) = dfloat(j-1)*dx
        enddo
          xl(npts+1) = xlen

       if(cfl.gt.0.d0)then

          dt = cfl*dx**2

       endif


       write(6,*)"dt,dx = ",dt,dx


c   set initial conditions

        do j=1,npts

           u(j) = uinit

           if(j.ge.jmidl .and. j.le.jmidr)then
              u(j) = u(j)*midfact
           endif

        enddo

            totmass = 0.d0

            do j=1,npts

                totmass = totmass+u(j)

            enddo

            totmass = totmass * dx

           write(6,*)"init totmass before = ",totmass


        if(is_hybrid.eq.1)then

           call init_part(u,part,xl,jpartl,jpartr,num_part,dx,ndim,pdim)
           call avg_down(u,part,xl,jpartl,jpartr,num_part,dx,
     1         ndim,pdim)

        endif

            totmass = 0.d0

            do j=1,npts

                totmass = totmass+u(j)

            enddo

            totmass = totmass * dx

           write(6,*)"init totmass after = ",totmass


        call blinitrand(seed)


c   set up stuff to gather statistics

      nstep = nstep + ntherm
      istat = 0


      call output(x,xl,u,npts,0,0.d0,ndim)

      do 200 n=1,nstep


         time = dfloat(n)*dt

c        if(jpartl.gt.1 .and.jpartr .lt.npts)then
         if(is_hybrid.eq.0 .or .is_pure_part.eq.0)then

         do j=1,npts+1
c        do l=1,2

c           call blinvnormdist(ranflux(j,l))
            call blinvnormdist(ranflux(j,1))

c        enddo
         enddo


         weight =0.d0
c        weight =(2.d0*sqrt(2.d0)+sqrt(3.d0))/5.d0
         call compflux(flux,ranflux,u,npts,weight,
     1        dx,dt,dorand,ndim,iper,ires)



      do  j=1,npts

           unew(j) = u(j)+dt*(flux(j+1)-flux(j))/dx

      enddo

           fluxregl = flux(jpartl)   * dt
           fluxregr = flux(jpartr+1) * dt

      do  j=1,npts

           u(j) = unew(j)

      enddo

      endif

      if(is_hybrid .eq. 1)then

        call fill_part(u,part,xl,jpartl,jpartr,num_part,dx,
     1         npts,ndim,pdim,npghost)

        call adv_part(part,part_new,xl,jpartl,jpartr,crossl,crossr,
     1    num_part,dx,dt,npts,ndim,pdim)

        call resort(part_new,part,xl,jpartl,jpartr,
     1    num_part,num_part_new,dx,dt,npts,iper,ndim,pdim)

        num_part = num_part_new

        call  avg_down(u,part,xl,jpartl,jpartr,num_part,dx,
     1         ndim,pdim)

c       totmass = 0.d0
c       do j=jpartl,jpartr
c           totmass = totmass + u(j)
c       enddo
c       write(6,*)" mass in particle region ",totmass*dx

        call reflux(u,fluxregl,fluxregr,crossl,crossr,
     1          jpartl,jpartr,npts,dx, ndim)

      endif

c         if(iper.eq.0)then
c            conpre(npts+1,2) = 0.d0
c            conpre(0,2) = 0.d0
c         else
c            conpre(npts+1,2) = conpre(1,2)
c            conpre(0,2) = conpre(npts+1,2)
c         endif

c         weight = (-4.d0*sqrt(2.d0)+3.d0*sqrt(3.d0))/5.d0
c         call flux(fluxr,fluxrm,fluxrt1,fluxrt2,
c     1          fluxe,fluxC,conpre,npts,
c     1          rmom,rmomt1,rmomt2,rtemp,rspec,weight,
c     1        dx,dt,dorand,vcell,ndim,iper,ires,iconc,islip,iadiab)
c
c      do  j=1,npts
c
c           conpre2(j,1) = 0.25*(3.d0*con(j,1)+conpre(j,1)
c     1          -dt*(fluxr(j+1)-fluxr(j))/dx)
c      enddo
c

c        weight = (sqrt(2.d0)-2.d0*sqrt(3.d0))/10.d0
c        call flux(fluxr,fluxrm,fluxrt1,fluxrt2,
c    1          fluxe,fluxC,conpre2,npts,
c    1          rmom,rmomt1,rmomt2,rtemp,rspec,weight,
c    1        dx,dt,dorand,vcell,ndim,iper,ires,iconc,islip,iadiab)

c     do  j=1,npts

c          con(j,1) = 1.d0/3.d0*(con(j,1)+2.d0*conpre2(j,1)
c    1          -2.d0*dt*(fluxr(j+1)-fluxr(j))/dx)
c      enddo

         if(mod(n,nout).eq.0)then

            totmass = 0.d0
            minmass = 1.d20

            do j=1,npts

                totmass = totmass+u(j)
                minmass = min(minmass,u(j))

            enddo

            totmass = totmass * dx

           write(6,*)n, "time = ",time, " totmass = ",totmass,
     1            " minmass = ",minmass
        endif

        if(n.gt.ntherm .and. mod(n,nstat).eq.0 .and. nstat .gt.0)then

c     gather statistics

          istat = istat+1

          iprint = 0

        endif


        if(mod(n,nout).eq.0)then

           call output(x,xl,u,npts,n,time,ndim)
           iprint = 1

        if(n.gt.ntherm .and. istat .gt.0)then
           factor = dfloat(istat)
1001    format(7e15.7)
           write(label,'(i10.10)') n
           uvarfile  = "uvar_" // label

           open(10,file=uvarfile,form='formatted')
  102 format(1p6e20.12)

c          write statistics

        close(10)

         endif

        endif

9001    format(5e15.7)



  200 continue

1000    format(1p37e17.7)

        if(iprint.eq.0)then

           time=nstep*dt
           call output(x,xl,u,npts,n,time,ndim)

        endif



      stop
      end

      subroutine output(x,xl,u,jmax,n,time,ndim)
      implicit none
      integer ndim,n
      double precision time
      double precision u(0:ndim+1)
      double precision x(1:jmax)
      double precision xl(1:jmax+1)

      integer jmax,j

      character*12 ufile
      character*10 step


      write(step,'(i10.10)') n
      ufile = "u_" // step

      open(20,file=ufile,form='formatted')

      write(6,100)n,time
  100 format(//" step,time = ",i10,1e15.7/)

      do j=1,jmax

      write(6,101)x(j),u(j)
      write(20,102)x(j),u(j)

      enddo

      close(20)

  101 format(1p8e15.7)
  102 format(1p6e20.12)
      return
      end
