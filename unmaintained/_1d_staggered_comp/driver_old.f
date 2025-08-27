program burger
c
c************************************************************
c
c     this program implements the muscl scheme for the buckley
c
c     by setting the flag igud to 1 the program sets
c     all slope to zero and implements godunov's
c     method
      implicit none

      integer ndim,ncomp
      parameter (ndim = 2006)
      parameter (ncomp=1)
c
       double precision fluxr(-5:ndim)   ,fluxrm(-5:ndim)  ,
     1 fluxe(-5:ndim)   ,fluxC(-5:ndim)

       double precision rmom(-5:ndim,2),rtemp(-5:ndim,2),
     1        rspec(-5:ndim,2)



      double precision con(-5:ndim,4)
      double precision conold(-5:ndim,4)
      double precision conpre(-5:ndim,4)
      double precision conpre2(-5:ndim,4)
      double precision x(1:ndim),xl(1:ndim+1)


      double precision vmean(1:ndim,4),vcor(1:ndim,4)
      double precision vmeant(1:ndim,4),vcort(1:ndim,4)
c       vmsq(1:ndim,4,4)
c     double precision vmeant(1:ndim,4),vcort(1:ndim,4,4),
c    1            vmsqt(1:ndim,4,4)

      double precision ustar,Tstar,ustarbar,Tstarbar

      double precision dt,dx,across,T0,rhoamb,camb,vcell,xlen,pamb
      double precision weight,dorand,gamma,factor,spdmx,time
      double precision toteng,totmass,totmom,totspe

      double precision conc,cvmix,eint,keng,rho,T

      character*15 rhoustarfile,rhotstarfile,ttstarfile
      character*10 label


      integer seed,ntherm,nstep,i,j,k,l,iper,iprint,istat,jstart
      integer icor,n,nout,npts,nstat,nsamples

      double precision  mass0,mass1,d0,d1,cv0,cv1,rgas0,rgas1,kboltz,
     1   etaC0,etaC1,kappaC0, kappaC1,Tleft,Tright,pi

      common /params/ mass0,mass1,d0,d1,cv0,cv1,rgas0,rgas1,kboltz,
     1   etaC0,etaC1,kappaC0, kappaC1,Tleft,Tright,pi

      namelist /input_param/ mass0,mass1,d0,d1,Tleft,Tright,npts,xlen,
     1      iper,dorand,seed,T0,rhoamb,camb,across,ntherm,nstep,dt,nout,
     2      icor,nstat

c
       mass1 = 6.63d-23
       mass0 = mass1*3.d0
       d0=3.66d-8
       d1=3.66d-8
       seed = 632489
       ustarbar=0.d0
       Tstarbar=0.d0
       nstat = -1

       read(5,input_param)
       write(6,input_param)

       kboltz = 1.3806d-16
       rgas1 = kboltz/mass1
       rgas0 = kboltz/mass0

       pi = 4.d0*atan2(1.d0,1.d0)

       cv0 = 1.5*rgas0
       cv1 = 1.5*rgas1

       pamb = .5*rhoamb*T0*(rgas0+rgas1)

       etaC0 = 5.d0/(16.d0*d0**2)*sqrt(mass0*kboltz/pi)
       kappaC0 = 15.d0*kboltz/(4.d0*mass0)*etaC0
       etaC1 = 5.d0/(16.d0*d1**2)*sqrt(mass1*kboltz/pi)
       kappaC1 = 15.d0*kboltz/(4.d0*mass1)*etaC1

       time =0.0
       dx = xlen/dfloat(npts)

       do j=1,npts
         x(j) = (dfloat(j)-0.5d0)*dx
         xl(j) = dfloat(j-1)*dx
       enddo
       xl(npts+1) = xlen

       vcell = across*dx
       write(6,*)"vcell = ",vcell
       gamma = 5.d0/3.d0
       spdmx = sqrt(gamma*pamb/rhoamb)

c   set initial conditions




c  estimate diffusive and cfl time steps


c     eps = max(4.d0*eta0/(3.d0*rhoamb),akappa0/(rhoamb*cv))
c     dthyp = dx/spdmx
c     dtvisc = 0.5d0*dx**2/eps

       call blinitrand(seed)


c   set up stuff to gather statistics
       do k=1,4

         do j=1,npts

           vmean(j,k) = 0.d0
           vcor(j,k) = 0.d0

         enddo
       enddo

       nstep = nstep + ntherm
       istat = 0


       do i=1,npts
         con(i,1) = rhoamb
         con(i,2) = 0.d0
         con(i,3) = rhoamb*T0*(cv1*camb+cv0*(1.d0-camb))
         con(i,4) = rhoamb*camb
       enddo
       con(npts+1,2) = 0.d0

       call output(x,xl,con,npts,0,0.d0,ndim)

       do 200 n=1,nstep

c     if(n.gt.ntherm)dorand = 0.d0

         if(iper.eq.0)then
           jstart = 2
         else
           jstart = 1
         endif

         time = dfloat(n)*dt

         do j=0,npts+2
           do l=1,2

             call blinvnormdist(rmom(j,l))
             call blinvnormdist(rtemp(j,l))
             call blinvnormdist(rspec(j,l))

           enddo
         enddo


         weight =(2.d0*sqrt(2.d0)+sqrt(3.d0))/5.d0
         call flux(fluxr,fluxrm,
     1          fluxe,fluxC,con,npts,
     1          rmom,rtemp,rspec,weight,
     1        dx,dt,dorand,vcell,ndim,iper)



         do j=1,npts

           conpre(j,1) = con(j,1)-dt*(fluxr(j+1)-fluxr(j))/dx
           conpre(j,3) = con(j,3)-dt*(fluxe(j+1)-fluxe(j))/dx
           conpre(j,4) = con(j,4)-dt*(fluxC(j+1)-fluxC(j))/dx

         enddo
         do j=jstart,npts
           conpre(j,2) = con(j,2)-dt*(fluxrm(j)-fluxrm(j-1))/dx
         enddo

         if(iper.eq.0)then
           conpre(npts+1,2) = 0.d0
           conpre(1,2) = 0.d0
         else
           conpre(npts+1,2) = conpre(1,2)
         endif

         weight = (-4.d0*sqrt(2.d0)+3.d0*sqrt(3.d0))/5.d0
         call flux(fluxr,fluxrm,
     1          fluxe,fluxC,conpre,npts,
     1          rmom,rtemp,rspec,weight,
     1        dx,dt,dorand,vcell,ndim,iper)

         do j=1,npts

           conpre2(j,1) = 0.25*(3.d0*con(j,1)+conpre(j,1)
     1          -dt*(fluxr(j+1)-fluxr(j))/dx)
           conpre2(j,3) = 0.25*(3.d0*con(j,3)+conpre(j,3)
     1         -dt*(fluxe(j+1)-fluxe(j))/dx)
           conpre2(j,4) = 0.25*(3.d0*con(j,4)+conpre(j,4)
     1         -dt*(fluxC(j+1)-fluxC(j))/dx)
         enddo
         do j=jstart,npts
           conpre2(j,2) = 0.25*(3.d0*con(j,2)+conpre(j,2)
     1             -dt*(fluxrm(j)-fluxrm(j-1))/dx)
         enddo

         if(iper.eq.0)then
           conpre2(npts+1,2) = 0.d0
           conpre2(1,2) = 0.d0
         else
           conpre2(npts+1,2) = conpre2(1,2)
         endif

         weight = (sqrt(2.d0)-2.d0*sqrt(3.d0))/10.d0
         call flux(fluxr,fluxrm,
     1          fluxe,fluxC,conpre2,npts,
     1          rmom,rtemp,rspec,weight,
     1        dx,dt,dorand,vcell,ndim,iper)

         do j=1,npts

           con(j,1) = 1.d0/3.d0*(con(j,1)+2.d0*conpre2(j,1)
     1          -2.d0*dt*(fluxr(j+1)-fluxr(j))/dx)
           con(j,3) = 1.d0/3.d0*(con(j,3)+2.d0*conpre2(j,3)
     1         -2.d0*dt*(fluxe(j+1)-fluxe(j))/dx)
           con(j,4) = 1.d0/3.d0*(con(j,4)+2.d0*conpre2(j,4)
     1         -2.d0*dt*(fluxC(j+1)-fluxC(j))/dx)
         enddo
         do j=jstart,npts
           con(j,2) = 1.d0/3.d0*(con(j,2)+2.d0*conpre2(j,2)
     1             -2.d0*dt*(fluxrm(j)-fluxrm(j-1))/dx)
         enddo

         if(iper.eq.0)then
           con(npts+1,2) = 0.d0
           con(1,2) = 0.d0
         else
           con(npts+1,2) = con(1,2)
         endif

         if(mod(n,nout).eq.0)then

           totmass = 0.d0
           totmom = 0.d0
           toteng = 0.d0
           totspe = 0.d0

           do j=1,npts

             totmass = totmass+con(j,1)
             totmom = totmom+con(j,2)
             toteng = toteng+con(j,3)
             totspe = totspe+con(j,4)

           enddo

           write(6,*)n, "time = ",time, " totmass = ",totmass," totmom = ",
     1        totmom
           write(6,*)n,"time = ",time," toteng = ",toteng," totspe = ",
     1       totspe
         endif

         if(n.gt.ntherm .and. mod(n,nstat).eq.0 .and. nstat .gt.0)then

           istat = istat+1
           rho = con(icor,1)
           conc = con(icor,4)/con(icor,1)
           keng = 0.125d0*(con(icor,2)+con(icor+1,2))**2/con(icor,1)
           eint = (con(icor,3)-keng)/con(icor,1)
           cvmix = cv1*conc+cv0*(1.d0-conc)
           Tstar = eint/cvmix
           ustar = 0.5d0*(con(icor,2)+con(icor+1,2))/rho

           Tstarbar = Tstarbar + Tstar
           ustarbar = ustarbar + ustar


           do j=1,npts

             rho = con(j,1)
             conc = con(j,4)/con(j,1)
             keng = 0.125d0*(con(j,2)+con(j+1,2))**2/con(j,1)
             eint = (con(j,3)-keng)/con(j,1)
             cvmix = cv1*conc+cv0*(1.d0-conc)
             T = eint/cvmix

             vmean(j,1) = vmean(j,1)+rho
             vmean(j,2) = vmean(j,2)+T
             vcor(j,1) = vcor(j,1)+rho*ustar
             vcor(j,2) = vcor(j,2)+rho*Tstar
             vcor(j,3) = vcor(j,3)+T*Tstar

           enddo

         endif

         iprint = 0


         if(mod(n,nout).eq.0)then

           call output(x,xl,con,npts,n,time,ndim)
           iprint = 1

           if(n.gt.ntherm .and. istat .gt.0)then
             factor = dfloat(istat)
             do k=1,2
               do j=1,npts
                 vmeant(j,k) = vmean(j,k)/factor
               enddo
             enddo
             ustar = ustarbar/factor
             Tstar = Tstarbar/factor
             write(label,'(i10.10)') n
             rhoustarfile = "r_us_" // label
             rhotstarfile = "r_Ts_" // label
             ttstarfile = "T_Ts_" // label

             open(10,file=rhoustarfile,form='formatted')
             open(11,file=rhotstarfile,form='formatted')
             open(12,file=ttstarfile,form='formatted')
             do j=1,npts
               vcort(j,1) = vcor(j,1)/factor
     1            - vmeant(j,1)*ustar
               vcort(j,2) = vcor(j,2)/factor
     1            - vmeant(j,1)*Tstar
               vcort(j,3) = vcor(j,3)/factor
     1            - vmeant(j,2)*Tstar
               write(10,102) x(j),vcort(j,1)
               write(11,102) x(j),vcort(j,2)
               write(12,102) x(j),vcort(j,3)
  102          format(1p6e20.12)
             enddo

             close(10)
             close(11)
             close(12)

           endif

         endif

       9001    format(5e15.7)

  200 continue

       write(6,*)"end of loop"


!       factor = dfloat(istat)
!       do k=1,4
!       do j=1,npts
!           vmean(j,k) = vmean(j,k)/factor
!       enddo
!       enddo
!       do l=1,4
!       do k=1,4
!       do j=1,npts
!           vmsq(j,k,l) = vmsq(j,k,l)/factor
!    1            - vmean(j,k)*vmean(j,l)
!           vcor(j,k,l) = vcor(j,k,l)/factor
!    1            - vmean(j,k)*vmean(icor,l)
!       enddo
!       enddo
!       enddo
!       do j=1,npts
!           write(9,1000)x(j),vmean(j,1),vmean(j,2),vmean(j,3)
!    1        ,vmean(j,4)
!           write(10,1000)x(j),((vmsq(j,k,l),l=1,4),k=1,4)
!           write(11,1000)x(j),((vcor(j,k,l),l=1,4),k=1,4)
!       enddo
 1000    format(1p37e17.7)

       if(iprint.eq.0)then

         time=nstep*dt
         call output(x,xl,con,npts,n,time,ndim)

       endif



       stop
       end

       subroutine output(x,xl,con,jmax,n,time,ndim)
       integer ndim,n
       double precision time
       double precision con (-5:ndim,4)
       double precision x(1:jmax)
       double precision xl(1:jmax+1)
       double precision rho,comp,keng,eing,T,p


       double precision  mass0,mass1,d0,d1,cv0,cv1,rgas0,rgas1,kboltz,
     1   etaC0,etaC1,kappaC0, kappaC1,Tleft,Tright,pi

       common /params/ mass0,mass1,d0,d1,cv0,cv1,rgas0,rgas1,kboltz,
     1   etaC0,etaC1,kappaC0, kappaC1,Tleft,Tright,pi

       character*14 rhofile
       character*15 tempfile
       character*15 concfile
       character*15 presfile
       character*15 xmomfile
       character*10 step

       write(step,'(i10.10)') n
       rhofile = "rho_" // step
       tempfile = "temp_" // step
       concfile = "conc_" // step
       presfile = "pres_" // step
       xmomfile = "xmom_" // step

       open(20,file=rhofile,form='formatted')
       open(21,file=tempfile,form='formatted')
       open(22,file=concfile,form='formatted')
       open(23,file=presfile,form='formatted')
       open(24,file=xmomfile,form='formatted')

       write(6,100)n,time
  100 format(//" step,time = ",i10,1e15.7/)

       do j=1,jmax

         rho = con(j,1)
         conc = con(j,4)/con(j,1)
         keng = 0.125d0*(con(j,2)+con(j+1,2))**2/con(j,1)
         eint = (con(j,3)-keng)/con(j,1)
         cvmix = cv1*conc+cv0*(1.d0-conc)
         T = eint/cvmix
         p = rho*T*(rgas1*conc+rgas0*(1.d0-conc))

         write(6,101)x(j),rho,p,T,conc,con(j,2)
         write(20,102)x(j),rho
         write(21,102)x(j),T
         write(22,102)x(j),conc
         write(23,102)x(j),p
       enddo

       close(20)
       close(21)
       close(22)
       close(23)

       write(6,101)xl(jmax+1),con(jmax+1,2)
       write(6,*)" "
       do j=1,jmax+1
c        write(6,101)xl(j),con(j,2)
         write(24,102)xl(j),con(j,2)
       enddo
       close(24)

  101 format(1p8e15.7)
  102 format(1p6e20.12)
       return
       end