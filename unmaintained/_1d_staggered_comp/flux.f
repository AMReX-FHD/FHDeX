      subroutine flux(fluxr,fluxrm,fluxe,fluxC,con,npts,
     1       rmom, rtemp, rspec, weight,
     1       dx,dt,dorand,vcell,ndim,iper)

      implicit none
      integer nloc,ndim,iper,npts

      parameter (nloc = 1006)

       double precision con(-5:ndim,4),
     1     fluxr(-5:ndim)   ,fluxrm(-5:ndim)  ,fluxe(-5:ndim),
     1       fluxC(-5:ndim)
       double precision rmom(-5:ndim,2),rtemp(-5:ndim,2),
     1        rspec(-5:ndim,2)
      double precision dx,dt,dorand,vcell,weight


      double precision x(-5:nloc)
      double precision viscmomf(-5:nloc)
      double precision viscenef(-5:nloc)
      double precision viscspec(-5:nloc)
      double precision ranmomf(-5:nloc)
      double precision ranenef(-5:nloc)
      double precision ranspec(-5:nloc)
      double precision tau(-5:nloc)
      double precision q(-5:nloc)
      double precision T(-5:nloc)
      double precision p(-5:nloc)
      double precision rho(-5:nloc)
      double precision u(-5:nloc)
      double precision eint(-5:nloc)
      double precision comp(-5:nloc)

      double precision eta(-5:nloc)
      double precision prim(-5:nloc,6)
      double precision kappa(-5:nloc)
      double precision dfic(-5:nloc)
      double precision kT(-5:nloc)
      double precision kP(-5:nloc)
      double precision kTdmudc(-5:nloc)
      double precision kap0,kap1,kapx,xkap,ykap,zkap
      double precision cu0, cu1,cuy,cuz, conc,temp
      double precision diamx,mx,omc,ctemp
      double precision  s0, s1, dmudcm1,Jx

      double precision prom,ratm,difm,summ,omcsq,comc,csq,cvmix

      double precision kxp
      double precision mxp,kTxp,kpxp
      real*8 cont1,cont2,cont3,cont4,rhot,ut,pt,compt,eta1,eta0
      real*8 cont5,cont6,vt,wt,  keng, concfl

      real*8 etax,xeta,yeta,zeta,dxp,soret,vloc
      real*8 varflv,ctempp,concp,omcp,gradc,gradt

      integer i,j
      integer jend,jstart

      double precision  mass0,mass1,d0,d1,cv0,cv1,rgas0,rgas1,kboltz,
     1   etaC0,etaC1,kappaC0, kappaC1,Tleft,Tright,pi

      common /params/ mass0,mass1,d0,d1,cv0,cv1,rgas0,rgas1,kboltz,
     1   etaC0,etaC1,kappaC0, kappaC1,Tleft,Tright,pi




       ratm = mass0/mass1

       do j=1,npts

          rho(j) = con(j,1)
          comp(j) = con(j,4)/con(j,1)
          keng = 0.125d0*(con(j,2)+con(j+1,2))**2/con(j,1)
          eint(j) = (con(j,3)-keng)/con(j,1)
          concfl = max(0.d0,min(1.d0,comp(j)))
          cvmix = cv1*concfl+cv0*(1.d0-concfl)
          T(j) = eint(j)/cvmix
          p(j) = rho(j)*T(j)*(rgas1*comp(j)+rgas0*(1.d0-comp(j)))

        enddo

        do j=2,npts

           u(j) = con(j,2)*2.d0/(rho(j-1)+rho(j))

        enddo

        if(iper.eq.1)then

           rho(0) = rho(npts)
           u(0) = u(npts)
           u(1) = con(1,2)*2.d0/(rho(1)+rho(0))
           u(npts+1) = u(1)
           u(npts+2) = u(2)

           comp(0) = comp(npts)
           comp(npts+1) = comp(1)
           eint(0) = eint(npts)
           eint(npts+1) = eint(1)
           p(0) = p(npts)
           p(npts+1) = p(1)
           T(0) = T(npts)
           T(npts+1) = T(1)

         else

             u(1) = 0.d0
             u(npts+1) = 0.d0
             T(0) = Tleft
             T(npts+1) = Tright

             comp(0) = comp(1)
             comp(npts+1) = comp(npts)
             rho(0) = rho(1)
             rho(npts+1) = rho(npts)
             p(0) = p(1)
             p(npts+1) = p(npts)




         endif



1000   format(4e15.7)


      prom = mass0*mass1
      ratm = mass0/mass1
      difm = (mass0-mass1)**2/prom
      summ = (mass0+mass1)**2/prom

      mx = 2.d0*mass0*mass1/(mass0+mass1)
      diamx = 0.5d0*(d0+d1)

         do j=0,npts+1


        ctemp = comp(j)
        ctemp = max(0.d0,min(1.d0,ctemp))
        conc = ratm*ctemp/(1.d0-(1.d0-ratm)*ctemp)

        omc = 1.d0-conc
        temp = T(j)
        eta1 = etaC1*sqrt(temp)
        eta0 = etaC0*sqrt(temp)

        etax = 5.d0*sqrt(mx*kboltz*temp/pi)/(16.d0*diamx**2)


        csq = conc*conc
        omcsq = omc*omc
        comc = conc*omc

        xeta = omcsq/eta0+2.d0*comc/etax+csq/eta1
        yeta = .6d0*( omcsq*ratm/eta0+
     1       .5d0*comc*etax*summ/(eta0*eta1) + csq/(eta1*ratm))
        zeta = .6*(omcsq*ratm+2.d0*comc*(0.25d0*summ*(etax/eta0
     2      +etax/eta1)-1.d0)+csq/ratm)

        eta(j) = (1.d0+zeta)/(xeta+yeta)
        dfic(j) = 6.d0*etax/(5.d0*rho(j))


        kap0 = 15*kboltz*eta0/(4.d0*mass0)
        kap1 = 15*kboltz*eta1/(4.d0*mass1)
        kapx = 15*kboltz*etax/(4.d0*mx)

        cu0 = 4.d0/15.d0-17.d0*ratm/60.d0+ 0.5d0*  difm
        cu1 = 4.d0/15.d0-17.d0/(60.d0*ratm)+ 0.5d0*  difm
        cuy = summ*kapx**2/
     1    (15.d0*kap0*kap1) -17.d0/60.d0 +
     2   13.d0/32.d0*difm
        cuz = (summ*(kapx/kap0+kapx/kap1)-4.d0)/15.d0-17.d0/60.d0

        xkap = omcsq/kap0+2.d0*comc/kapx+csq/kap1
        ykap = omcsq*cu0/kap0+2*comc*cuy/kapx+csq*cu1/kap1
        zkap = omcsq*cu0+2*comc*cuz+csq*cu1


        kappa(j) = (1.d0+zkap)/(xkap+ykap)

        s0 = 0.5d0*(mass0+mass1)*kapx/(mass0*kap0) -
     1        15d0*(mass1-mass0)/(8.d0*mass0)-1.d0
        s1 = 0.5d0*(mass0+mass1)*kapx/(mass1*kap1) +
     1        15d0*(mass1-mass0)/(8.d0*mass1)-1.d0

        dmudcm1 =  ratm/(ratm + (1.d0-ratm)*conc)**2*
     &    mass0*mass1*conc*omc/(mass0*omc+mass1*conc)


        kT(j) = conc*omc*(s1*conc-s0*omc)/(6.d0*kapx*(xkap+ykap))
        kP(j) = (mass0-mass1)*ctemp*(1.d0-ctemp)*
     &       (ctemp/mass1+(1.d0-ctemp)/mass0)
        kTdmudc(j) = kboltz*T(j)*
     &         (s1*conc-s0*omc)/(6.d0*kapx*(xkap+ykap))
     &    *(ratm + (1.d0-ratm)*conc)**2*(omc/mass1+conc/mass0)/
     &         ratm


          enddo


         if(iper.eq.1)then
            jstart = 1
         else
            jstart = 2
         endif

         do j=jstart,npts

            fluxr(j) = con(j,2)
            fluxe(j) = 0.5d0*u(j)*(con(j,3)+con(j-1,3) + p(j) + p(j-1))
            fluxC(j) = 0.5d0*u(j)*(con(j,4)+con(j-1,4))

         enddo

         if(iper.eq.0)then

              fluxr(1) = 0.d0
              fluxe(1) = 0.d0
              fluxC(1) = 0.d0
              fluxr(npts+1) = 0.d0
              fluxe(npts+1) = 0.d0
              fluxC(npts+1) = 0.d0

          else

               fluxr(npts+1) = fluxr(1)
               fluxe(npts+1) = fluxe(1)
               fluxC(npts+1) = fluxC(1)

          endif

          do j=1,npts

                fluxrm(j) = p(j) + 0.25d0*(u(j)+u(j+1))*
     1               (con(j,2)+con(j+1,2))

          enddo

          if(iper.eq.1)then

                 fluxrm(0) = fluxrm(npts)
                 fluxrm(npts+1) = fluxrm(1)

          endif


          do j=jstart,npts

           dxp = 0.5d0*(dfic(j-1)*rho(j-1)+dfic(j)*rho(j))
           kTxp = 0.5d0*(kT(j-1)/T(j-1)+kT(j)/T(j))
           kpxp = 0.5d0*(kp(j-1)/p(j-1)+kp(j)/p(j))


        Jx = dxp*(comp(j)-comp(j-1))/dx
     &      + dxp*(kTxp*(T(j)-T(j-1))/dx +
     &          kpxp*(p(j)-p(j-1))/dx)


        soret = 0.5d0*(kTdmudc(j-1)+kTdmudc(j))
        ctempp = comp(j-1)
        ctemp = comp(j)

        mxp = ctempp*(1.d0-ctempp)*(mass1*(1.d0-ctempp)
     &           + mass0*ctempp)
        mx = ctemp*(1.d0-ctemp)*(mass1*(1.d0-ctemp)
     &           + mass0*ctemp)




              q(j) = 0.5d0*(kappa(j)+kappa(j-1))*
     1               (T(j)-T(j-1))/dx

             viscenef(j) = q(j)
     &     +(soret+0.5d0* (cv1+rgas1-cv0-rgas0)*
     &         (T(j-1)+T(j)))*Jx

             viscspec(j) =  Jx


          enddo

          if(iper.eq.1)then
              viscenef(npts+1) = viscenef(1)
              viscspec(npts+1) = viscspec(1)
          else

             Jx = 0.d0
             q(1) = 2.d0*kappa(0)*(T(1)-T(0))/dx
            viscenef(1) = q(1)
            viscspec(1) =  Jx
            viscenef(npts+1) = 2.d0*kappa(npts+1)*(T(npts+1)-T(npts))/dx
            viscspec(npts+1) =  Jx
          endif


          do j=1,npts

             tau(j) = 4.d0*eta(j)/3.d0*
     1            (u(j+1)-u(j))/dx
             viscmomf(j) =tau(j)
           enddo

           if(iper.eq.1)then
              viscmomf(0) = viscmomf(npts)
              viscmomf(npts+1) = viscmomf(1)
           endif

      do j=jstart,npts
           viscenef(j) = viscenef(j)+ 0.5d0*u(j)*(tau(j)+tau(j-1))
      enddo

           if (iper.eq.1)then
              viscenef(npts+1) = viscenef(1)
           endif

      do j=1,npts

         varflv = 8.d0*kboltz/(3.d0*dt*vcell)*eta(j)*T(j)
         ranmomf(j) = sqrt(varflv)*(rmom(j,1)+weight*rmom(j,2))
     1             *dorand

      enddo

        if(iper.eq.1)then

             ranmomf(0) = ranmomf(npts)
             ranmomf(npts+1) = ranmomf(1)

        endif

      do j=1,npts+1

        gradc = rspec(j,1)+weight*rspec(j,2)
        gradt = rtemp(j,1)+weight*rtemp(j,2)

        if(j.eq.1 .and. iper.eq.0)then

        gradc = 0.d0
        gradt = gradt*sqrt(2.d0)

        kxp = 2.d0* kappa(j-1)*kboltz*T(j-1)**2

          soret = kTdmudc(j-1)


        ctempp = comp(j-1)
        ctempp = max(0.d0,min(1.d0,ctemp))
        concp = ratm*ctempp/(1.d0-(1.d0-ratm)*ctempp)
        omcp = 1.d0-concp

        mxp = ctempp*(1.d0-ctempp)*(mass1*(1.d0-ctempp)
     &           + mass0*ctempp)

        dxp = 2.d0*dfic(j-1)   *rho(j-1)*mxp

        elseif(j.eq.npts+1 .and. iper.eq.0)then
        gradc = 0.d0
        gradt = gradt*sqrt(2.d0)

        kxp = 2.d0*  kappa(j)*kboltz*T(j)**2

          soret = kTdmudc(j)


        ctemp = comp(j)
        ctemp = max(0.d0,min(1.d0,ctemp))
        conc = ratm*ctemp/(1.d0-(1.d0-ratm)*ctemp)
        omc = 1.d0-conc

        mx = ctemp*(1.d0-ctemp)*(mass1*(1.d0-ctemp)
     &           + mass0*ctemp)

        dxp = 2.d0* dfic(j)  *rho(j)*mx


        else

        kxp = kappa(j-1)*kboltz*T(j-1)**2 +
     1   kappa(j)*kboltz*T(j)**2

          soret = 0.5d0*(kTdmudc(j-1)+kTdmudc(j))


        ctemp = comp(j)
        ctemp = max(0.d0,min(1.d0,ctemp))
        conc = ratm*ctemp/(1.d0-(1.d0-ratm)*ctemp)
        ctempp = comp(j-1)
        ctempp = max(0.d0,min(1.d0,ctemp))
        concp = ratm*ctempp/(1.d0-(1.d0-ratm)*ctempp)
        omcp = 1.d0-concp
        omc = 1.d0-conc

        mxp = ctempp*(1.d0-ctempp)*(mass1*(1.d0-ctempp)
     &           + mass0*ctempp)
        mx = ctemp*(1.d0-ctemp)*(mass1*(1.d0-ctemp)
     &           + mass0*ctemp)

        dxp = dfic(j-1)   *rho(j-1)*mxp
     &   + dfic(j)  *rho(j)*mx

        endif


        gradc = dorand * gradc*sqrt(dxp/(vcell*dt))
        gradt = dorand * gradt*sqrt(kxp/(vcell*dt))

         ranenef(j) = gradt
     &     +(0.5*(cv1+rgas1-cv0-rgas0)*(T(j)+T(j+1)) + soret)
     &     *gradc
         ranspec(j) = gradc


      enddo

      if(iper.eq.1)then
         jstart = 1
      else
         jstart=2
      endif

      do j=jstart,npts


         vloc = 0.5d0*u(j)*(ranmomf(j)+ranmomf(j-1))
         ranenef(j) = ranenef(j)+vloc

      enddo

      if(iper.eq.1)then
         viscenef(npts+1) = viscenef(1)
         viscspec(npts+1) = viscspec(1)
         ranenef(npts+1) = ranenef(1)
         ranspec(npts+1) = ranspec(1)
      else
         ranspec(1) = 0.d0
         ranspec(npts+1) = 0.d0
         viscspec(1) = 0.d0
         viscspec(npts+1) = 0.d0

      endif

      do  j=1,npts+1

           fluxrm(j) = fluxrm(j)- viscmomf(j)
     1          - ranmomf(j)
           fluxe(j) = fluxe(j)- viscenef(j)
     1          - ranenef(j)
           fluxC(j) = fluxC(j)- viscspec(j)
     1          - ranspec(j)
      enddo
        if(iper.eq.1)then
           fluxrm(0) = fluxrm(npts)
        endif

      return
      end
