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

    integer ndim,pdim,edim,tdim
    parameter (ndim = 200)
    parameter (pdim = 20000)
    parameter (edim = 10000)
    parameter (tdim = 20)

    double precision part(pdim),part_new(pdim)
c
    double precision flux(1:ndim+1,2)
    double precision ranflux(1:ndim+1,2)

    double precision u(0:ndim+1,2),unew(0:ndim+1,2)
    double precision x(1:ndim),xl(1:ndim+1)
    double precision mean(ndim),var(ndim),skew(ndim),kur(ndim)
    double precision sqr(ndim),cub(ndim),for(ndim)
    double precision mean_acc(ndim),sqr_acc(ndim),cub_acc(ndim),
   1                 for_acc(ndim)
    double precision skew2(ndim),kur2(ndim)
    double precision pdf(-1:201),binlo,binhi,dbin
    double precision av_mean, av_var, av_skew, av_kur, meanloc
    double precision var_mean, var_var, var_skew, var_kur
    integer nbins,ipdf,ibin,totpts


    character*17 uvarfile
    character*16 label
    character*3  tag

    character*28 meanfile
    character*27 varfile
    character*28 skewfile
    character*28 kurfile
    character*27 pdffile



    integer seed,ntherm,nstep,i,j,k,l,iper,iprint,istat,jstart
    integer icor,n,nout,npts,nstat,nsamples,ires
    integer jmidl,jmidr,jpartl,jpartr,is_hybrid
    integer num_part, npart, npghost, num_part_new, is_pure_part
    integer ensemble,ens,ensout,tout,is_gaussian,ncoef, irestart,icnt

    double precision xlen,dt,dx,factor,time,totmass,dorand,weight,cfl
    double precision minmass,midfact

    double precision crossl,crossr,fluxregl, fluxregr

    double precision  uleft, uright, pi, uinit

    common /params/  uleft, uright, pi, uinit

    namelist /input_param/ npts,xlen,iper,dorand,seed,ntherm,nstep,
   1   dt,nout,ires,nstat,icor,ires,uleft,uright, uinit, cfl,
   1   jmidl,jmidr,midfact,jpartl,jpartr,is_hybrid,ensemble,
   1   is_gaussian,irestart, ipdf, binlo, nbins, dbin, ensout


    xlen = 1.d0
    iper = 1
    irestart = 0
    ntherm = 0
    dorand = 1.d0
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
    is_gaussian = 0
    ncoef = 1
    icnt = 0
    nbins = 150
    ipdf = 0
    ensout = 1


    read(5,input_param)
    write(6,input_param)

    if(ipdf.eq. 1)then
        binhi = binlo+dbin*dfloat(nbins)
        do ibin = -1,nbins+1
            pdf(ibin) = 0.d0
        enddo
    endif
    tag = "fv_"
    if(is_hybrid.eq.1)then
        tag = "hy_"
        if(jpartl.eq.1 .and. jpartr.eq.npts)then
            is_pure_part = 1
            tag = "pp_"
            write(6,*)" pure particle simulation"
        elseif(jpartl.le.npghost .or. jpartr.gt.npts-npghost)then
            write(6,*)" patch close to boundary"
            if(iper.eq.1)then
                write(6,*)" can't be periodic with particles",
   1               " patch touching boundary"
                stop
            endif
        endif
    endif

    if(is_gaussian .eq. 1) then
        ncoef = 2
        tag = "ga_"
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

    do j=1,npts
        mean_acc(j) = 0.d0
        sqr_acc(j) = 0.d0
        cub_acc(j) = 0.d0
        for_acc(j) = 0.d0
    enddo

    istat = 0

    write(6,*)"delreg,delnreg",4.d0*dx**2, 4.d0*dx**2*uinit


    call blinitrand(seed)

    nstep = nstep + ntherm
c   set initial conditions
    do ens = 1,ensemble

        tout = 0

        do k=1,ncoef
            do j=1,npts
                u(j,k) = uinit
                if(j.ge.jmidl .and. j.le.jmidr)then
                    u(j,k) = u(j,k)*midfact
                endif
            enddo
        enddo

        totmass = 0.d0
        do j=1,npts
            totmass = totmass+u(j,1)
        enddo
        totmass = totmass * dx

c          write(6,*)"init totmass before = ",totmass


        if(is_hybrid.eq.1)then
            call init_part(u,part,xl,jpartl,jpartr,num_part,dx,ndim,pdim)
            call avg_down(u,part,xl,jpartl,jpartr,num_part,dx,
   1         ndim,pdim)
        endif

        totmass = 0.d0
        do j=1,npts
            totmass = totmass+u(j,1)
        enddo
        totmass = totmass * dx

c          write(6,*)"init totmass after = ",totmass



c   set up stuff to gather statistics



        call output(x,xl,u,npts,0,0.d0,ens,ndim)

        do 200 n=1,nstep

            time = dfloat(n)*dt

c        if(jpartl.gt.1 .and.jpartr .lt.npts)then
c        if(is_hybrid.eq.0 .or .is_pure_part.eq.0)then
            if(is_pure_part.eq.0)then

                do j=1,npts+1
c           call blinvnormdist(ranflux(j,l))
                    call blinvnormdist(ranflux(j,1))
                enddo

                weight =0.d0
c        weight =(2.d0*sqrt(2.d0)+sqrt(3.d0))/5.d0
                call compflux(flux,ranflux,u,npts,weight,
   1            dx,dt,dorand,ndim,iper,ires,ncoef)

                do k=1,ncoef
                    do j=1,npts
                        unew(j,k) = u(j,k)+dt*(flux(j+1,k)-flux(j,k))/dx
                    enddo
                enddo

                fluxregl = flux(jpartl,1)   * dt
                fluxregr = flux(jpartr+1,1) * dt

                do k=1,ncoef
                    do j=1,npts
                        u(j,k) = unew(j,k)
                    enddo
                enddo

            endif

            if(is_hybrid .eq. 1)then
                call fill_part(u,part,xl,jpartl,jpartr,num_part,dx,
   1         npts,ndim,pdim,npghost)

                call adv_part(part,part_new,xl,jpartl,jpartr,crossl,crossr,
   1        num_part,dx,dt,npts,ndim,pdim)

                call resort(part_new,part,xl,jpartl,jpartr,
   1        num_part,num_part_new,dx,dt,npts,iper,ndim,pdim)

                num_part = num_part_new

c           write(6,*)" before avg ",u(1,1)

                call  avg_down(u,part,xl,jpartl,jpartr,num_part,dx,
   1         ndim,pdim)

c           write(6,*)" after avg ",u(1,1)

c           totmass = 0.d0
c           do j=jpartl,jpartr
c               totmass = totmass + u(j)
c           enddo
c           write(6,*)" mass in particle region ",totmass*dx

                call reflux(u,fluxregl,fluxregr,crossl,crossr,
   1          jpartl,jpartr,npts,dx, ndim)

c           write(6,*)" after rflux ",u(1,1)

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

            if(n.gt.ntherm .and. mod(n,nout).eq.0)then
                totmass = 0.d0
                minmass = 1.d20

                do j=1,npts
                    totmass = totmass+u(j,1)
                    minmass = min(minmass,u(j,1))
                enddo

                totmass = totmass * dx

c              write(6,*)n, "time = ",time, " totmass = ",totmass,
c    1            " minmass = ",minmass
            endif

            if(n.gt.ntherm .and. mod(n,nstat).eq.0 .and. nstat .gt.0)then
c     gather statistics
                istat = istat+1

                do j = 1,npts
                    mean_acc(j) = mean_acc(j) + u(j,1)
                    sqr_acc(j) = sqr_acc(j) + u(j,1)**2
                    cub_acc(j) = cub_acc(j) + u(j,1)**3
                    for_acc(j) = for_acc(j) + u(j,1)**4
                enddo

                if(ipdf .eq. 1)then
                    do j=1,npts
                        if(u(j,1) .lt. binlo)then
                            ibin = -1
                        elseif (u(j,1).gt.binhi) then
                            ibin = nbins
                        else
                            ibin = (u(j,1)-binlo)/dbin
                        endif
                        pdf(ibin) = pdf(ibin)+1
                        totpts = totpts+1
                    enddo
                endif

                iprint = 0
            endif


            if(n.gt.ntherm .and. mod(n,nout).eq.0 .and.
   1         mod(ens,ensout).eq.0)then

                call output(x,xl,u,npts,n,time,ens,ndim)

                do j=1,npts
                    mean(j) = mean_acc(j) / dfloat(istat)
                    var(j) = sqr_acc(j) / dfloat(istat)- mean(j)**2
                    sqr(j) = sqr_acc(j) / dfloat(istat)
                    cub(j) = cub_acc(j) / dfloat(istat)
                    for(j) = for_acc(j) / dfloat(istat)
                enddo

                av_mean = 0.d0
                av_var  = 0.d0
                av_skew = 0.d0
                av_kur  = 0.d0
                var_mean = 0.d0
                var_var  = 0.d0
                var_skew = 0.d0
                var_kur  = 0.d0

                do j=1,npts
                    skew(j) = cub(j) - 3.d0*sqr(j)*mean(j) + 2.d0*mean(j)**3
                    kur(j) = for(j) - 4.d0*cub(j)*mean(j) +6.d0*sqr(j)*mean(j)**2
   1                - 3.d0*mean(j)**4

                    if(var(j).gt.0.d0)then
c                       skew(j) = skew(j)/(dfloat(ensemble)*var(j)**1.5)
c                       kur(j) = kur(j)/(dfloat(ensemble)*var(j)**2)
                        skew(j) = skew(j)/(var(j)**1.5)
                        kur(j) = kur(j)/(var(j)**2)
c                    else
c                       skew(j) = skew(j)/(dfloat(ensemble))
c                       kur(j) = kur(j)/(dfloat(ensemble))
                    endif

                    av_mean = av_mean + mean(j)
                    av_var  = av_var  + var(j)
                    av_skew = av_skew + skew(j)
                    av_kur  = av_kur  + kur(j)
                    var_mean = var_mean + mean(j)**2
                    var_var  = var_var  + var(j)**2
                    var_skew = var_skew + skew(j)**2
                    var_kur  = var_kur  + kur(j)**2
                enddo

                write(6,*)"istats",istat,npts

                av_mean = av_mean / dfloat(npts)
                av_var  = av_var  / dfloat(npts)
                av_skew = av_skew / dfloat(npts)
                av_kur  = av_kur  / dfloat(npts)
                var_mean = var_mean / dfloat(npts) - av_mean**2
                var_var  = var_var  / dfloat(npts) - av_var**2
                var_skew = var_skew / dfloat(npts) - av_skew**2
                var_kur  = var_kur  / dfloat(npts) - av_kur**2

                write(6,*)"aver stats ", av_mean, av_var, av_skew,av_kur
                write(6,*)"error bars ",sqrt(var_mean/dfloat(npts)),
   1                            sqrt(var_var/dfloat(npts)),
   1                            sqrt(var_skew/dfloat(npts)),
   1                            sqrt(var_kur/dfloat(npts))

                write(label,'(i6.6,"_",i9.9)')ens,n

                meanfile = "mean_" // tag // label // ".dat"
                varfile  =  "var_" // tag // label // ".dat"
                skewfile = "skew_" // tag // label // ".dat"
                kurfile = "kurt_" // tag // label // ".dat"

                open(20,file=meanfile,form='formatted')
                open(21,file=varfile,form='formatted')
                open(22,file=skewfile,form='formatted')
                open(23,file=kurfile,form='formatted')

                do j=1,npts
                    write(20,102)x(j),mean(j)
                    write(21,102)x(j),var(j)
                    write(22,102)x(j),skew(j)
                    write(23,102)x(j),kur(j)
                enddo

                close(20)
                close(21)
                close(22)
                close(23)

                if(ipdf .eq.1)then
                    pdffile = "pdf_" // tag // label // ".dat"
                    open(24,file=pdffile,form='formatted')
                    do ibin = -1,nbins
                        write(24,102)binlo+(dfloat(ibin)+.5)*dbin,
   1                   pdf(ibin)/dfloat(totpts)
                    enddo
                    close(24)
                endif

            endif

  102 format(1p6e20.12)


9001    format(5e15.7)



  200 continue

    write(6,*)"completed ",ens

1000    format(1p37e17.7)

    if(iprint.eq.0 .and. mod(ens,ensout).eq.0)then
        time=nstep*dt
        call output(x,xl,u,npts,n,time,ens,ndim)
    endif

    enddo

    stop
end

subroutine output(x,xl,u,jmax,n,time,ens,ndim)
    implicit none
    integer ndim,n,ens
    double precision time
    double precision u(0:ndim+1)
    double precision x(1:jmax)
    double precision xl(1:jmax+1)

    integer jmax,j

    character*15 ufile
    character*5 step
    character*7  ensemble


    write(ensemble,'(i7.7)') ens
    write(step,'(i5.5)') n
    ufile = "u_" // ensemble // "_" // step

c   open(20,file=ufile,form='formatted')

    write(6,100)n,time,ens
100 format(//" step,time = ",i10,1e15.7/)

    do j=1,jmax
        write(6,101)x(j),u(j)
c   write(20,102)x(j),u(j)
    enddo

c   close(20)

101 format(1p8e15.7)
102 format(1p6e20.12)
    return
end