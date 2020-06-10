module multispec_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : n_cells, ngc, k_b, Runiv, nprimvars, nspecies, molmass, diameter, hcp, hcv
  use conv_module, only : get_molfrac

  implicit none

  private

  public :: cholesky_decomp, makecoef, visc_lin, lambda_lin, D_GIOVANGIGLI, thermalDiff

contains

  subroutine cholesky_decomp(a,np,sqda) bind(C,name="cholesky_decomp")  

    integer :: np

    real(amrex_real), intent(inout) :: a(np,np)
    real(amrex_real), intent(inout) :: sqda(np,np)

    real(amrex_real) :: p(np), sqda2(np,np), dij(np,np)
    real(amrex_real) :: dd(np,np)
    real(amrex_real) :: yy(np), mwmix 

    integer :: i, j, k, ii, jj
    real(amrex_real) :: sum1
    real(amrex_real) :: small_number
    real(amrex_real) :: sum

    integer :: idiag,ising

    small_number = 0.d0
    sum1 = 0.d0

    ! NOTE: For idiag=1, please refer to original LLNS code
    idiag = 0

    do i = 1, np
       ising = 0

       do j = i, np
          sum1 = a(i,j)

          do k = i-1, 1, -1
             sum1 = sum1 - a(i,k)*a(j,k)
          enddo

          if(i.eq.j) then
             if(sum1.le.small_number) then
                p(i) = 0.d0
                ising = 1
             else
                p(i) = sqrt(sum1)
             endif
          else
             if(ising.eq.0)then
                a(j,i) = sum1/p(i)
             else
                a(j,i) = 0.d0
             endif
          endif

       enddo
    enddo

    sqda = 0.0d0

    do i = 1, np
       do j = i-1, 1, -1
          sqdA(i,j) = a(i,j)
       enddo
       sqdA(i,i) = p(i)
    enddo

  end subroutine cholesky_decomp

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Compute Coefficients

  subroutine makecoef(lo,hi,prim,eta,zeta,kappa,chi,Dij) bind(C,name="makecoef")  

    integer         , intent(in   ) :: lo(3),hi(3)

    real(amrex_real), intent(in   ) :: prim(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)

    real(amrex_real), intent(inout) :: eta(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
    real(amrex_real), intent(inout) :: zeta(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
    real(amrex_real), intent(inout) :: kappa(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))

    real(amrex_real), intent(inout) :: chi(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nspecies)
    real(amrex_real), intent(inout) :: Dij(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nspecies,nspecies)
    !------------------------------------
    integer i, j, k, ll, kk, ns, d
    integer kbounds
    integer iwrk
    real(amrex_real) rwrk, sumYk
    real(amrex_real) :: Yk_fixed(1:nspecies)
    real(amrex_real) :: Xk_fixed(1:nspecies)

    real(amrex_real) etaf,kappaf,zetaf,DF(1:nspecies,1:nspecies),chif(1:nspecies),Yf(nspecies),Xf(nspecies)
    
    if(n_cells(3).eq.1)then
       kbounds = 0
    else
       kbounds = ngc(3)
    endif
    
    do k = lo(3)-kbounds,hi(3)+kbounds
       do j = lo(2)-ngc(2),hi(2)+ngc(2)
          do i = lo(1)-ngc(1),hi(1)+ngc(1)

             sumYk = 0.d0
             do ns = 1, nspecies
                Yk_fixed(ns) = max(0.d0,min(1.d0,prim(i,j,k,6+ns)))
                sumYk = sumYk + Yk_fixed(ns)

             enddo

             Yk_fixed(:) = Yk_fixed(:)/sumYk

             call get_molfrac(Yk_fixed, Xk_fixed)

             call ideal_mixture_transport(prim(i,j,k,1), prim(i,j,k,5), prim(i,j,k,6), &
                  Yk_fixed,Xk_fixed, &
                  eta(i,j,k), kappa(i,j,k), zeta(i,j,k), Dij(i,j,k,1:nspecies,1:nspecies),chi(i,j,k,1:nspecies))

             !   want this multiplied by rho for all times
             do kk = 1,nspecies
                do ll = 1,nspecies
                   Dij(i,j,k,kk,ll) = Dij(i,j,k,kk,ll)*prim(i,j,k,1)

                   if ( isnan(Dij(i,j,k,kk,ll)) ) then
                      print*, "Hack 1, (makecoef) ", i,j,k, Dij(i,j,k,:,:)
                      print*, prim(i,j,k,1:6), Yk_fixed,Xk_fixed, &
                           eta(i,j,k), kappa(i,j,k), zeta(i,j,k)
                      stop
                   end if

                enddo
             enddo

          enddo
       enddo
    enddo

  end subroutine makecoef
  !-------------------------------------------------

  subroutine ideal_mixture_transport(density,temperature,pressure,Yk,Xk,eta,kappa,zeta,diff_ij,chitil) bind(C,name="ideal_mixture_transport")  

    real(amrex_real), intent(in   ) :: density,temperature,pressure
    real(amrex_real), intent(in   ) :: Yk(nspecies),Xk(nspecies)
    real(amrex_real), intent(inout) :: eta,kappa,zeta
    real(amrex_real), intent(inout) :: diff_ij(nspecies,nspecies)
    real(amrex_real), intent(inout) :: chitil(nspecies)

    ! Local variables
    real(amrex_real) :: Dbin(nspecies,nspecies)
    real(amrex_real) :: omega11(nspecies,nspecies)
    real(amrex_real) :: sigma11(nspecies,nspecies)
    real(amrex_real) :: a_ij1(nspecies,nspecies)
    real(amrex_real) :: a_ij2(nspecies,nspecies)
    real(amrex_real) :: alphabar(nspecies,nspecies)

    real(amrex_real) :: xxtr(nspecies), yytr(nspecies), molecular_mass(nspecies)
    
    real(amrex_real) :: Mwmix, sqrtT

    integer :: ii, i, j

    real(amrex_real) :: sigma11bar,diamat

    real(amrex_real) :: mu, Fijstar, Fij

    real(amrex_real), PARAMETER :: pi= 3.1415926535897932d0

    !====================================================
    
    ! compute molecular masses by dividing by Avogadro's
    do ii = 1, nspecies
       molecular_mass(ii) = molmass(ii)*(k_B/Runiv)
    enddo

    ! mole fractions correction - EGLIB
    do ii = 1, nspecies
       ! GM: Why this factor of 1E-15???
       xxtr(ii) = Xk(ii) + (1.0d-15)*(sum(Xk(:))/dble(nspecies)-Xk(ii))
    enddo

    ! molecular weight of mixture - EGLIB
    Mwmix = 0.0d0
    do ii = 1, nspecies
       MWmix = MWmix + xxtr(ii)*molmass(ii)
    enddo

    ! mass fractions correction - EGLIB
    do ii = 1, nspecies
       yytr(ii) = molmass(ii)/MWmix*xxtr(ii)
    enddo

    Fijstar = -27.0d0/5.0d0

    ! find binary diffusion coefficients  
    ! HCB 8.2-9   
    sqrtT = dsqrt(temperature)

    ! print*, "Hack (imt): ", xxtr, yytr
    ! stop

    do i = 1, nspecies
       do j = 1, nspecies

          ! These matrices are computed in init_chemistry in FluctHydro code

          diamat = 0.5d0*(diameter(i) + diameter(j))
          mu = molecular_mass(i)*molecular_mass(j)/(molecular_mass(i) + molecular_mass(j))
          Fij = (6.0d0*molecular_mass(i)*molecular_mass(i) + 13.0d0/5.0d0*molecular_mass(j)*molecular_mass(j) +   &
               16.0d0/5.0d0*molecular_mass(i)*molecular_mass(j))/((molecular_mass(i)+molecular_mass(j))**2.0d0)

          !!

          sigma11bar = sqrt( k_b/(2.0d0*pi*mu) )*pi*(diamat**2.0d0)
          
          !!

          alphabar(i,j) = 8.0d0/(3.0d0*k_b)*mu*mu*(-.5d0*sigma11bar)

          !!!!

          Dbin(i,j) = 3.0d0/16.0d0*sqrt(2.0d0*pi*k_b**3.d00    &   
               *(molecular_mass(i)+molecular_mass(j))/molecular_mass(i)/molecular_mass(j))/(pi*diamat**2.0d0) &
               *temperature*sqrtT/pressure
          omega11(i,j) = sqrt(pi*k_b/(2.0d0*mu))*diamat**2.0d0*sqrtT
          sigma11(i,j) = sigma11bar*sqrtT
          a_ij1(i,j) = (5.0d0/(k_b)*molecular_mass(i)*molecular_mass(j)/    &
               (molecular_mass(i)+molecular_mass(j))*Fij*sigma11bar) /sqrtT
          a_ij2(i,j) = (5.0d0/(k_b)*molecular_mass(i)*molecular_mass(j)*molecular_mass(i)*molecular_mass(j)/    &
               ((molecular_mass(i)+molecular_mass(j))**3.0d0)*Fijstar*sigma11bar) /sqrtT

       enddo
    enddo

    call visc_lin(omega11,yytr,temperature,density,molecular_mass,eta)
    ! eta = 0.0d0
   
    ! ! GCM: Why hard-coded in? from original LLNS
    zeta = 0.d0
    
    call lambda_lin(Dbin,omega11,yytr,temperature,density,molecular_mass,kappa)
    ! kappa = 0.d0

    call D_GIOVANGIGLI(Dbin,yytr,xxtr,diff_ij)
    ! diff_ij = 0.0d0

    call thermalDiff(sigma11,a_ij1,a_ij2,alphabar,xxtr,sqrtT,molecular_mass,chitil)
    ! chitil = 0.0d0

  end subroutine ideal_mixture_transport

  !-------------------------------------------------------------------------------------

  subroutine visc_lin(omega11,Ykp,T,rho,mk,etaMix) bind(C,name="visc_lin")  

    implicit none 

    real(amrex_real), intent(in   ) :: omega11(nspecies,nspecies), mk(nspecies), Ykp(nspecies)
    real(amrex_real), intent(in   ) :: T, rho
    real(amrex_real), intent(inout) :: etaMix

    integer :: ii, jj, kk
    real(amrex_real) :: nk(nspecies), QoR(nspecies,nspecies), bSonine(nspecies), diag(nspecies)
    real(amrex_real) :: sum1
    integer :: ip(nspecies)


    ! find number density
    do ii = 1, nspecies
       nk(ii) = rho*Ykp(ii)/mk(ii)
    enddo

    ! print*, "Hack (visc_lin): T = ", T
    ! stop

    do  ii = 1, nspecies
       diag(ii) = 0.d0
       do kk = 1, nspecies
          diag(ii) = diag(ii) + nk(kk)*mk(kk)/(mk(ii)+mk(kk))**2.0d0 *      &
               (5.0d0*mk(ii)*omega11(ii,kk) +  &
               3.0d0*mk(kk)*omega11(ii,kk))

       enddo
    enddo

    ! print*, "Hack (visc_lin): ", mk, nk, Ykp
    ! stop

    do  ii = 1, nspecies
       do jj = 1, nspecies
          sum1 = -nk(jj)*mk(jj)/(mk(ii)+mk(jj))**2.0d0 *      &
               2.0d0*mk(jj)*omega11(ii,jj) 

          if(ii.eq.jj)then
             sum1 = sum1 + diag(ii)
          endif

          QoR(ii,jj) = -(16.0d0/15.0d0)*(mk(ii)/mk(jj))*sum1  
          ! HCB 7.4-62,63
       enddo
       bsonine(ii) = -1.0d0
    enddo

    ! print*, "Hack (visc_lin): original QoR = ", QoR

    ! QoR = 1.0d24*QoR
    ! bsonine = 1.0d24*bsonine

    ! print*, "Hack (visc_lin): scaled QoR = ", QoR
        
    call decomp(nspecies,nspecies,QoR,ip)
    call solve(nspecies,nspecies,QoR,bsonine,ip)

    ! print*, "Hack (visc_lin): factored QoR = ", QoR
    ! stop

    sum1 = 0.0d0
    do ii = 1, nspecies
       sum1 = sum1 + nk(ii)*bSonine(ii)
       ! HCB 7.4-56
    enddo


    etaMix = 0.5d0*k_b*T*sum1

    ! print*, "Hack (visc_lin): etaMix = ", etaMix, T
    ! stop

  end subroutine visc_lin

  !-------------------------------------------------------------------------

  subroutine lambda_lin(Dbin,omega11,Ykp,T,rho,mk,lammix) bind(C,name="lambda_lin")  

    implicit none 

    real(amrex_real), intent(in   ) :: T, rho, mk(nspecies)
    real(amrex_real), intent(inout) :: lammix
    real(amrex_real), intent(in   ) :: Ykp(nspecies)
    real(amrex_real), intent(in   ) :: omega11(nspecies,nspecies)
    real(amrex_real), intent(in   ) :: Dbin(nspecies,nspecies)
    
    real(amrex_real) :: mw
    real(amrex_real) :: nk(nspecies), beta(nspecies)
    integer :: i, j, k
    real(amrex_real) :: sum00, sum01, sum11, mu
    real(amrex_real) :: sum1, lamdaprime, ntotal
    
    real(amrex_real) :: QQ(2*nspecies,2*nspecies)
    real(amrex_real) :: aSonine(2*nspecies), D_T(nspecies)

    real(amrex_real) :: diag1(nspecies),diag2(nspecies),diag3(nspecies),diag4(nspecies)
    real(amrex_real) :: ratm,sqratm,scr
    integer :: ip(2*nspecies)


    ! find number density
    do i = 1, nspecies
       nk(i) = rho*Ykp(i)/mk(i)
    enddo


    diag1=0.d0
    diag2=0.d0
    diag3=0.d0
    diag4=0.d0

    do i=1,nspecies
       do k=1,nspecies

          diag1(i) =diag1(i)+ nk(k)*mk(k)/(mk(i)+mk(k)) *   &
               nk(i)*mk(i)*omega11(i,k)
          if(k.ne.i)then
             diag2(i) = diag2(i)+nk(k)*mk(k)/(mk(i)+mk(k))*omega11(i,k)
          endif
          diag3(i) = diag3(i) + nk(i)*nk(k)*mk(k)**2.0d0/(mk(i)+mk(k))**2.0d0  &
               *  .5d0*omega11(i,k)

          diag4(i) = diag4(i)+ nk(i)*nk(k)*mk(k)/(mk(i)+mk(k))**3.0d0 *    &
               ((5.0d0/4.0d0*(6.0d0*mk(i)**2.0d0+5.0d0*mk(k)**2.0d0)*omega11(i,k) &
               - 3.0d0*mk(k)**2.0d0*omega11(i,k) ) & 
               +  4.0d0*mk(i)*mk(k)*omega11(i,k) )

       enddo
    enddo

    do i = 1, nspecies
       do j = 1, nspecies

          sum00 = -nk(j)*mk(j)*diag2(i) - nk(j)*mk(j)/(mk(i)+mk(j))*nk(i)*mk(i) &
               *omega11(i,j)
          sum01 = - nk(i)*nk(j)*mk(j)**2.0d0/(mk(i)+mk(j))**2.0d0 *  &
               .5d0*omega11(i,j)

          sum11 =  nk(i)*nk(j)*mk(j)/(mk(i)+mk(j))**3.0d0 *         &
               ( -(5.0d0/4.0d0*(6.0d0*mk(j)**2.0d0+5.0d0*mk(j)**2.0d0)*omega11(i,j) &
               - 3.0d0*mk(j)**2.0d0*omega11(i,j) ) & 
               +  4.0d0*mk(j)*mk(j)*omega11(i,j) )

          if(i.eq.j)then
             sum00 = sum00+diag1(i)
             sum01 = sum01+diag3(i)
             sum11 = sum11+diag4(i)
          endif

          ratm = mk(i)/mk(j)
          sqratm = sqrt(ratm)
          scr = ratm*sqratm
          QQ(i,j) =  8.0d0*sqratm/mk(i) * sum00         ! HCB 7.4-50 
          QQ(i,j+nspecies) = -8.0d0*scr * sum01    ! HCB 7.4-51
          QQ(i+nspecies,j) = -8.d0*sqratm*sum01
          QQ(i+nspecies,j+nspecies) = 8.0d0*ratm*sqratm* sum11  ! HCB 7.4-53

       enddo
    enddo

    ! Build vector rhs ; see HCB 7.4-54
    do i = 1, nspecies
       aSonine(i) = 0.0d0
       beta(i)= sqrt(2.0d0*k_b*T/mk(i))
       aSonine(i+nspecies) = -(15.0d0/4.0d0)*nk(i)*beta(i)
    enddo

    ! print*, "Hack (lambda_lin) predecomp: ", QQ

    ! NOTE: the minus sign below; see HCB p 488
    call decomp(2*nspecies,2*nspecies,QQ,ip)
    call solve(2*nspecies,2*nspecies,QQ,aSonine,ip)

    ! print*, "Hack (lambda_lin) postdecomp: ", QQ
    ! stop
    
    do i = 1, nspecies           
       ! HCB 7.4-9
       D_T(i) = 0.5d0*nk(i)*beta(i)*mk(i)*aSonine(i)
    enddo

    sum1=0              
    ! HCB 7.4-33
    do i = 1, nspecies
       sum1 = sum1 + nk(i)*beta(i)* aSonine(nspecies+i)
    enddo
    lamdaprime = -5.0d0/4.0d0 * k_b * sum1


    sum1 = 0
    ! HCB 7.4-65
    do i = 1, nspecies
       do j = 1, nspecies 
          sum1 = sum1 + nk(i)*nk(j)/Dbin(i,j) *     & 
               (D_T(i)/(nk(i)*mk(i)) - D_T(j)/(nk(j)*mk(j)))**2.0d0
       enddo
    enddo

    ntotal = sum(nk(:))


    lammix = lamdaprime - 0.5d0*k_b/nTotal * sum1  ! HCB 7.4-65

    ! print*, "Hack (lambda_lin): ", D_T, nk, beta, mk, aSonine
    ! print*, "Hack (lambda_lin): ", sum1, nk, Dbin, D_T, mk
    ! print*, "Hack (lambda_lin): ", lammix, lamdaprime, nTotal, sum1
    ! print*, "Hack (lambda_lin) inputs: ", Dbin,omega11,Ykp,T,rho,mk,lammix
    ! stop

  end subroutine lambda_lin
  !-------------------------------------------------------------------------

  !----------------------------------------------------------------------

  subroutine D_GIOVANGIGLI(Dbin,Ykp,Xkp,D_tilde) bind(C,name="D_GIOVANGIGLI")  

    implicit none 

    real(amrex_real), intent(in   ) :: Ykp(nspecies), Xkp(nspecies)
    real(amrex_real), intent(in   ) :: Dbin(nspecies,nspecies) 
    real(amrex_real), intent(inout) :: D_tilde(nspecies,nspecies) 

    integer :: i, j, k, jj
    real(amrex_real) :: term1, term2, Deltamat
    real(amrex_real) :: Di(nspecies), Diff_ij(nspecies,nspecies)
    real(amrex_real) :: Zmat(nspecies,nspecies)
    real(amrex_real) :: Pmat(nspecies,nspecies), Jmat(nspecies,nspecies)
    real(amrex_real) :: Minv(nspecies), Mmat(nspecies)
    real(amrex_real) :: PJ(nspecies,nspecies), matrix1(nspecies,nspecies), matrix2(nspecies,nspecies)
    real(amrex_real) :: scr

    integer :: jmax

    jmax = 3 

    ! Find Di matrix 
    do i = 1, nspecies
       term2 = 0.0d0
       do j = 1, nspecies
          if(j.ne.i) then
             term2 = term2 + Xkp(j)/Dbin(i,j)
          endif
       enddo
       Di(i) = (1.d0-Ykp(i))/term2 
    enddo


    ! Compute Mmat and Minv
    do i = 1, nspecies
       Mmat(i) = Xkp(i)/Di(i)
       Minv(i) = Di(i)/Xkp(i)
    enddo


    ! Compute P matrix
    Pmat = 0.0d0
    do i = 1, nspecies
       do j = 1, nspecies
          Pmat(i,j) = - Ykp(j) 
          if(i.eq.j) then
             Pmat(i,j) =  Pmat(i,j) + 1.0d0  
          endif
       enddo
    enddo


    ! Compute Deltamat
    Deltamat = 0.0d0 
    do i = 1, nspecies
       do j = 1, nspecies
          if(i.eq.j) then
             term1 = 0.0d0
             do k = 1, nspecies
                if(k.ne.i) then
                   term1 = term1 + Xkp(i)*Xkp(k)/Dbin(i,k)
                endif
             enddo
             Deltamat = term1
          else
             Deltamat = -Xkp(i)*Xkp(j)/Dbin(i,j) 
          endif
          Zmat(i,j) = -Deltamat
       enddo
    enddo


    ! Compute Zmat
    do i = 1, nspecies
       Zmat(i,i) = Zmat(i,i) + Mmat(i)
    enddo


    ! Compute Jmat
    do i = 1, nspecies
       do j = 1, nspecies
          Jmat(i,j) = Minv(i)*Zmat(i,j)
       enddo
    enddo

    ! Compute PJ
    PJ = 0.0d0
    do i = 1, nspecies
       do j = 1, nspecies
          do k = 1, nspecies
             PJ(i,j) = PJ(i,j) + Pmat(i,k)*Jmat(k,j)
          enddo
       enddo
    enddo



    ! Compute P M^-1 Pt; store it in matrix2
    do i = 1, nspecies
       do j = 1, nspecies
          scr = 0.d0
          do k = 1, nspecies
             scr = scr + Pmat(i,k)*Minv(k)*Pmat(j,k) 
             ! notice the change in indices for Pmat to represent Pmat^t
          enddo
          matrix2(i,j) = scr
          Diff_ij(i,j) = scr
       enddo
    enddo

    if(jmax.gt.0)then

       do jj = 1,jmax

          !         matrix1=0
          do i = 1, nspecies
             do j = 1, nspecies
                scr = 0.d0
                do k = 1, nspecies
                   scr = scr + PJ(i,k)*Diff_ij(k,j)
                enddo
                matrix1(i,j) = scr+matrix2(i,j)
             enddo
          enddo

          Diff_ij=matrix1

       enddo

    endif



    ! Compute D_tilde
    do i = 1, nspecies
       do j = 1, nspecies
          D_tilde(i,j) = Diff_ij(i,j)*Ykp(i)
       enddo
    enddo

  end subroutine D_GIOVANGIGLI


  !----------------------------------------------------------------------

  subroutine thermalDiff(sigma11,a_ij1,a_ij2,alphabar,Xkp,sqrtT,mk,kT) bind(C,name="thermalDiff")  

    implicit none

    real(amrex_real), intent(in   ) :: mk(nspecies),  Xkp(nspecies)
    real(amrex_real), intent(in   ) :: alphabar(nspecies,nspecies)
    real(amrex_real), intent(in   ) :: sqrtT
    real(amrex_real), intent(inout) :: kT(nspecies)
    real(amrex_real), intent(in   ) :: a_ij1(nspecies,nspecies), a_ij2(nspecies,nspecies), sigma11(nspecies,nspecies)
    
    integer :: i, j, k
    real(amrex_real) :: pi

    real(amrex_real) :: Aij(nspecies,nspecies)
    real(amrex_real) :: AA(nspecies)
    real(amrex_real) :: fact1, sumTemp, sum1
    real(amrex_real) :: alphaij(nspecies,nspecies)
    integer :: ip(nspecies)

    ! Based on Valk 1963 (Waldmann)

    do i = 1, nspecies
       do j = 1, nspecies
          Aij(i,j) = Xkp(j)*a_ij2(i,j)
       enddo
       AA(i) = 15.d0/4.d0
    enddo

    do i = 1, nspecies
       sumTemp = 0.0d0
       do k = 1, nspecies
          sumTemp = sumTemp + Xkp(k)*a_ij1(i,k)
       enddo
       Aij(i,i) = Aij(i,i) + sumTemp
    enddo

    !          call decomp(nspecies,nspecies,Aij,ip)
    !          call solve(nspecies,nspecies,Aij,AA,ip)
    call decompnp(nspecies,nspecies,Aij)
    call solvenp(nspecies,nspecies,Aij,AA)

    do i = 1, nspecies
       do j = 1, nspecies

          alphaij(i,j) = alphabar(i,j)*(AA(i)/mk(i) - AA(j)/mk(j))/sqrtT

       enddo
    enddo

    do k = 1, nspecies
       sumTemp = 0.0d0
       do j = 1, nspecies
          sumTemp = sumTemp + Xkp(j)*alphaij(k,j)
       enddo
       kT(k) = sumTemp
    enddo

  end subroutine thermalDiff


  !-------------------------------------------------


end module multispec_module
