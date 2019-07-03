module flux_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, molmass, cell_depth, k_b, runiv, bc_lo, bc_hi, n_cells, membrane_cell, visc_type, algorithm_type
  use conv_module, only : get_temperature, get_pressure_gas, get_energy, get_enthalpies, get_temperature_gas, get_density_gas, get_energy_gas, get_hc_gas
  use multispec_module, only : cholesky_decomp

  implicit none

  private

  public :: diff_flux, stoch_flux_BC

contains


  subroutine stoch_flux(lo,hi, cons, fluxx, fluxy, &
#if (AMREX_SPACEDIM == 3)
       fluxz, &
#endif
       ranfluxx, ranfluxy, &
#if (AMREX_SPACEDIM == 3)
       ranfluxz, &
#endif
       dx, dt) bind(C,name="stoch_flux")

    integer         , intent(in   ) :: lo(3),hi(3)
    real(amrex_real), intent(in   ) :: dx(3), dt
    real(amrex_real), intent(inout) :: fluxx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
    real(amrex_real), intent(inout) :: fluxy(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
    real(amrex_real), intent(inout) :: fluxz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif

    real(amrex_real), intent(inout) :: ranfluxx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
    real(amrex_real), intent(inout) :: ranfluxy(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
    real(amrex_real), intent(inout) :: ranfluxz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif

    real(amrex_real), intent(in   ) :: cons(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)


    real(amrex_real) ::etatF, kappattF, dtinv, volinv, sFac, qFac, velu, velv, velw, wgt1, wgt2, weiner(5+nspecies), fweights(5+nspecies), nweight, muzepp, muzemp, muzepm, muzemm, phiflx, muxp, muyp, muzp, kxp, kyp, kzp, meanT

    real(amrex_real) :: hk(nspecies), yy(nspecies), yyp(nspecies), sumy, sumyp, DijY_edge(nspecies,nspecies), sqD(nspecies,nspecies), soret, MWmix

    integer :: i,j,k,l
    integer :: ll, ns

    dtinv = 1d0/dt
#if (AMREX_SPACEDIM == 3)
    volinv = 1d0/(dx(1)*dx(2)*dx(3))
#endif

#if (AMREX_SPACEDIM == 2)
    volinv = 1d0/(dx(1)*dx(2)*cell_depth)
#endif

    sFac = 2d0*4d0*k_b*volinv*dtinv/3d0
    qFac = 2d0*k_b*volinv*dtinv

    if (abs(visc_type) .gt. 1) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! JB's tensor form !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!! x-flux !!!!!!!!!!!!!!!!!!!

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)+1




             end do
          end do
       end do

       ! print*, "Hack: got here (end) stochflux"

!!!!!!!!!!!!!!!!!!! y-flux !!!!!!!!!!!!!!!!!!!

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)

 
             end do
          end do
       end do

!!!!!!!!!!!!!!!!!!! z-flux !!!!!!!!!!!!!!!!!!!

       do k = lo(3),hi(3)+1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

 

             end do
          end do
       end do

    else

    endif

!!!!!!!!!!!!! Enforce flux boundary conditions !!!!!!!!!!!!!

    call stoch_flux_BC(lo,hi, fluxx, fluxy, &
#if (AMREX_SPACEDIM == 3)
         fluxz &
#endif
         )

  end subroutine stoch_flux

  subroutine stoch_flux_BC(lo,hi, xflux, yflux, &
#if (AMREX_SPACEDIM == 3)
       zflux &
#endif
       ) bind(C,name="stoch_flux_bounds")

    integer         , intent(in   ) :: lo(3),hi(3)
    real(amrex_real), intent(inout) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
    real(amrex_real), intent(inout) :: yflux(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
    real(amrex_real), intent(inout) :: zflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif

    real(amrex_real) :: sqrtTwo

    integer :: i,j,k,l

    sqrtTwo = sqrt(2.0)

!!!!!!!!!!!!!! x-flux BCs !!!!!!!!!!!!!!

    !if on lower bound and specular
    if((lo(1) .eq. 0) .and. (bc_lo(1) .eq. 1)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)

                xflux(0,j,k,l) = 0        

             end do
          end do
       end do
    endif

    !if on upper bound and specular
    if((hi(1) .eq. n_cells(1)-1) .and. (bc_hi(1) .eq. 1)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)

                xflux(hi(1)+1,j,k,l) = 0           

             end do
          end do
       end do
    endif

    !if on lower bound and diff
    if((lo(1) .eq. 0) .and. (bc_lo(1) .eq. 2)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)

                xflux(0,j,k,l) = sqrtTwo*xflux(0,j,k,l)

             end do
          end do
       end do
    endif

    !if on upper bound and diff
    if((hi(1) .eq. n_cells(1)-1) .and. (bc_hi(1) .eq. 2)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)

                xflux(hi(1)+1,j,k,l) = sqrtTwo*xflux(hi(1)+1,j,k,l)        

             end do
          end do
       end do
    endif

!!!!!!!!!!!!!! y-flux BCs !!!!!!!!!!!!!!

    !if on lower bound and specular
    if((lo(2) .eq. 0) .and. (bc_lo(2) .eq. 1)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)

                yflux(i,0,k,l) = 0        

             end do
          end do
       end do
    endif

    !if on upper bound and specular
    if((hi(2) .eq. n_cells(2)-1) .and. (bc_hi(2) .eq. 1)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)

                yflux(i,hi(2)+1,k,l) = 0           

             end do
          end do
       end do
    endif

    !if on lower bound and diff
    if((lo(2) .eq. 0) .and. (bc_lo(2) .eq. 2)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)

                yflux(i,0,k,l) = sqrtTwo*yflux(i,0,k,l)

             end do
          end do
       end do
    endif

    !if on upper bound and diff
    if((hi(2) .eq. n_cells(2)-1) .and. (bc_hi(2) .eq. 2)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)

                yflux(i,hi(2)+1,k,l) = sqrtTwo*yflux(i,hi(2)+1,k,l)        

             end do
          end do
       end do
    endif

!!!!!!!!!!!!!! z-flux BCs !!!!!!!!!!!!!!

    !if on lower bound and specular
    if((lo(3) .eq. 0) .and. (bc_lo(3) .eq. 1)) then
       do l = 2,nvars
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                zflux(i,j,0,l) = 0

             end do
          end do
       end do
    endif

    !if on upper bound and specular
    if((hi(3) .eq. n_cells(3)-1) .and. (bc_hi(3) .eq. 1)) then
       do l = 2,nvars
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                zflux(i,j,hi(3)+1,l) = 0           

             end do
          end do
       end do
    endif

    !if on lower bound and diff
    if((lo(3) .eq. 0) .and. (bc_lo(3) .eq. 2)) then
       do l = 2,nvars
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                zflux(i,j,0,l) = sqrtTwo*zflux(i,j,0,l)

             end do
          end do
       end do
    endif

    !if on upper bound and diff
    if((hi(3) .eq. n_cells(3)-1) .and. (bc_hi(3) .eq. 2)) then
       do l = 2,nvars
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                zflux(i,j,hi(3)+1,l) = sqrtTwo*zflux(i,j,hi(3)+1,l)        

             end do
          end do
       end do
    endif


  end subroutine stoch_flux_BC

  subroutine diff_flux(lo,hi, cons, prim, & 
       eta, zeta, kappa, & 
       chi, Dij, &
       fluxx, fluxy, &
#if (AMREX_SPACEDIM == 3)
       fluxz, &
#endif
       cornux, cornvx, cornwx, & 
       cornuy, cornvy, cornwy, & 
       cornuz, cornvz, cornwz, & 
       visccorn, dx) bind(C,name="diff_flux")

    integer         , intent(in   ) :: lo(3),hi(3)
    real(amrex_real), intent(in   ) :: dx(3)
    real(amrex_real), intent(inout) :: fluxx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),nvars)
    real(amrex_real), intent(inout) :: fluxy(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),nvars)
#if (AMREX_SPACEDIM == 3)
    real(amrex_real), intent(inout) :: fluxz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,nvars)
#endif

    real(amrex_real), intent(inout) :: cornux(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
    real(amrex_real), intent(inout) :: cornvx(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
    real(amrex_real), intent(inout) :: cornwx(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)

    real(amrex_real), intent(inout) :: cornuy(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
    real(amrex_real), intent(inout) :: cornvy(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
    real(amrex_real), intent(inout) :: cornwy(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)

    real(amrex_real), intent(inout) :: cornuz(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
    real(amrex_real), intent(inout) :: cornvz(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
    real(amrex_real), intent(inout) :: cornwz(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)

    real(amrex_real), intent(inout) :: visccorn(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)

    real(amrex_real), intent(in   ) :: cons(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nvars)
    real(amrex_real), intent(in   ) :: prim(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nprimvars)

    real(amrex_real), intent(in   ) :: eta(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
    real(amrex_real), intent(in   ) :: zeta(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
    real(amrex_real), intent(in   ) :: kappa(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))

    real(amrex_real), intent(in   ) :: chi(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nspecies)
    real(amrex_real), intent(in   ) :: Dij(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nspecies,nspecies)

    integer :: i,j,k 
    real(amrex_real) :: half, two, dxinv(3), twothirds, muxp, kxp, tauxxp, tauyxp, tauzxp, divxp, phiflx, muyp, kyp, tauxyp, tauyyp, tauzyp, divyp, zetaxp, onetwelfth
#if (AMREX_SPACEDIM == 3)
    real(amrex_real) :: muzp, kzp, tauxzp, tauyzp, tauzzp, divzp
#endif

    ! Multispecies local
    real(amrex_real) :: term1, term2, Q5
    real(amrex_real) :: dk(nspecies), Fk(nspecies), hk(nspecies), soret(nspecies), meanXk(nspecies), meanYk(nspecies)
    real(amrex_real) :: meanT, meanP
    integer :: ns, kk, ll

    dxinv = 1d0/dx

    two = 2.d0
    half = 0.5d0
    twothirds = 2d0/3d0
    onetwelfth = 1d0/12d0

    !x flux
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1

             muxp = half*(eta(i,j,k) + eta(i-1,j,k))
             kxp = half*(kappa(i,j,k) + kappa(i-1,j,k))

             tauxxp = muxp*(prim(i,j,k,2) - prim(i-1,j,k,2))/dx(1)
             tauyxp = muxp*(prim(i,j,k,3) - prim(i-1,j,k,3))/dx(1)
             tauzxp = muxp*(prim(i,j,k,4) - prim(i-1,j,k,4))/dx(1)
             divxp = 0

             phiflx =  tauxxp*(prim(i-1,j,k,2)+prim(i,j,k,2))                &
                  &      +  divxp*(prim(i-1,j,k,2)+prim(i,j,k,2))                    &
                  &      +  tauyxp*(prim(i-1,j,k,3)+prim(i,j,k,3))                   &
                  &      +  tauzxp*(prim(i-1,j,k,4)+prim(i,j,k,4))

             fluxx(i,j,k,2) = fluxx(i,j,k,2) - (tauxxp+divxp)
             fluxx(i,j,k,3) = fluxx(i,j,k,3) - tauyxp
             fluxx(i,j,k,4) = fluxx(i,j,k,4) - tauzxp
             fluxx(i,j,k,5) = fluxx(i,j,k,5) - (half*phiflx                  &
                  &                      + kxp*(prim(i,j,k,5)-prim(i-1,j,k,5))/dx(1))

             meanT = 0.5d0*(prim(i-1,j,k,5)+prim(i,j,k,5))
             meanP = 0.5d0*(prim(i-1,j,k,6)+prim(i,j,k,6))

             if(algorithm_type.eq.2) then

                ! compute dk
                do ns = 1, nspecies
                   term1 = (prim(i,j,k,6+nspecies+ns)-prim(i-1,j,k,6+nspecies+ns))/dx(1)
                   meanXk(ns) = 0.5d0*(prim(i-1,j,k,6+nspecies+ns)+prim(i,j,k,6+nspecies+ns))
                   meanYk(ns) = 0.5d0*(prim(i-1,j,k,6+ns)+prim(i,j,k,6+ns))
                   term2 = (meanXk(ns)-meanYk(ns))*(prim(i,j,k,6)-prim(i-1,j,k,6))/dx(1)/meanP
                   dk(ns) = term1 + term2 
                   soret(ns) = 0.5d0*(chi(i-1,j,k,ns)*prim(i-1,j,k,6+nspecies+ns)+chi(i,j,k,ns)*prim(i,j,k,6+nspecies+ns))  &
                        *(prim(i,j,k,5)-prim(i-1,j,k,5))/dx(1)/meanT
                enddo

                ! compute Fk (based on Eqn. 2.5.24, Giovangigli's book)
                Fk = 0.0d0 
                do kk = 1, nspecies
                   do ll = 1, nspecies
                      Fk(kk) = Fk(kk) - half*(Dij(i-1,j,k,kk,ll)+Dij(i,j,k,kk,ll))*( dk(ll) +soret(ll))
                   enddo
                enddo

                ! compute Q (based on Eqn. 2.5.25, Giovangigli's book)
                call get_enthalpies(meanT, hk)

                Q5 = 0.0d0
                do ns = 1, nspecies
                   Q5 = Q5 + (hk(ns) + 0.5d0 * Runiv*meanT*(chi(i-1,j,k,ns)+chi(i,j,k,ns))/molmass(ns))*Fk(ns)  
                enddo
                ! heat conduction already included in flux(5)       


                fluxx(i,j,k,5) = fluxx(i,j,k,5) + Q5

                do ns = 1, nspecies  
                   fluxx(i,j,k,5+ns) = fluxx(i,j,k,5+ns) + Fk(ns)
                enddo

                do ll = 1, nvars
                   if ( isnan(fluxx(i,j,k,ll)) ) then
                      print*, "Hack 1, (diff_flux) in x = ", i,j,k, fluxx(i,j,k,:)
                      stop
                   end if
                end do

                ! print*, "Hack: loc = ", i,j,k
                ! print*, "Hack: chi = ", chi(i,j,k,:)
                ! print*, "Hack: chim = ", chi(i-1,j,k,:)
                ! print*, "Hack: Dij = ", Dij(i,j,k,:,:)
                ! print*, "Hack: Dijm = ", Dij(i-1,j,k,:,:)
                ! print*, "Hack: flux = ", fluxx(i,j,k,:)
                ! stop

             end if

          end do
       end do
    end do

    !y flux
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)

             muyp = half*(eta(i,j,k) + eta(i,j-1,k))
             kyp = half*(kappa(i,j,k) + kappa(i,j-1,k))

             tauxyp =  muyp*(prim(i,j,k,2) - prim(i,j-1,k,2))/dx(2)
             tauyyp =  muyp*(prim(i,j,k,3) - prim(i,j-1,k,3))/dx(2)
             tauzyp =  muyp*(prim(i,j,k,4) - prim(i,j-1,k,4))/dx(2)
             divyp = 0

             phiflx =                                                        &
                  &         tauxyp*(prim(i,j,k,2)+prim(i,j-1,k,2))                   &
                  &      +  tauyyp*(prim(i,j,k,3)+prim(i,j-1,k,3))                   &
                  &      +  divyp*(prim(i,j,k,3)+prim(i,j-1,k,3))                    &
                  &      +  tauzyp*(prim(i,j,k,4)+prim(i,j-1,k,4))

             fluxy(i,j,k,2) = fluxy(i,j,k,2) - tauxyp
             fluxy(i,j,k,3) = fluxy(i,j,k,3) - (tauyyp+divyp)
             fluxy(i,j,k,4) = fluxy(i,j,k,4) - tauzyp
             fluxy(i,j,k,5) = fluxy(i,j,k,5) - (half*phiflx                  &
                  &                      + kyp*(prim(i,j,k,5)-prim(i,j-1,k,5))/dx(2))

             meanT = 0.5d0*(prim(i,j-1,k,5)+prim(i,j,k,5))
             meanP = 0.5d0*(prim(i,j-1,k,6)+prim(i,j,k,6))

             if(algorithm_type.eq.2) then
                ! compute dk  

                do ns = 1, nspecies
                   term1 = (prim(i,j,k,6+nspecies+ns)-prim(i,j-1,k,6+nspecies+ns))/dx(2)
                   meanXk(ns) = 0.5d0*(prim(i,j-1,k,6+nspecies+ns)+prim(i,j,k,6+nspecies+ns))
                   meanYk(ns) = 0.5d0*(prim(i,j-1,k,6+ns)+prim(i,j,k,6+ns))
                   term2 = (meanXk(ns)-meanYk(ns))*(prim(i,j,k,6)-prim(i,j-1,k,6))/dx(2)/meanP
                   dk(ns) = term1 + term2 
                   soret(ns) = 0.5d0*(chi(i,j-1,k,ns)*prim(i,j-1,k,6+nspecies+ns)+chi(i,j,k,ns)*prim(i,j,k,6+nspecies+ns))  &
                        *(prim(i,j,k,5)-prim(i,j-1,k,5))/dx(2)/meanT
                enddo

                ! compute Fk (based on Eqn. 2.5.24, Giovangigli's book)
                Fk = 0.0d0
                do kk = 1, nspecies
                   do ll = 1, nspecies
                      Fk(kk) = Fk(kk) - half*(Dij(i,j-1,k,kk,ll)+Dij(i,j,k,kk,ll))*( dk(ll) +soret(ll))
                   enddo
                enddo

                ! compute Q (based on Eqn. 2.5.25, Giovangigli's book)
                call get_enthalpies(meanT, hk)

                Q5 = 0.0d0
                do ns = 1, nspecies
                   Q5 = Q5 + (hk(ns) + 0.5d0 * Runiv*meanT*(chi(i,j-1,k,ns)+chi(i,j,k,ns))/molmass(ns))*Fk(ns)  
                enddo

                ! heat conduction already included in flux(5)

                fluxy(i,j,k,5) = fluxy(i,j,k,5) + Q5

                do ns = 1, nspecies
                   fluxy(i,j,k,5+ns) = fluxy(i,j,k,5+ns) + Fk(ns)
                enddo

             end if

          end do
       end do
    end do

#if (AMREX_SPACEDIM == 3)
    !z flux
    do k = lo(3),hi(3)+1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             muzp = half*(eta(i,j,k) + eta(i,j,k-1))
             kzp = half*(kappa(i,j,k) + kappa(i,j,k-1))

             tauxzp =  muzp*(prim(i,j,k,2) - prim(i,j,k-1,2))/dx(3)
             tauyzp =  muzp*(prim(i,j,k,3) - prim(i,j,k-1,3))/dx(3)
             tauzzp =  muzp*(prim(i,j,k,4) - prim(i,j,k-1,4))/dx(3)
             divzp = 0

             phiflx =                                                        &
                  &      +  tauxzp*(prim(i,j,k-1,2)+prim(i,j,k,2))                   &
                  &      +  tauyzp*(prim(i,j,k-1,3)+prim(i,j,k,3))                   &
                  &      +  tauzzp*(prim(i,j,k-1,4)+prim(i,j,k,4))                   &
                  &      +  divzp*(prim(i,j,k-1,4)+prim(i,j,k,4))

             fluxz(i,j,k,2) = fluxz(i,j,k,2) - tauxzp
             fluxz(i,j,k,3) = fluxz(i,j,k,3) - tauyzp
             fluxz(i,j,k,4) = fluxz(i,j,k,4) - (tauzzp+divzp)
             fluxz(i,j,k,5) = fluxz(i,j,k,5) - (half*phiflx                  &
                  &                      +kzp*(prim(i,j,k,5)-prim(i,j,k-1,5))/dx(3))

             meanT = 0.5d0*(prim(i,j,k-1,5)+prim(i,j,k,5))
             meanP = 0.5d0*(prim(i,j,k-1,6)+prim(i,j,k,6))

             if(algorithm_type.eq.2) then
                ! compute dk  

                do ns = 1, nspecies
                   term1 = (prim(i,j,k,6+nspecies+ns)-prim(i,j,k-1,6+nspecies+ns))/dx(3)
                   meanXk(ns) = 0.5d0*(prim(i,j,k-1,6+nspecies+ns)+prim(i,j,k,6+nspecies+ns))
                   meanYk(ns) = 0.5d0*(prim(i,j,k-1,6+ns)+prim(i,j,k,6+ns))
                   term2 = (meanXk(ns)-meanYk(ns))*(prim(i,j,k,6)-prim(i,j,k-1,6))/dx(3)/meanP
                   dk(ns) = term1 + term2 
                   soret(ns) = 0.5d0*(chi(i,j,k,ns)*prim(i,j,k-1,6+nspecies+ns)+chi(i,j,k+1,ns)*prim(i,j,k,6+nspecies+ns))  &
                        *(prim(i,j,k,5)-prim(i,j,k-1,5))/dx(3)/meanT
                enddo

                ! compute Fk (based on Eqn. 2.5.24, Giovangigli's book)
                Fk = 0.0d0
                do kk = 1, nspecies
                   do ll = 1, nspecies
                      Fk(kk) = Fk(kk) - half*(Dij(i,j,k-1,kk,ll)+Dij(i,j,k,kk,ll))*( dk(ll) +soret(ll))
                   enddo
                enddo

                ! compute Q (based on Eqn. 2.5.25, Giovangigli's book)
                call get_enthalpies(meanT, hk)

                Q5 = 0.0d0
                do ns = 1, nspecies
                   Q5 = Q5 + (hk(ns) + 0.5d0 * Runiv*meanT*(chi(i,j,k,ns)+chi(i,j,k,ns))/molmass(ns))*Fk(ns)  
                enddo

                ! heat conduction already included in flux(5)

                fluxz(i,j,k,5) = fluxz(i,j,k,5) + Q5

                do ns = 1, nspecies
                   fluxz(i,j,k,5+ns) = fluxz(i,j,k,5+ns) + Fk(ns)
                enddo

             end if

          end do
       end do
    end do
#endif


    do k = lo(3),hi(3)+1
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)+1

             ! Corner viscosity
             muxp = 0.125d0*(eta(i,j-1,k-1) + eta(i-1,j-1,k-1) + eta(i,j,k-1) + eta(i-1,j,k-1)+ &
                  eta(i,j-1,k) + eta(i-1,j-1,k) + eta(i,j,k) + eta(i-1,j,k))
             if (abs(visc_type) .eq. 3) then
                zetaxp = 0.125d0*(zeta(i,j-1,k-1) + zeta(i-1,j-1,k-1) + zeta(i,j,k-1) + zeta(i-1,j,k-1)+ &
                     zeta(i,j-1,k) + zeta(i-1,j-1,k) + zeta(i,j,k) + zeta(i-1,j,k))
             else
                zetaxp = 0.0
             endif

             cornux(i,j,k) = 0.25d0*muxp*(prim(i,j-1,k-1,2)-prim(i-1,j-1,k-1,2) + prim(i,j,k-1,2)-prim(i-1,j,k-1,2)+ &
                  prim(i,j-1,k,2)-prim(i-1,j-1,k,2) + prim(i,j,k,2)-prim(i-1,j,k,2))/dx(1)
             cornvx(i,j,k) = 0.25d0*muxp*(prim(i,j-1,k-1,3)-prim(i-1,j-1,k-1,3) + prim(i,j,k-1,3)-prim(i-1,j,k-1,3)+ &
                  prim(i,j-1,k,3)-prim(i-1,j-1,k,3) + prim(i,j,k,3)-prim(i-1,j,k,3))/dx(1)
             cornwx(i,j,k) = 0.25d0*muxp*(prim(i,j-1,k-1,4)-prim(i-1,j-1,k-1,4) + prim(i,j,k-1,4)-prim(i-1,j,k-1,4)+ &
                  prim(i,j-1,k,4)-prim(i-1,j-1,k,4) + prim(i,j,k,4)-prim(i-1,j,k,4))/dx(1)

             cornuy(i,j,k) = 0.25d0*muxp* (prim(i-1,j,k-1,2)-prim(i-1,j-1,k-1,2) + prim(i,j,k-1,2)-prim(i,j-1,k-1,2) + &
                  prim(i-1,j,k,2)-prim(i-1,j-1,k,2) + prim(i,j,k,2)-prim(i,j-1,k,2))/dx(2)
             cornvy(i,j,k) = 0.25d0*muxp* (prim(i-1,j,k-1,3)-prim(i-1,j-1,k-1,3) + prim(i,j,k-1,3)-prim(i,j-1,k-1,3) + &
                  prim(i-1,j,k,3)-prim(i-1,j-1,k,3) + prim(i,j,k,3)-prim(i,j-1,k,3))/dx(2)
             cornwy(i,j,k) = 0.25d0*muxp* (prim(i-1,j,k-1,4)-prim(i-1,j-1,k-1,4) + prim(i,j,k-1,4)-prim(i,j-1,k-1,4) + &
                  prim(i-1,j,k,4)-prim(i-1,j-1,k,4) + prim(i,j,k,4)-prim(i,j-1,k,4))/dx(2)

             cornuz(i,j,k) = 0.25d0*muxp*(prim(i-1,j-1,k,2)-prim(i-1,j-1,k-1,2) + prim(i,j-1,k,2)-prim(i,j-1,k-1,2) + &
                  prim(i-1,j,k,2)-prim(i-1,j,k-1,2) + prim(i,j,k,2)-prim(i,j,k-1,2))/dx(3)
             cornvz(i,j,k) = 0.25d0*muxp*(prim(i-1,j-1,k,3)-prim(i-1,j-1,k-1,3) + prim(i,j-1,k,3)-prim(i,j-1,k-1,3) + &
                  prim(i-1,j,k,3)-prim(i-1,j,k-1,3) + prim(i,j,k,3)-prim(i,j,k-1,3))/dx(3)
             cornwz(i,j,k) = 0.25d0*muxp*(prim(i-1,j-1,k,4)-prim(i-1,j-1,k-1,4) + prim(i,j-1,k,4)-prim(i,j-1,k-1,4) + &
                  prim(i-1,j,k,4)-prim(i-1,j,k-1,4) + prim(i,j,k,4)-prim(i,j,k-1,4))/dx(3)

             visccorn(i,j,k) =  (muxp/12d0+zetaxp/4d0)*( & ! Divergence stress
                  (prim(i,j-1,k-1,2)-prim(i-1,j-1,k-1,2))/dx(1) + (prim(i,j,k-1,2)-prim(i-1,j,k-1,2))/dx(1) + &
                  (prim(i,j-1,k,2)-prim(i-1,j-1,k,2))/dx(1) + (prim(i,j,k,2)-prim(i-1,j,k,2))/dx(1) + &
                  (prim(i-1,j,k-1,3)-prim(i-1,j-1,k-1,3))/dx(2) + (prim(i,j,k-1,3)-prim(i,j-1,k-1,3))/dx(2) + &
                  (prim(i-1,j,k,3)-prim(i-1,j-1,k,3))/dx(2) + (prim(i,j,k,3)-prim(i,j-1,k,3))/dx(2) + &
                  (prim(i-1,j-1,k,4)-prim(i-1,j-1,k-1,4))/dx(3) + (prim(i,j-1,k,4)-prim(i,j-1,k-1,4))/dx(3) + &
                  (prim(i-1,j,k,4)-prim(i-1,j,k-1,4))/dx(3) + (prim(i,j,k,4)-prim(i,j,k-1,4))/dx(3))

          end do
       end do
    end do

    !x flux
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1

             fluxx(i,j,k,2) = fluxx(i,j,k,2) - 0.25d0*(visccorn(i,j+1,k+1)+visccorn(i,j,k+1) +&
                  visccorn(i,j+1,k)+visccorn(i,j,k)) ! Viscous "divergence" stress

             fluxx(i,j,k,2) = fluxx(i,j,k,2) + .25d0*   &
                  (cornvy(i,j+1,k+1)+cornvy(i,j,k+1)+cornvy(i,j+1,k)+cornvy(i,j,k)  + &
                  cornwz(i,j+1,k+1)+cornwz(i,j,k+1)+cornwz(i,j+1,k)+cornwz(i,j,k))

             fluxx(i,j,k,3) = fluxx(i,j,k,3) - .25d0*   &
                  (cornuy(i,j+1,k+1)+cornuy(i,j,k+1)+cornuy(i,j+1,k)+cornuy(i,j,k))

             fluxx(i,j,k,4) = fluxx(i,j,k,4) - .25d0*   &
                  (cornuz(i,j+1,k+1)+cornuz(i,j,k+1)+cornuz(i,j+1,k)+cornuz(i,j,k))

             phiflx =  0.25d0*(visccorn(i,j+1,k+1)+visccorn(i,j,k+1) + &
                  visccorn(i,j+1,k)+visccorn(i,j,k) &
                  -(cornvy(i,j+1,k+1)+cornvy(i,j,k+1)+cornvy(i,j+1,k)+cornvy(i,j,k)  + &
                  cornwz(i,j+1,k+1)+cornwz(i,j,k+1)+cornwz(i,j+1,k)+cornwz(i,j,k))) * &
                  (prim(i-1,j,k,2)+prim(i,j,k,2))

             phiflx = phiflx + .25d0*   &
                  (cornuy(i,j+1,k+1)+cornuy(i,j,k+1)+cornuy(i,j+1,k)+cornuy(i,j,k)) * &
                  (prim(i-1,j,k,3)+prim(i,j,k,3))


             phiflx = phiflx + .25d0*   &
                  (cornuz(i,j+1,k+1)+cornuz(i,j,k+1)+cornuz(i,j+1,k)+cornuz(i,j,k)) * &
                  (prim(i-1,j,k,4)+prim(i,j,k,4))

             fluxx(i,j,k,5) = fluxx(i,j,k,5)-0.5d0*phiflx

          end do
       end do
    end do

    !y flux
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)

             fluxy(i,j,k,3) = fluxy(i,j,k,3) - &
                  0.25d0*(visccorn(i+1,j,k+1)+visccorn(i,j,k+1)+visccorn(i+1,j,k)+visccorn(i,j,k))

             fluxy(i,j,k,3) = fluxy(i,j,k,3) + .25d0*   &
                  (cornux(i+1,j,k+1)+cornux(i,j,k+1)+cornux(i+1,j,k)+cornux(i,j,k)  + &
                  cornwz(i+1,j,k+1)+cornwz(i,j,k+1)+cornwz(i+1,j,k)+cornwz(i,j,k))

             fluxy(i,j,k,2) = fluxy(i,j,k,2) - .25d0*   &
                  (cornvx(i+1,j,k+1)+cornvx(i,j,k+1)+cornvx(i+1,j,k)+cornvx(i,j,k))

             fluxy(i,j,k,4) = fluxy(i,j,k,4) - .25d0*   &
                  (cornvz(i+1,j,k+1)+cornvz(i,j,k+1)+cornvz(i+1,j,k)+cornvz(i,j,k))

             phiflx = 0.25d0*(visccorn(i+1,j,k+1)+visccorn(i,j,k+1)+visccorn(i+1,j,k)+visccorn(i,j,k) &
                  -(cornux(i+1,j,k+1)+cornux(i,j,k+1)+cornux(i+1,j,k)+cornux(i,j,k)  + &
                  cornwz(i+1,j,k+1)+cornwz(i,j,k+1)+cornwz(i+1,j,k)+cornwz(i,j,k))) * &
                  (prim(i,j-1,k,3)+prim(i,j,k,3))

             phiflx = phiflx + .25d0*   &
                  (cornvx(i+1,j,k+1)+cornvx(i,j,k+1)+cornvx(i+1,j,k)+cornvx(i,j,k)) * &
                  (prim(i,j-1,k,2)+prim(i,j,k,2))

             phiflx = phiflx + .25d0*   &
                  (cornvz(i+1,j,k+1)+cornvz(i,j,k+1)+cornvz(i+1,j,k)+cornvz(i,j,k)) * &
                  (prim(i,j-1,k,4)+prim(i,j,k,4))

             fluxy(i,j,k,5) = fluxy(i,j,k,5)-0.5d0*phiflx

          end do
       end do
    end do

    !z flux
    do k = lo(3),hi(3)+1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             fluxz(i,j,k,4) = fluxz(i,j,k,4) - &
                  0.25d0*(visccorn(i+1,j+1,k)+visccorn(i,j+1,k)+visccorn(i+1,j,k)+visccorn(i,j,k))

             fluxz(i,j,k,4) = fluxz(i,j,k,4) + .25d0*   &
                  (cornvy(i+1,j+1,k)+cornvy(i+1,j,k)+cornvy(i,j+1,k)+cornvy(i,j,k)  + &
                  cornux(i+1,j+1,k)+cornux(i+1,j,k)+cornux(i,j+1,k)+cornux(i,j,k))

             fluxz(i,j,k,2) = fluxz(i,j,k,2) - .25d0*   &
                  (cornwx(i+1,j+1,k)+cornwx(i+1,j,k)+cornwx(i,j+1,k)+cornwx(i,j,k))

             fluxz(i,j,k,3) = fluxz(i,j,k,3) - .25d0*   &
                  (cornwy(i+1,j+1,k)+cornwy(i+1,j,k)+cornwy(i,j+1,k)+cornwy(i,j,k)) 

             phiflx = 0.25d0*(visccorn(i+1,j+1,k)+visccorn(i,j+1,k)+visccorn(i+1,j,k)+visccorn(i,j,k) &
                  -(cornvy(i+1,j+1,k)+cornvy(i+1,j,k)+cornvy(i,j+1,k)+cornvy(i,j,k)  + &
                  cornux(i+1,j+1,k)+cornux(i+1,j,k)+cornux(i,j+1,k)+cornux(i,j,k))) *  &
                  (prim(i,j,k-1,4)+prim(i,j,k,4))

             phiflx = phiflx + .25d0*   &
                  (cornwx(i+1,j+1,k)+cornwx(i+1,j,k)+cornwx(i,j+1,k)+cornwx(i,j,k))* &
                  (prim(i,j,k-1,2)+prim(i,j,k,2))

             phiflx = phiflx + .25d0*   &
                  (cornwy(i+1,j+1,k)+cornwy(i+1,j,k)+cornwy(i,j+1,k)+cornwy(i,j,k)) * &
                  (prim(i,j,k-1,3)+prim(i,j,k,3))

             fluxz(i,j,k,5) = fluxz(i,j,k,5)-0.5d0*phiflx

          end do
       end do
    end do

  end subroutine diff_flux

end module flux_module

































