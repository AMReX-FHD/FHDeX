module flux_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, molmass, &
                                     cell_depth, k_b, runiv, n_cells, membrane_cell, &
                                     visc_type, algorithm_type, &
                                     bc_mass_lo, bc_mass_hi, bc_therm_lo, bc_therm_hi, &
                                     bc_vel_lo, bc_vel_hi
  use conv_module, only : get_temperature, get_pressure_gas, get_energy, &
                          get_enthalpies, get_temperature_gas, get_density_gas, &
                          get_energy_gas, get_hc_gas
  use multispec_module, only : cholesky_decomp
  implicit none

  private

  public :: stoch_flux_BC

contains

  subroutine stoch_flux(lo,hi, cons, prim, fluxx, fluxy, &
#if (AMREX_SPACEDIM == 3)
                        fluxz, &
#endif
                        ranfluxx, ranfluxy, &
#if (AMREX_SPACEDIM == 3)
                        ranfluxz, &
#endif
                        rancorn, & 
                        eta, zeta, kappa, &
                        chi, Dij, & 
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
    real(amrex_real), intent(in   ) :: cons(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nvars)
    real(amrex_real), intent(in   ) :: prim(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nprimvars)
    real(amrex_real), intent(in   ) :: rancorn(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
    real(amrex_real), intent(in   ) :: eta  (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
    real(amrex_real), intent(in   ) :: zeta (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
    real(amrex_real), intent(in   ) :: kappa(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
    real(amrex_real), intent(in   ) :: chi  (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nspecies)
    real(amrex_real), intent(in   ) :: Dij  (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nspecies,nspecies)

    real(amrex_real) ::etatF, kappattF, dtinv, volinv,  velu, velv, velw, wgt1, wgt2
    real(amrex_real) :: weiner(5+nspecies), fweights(5+nspecies), nweight, muzepp, muzemp, muzepm, muzemm
    real(amrex_real) :: phiflx, muxp, muyp, muzp, kxp, kyp, kzp, meanT

    real(amrex_real) :: hk(nspecies), yy(nspecies), yyp(nspecies), sumy, sumyp, DijY_edge(nspecies,nspecies), sqD(nspecies,nspecies)
    real(amrex_real) :: soret, MWmix

    integer :: i,j,k,l
    integer :: ll, ns

    dtinv = 1.d0/dt
#if (AMREX_SPACEDIM == 3)
    volinv = 1.d0/(dx(1)*dx(2)*dx(3))
#endif

#if (AMREX_SPACEDIM == 2)
    volinv = 1.d0/(dx(1)*dx(2)*cell_depth)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! JB's tensor form !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!!!!!!!!!!!!!!!!!!! x-flux !!!!!!!!!!!!!!!!!!!

#if 0
    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
    do i = lo(1),hi(1)+1

       muxp = (eta(i,j,k)*prim(i,j,k,5) + eta(i-1,j,k)*prim(i-1,j,k,5))
       kxp = (kappa(i,j,k)*prim(i,j,k,5)**2 + kappa(i-1,j,k)*prim(i-1,j,k,5)**2)

       !! Look into: zeta*p is not used in this function (original FluctHydro)
       ! if (abs(visc_type) .eq. 3) then 
       !    zetaxp = (zeta(i,j,k)*prim(i,j,k,5) + zeta(i-1,j,k)*prim(i-1,j,k,5))
       ! else
       !    zetaxp = 0.0
       ! endif

       meanT = 0.5d0*(prim(i,j,k,5)+prim(i-1,j,k,5))

       ! Weights for facial fluxes:
       fweights(1) = 0 ! No mass flux
       fweights(2:4)=sqrt(k_b*muxp*volinv*dtinv)
       fweights(5)=sqrt(k_b*kxp*volinv*dtinv)

       ! Construct the random increments
       weiner(1:5) = fweights(1:5)*ranfluxx(i,j,k,1:5)

       nweight=sqrt(k_b*volinv*dtinv)

       if(n_cells(3).gt.1) then

          ! Corner viscosity coefficients in 3D
          muzepp = 0.25d0*(eta(i,j,k)*prim(i,j,k,5) + eta(i-1,j,k)*prim(i-1,j,k,5) + &
               eta(i,j+1,k)*prim(i,j+1,k,5) + eta(i-1,j+1,k)*prim(i-1,j+1,k,5) + &
               eta(i,j,k+1)*prim(i,j,k+1,5) + eta(i-1,j,k+1)*prim(i-1,j,k+1,5) + &
               eta(i,j+1,k+1)*prim(i,j+1,k+1,5) + eta(i-1,j+1,k+1)*prim(i-1,j+1,k+1,5) )/3.d0
          muzemp = 0.25d0*(eta(i,j-1,k)*prim(i,j-1,k,5) + eta(i-1,j-1,k)*prim(i-1,j-1,k,5) + &
               eta(i,j,k)*prim(i,j,k,5) + eta(i-1,j,k)*prim(i-1,j,k,5) + &
               eta(i,j-1,k+1)*prim(i,j-1,k+1,5) + eta(i-1,j-1,k+1)*prim(i-1,j-1,k+1,5) + &
               eta(i,j,k+1)*prim(i,j,k+1,5) + eta(i-1,j,k+1)*prim(i-1,j,k+1,5) )/3.d0
          muzepm = 0.25d0*(eta(i,j,k-1)*prim(i,j,k-1,5) + eta(i-1,j,k-1)*prim(i-1,j,k-1,5) + &
               eta(i,j+1,k-1)*prim(i,j+1,k-1,5) + eta(i-1,j+1,k-1)*prim(i-1,j+1,k-1,5) + &
               eta(i,j,k)*prim(i,j,k,5) + eta(i-1,j,k)*prim(i-1,j,k,5) + &
               eta(i,j+1,k)*prim(i,j+1,k,5) + eta(i-1,j+1,k)*prim(i-1,j+1,k,5) )/3.d0
          muzemm = 0.25d0*(eta(i,j-1,k-1)*prim(i,j-1,k-1,5) + eta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,5) + &
               eta(i,j,k-1)*prim(i,j,k-1,5) + eta(i-1,j,k-1)*prim(i-1,j,k-1,5) + &
               eta(i,j-1,k)*prim(i,j-1,k,5) + eta(i-1,j-1,k)*prim(i-1,j-1,k,5) + &
               eta(i,j,k)*prim(i,j,k,5) + eta(i-1,j,k)*prim(i-1,j,k,5) )/3.d0

          if (abs(visc_type) .eq. 3) then

             muzepp = muzepp + 0.25d0*(zeta(i,j,k)*prim(i,j,k,5) + zeta(i-1,j,k)*prim(i-1,j,k,5) + &
                  zeta(i,j+1,k)*prim(i,j+1,k,5) + zeta(i-1,j+1,k)*prim(i-1,j+1,k,5) + &
                  zeta(i,j,k+1)*prim(i,j,k+1,5) + zeta(i-1,j,k+1)*prim(i-1,j,k+1,5) + &
                  zeta(i,j+1,k+1)*prim(i,j+1,k+1,5) + zeta(i-1,j+1,k+1)*prim(i-1,j+1,k+1,5) )
             muzemp = muzemp + 0.25d0*(zeta(i,j-1,k)*prim(i,j-1,k,5) + zeta(i-1,j-1,k)*prim(i-1,j-1,k,5) + &
                  zeta(i,j,k)*prim(i,j,k,5) + zeta(i-1,j,k)*prim(i-1,j,k,5) + &
                  zeta(i,j-1,k+1)*prim(i,j-1,k+1,5) + zeta(i-1,j-1,k+1)*prim(i-1,j-1,k+1,5) + &
                  zeta(i,j,k+1)*prim(i,j,k+1,5) + zeta(i-1,j,k+1)*prim(i-1,j,k+1,5) )
             muzepm = muzepm + 0.25d0*(zeta(i,j,k-1)*prim(i,j,k-1,5) + zeta(i-1,j,k-1)*prim(i-1,j,k-1,5) + &
                  zeta(i,j+1,k-1)*prim(i,j+1,k-1,5) + zeta(i-1,j+1,k-1)*prim(i-1,j+1,k-1,5) + &
                  zeta(i,j,k)*prim(i,j,k,5) + zeta(i-1,j,k)*prim(i-1,j,k,5) + &
                  zeta(i,j+1,k)*prim(i,j+1,k,5) + zeta(i-1,j+1,k)*prim(i-1,j+1,k,5) )
             muzemm = muzemm + 0.25d0*(zeta(i,j-1,k-1)*prim(i,j-1,k-1,5) + zeta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,5) + &
                  zeta(i,j,k-1)*prim(i,j,k-1,5) + zeta(i-1,j,k-1)*prim(i-1,j,k-1,5) + &
                  zeta(i,j-1,k)*prim(i,j-1,k,5) + zeta(i-1,j-1,k)*prim(i-1,j-1,k,5) + &
                  zeta(i,j,k)*prim(i,j,k,5) + zeta(i-1,j,k)*prim(i-1,j,k,5) )

          endif

          weiner(2) = weiner(2) + 0.25d0*nweight*(sqrt(muzepp)*rancorn(i,j+1,k+1)+ &
               sqrt(muzemp)*rancorn(i,j,k+1) + sqrt(muzepm)* rancorn(i,j+1,k)+  &
               sqrt(muzemm)*rancorn(i,j,k)) ! Random "divergence" stress

       elseif(n_cells(3).eq.1) then

          ! Corner viscosity coefficients in 2D
          muzepp = 0.5d0*(eta(i,j,k)*prim(i,j,k,5) + eta(i-1,j,k)*prim(i-1,j,k,5) + &
               eta(i,j+1,k)*prim(i,j+1,k,5) + eta(i-1,j+1,k)*prim(i-1,j+1,k,5) )/3.d0
          muzemp = 0.5d0*(eta(i,j-1,k)*prim(i,j-1,k,5) + eta(i-1,j-1,k)*prim(i-1,j-1,k,5) + &
               eta(i,j,k)*prim(i,j,k,5) + eta(i-1,j,k)*prim(i-1,j,k,5) )/3.d0

          if (abs(visc_type) .eq. 3) then

             muzepp = muzepp + 0.25d0*(zeta(i,j,k)*prim(i,j,k,5) + zeta(i-1,j,k)*prim(i-1,j,k,5) + &
                  zeta(i,j+1,k)*prim(i,j+1,k,5) + zeta(i-1,j+1,k)*prim(i-1,j+1,k,5) )
             muzemp = muzemp + 0.25d0*(zeta(i,j-1,k)*prim(i,j-1,k,5) + zeta(i-1,j-1,k)*prim(i-1,j-1,k,5) + &
                  zeta(i,j,k)*prim(i,j,k,5) + zeta(i-1,j,k)*prim(i-1,j,k,5) )

          endif

          weiner(2) = weiner(2) + 0.5d0*nweight*(sqrt(muzepp)*rancorn(i,j+1,k)+ &
               sqrt(muzemp)*rancorn(i,j,k)) ! Random "divergence" stress

       endif

       fluxx(i,j,k,2:5) = fluxx(i,j,k,2:5) + weiner(2:5)

       ! Viscous heating:
       phiflx =  weiner(2)*(prim(i-1,j,k,2)+prim(i,j,k,2)) +           &
            weiner(3)*(prim(i-1,j,k,3)+prim(i,j,k,3)) + &
            weiner(4)*(prim(i-1,j,k,4)+prim(i,j,k,4))
       phiflx =  - 0.5d0*phiflx

       fluxx(i,j,k,5) = fluxx(i,j,k,5) - phiflx

       if(algorithm_type.eq.2) then

          weiner(6:5+nspecies) = 0.0d0

          do ns = 1,nspecies

             yy(ns) = max(0.d0,min(1.d0,prim(i-1,j,k,6+ns)))
             yyp(ns) = max(0.d0,min(1.d0,prim(i,j,k,6+ns)))

          enddo
          sumy = sum(yy(:))
          sumyp = sum(yyp(:))
          yy(:) = yy(:)/sumy
          yyp(:) = yyp(:)/sumyp

          MWmix = 0.d0

          do ns = 1, nspecies

             MWmix = MWmix + 0.5d0*(yy(ns)+yyp(ns))/molmass(ns)

             do ll = 1, nspecies

                DijY_edge(ns,ll) = 0.5d0*(Dij(i-1,j,k,ns,ll)*yy(ll) + &
                     Dij(i,j,k,ns,ll)*yyp(ll) &
                     + (Dij(i-1,j,k,ll,ns)*yy(ns) + &
                     Dij(i,j,k,ll,ns)*yyp(ns) ))

             enddo

          enddo

          do ns=1,nspecies
             if(abs(yy(ns)) + abs(yyp(ns)) .le. 1.d-12)then
                DijY_edge(ns,1:nspecies)=0.d0
                DijY_edge(1:nspecies,ns)=0.d0
             endif
          enddo

          MWmix = 1.d0 / MWmix

          call cholesky_decomp(DijY_edge,nspecies,sqD)

          do ns = 1, nspecies
             do ll = 1, ns
                fweights(5+ll)=sqrt(k_b*MWmix*volinv/(Runiv*dt))*sqD(ns,ll)
                weiner(5+ns) = weiner(5+ns) + fweights(5+ll)*ranfluxx(i,j,k,5+ll)
             enddo
             fluxx(i,j,k,5+ns) = weiner(5+ns)
          enddo

          call get_enthalpies(meanT, hk)

          soret = 0.d0

          do ns = 1, nspecies
             soret = soret + (hk(ns) + Runiv*meanT/molmass(ns) & 
                  *0.5d0*(chi(i-1,j,k,ns)+chi(i,j,k,ns)))*weiner(5+ns)
          enddo
          fluxx(i,j,k,5) = fluxx(i,j,k,5) +  soret

       end if

    end do
    end do
 end do
#endif

!!!!!!!!!!!!!!!!!!! y-flux !!!!!!!!!!!!!!!!!!!
#if 0
    do k = lo(3),hi(3)
    do j = lo(2),hi(2)+1
    do i = lo(1),hi(1)

       muyp = eta(i,j,k)*prim(i,j,k,5) + eta(i,j-1,k)*prim(i,j-1,k,5)
       kyp = kappa(i,j,k)*prim(i,j,k,5)**2 + kappa(i,j-1,k)*prim(i,j-1,k,5)**2

       !! Look into: zeta*p is not used in this function (original FluctHydro)
       ! if (abs(visc_type) .eq. 3) then 
       !    zetayp = (zeta(i,j,k)*prim(i,j,k,5) + zeta(i,j-1,k)*prim(i,j-1,k,5))
       ! else
       !    zetayp = 0.0
       ! endif

       meanT = 0.5d0*(prim(i,j,k,5)+prim(i,j-1,k,5))

       ! Weights for facial fluxes:
       fweights(1)=0 ! No mass flux
       fweights(2:4)=sqrt(k_b*muyp*volinv*dtinv)
       fweights(5)=sqrt(k_b*kyp*volinv*dtinv)

       ! Construct the random increments
       weiner(1:5) = fweights(1:5)*ranfluxy(i,j,k,1:5)

       nweight=sqrt(k_b*volinv*dtinv)

       if(n_cells(3).gt.1) then

          ! Corner viscosity coefficients 3D
          muzepp = 0.25d0*(eta(i+1,j-1,k)*prim(i+1,j-1,k,5) + eta(i,j-1,k)*prim(i,j-1,k,5) + &
               eta(i+1,j,k)*prim(i+1,j,k,5) + eta(i,j,k)*prim(i,j,k,5) + &
               eta(i+1,j-1,k+1)*prim(i+1,j-1,k+1,5) + eta(i,j-1,k+1)*prim(i,j-1,k+1,5) + &
               eta(i+1,j,k+1)*prim(i+1,j,k+1,5) + eta(i,j,k+1)*prim(i,j,k+1,5) )/3.d0

          muzemp = 0.25d0*(eta(i-1,j,k)*prim(i-1,j,k,5) + eta(i,j,k)*prim(i,j,k,5) + &
               eta(i-1,j-1,k)*prim(i-1,j-1,k,5) + eta(i,j-1,k)*prim(i,j-1,k,5) + &
               eta(i-1,j,k+1)*prim(i-1,j,k+1,5) + eta(i,j,k+1)*prim(i,j,k+1,5) + &
               eta(i-1,j-1,k+1)*prim(i-1,j-1,k+1,5) + eta(i,j-1,k+1)*prim(i,j-1,k+1,5) )/3.d0

          muzepm = 0.25d0*(eta(i+1,j-1,k-1)*prim(i+1,j-1,k-1,5) + eta(i,j-1,k-1)*prim(i,j-1,k-1,5) + &
               eta(i+1,j,k-1)*prim(i+1,j,k-1,5) + eta(i,j,k-1)*prim(i,j,k-1,5) + &
               eta(i+1,j-1,k)*prim(i+1,j-1,k,5) + eta(i,j-1,k)*prim(i,j-1,k,5) + &
               eta(i+1,j,k)*prim(i+1,j,k,5) + eta(i,j,k)*prim(i,j,k,5) )/3.d0

          muzemm = 0.25d0*(eta(i-1,j,k-1)*prim(i-1,j,k-1,5) + eta(i,j,k-1)*prim(i,j,k-1,5) + &
               eta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,5) + eta(i,j-1,k-1)*prim(i,j-1,k-1,5) + &
               eta(i-1,j,k)*prim(i-1,j,k,5) + eta(i,j,k)*prim(i,j,k,5) + &
               eta(i-1,j-1,k)*prim(i-1,j-1,k,5) + eta(i,j-1,k)*prim(i,j-1,k,5) )/3.d0

          if (abs(visc_type) .eq. 3) then

             muzepp = muzepp + 0.25d0*(zeta(i+1,j-1,k)*prim(i+1,j-1,k,5) + zeta(i,j-1,k)*prim(i,j-1,k,5) + &
                  zeta(i+1,j,k)*prim(i+1,j,k,5) + zeta(i,j,k)*prim(i,j,k,5) + &
                  zeta(i+1,j-1,k+1)*prim(i+1,j-1,k+1,5) + zeta(i,j-1,k+1)*prim(i,j-1,k+1,5) + &
                  zeta(i+1,j,k+1)*prim(i+1,j,k+1,5) + zeta(i,j,k+1)*prim(i,j,k+1,5) )

             muzemp = muzemp + 0.25d0*(zeta(i-1,j,k)*prim(i-1,j,k,5) + zeta(i,j,k)*prim(i,j,k,5) + &
                  zeta(i-1,j-1,k)*prim(i-1,j-1,k,5) + zeta(i,j-1,k)*prim(i,j-1,k,5) + &
                  zeta(i-1,j,k+1)*prim(i-1,j,k+1,5) + zeta(i,j,k+1)*prim(i,j,k+1,5) + &
                  zeta(i-1,j-1,k+1)*prim(i-1,j-1,k+1,5) + zeta(i,j-1,k+1)*prim(i,j-1,k+1,5) )

             muzepm =  muzepm +0.25d0*(zeta(i+1,j-1,k-1)*prim(i+1,j-1,k-1,5) + zeta(i,j-1,k-1)*prim(i,j-1,k-1,5) + &
                  zeta(i+1,j,k-1)*prim(i+1,j,k-1,5) + zeta(i,j,k-1)*prim(i,j,k-1,5) + &
                  zeta(i+1,j-1,k)*prim(i+1,j-1,k,5) + zeta(i,j-1,k)*prim(i,j-1,k,5) + &
                  zeta(i+1,j,k)*prim(i+1,j,k,5) + zeta(i,j,k)*prim(i,j,k,5) )


             muzemm = muzemm + 0.25d0*(zeta(i-1,j,k-1)*prim(i-1,j,k-1,5) + zeta(i,j,k-1)*prim(i,j,k-1,5) + &
                  zeta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,5) + zeta(i,j-1,k-1)*prim(i,j-1,k-1,5) + &
                  zeta(i-1,j,k)*prim(i-1,j,k,5) + zeta(i,j,k)*prim(i,j,k,5) + &
                  zeta(i-1,j-1,k)*prim(i-1,j-1,k,5) + zeta(i,j-1,k)*prim(i,j-1,k,5) )

          endif

          weiner(3) = weiner(3) + 0.25d0*nweight*    &
               (sqrt(muzepp)*rancorn(i+1,j,k+1)+ sqrt(muzemp)*rancorn(i,j,k+1) +  &
               sqrt(muzepm)* rancorn(i+1,j,k)+ sqrt(muzemm)*rancorn(i,j,k)) ! Random "divergence" stress

       elseif(n_cells(3).eq.1) then

          ! Corner viscosity coefficients 2D
          muzepp = 0.5d0*(eta(i+1,j-1,k)*prim(i+1,j-1,k,5) + eta(i,j-1,k)*prim(i,j-1,k,5) + &
               eta(i+1,j,k)*prim(i+1,j,k,5) + eta(i,j,k)*prim(i,j,k,5) )/3.d0

          muzemp = 0.5d0*(eta(i-1,j,k)*prim(i-1,j,k,5) + eta(i,j,k)*prim(i,j,k,5) + &
               eta(i-1,j-1,k)*prim(i-1,j-1,k,5) + eta(i,j-1,k)*prim(i,j-1,k,5) )/3.d0

          if (abs(visc_type) .eq. 3) then

             muzepp = muzepp + 0.25d0*(zeta(i+1,j-1,k)*prim(i+1,j-1,k,5) + zeta(i,j-1,k)*prim(i,j-1,k,5) + &
                  zeta(i+1,j,k)*prim(i+1,j,k,5) + zeta(i,j,k)*prim(i,j,k,5) )

             muzemp = muzemp + 0.25d0*(zeta(i-1,j,k)*prim(i-1,j,k,5) + zeta(i,j,k)*prim(i,j,k,5) + &
                  zeta(i-1,j-1,k)*prim(i-1,j-1,k,5) + zeta(i,j-1,k)*prim(i,j-1,k,5) )

          endif

          weiner(3) = weiner(3) + 0.5d0*nweight*    &
               (sqrt(muzepp)*rancorn(i+1,j,k) + sqrt(muzemp)*rancorn(i,j,k)) ! Random "divergence" stress

       endif

       fluxy(i,j,k,2:5) = fluxy(i,j,k,2:5) + weiner(2:5)

       ! Viscous heating:
       phiflx =  weiner(2)*(prim(i,j-1,k,2)+prim(i,j,k,2)) +          &
            weiner(3)*(prim(i,j-1,k,3)+prim(i,j,k,3)) + &
            weiner(4)*(prim(i,j-1,k,4)+prim(i,j,k,4))
       phiflx =  - 0.5d0*phiflx

       fluxy(i,j,k,5) = fluxy(i,j,k,5) - phiflx

       if(algorithm_type.eq.2) then

          weiner(6:5+nspecies) = 0.0d0

          do ns = 1,nspecies

             yy(ns) = max(0.d0,min(1.d0,prim(i,j-1,k,6+ns)))
             yyp(ns) = max(0.d0,min(1.d0,prim(i,j,k,6+ns)))

          enddo
          sumy = sum(yy(:))
          sumyp = sum(yyp(:))
          yy(:) = yy(:)/sumy
          yyp(:) = yyp(:)/sumyp

          MWmix = 0.d0

          do ns = 1, nspecies

             MWmix = MWmix + 0.5d0*(yy(ns)+yyp(ns))/molmass(ns)

             do ll = 1, nspecies

                DijY_edge(ns,ll) = 0.5d0*(Dij(i,j-1,k,ns,ll)*yy(ll) + &
                     Dij(i,j,k,ns,ll)*yyp(ll) &
                     + (Dij(i,j-1,k,ll,ns)*yy(ns) + &
                     Dij(i,j,k,ll,ns)*yyp(ns) ))

             enddo

          enddo

          do ns=1,nspecies
             if(abs(yy(ns)) + abs(yyp(ns)) .le. 1.d-12)then
                DijY_edge(ns,1:nspecies)=0.d0
                DijY_edge(1:nspecies,ns)=0.d0
             endif
          enddo

          MWmix = 1.d0 / MWmix

          call cholesky_decomp(DijY_edge,nspecies,sqD)

          do ns = 1, nspecies

             do ll = 1, ns

                fweights(5+ll)=sqrt(k_b*MWmix*volinv/(Runiv*dt))*sqD(ns,ll)
                weiner(5+ns) = weiner(5+ns) + fweights(5+ll)*ranfluxy(i,j,k,5+ll)

             enddo

             fluxy(i,j,k,5+ns) = weiner(5+ns)

          enddo

          call get_enthalpies(meanT, hk)

          soret = 0.d0

          do ns = 1, nspecies
             soret = soret + (hk(ns) + Runiv*meanT/molmass(ns) & 
                  *0.5d0*(chi(i,j-1,k,ns)+chi(i,j,k,ns)))*weiner(5+ns)
          enddo
          fluxy(i,j,k,5) = fluxy(i,j,k,5) +  soret

       end if

    end do
    end do
    end do
#endif
    if(n_cells(3).gt.1) then

!!!!!!!!!!!!!!!!!!! z-flux !!!!!!!!!!!!!!!!!!!

       do k = lo(3),hi(3)+1
       do j = lo(2),hi(2)
       do i = lo(1),hi(1)

          muzp = eta(i,j,k)*prim(i,j,k,5) + eta(i,j,k-1)*prim(i,j,k-1,5)
          kzp = kappa(i,j,k)*prim(i,j,k,5)**2 + kappa(i,j,k-1)*prim(i,j,k-1,5)**2

          !! Look into: zeta*p is not used in this function (original FluctHydro)
          ! if (abs(visc_type) .eq. 3) then 
          !    zetazp = (zeta(i,j,k)*prim(i,j,k,5) + zeta(i,j,k-1)*prim(i,j,k-1,5))
          ! else
          !    zetazp = 0.0
          ! endif

          meanT = 0.5d0*(prim(i,j,k,5)+prim(i,j,k-1,5))

          ! Weights for facial fluxes:
          fweights(1)=0 ! No mass flux
          fweights(2:4)=sqrt(k_b*muzp*volinv*dtinv)
          fweights(5)=sqrt(k_b*kzp*volinv*dtinv)

          ! Construct the random increments
          weiner(1:5) = fweights(1:5)*ranfluxz(i,j,k,1:5)

          nweight=sqrt(k_b*volinv*dtinv)

          ! Corner viscosity coefficients
          muzepp = 0.25d0*(eta(i+1,j,k-1)*prim(i+1,j,k-1,5) + eta(i,j,k-1)*prim(i,j,k-1,5) + &
               eta(i+1,j+1,k-1)*prim(i+1,j+1,k-1,5) + eta(i,j+1,k-1)*prim(i,j+1,k-1,5) + &
               eta(i+1,j,k)*prim(i+1,j,k,5) + eta(i,j,k)*prim(i,j,k,5) + &
               eta(i+1,j+1,k)*prim(i+1,j+1,k,5) + eta(i,j+1,k)*prim(i,j+1,k,5) )/3.d0

          muzemp = 0.25d0*(eta(i-1,j+1,k-1)*prim(i-1,j+1,k-1,5) + eta(i,j+1,k-1)*prim(i,j+1,k-1,5) + &
               eta(i-1,j,k-1)*prim(i-1,j,k-1,5) + eta(i,j,k-1)*prim(i,j,k-1,5) + &
               eta(i-1,j+1,k)*prim(i-1,j+1,k,5) + eta(i,j+1,k)*prim(i,j+1,k,5) + &
               eta(i-1,j,k)*prim(i-1,j,k,5) + eta(i,j,k)*prim(i,j,k,5) )/3.d0

          muzepm = 0.25d0*(eta(i+1,j,k)*prim(i+1,j,k,5) + eta(i,j,k-2)*prim(i,j,k,5) + &
               eta(i+1,j-1,k)*prim(i+1,j-1,k,5) + eta(i,j-1,k-2)*prim(i,j-1,k,5) + &
               eta(i+1,j,k-1)*prim(i+1,j,k-1,5) + eta(i,j,k-1)*prim(i,j,k-1,5) + &
               eta(i+1,j-1,k-1)*prim(i+1,j-1,k-1,5) + eta(i,j-1,k-1)*prim(i,j-1,k-1,5) )/3.d0

          muzemm = 0.25d0*(eta(i-1,j-1,k)*prim(i-1,j-1,k,5) + eta(i,j-1,k)*prim(i,j-1,k,5) + &
               eta(i-1,j,k)*prim(i-1,j,k,5) + eta(i,j,k)*prim(i,j,k,5) + &
               eta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,5) + eta(i,j-1,k-1)*prim(i,j-1,k-1,5) + &
               eta(i-1,j,k-1)*prim(i-1,j,k-1,5) + eta(i,j,k-1)*prim(i,j,k-1,5) )/3.d0

          if (abs(visc_type) .eq. 3) then

             muzepp = muzepp+ 0.25d0*(zeta(i+1,j,k-1)*prim(i+1,j,k-1,5) + zeta(i,j,k-1)*prim(i,j,k-1,5) + &
                  zeta(i+1,j+1,k-1)*prim(i+1,j+1,k-1,5) + zeta(i,j+1,k-1)*prim(i,j+1,k-1,5) + &
                  zeta(i+1,j,k)*prim(i+1,j,k,5) + zeta(i,j,k)*prim(i,j,k,5) + &
                  zeta(i+1,j+1,k)*prim(i+1,j+1,k,5) + zeta(i,j+1,k)*prim(i,j+1,k,5) )

             muzemp = muzemp + 0.25d0*(zeta(i-1,j+1,k-1)*prim(i-1,j+1,k-1,5) + zeta(i,j+1,k-1)*prim(i,j+1,k-1,5) + &
                  zeta(i-1,j,k-1)*prim(i-1,j,k-1,5) + zeta(i,j,k-1)*prim(i,j,k-1,5) + &
                  zeta(i-1,j+1,k)*prim(i-1,j+1,k,5) + zeta(i,j+1,k)*prim(i,j+1,k,5) + &
                  zeta(i-1,j,k)*prim(i-1,j,k,5) + zeta(i,j,k)*prim(i,j,k,5) )

             muzepm = muzepm + 0.25d0*(zeta(i+1,j,k)*prim(i+1,j,k,5) + zeta(i,j,k-2)*prim(i,j,k,5) + &
                  zeta(i+1,j-1,k)*prim(i+1,j-1,k,5) + zeta(i,j-1,k-2)*prim(i,j-1,k,5) + &
                  zeta(i+1,j,k-1)*prim(i+1,j,k-1,5) + zeta(i,j,k-1)*prim(i,j,k-1,5) + &
                  zeta(i+1,j-1,k-1)*prim(i+1,j-1,k-1,5) + zeta(i,j-1,k-1)*prim(i,j-1,k-1,5) )

             muzemm = muzemm + 0.25d0*(zeta(i-1,j-1,k)*prim(i-1,j-1,k,5) + zeta(i,j-1,k)*prim(i,j-1,k,5) + &
                  zeta(i-1,j,k)*prim(i-1,j,k,5) + zeta(i,j,k)*prim(i,j,k,5) + &
                  zeta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,5) + zeta(i,j-1,k-1)*prim(i,j-1,k-1,5) + &
                  zeta(i-1,j,k-1)*prim(i-1,j,k-1,5) + zeta(i,j,k-1)*prim(i,j,k-1,5) )

          endif

          weiner(4) = weiner(4) + 0.25d0*nweight*    &
               (sqrt(muzepp)*rancorn(i+1,j+1,k)+ sqrt(muzemp)*rancorn(i,j+1,k) +  &
               sqrt(muzepm)* rancorn(i+1,j,k)+ sqrt(muzemm)*rancorn(i,j,k)) ! Random "divergence" stress

          fluxz(i,j,k,2:5) = fluxz(i,j,k,2:5) + weiner(2:5)

          ! Viscous heating:
          phiflx =  weiner(2)*(prim(i,j,k-1,2)+prim(i,j,k,2)) +          &
               weiner(3)*(prim(i,j,k-1,3)+prim(i,j,k,3)) + &
               weiner(4)*(prim(i,j,k-1,4)+prim(i,j,k,4))
          phiflx =  - 0.5d0*phiflx

          fluxz(i,j,k,5) = fluxz(i,j,k,5) - phiflx

          if(algorithm_type.eq.2) then

             weiner(6:5+nspecies) = 0.0d0

             do ns = 1,nspecies

                yy(ns) = max(0.d0,min(1.d0,prim(i,j,k-1,6+ns)))
                yyp(ns) = max(0.d0,min(1.d0,prim(i,j,k,6+ns)))

             enddo
             sumy = sum(yy(:))
             sumyp = sum(yyp(:))
             yy(:) = yy(:)/sumy
             yyp(:) = yyp(:)/sumyp

             MWmix = 0.d0

             do ns = 1, nspecies

                MWmix = MWmix + 0.5d0*(yy(ns)+yyp(ns))/molmass(ns)

                do ll = 1, nspecies

                   DijY_edge(ns,ll) = 0.5d0*(Dij(i,j,k-1,ns,ll)*yy(ll) + &
                        Dij(i,j,k,ns,ll)*yyp(ll) &
                        + (Dij(i,j,k-1,ll,ns)*yy(ns) + &
                        Dij(i,j,k,ll,ns)*yyp(ns) ))

                enddo

             enddo

             do ns=1,nspecies
                if(abs(yy(ns)) + abs(yyp(ns)) .le. 1.d-12)then
                   DijY_edge(ns,1:nspecies)=0.d0
                   DijY_edge(1:nspecies,ns)=0.d0
                endif
             enddo

             MWmix = 1.d0 / MWmix

             call cholesky_decomp(DijY_edge,nspecies,sqD)

             do ns = 1, nspecies

                do ll = 1, ns

                   fweights(5+ll)=sqrt(k_b*MWmix*volinv/(Runiv*dt))*sqD(ns,ll)
                   weiner(5+ns) = weiner(5+ns) + fweights(5+ll)*ranfluxz(i,j,k,5+ll)

                enddo

                fluxz(i,j,k,5+ns) = weiner(5+ns)

             enddo

             call get_enthalpies(meanT, hk)

             soret = 0.d0

             do ns = 1, nspecies
                soret = soret + (hk(ns) + Runiv*meanT/molmass(ns) & 
                     *0.5d0*(chi(i,j,k-1,ns)+chi(i,j,k,ns)))*weiner(5+ns)
             enddo
             fluxz(i,j,k,5) = fluxz(i,j,k,5) +  soret

          end if

       end do
       end do
       end do

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
    real(amrex_real), intent(inout) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),nvars)
    real(amrex_real), intent(inout) :: yflux(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),nvars)
#if (AMREX_SPACEDIM == 3)
    real(amrex_real), intent(inout) :: zflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,nvars)
#endif

    real(amrex_real) :: sqrtTwo

    integer :: bc_iter, bc_tmp, indx_lo, indx_hi
    integer :: mass_ind_lo,mass_ind_hi, therm_ind_lo,therm_ind_hi, vel_ind_lo,vel_ind_hi

    integer :: i,j,k,l

    sqrtTwo = sqrt(2.0)
    
    ! Set index ranges

    mass_ind_lo = 3 + AMREX_SPACEDIM
    mass_ind_hi = 2 + AMREX_SPACEDIM + nspecies

    therm_ind_lo = 2 + AMREX_SPACEDIM
    therm_ind_hi = 2 + AMREX_SPACEDIM

    vel_ind_lo = 2
    vel_ind_hi = 1 + AMREX_SPACEDIM

!!!!!!!!

    !! Template:

    ! do bc_iter = 1,3 ! iterate over 1) mass, 2) temperature, and 3) velocity BCs

    !    SELECT CASE (bc_iter)
    !    CASE (1) ! mass boundary conditions
    !       bc_tmp = bc_mass_()
    !       indx_lo = mass_ind_lo
    !       indx_hi = mass_ind_hi
    !    CASE (2) ! temperature boundary conditions
    !       bc_tmp = bc_therm_()
    !       indx_lo = therm_ind_lo
    !       indx_hi = therm_ind_hi
    !    CASE (3) ! velocity boundary conditions
    !       bc_tmp = bc_vel_()
    !       indx_lo = vel_ind_lo
    !       indx_hi = vel_ind_hi
    !    END SELECT

    !    if(bc_tmp .eq. 1) then ! neumann (0 scaling)
    !       do l = indx_lo,indx_hi
    !          if(l.ne.) then ! neumann if not normal velocity
    !          else           ! dirichlet if normal velocity
    !          endif
    !       enddo
    !    elseif(bc_tmp .eq. 2) then ! dirichlet (root 2 scaling) 
    !       do l = indx_lo,indx_hi
    !       enddo
    !    endif

    ! enddo

!!!!!!!!!!!!!! x-flux BCs !!!!!!!!!!!!!!

    !if on lower bound
    if(lo(1) .eq. 0) then !lower x bound
       
       do bc_iter = 1,3 ! iterate over 1) mass, 2) temperature, and 3) velocity BCs
          
          SELECT CASE (bc_iter)
          CASE (1) ! mass boundary conditions
             bc_tmp = bc_mass_lo(1)
             indx_lo = mass_ind_lo
             indx_hi = mass_ind_hi
          CASE (2) ! temperature boundary conditions
             bc_tmp = bc_therm_lo(1)
             indx_lo = therm_ind_lo
             indx_hi = therm_ind_hi
          CASE (3) ! velocity boundary conditions
             bc_tmp = bc_vel_lo(1)
             indx_lo = vel_ind_lo
             indx_hi = vel_ind_hi
          END SELECT

          if(bc_tmp .eq. 1) then ! neumann (0 scaling)

             do l = indx_lo,indx_hi

                if(l.ne.2) then ! neumann if not normal velocity
                   do k = lo(3),hi(3)
                      do j = lo(2),hi(2)
                         xflux(0,j,k,l) = 0
                      end do
                   end do
                else            ! dirichlet if normal velocity
                   do k = lo(3),hi(3)
                      do j = lo(2),hi(2)
                         xflux(0,j,k,l) = sqrtTwo*xflux(0,j,k,l)
                      end do
                   end do
                endif

             end do

          elseif(bc_tmp .eq. 2) then ! dirichlet (root 2 scaling) 

             do l = indx_lo,indx_hi

                do k = lo(3),hi(3)
                   do j = lo(2),hi(2)
                      xflux(0,j,k,l) = sqrtTwo*xflux(0,j,k,l)
                   end do
                end do

             end do

          endif

       enddo

    endif

    !if on upper bound
    if(hi(1) .eq. n_cells(1)-1) then

       do bc_iter = 1,3 ! iterate over 1) mass, 2) temperature, and 3) velocity BCs
          
          SELECT CASE (bc_iter)
          CASE (1) ! mass boundary conditions
             bc_tmp = bc_mass_hi(1)
             indx_lo = mass_ind_lo
             indx_hi = mass_ind_hi
          CASE (2) ! temperature boundary conditions
             bc_tmp = bc_therm_hi(1)
             indx_lo = therm_ind_lo
             indx_hi = therm_ind_hi
          CASE (3) ! velocity boundary conditions
             bc_tmp = bc_vel_hi(1)
             indx_lo = vel_ind_lo
             indx_hi = vel_ind_hi
          END SELECT

          if(bc_tmp .eq. 1) then ! neumann (0 scaling)

             do l = indx_lo,indx_hi

                if(l.ne.2) then ! neumann if not normal velocity
                   do k = lo(3),hi(3)
                      do j = lo(2),hi(2)
                         xflux(hi(1)+1,j,k,l) = 0
                      end do
                   end do
                else           ! dirichlet if normal velocity
                   do k = lo(3),hi(3)
                      do j = lo(2),hi(2)
                         xflux(hi(1)+1,j,k,l) = sqrtTwo*xflux(hi(1)+1,j,k,l)        
                      end do
                   end do
                endif

             enddo

          elseif(bc_tmp .eq. 2) then ! dirichlet (root 2 scaling) 

             do l = indx_lo,indx_hi

                do k = lo(3),hi(3)
                   do j = lo(2),hi(2)
                      xflux(hi(1)+1,j,k,l) = sqrtTwo*xflux(hi(1)+1,j,k,l)        
                   end do
                end do

             enddo

          endif

       enddo

    endif

      !wall cell - hard wired for specular adiabatic for now
      if(lo(1) .eq. membrane_cell) then
        do k = lo(3),hi(3)
          do j = lo(2),hi(2)

              xflux(membrane_cell,j,k,2) = 0
              xflux(membrane_cell,j,k,3) = 0
              xflux(membrane_cell,j,k,4) = 0
              xflux(membrane_cell,j,k,5) = 0

!              xflux(membrane_cell,j,k,2) = 1.4142*xflux(membrane_cell,j,k,2)        
!              xflux(membrane_cell,j,k,5) = 1.4142*xflux(membrane_cell,j,k,5)        

          end do
        end do
      endif

      !wall cell
      if(hi(1) .eq. membrane_cell-1) then
        do k = lo(3),hi(3)
          do j = lo(2),hi(2)

              xflux(membrane_cell,j,k,2) = 0
              xflux(membrane_cell,j,k,3) = 0
              xflux(membrane_cell,j,k,4) = 0
              xflux(membrane_cell,j,k,5) = 0

!              xflux(membrane_cell,j,k,2) = 1.4142*xflux(membrane_cell,j,k,2)        
!              xflux(membrane_cell,j,k,5) = 1.4142*xflux(membrane_cell,j,k,5)        

          end do
        end do
      endif


!!!!!!!!!!!!!! y-flux BCs !!!!!!!!!!!!!!

    !if on lower bound
    if(lo(2) .eq. 0) then
       
       do bc_iter = 1,3 ! iterate over 1) mass, 2) temperature, and 3) velocity BCs

          SELECT CASE (bc_iter)
          CASE (1) ! mass boundary conditions
             bc_tmp = bc_mass_lo(2)
             indx_lo = mass_ind_lo
             indx_hi = mass_ind_hi
          CASE (2) ! temperature boundary conditions
             bc_tmp = bc_therm_lo(2)
             indx_lo = therm_ind_lo
             indx_hi = therm_ind_hi
          CASE (3) ! velocity boundary conditions
             bc_tmp = bc_vel_lo(2)
             indx_lo = vel_ind_lo
             indx_hi = vel_ind_hi
          END SELECT

          if(bc_tmp .eq. 1) then ! neumann (0 scaling)

             do l = indx_lo,indx_hi

                if(l.ne.3) then ! neumann if not normal velocity
                   do k = lo(3),hi(3)
                      do i = lo(1),hi(1)
                         yflux(i,0,k,l) = 0
                      end do
                   end do
                else           ! dirichlet if normal velocity
                   do k = lo(3),hi(3)
                      do i = lo(1),hi(1)
                         yflux(i,0,k,l) = sqrtTwo*yflux(i,0,k,l)
                      end do
                   end do
                endif

             enddo

          elseif(bc_tmp .eq. 2) then ! dirichlet (root 2 scaling) 

             do l = indx_lo,indx_hi

                do k = lo(3),hi(3)
                   do i = lo(1),hi(1)
                      yflux(i,0,k,l) = sqrtTwo*yflux(i,0,k,l)
                   end do
                end do

             enddo

          endif

       enddo

    endif

    !if on upper bound
    if(hi(2) .eq. n_cells(2)-1) then

       do bc_iter = 1,3 ! iterate over 1) mass, 2) temperature, and 3) velocity BCs

          SELECT CASE (bc_iter)
          CASE (1) ! mass boundary conditions
             bc_tmp = bc_mass_hi(2)
             indx_lo = mass_ind_lo
             indx_hi = mass_ind_hi
          CASE (2) ! temperature boundary conditions
             bc_tmp = bc_therm_hi(2)
             indx_lo = therm_ind_lo
             indx_hi = therm_ind_hi
          CASE (3) ! velocity boundary conditions
             bc_tmp = bc_vel_hi(2)
             indx_lo = vel_ind_lo
             indx_hi = vel_ind_hi
          END SELECT

          if(bc_tmp .eq. 1) then ! neumann (0 scaling)

             do l = indx_lo,indx_hi

                if(l.ne.3) then ! neumann if not normal velocity
                   do k = lo(3),hi(3)
                      do i = lo(1),hi(1)
                         yflux(i,hi(2)+1,k,l) = 0
                      end do
                   end do
                else           ! dirichlet if normal velocity
                   do k = lo(3),hi(3)
                      do i = lo(1),hi(1)
                         yflux(i,hi(2)+1,k,l) = sqrtTwo*yflux(i,hi(2)+1,k,l)        
                      end do
                   end do
                endif

             enddo

          elseif(bc_tmp .eq. 2) then ! dirichlet (root 2 scaling) 

             do l = indx_lo,indx_hi

                do k = lo(3),hi(3)
                   do i = lo(1),hi(1)
                      yflux(i,hi(2)+1,k,l) = sqrtTwo*yflux(i,hi(2)+1,k,l)        
                   end do
                end do

             enddo

          endif

       enddo

    endif

!!!!!!!!!!!!!! z-flux BCs !!!!!!!!!!!!!!

    !if on lower bound
    if(lo(3) .eq. 0) then

       do bc_iter = 1,3 ! iterate over 1) mass, 2) temperature, and 3) velocity BCs

          SELECT CASE (bc_iter)
          CASE (1) ! mass boundary conditions
             bc_tmp = bc_mass_lo(3)
             indx_lo = mass_ind_lo
             indx_hi = mass_ind_hi
          CASE (2) ! temperature boundary conditions
             bc_tmp = bc_therm_lo(3)
             indx_lo = therm_ind_lo
             indx_hi = therm_ind_hi
          CASE (3) ! velocity boundary conditions
             bc_tmp = bc_vel_lo(3)
             indx_lo = vel_ind_lo
             indx_hi = vel_ind_hi
          END SELECT

          if(bc_tmp .eq. 1) then ! neumann (0 scaling)

             do l = indx_lo,indx_hi

                if(l.ne.4) then ! neumann if not normal velocity
                   do j = lo(2),hi(2)
                      do i = lo(1),hi(1)
                         zflux(i,j,0,l) = 0
                      end do
                   end do
                else           ! dirichlet if normal velocity
                   do j = lo(2),hi(2)
                      do i = lo(1),hi(1)
                         zflux(i,j,0,l) = sqrtTwo*zflux(i,j,0,l)
                      end do
                   end do
                endif

             enddo

          elseif(bc_tmp .eq. 2) then ! dirichlet (root 2 scaling) 

             do l = indx_lo,indx_hi

                do j = lo(2),hi(2)
                   do i = lo(1),hi(1)
                      zflux(i,j,0,l) = sqrtTwo*zflux(i,j,0,l)
                   end do
                end do

             enddo

          endif

       enddo

    endif

    !if on upper bound
    if(hi(3) .eq. n_cells(3)-1) then

       do bc_iter = 1,3 ! iterate over 1) mass, 2) temperature, and 3) velocity BCs

          SELECT CASE (bc_iter)
          CASE (1) ! mass boundary conditions
             bc_tmp = bc_mass_hi(3)
             indx_lo = mass_ind_lo
             indx_hi = mass_ind_hi
          CASE (2) ! temperature boundary conditions
             bc_tmp = bc_therm_hi(3)
             indx_lo = therm_ind_lo
             indx_hi = therm_ind_hi
          CASE (3) ! velocity boundary conditions
             bc_tmp = bc_vel_hi(3)
             indx_lo = vel_ind_lo
             indx_hi = vel_ind_hi
          END SELECT

          if(bc_tmp .eq. 1) then ! neumann (0 scaling)

             do l = indx_lo,indx_hi

                if(l.ne.4) then ! neumann if not normal velocity
                   do j = lo(2),hi(2)
                      do i = lo(1),hi(1)
                         zflux(i,j,hi(3)+1,l) = 0
                      end do
                   end do
                else           ! dirichlet if normal velocity
                   do j = lo(2),hi(2)
                      do i = lo(1),hi(1)
                         zflux(i,j,hi(3)+1,l) = sqrtTwo*zflux(i,j,hi(3)+1,l)        
                      end do
                   end do
                endif

             enddo

          elseif(bc_tmp .eq. 2) then ! dirichlet (root 2 scaling) 

             do l = indx_lo,indx_hi

                do j = lo(2),hi(2)
                   do i = lo(1),hi(1)
                      zflux(i,j,hi(3)+1,l) = sqrtTwo*zflux(i,j,hi(3)+1,l)        
                   end do
                end do

             enddo

          endif

       enddo

    endif

  end subroutine stoch_flux_BC

end module flux_module

































