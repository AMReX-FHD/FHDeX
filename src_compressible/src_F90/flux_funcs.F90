module flux_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, molmass, cell_depth, k_b, runiv, bc_lo, bc_hi, n_cells, membrane_cell, visc_type
  use conv_module, only : get_temperature, get_pressure_gas, get_energy, get_enthalpies, get_temperature_gas, get_density_gas, get_energy_gas, get_hc_gas
  use multispec_module, only : cholesky_decomp

  implicit none

  private

  public :: diff_flux, stoch_flux_bounds

contains

  subroutine diff_flux(lo,hi, cons, prim, eta, zeta, kappa, xflux, yflux, &
#if (AMREX_SPACEDIM == 3)
                        zflux, &
#endif
                        dx) bind(C,name="diff_flux")

      integer         , intent(in   ) :: lo(3),hi(3)
      real(amrex_real), intent(in   ) :: dx(3)
      real(amrex_real), intent(inout) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
      real(amrex_real), intent(inout) :: yflux(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
      real(amrex_real), intent(inout) :: zflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif

      real(amrex_real), intent(in   ) :: cons(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      real(amrex_real), intent(in   ) :: prim(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)

      real(amrex_real), intent(in   ) :: eta(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
      real(amrex_real), intent(in   ) :: zeta(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
      real(amrex_real), intent(in   ) :: kappa(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))

      integer :: i,j,k 
      real(amrex_real) :: u, v, w, dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz, etaf, kappaf, zetaf, taux, tauy, tauz, dxinv(3), up, down, twothirds, div, dtx, dty, dtz

      dxinv = 1d0/dx
      twothirds = 2d0/3d0

      duz = 0
      dvz = 0
      dwz = 0
      dtz = 0

      !x flux
      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1)-1,hi(1)

            etaF   = 0.5*(eta(i+1,j,k) + eta(i,j,k))
            kappaF = 0.5*(kappa(i+1,j,k) + kappa(i,j,k))
            zetaF  = 0.5*(zeta(i+1,j,k) + zeta(i,j,k))

            ! x derivatives of u, v, and w on x face
            dux = (prim(i+1,j,k,2) - prim(i,j,k,2))*dxinv(1)
            dvx = (prim(i+1,j,k,3) - prim(i,j,k,3))*dxinv(1)
            dwx = (prim(i+1,j,k,4) - prim(i,j,k,4))*dxinv(1)

            ! y derivative of u on x face
            up    = 0.5*(prim(i,j+1,k,2)+prim(i+1,j+1,k,2))
            down  = 0.5*(prim(i,j-1,k,2)+prim(i+1,j-1,k,2))
            duy   = 0.5*(up - down)*dxinv(2)

#if (AMREX_SPACEDIM == 3)
            ! z derivative of u on x face
            up    = 0.5*(prim(i,j,k+1,2)+prim(i+1,j,k+1,2))
            down  = 0.5*(prim(i,j,k-1,2)+prim(i+1,j,k-1,2))
            duz   = 0.5*(up - down)*dxinv(3)
#endif
            ! y derivative of v on x face
            up    = 0.5*(prim(i,j+1,k,3)+prim(i+1,j+1,k,3))
            down  = 0.5*(prim(i,j-1,k,3)+prim(i+1,j-1,k,3))
            dvy   = 0.5*(up - down)*dxinv(2)
#if (AMREX_SPACEDIM == 3)
            ! z derivative of w on x face
            up    = 0.5*(prim(i,j,k+1,4)+prim(i+1,j,k+1,4))
            down  = 0.5*(prim(i,j,k-1,4)+prim(i+1,j,k-1,4))
            dwz   = 0.5*(up - down)*dxinv(3)
#endif
            !div = dux + dvy + dwz
            div = dux !1D
            
            taux = etaf*(2d0*dux - twothirds*div)
            if (abs(visc_type) .eq. 3) then
               taux = taux + zetaf*div
            endif
            tauy = etaf*(duy + dvx)
            tauz = etaf*(duz + dwx)

            !momentum fluxes
            xflux(i+1,j,k,2) = xflux(i+1,j,k,2) - taux
            xflux(i+1,j,k,3) = xflux(i+1,j,k,3) - tauy
            xflux(i+1,j,k,4) = xflux(i+1,j,k,4) - tauz
    
            !print *,  taux
            
            ! x derivative of T on x face
            dtx = (prim(i+1,j,k,5) - prim(i,j,k,5))*dxinv(1)

            ! y derivative of T on x face
            up    = 0.5*(prim(i,j+1,k,5)+prim(i+1,j+1,k,5))
            down  = 0.5*(prim(i,j-1,k,5)+prim(i+1,j-1,k,5))
            dty   = 0.5*(up - down)*dxinv(2)
#if (AMREX_SPACEDIM == 3)
            ! z derivative of T on x face
            up    = 0.5*(prim(i,j,k+1,5)+prim(i+1,j,k+1,5))
            down  = 0.5*(prim(i,j,k-1,5)+prim(i+1,j,k-1,5))
            dtz = 0.5*(up - down)*dxinv(3)
#endif
            !u, v, and w on x face
            u = 0.5*(prim(i+1,j,k,2) + prim(i,j,k,2))
            v = 0.5*(prim(i+1,j,k,3) + prim(i,j,k,3))
            w = 0.5*(prim(i+1,j,k,4) + prim(i,j,k,4))            

            !print *, "tau: ", taux, tauy, tauz
            !print *, "vel: ", u, v, w
            !energy flux
            xflux(i+1,j,k,5) = xflux(i+1,j,k,5) - (u*taux + v*tauy + w*tauz) - kappaf*(dtx + dty + dtz)
            
          end do
        end do
      end do

     !y flux
     do k = lo(3),hi(3)
       do j = lo(2)-1,hi(2)
         do i = lo(1),hi(1)

           etaF = 0.5*(eta(i,j+1,k) + eta(i,j,k))
           kappaF = 0.5*(kappa(i,j+1,k) + kappa(i,j,k))
           zetaF = 0.5*(zeta(i,j+1,k) + zeta(i,j,k))

           ! y derivatives of u, v, and w on y face
           duy = (prim(i,j+1,k,2) - prim(i,j,k,2))*dxinv(2)
           dvy = (prim(i,j+1,k,3) - prim(i,j,k,3))*dxinv(2)
           dwy = (prim(i,j+1,k,4) - prim(i,j,k,4))*dxinv(2)

           ! x derivative of v on y face
           up    = 0.5*(prim(i+1,j+1,k,3)+prim(i+1,j,k,3))
           down  = 0.5*(prim(i-1,j+1,k,3)+prim(i-1,j,k,3))
           dvx = 0.5*(up - down)*dxinv(1)

#if (AMREX_SPACEDIM == 3)
           ! z derivative of v on y face
           up    = 0.5*(prim(i,j+1,k+1,3)+prim(i,j,k+1,3))
           down  = 0.5*(prim(i,j+1,k-1,3)+prim(i,j,k-1,3))
           dvz = 0.5*(up - down)*dxinv(3)
#endif
           ! x derivative of u on y face
           up    = 0.5*(prim(i+1,j+1,k,2)+prim(i+1,j,k,2))
           down  = 0.5*(prim(i-1,j+1,k,2)+prim(i-1,j,k,2))
           dux = 0.5*(up - down)*dxinv(1)
#if (AMREX_SPACEDIM == 3)
           ! z derivative of w on y face
           up    = 0.5*(prim(i,j+1,k+1,4)+prim(i,j,k+1,4))
           down  = 0.5*(prim(i,j+1,k-1,4)+prim(i,j,k-1,4))
           dwz = 0.5*(up - down)*dxinv(3)
#endif
           div = dux + dvy + dwz

           taux = etaf*(duy +  dvx)
           tauy = etaf*(2d0*dvy - twothirds*div)                  
           if (abs(visc_type) .eq. 3) then
              tauy = tauy + zetaf*div                  
           endif
           tauz = etaf*(dvz + dwy)

           !momentum fluxes
           yflux(i,j+1,k,2) = yflux(i,j+1,k,2) - taux
           yflux(i,j+1,k,3) = yflux(i,j+1,k,3) - tauy
           yflux(i,j+1,k,4) = yflux(i,j+1,k,4) - tauz
           
           ! y derivative of T on y face
           dty = (prim(i,j+1,k,5) - prim(i,j,k,5))*dxinv(2)

           ! x derivative of T on y face
           up    = 0.5*(prim(i+1,j+1,k,5)+prim(i+1,j,k,5))
           down  = 0.5*(prim(i-1,j+1,k,5)+prim(i-1,j,k,5))
           dtx = 0.5*(up - down)*dxinv(1)
#if (AMREX_SPACEDIM == 3)
           ! z derivative of T on y face
           up    = 0.5*(prim(i,j+1,k+1,5)+prim(i,j,k+1,5))
           down  = 0.5*(prim(i,j+1,k-1,5)+prim(i,j,k-1,5))
           dtz = 0.5*(up - down)*dxinv(3)
#endif
           !u, v, and w on y face
           u = 0.5*(prim(i,j+1,k,2) - prim(i,j,k,2))
           v = 0.5*(prim(i,j+1,k,3) - prim(i,j,k,3))
           w = 0.5*(prim(i,j+1,k,4) - prim(i,j,k,4))            

           !energy flux
           yflux(i,j+1,k,5) = yflux(i,j+1,k,5) - (u*taux + v*tauy + w*tauz) - kappaf*(dtx + dty + dtz)

         end do
       end do
     end do

#if (AMREX_SPACEDIM == 3)
     !z flux
     do k = lo(3)-1,hi(3)
       do j = lo(2),hi(2)
         do i = lo(1),hi(1)

           etaF   = 0.5*(eta(i,j,k+1) + eta(i,j,k))
           kappaF = 0.5*(kappa(i,j,k+1) + kappa(i,j,k))
           zetaF  = 0.5*(zeta(i,j,k+1) + zeta(i,j,k))

           ! z derivatives of u, v, and w on z face
           duz = (prim(i,j,k+1,2) - prim(i,j,k,2))*dxinv(3)
           dvz = (prim(i,j,k+1,3) - prim(i,j,k,3))*dxinv(3)
           dwz = (prim(i,j,k+1,4) - prim(i,j,k,4))*dxinv(3)

           ! x derivative of w on z face
           up    = 0.5*(prim(i+1,j,k+1,4)+prim(i+1,j,k,4))
           down  = 0.5*(prim(i-1,j,k+1,4)+prim(i-1,j,k,4))
           dwx = 0.5*(up - down)*dxinv(1)

           ! y derivative of w on z face
           up    = 0.5*(prim(i,j+1,k,4)+prim(i,j+1,k+1,4))
           down  = 0.5*(prim(i,j-1,k,4)+prim(i,j-1,k+1,4))
           dwy   = 0.5*(up - down)*dxinv(2)

           ! x derivative of u on z face
           up    = 0.5*(prim(i+1,j,k+1,2)+prim(i+1,j,k,2))
           down  = 0.5*(prim(i-1,j,k+1,2)+prim(i-1,j,k,2))
           dux   = 0.5*(up - down)*dxinv(1)

           ! y derivative of v on z face
           up    = 0.5*(prim(i,j+1,k+1,3)+prim(i,j+1,k,3))
           down  = 0.5*(prim(i,j-1,k+1,3)+prim(i,j-1,k,3))
           dvy   = 0.5*(up - down)*dxinv(2)

           div = dux + dvy + dwz

           taux = etaf*(duz + dwx)
           tauy = etaf*(dvz + dwy)
           tauz = etaf*(2d0*dwz - twothirds*div)
           if (abs(visc_type) .eq. 3) then
              tauz = tauz + zetaf*div
           endif

           !momentum fluxes
           zflux(i,j,k+1,2) = zflux(i,j,k+1,2) - taux
           zflux(i,j,k+1,3) = zflux(i,j,k+1,3) - tauy
           zflux(i,j,k+1,4) = zflux(i,j,k+1,4) - tauz
           
           ! z derivative of T on z face
           dtz = (prim(i,j,k+1,5) - prim(i,j,k,5))*dxinv(3)

           ! y derivative of T on z face
           up    = 0.5*(prim(i,j+1,k,5)+prim(i,j+1,k+1,5))
           down  = 0.5*(prim(i,j-1,k,5)+prim(i,j-1,k+1,5))
           dty   = 0.5*(up - down)*dxinv(2)

           ! x derivative of T on z face
           up    = 0.5*(prim(i+1,j,k+1,5)+prim(i+1,j,k,5))
           down  = 0.5*(prim(i-1,j,k+1,5)+prim(i-1,j,k,5))
           dtx = 0.5*(up - down)*dxinv(1)

           !u, v, and w on z face
           u = 0.5*(prim(i,j,k+1,2) + prim(i,j,k,2))
           v = 0.5*(prim(i,j,k+1,3) + prim(i,j,k,3))
           w = 0.5*(prim(i,j,k+1,4) + prim(i,j,k,4))            

           !energy flux
           zflux(i,j,k+1,5) = zflux(i,j,k+1,5) - (u*taux + v*tauy + w*tauz) - kappaf*(dtx + dty + dtz)

         end do
       end do
     end do
#endif

  end subroutine diff_flux
  
  subroutine hyp_flux(lo,hi, cons, prim, xflux, yflux, &
#if (AMREX_SPACEDIM == 3)
                        zflux, &
#endif
                        dx) bind(C,name="hyp_flux")

      integer         , intent(in   ) :: lo(3),hi(3)
      real(amrex_real), intent(in   ) :: dx(3)
      real(amrex_real), intent(inout) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
      real(amrex_real), intent(inout) :: yflux(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
      real(amrex_real), intent(inout) :: zflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif

      real(amrex_real), intent(in   ) :: cons(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      real(amrex_real), intent(in   ) :: prim(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)

      real(amrex_real) :: conserved(nvars), primitive(nprimvars), wgt1, wgt2, vsqr, intenergy, specden(nspecies), Yk(nspecies), rho, temp, pt

      integer :: i,j,k,l,n

      wgt2 = 1.0/12.0 !fourth order interpolation
      wgt1 = 0.5 + wgt2 

      !wgt2 = 0
      !wgt1 = 0.5 + wgt2 !second order

      !wgt2 = (sqrt(7d0)-1d0)/4d0 !adjusted for correct variance fourth order interpolation - this apparently makes the overall spectrum worse
      !wgt1 = (sqrt(7d0)+1d0)/4d0


      !Interpolating conserved quantaties for conv term, apparently this has some advantge over interpolating primitives

      !x flux
      
      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1


            do l = 1,nprimvars             
               primitive(l) = wgt1*(prim(i,j,k,l)+prim(i-1,j,k,l)) -wgt2*(prim(i-2,j,k,l)+prim(i+1,j,k,l))  
            enddo

            temp = primitive(5)
            pt = primitive(6)
            rho = primitive(1)
            
            ! call get_density_gas(pt,rho, temp)
            ! conserved(1) = rho
            
            !  want sum of specden == rho
            do n=1,nspecies
               specden(n) = wgt1*(cons(i,j,k,5+n)+cons(i-1,j,k,5+n))                 &
                           -wgt2*(cons(i-2,j,k,5+n)+cons(i+1,j,k,5+n))

               Yk(n) = specden(n)/rho

            enddo

            call get_energy(intenergy, Yk, temp)
            ! call get_energy_gas(pt, intenergy)

            vsqr = primitive(2)**2 + primitive(3)**2 + primitive(4)**2

            conserved(5) = rho*intenergy + 0.5*rho*vsqr

            xflux(i,j,k,1) = xflux(i,j,k,1) + conserved(1)*primitive(2)
            xflux(i,j,k,2) = xflux(i,j,k,2) + conserved(1)*(primitive(2)**2)+primitive(6)
            xflux(i,j,k,3) = xflux(i,j,k,3) + conserved(1)*primitive(2)*primitive(3)
            xflux(i,j,k,4) = xflux(i,j,k,4) + conserved(1)*primitive(2)*primitive(4)

            ! print*, "Hack (hyp_flux): flux = ", xflux(i,j,k,5)
            ! stop

            xflux(i,j,k,5) = xflux(i,j,k,5) + primitive(2)*conserved(5) + primitive(6)*primitive(2)
 
            do n=1,nspecies
               xflux(i,j,k,5+n) = xflux(i,j,k,5+n) + specden(n)*primitive(2)
            enddo

          end do
        end do
      end do


     !y flux
    
     do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1
         do i = lo(1),hi(1)

           do l = 1,nprimvars 
             primitive(l) = wgt1*(prim(i,j,k,l)+prim(i,j-1,k,l)) -wgt2*(prim(i,j-2,k,l)+prim(i,j+1,k,l))
           enddo

           temp = primitive(5)
           pt = primitive(6)
           rho = primitive(1)

           ! call get_density_gas(pt,rho, temp)
           ! conserved(1) = rho
           
           !  want sum of specden == rho
           do n=1,nspecies
              specden(n) = wgt1*(cons(i,j,k,5+n)+cons(i,j-1,k,5+n))                 &
                          -wgt2*(cons(i,j-2,k,5+n)+cons(i,j+1,k,5+n))

              Yk(n) = specden(n)/rho

           enddo
           
           call get_energy(intenergy, Yk, temp)
           ! call get_energy_gas(pt, intenergy)

           vsqr = primitive(2)**2 + primitive(3)**2 + primitive(4)**2

           conserved(5) = intenergy + 0.5*rho*vsqr

           yflux(i,j,k,1) = yflux(i,j,k,1) + conserved(1)*primitive(3)
           yflux(i,j,k,2) = yflux(i,j,k,2) + conserved(1)*primitive(2)*primitive(3)
           yflux(i,j,k,3) = yflux(i,j,k,3) + conserved(1)*primitive(3)**2+primitive(6)
           yflux(i,j,k,4) = yflux(i,j,k,4) + conserved(1)*primitive(4)*primitive(3)

           yflux(i,j,k,5) = yflux(i,j,k,5) + primitive(3)*conserved(5) + primitive(6)*primitive(3)
           
           do n=1,nspecies
              yflux(i,j,k,5+n) = yflux(i,j,k,5+n) + specden(n)*primitive(3)
           enddo

         end do
       end do
     end do

     !z flux
     
     do k = lo(3),hi(3)+1
       do j = lo(2),hi(2)
         do i = lo(1),hi(1)

           do l = 1,nprimvars 
             primitive(l) = wgt1*(prim(i,j,k,l)+prim(i,j,k-1,l)) -wgt2*(prim(i,j,k-2,l)+prim(i,j,k+1,l))
           enddo

           temp = primitive(5)
           pt = primitive(6)
           rho = primitive(1)

           ! call get_density_gas(pt,rho, temp)
           ! conserved(1) = rho
           
           !  want sum of specden == rho
           do n=1,nspecies
              specden(n) = wgt1*(cons(i,j,k,5+n)+cons(i,j,k-1,5+n))                 &
                          -wgt2*(cons(i,j,k-2,5+n)+cons(i,j,k+1,5+n))

              Yk(n) = specden(n)/rho

           enddo
           
           call get_energy(intenergy, Yk, temp)
           ! call get_energy_gas(pt, intenergy)

           vsqr = primitive(2)**2 + primitive(3)**2 + primitive(4)**2

           conserved(5) = intenergy + 0.5*rho*vsqr

           zflux(i,j,k,1) = zflux(i,j,k,1) + conserved(1)*primitive(4)
           zflux(i,j,k,2) = zflux(i,j,k,2) + conserved(1)*primitive(2)*primitive(4)
           zflux(i,j,k,3) = zflux(i,j,k,3) + conserved(1)*primitive(3)*primitive(4)
           zflux(i,j,k,4) = zflux(i,j,k,4) + conserved(1)*primitive(4)**2+primitive(6)
           zflux(i,j,k,5) = zflux(i,j,k,5) + primitive(4)*conserved(5) + primitive(6)*primitive(4)
           
           do n=1,nspecies
              zflux(i,j,k,5+n) = zflux(i,j,k,5+n) + specden(n)*primitive(4)
           enddo

         end do
       end do
     end do


  end subroutine hyp_flux

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

      real(amrex_real), intent(in   ) :: cons(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      real(amrex_real), intent(in   ) :: prim(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)
  
      real(amrex_real), intent(in   ) :: rancorn(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)

      real(amrex_real), intent(in   ) :: eta(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
      real(amrex_real), intent(in   ) :: zeta(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
      real(amrex_real), intent(in   ) :: kappa(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))

      real(amrex_real), intent(in   ) :: chi(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nspecies)
      real(amrex_real), intent(in   ) :: Dij(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nspecies,nspecies)

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
                  
                  fluxx(i,j,k,2:5) = fluxx(i,j,k,2:5) + weiner(2:5)

                  ! Viscous heating:
                  phiflx =  weiner(2)*(prim(i-1,j,k,2)+prim(i,j,k,2)) +           &
                       weiner(3)*(prim(i-1,j,k,3)+prim(i,j,k,3)) + &
                       weiner(4)*(prim(i-1,j,k,4)+prim(i,j,k,4))
                  phiflx =  - 0.5d0*phiflx

                  fluxx(i,j,k,5) = fluxx(i,j,k,5) - phiflx

                  print*, "Hack 1, (stochflux) in x = ", fluxx(i,j,k,5), fweights(5), kxp
                  stop
                  
                  ! if(.not.single_component) then

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

                           ! if ((i.eq.0).and.(j.eq.0).and.(k.eq.0)) then
                           !    print*, "Hack = ", ll, ns, MWmix, volinv, Runiv, dt, sqD(ns,ll)
                           ! endif

                        enddo

                        ! if(lbc(1) == -1)then

                        !    if(i == lb(1)-1)then
                        !       weiner(5+ns)=0.d0
                        !    endif

                        ! elseif (lbc(1) == -2)then

                        !    if(i == lb(1)-1)then
                        !       weiner(5+ns)=weiner(5+ns)*sqrt(2.d0)
                        !    endif

                        ! endif

                        ! if(ubc(1) == -1)then

                        !    if(i == ub(1))then
                        !       weiner(5+ns)=0.d0
                        !    endif

                        ! elseif (ubc(1) == -2)then

                        !    if(i == ub(1))then
                        !       weiner(5+ns)=weiner(5+ns)*sqrt(2.d0)
                        !    endif

                        ! endif


                        fluxx(i,j,k,5+ns) = weiner(5+ns)

                        ! if ((i.eq.0).and.(j.eq.0).and.(k.eq.0)) then
                        !    print*, "Hack, massflux in x = ", weiner(5+ns), ns
                        ! endif

                     enddo

                     call get_enthalpies(meanT, hk)

                     soret = 0.d0

                     do ns = 1, nspecies
                        soret = soret + (hk(ns) + Runiv*meanT/molmass(ns) & 
                             *0.5d0*(chi(i-1,j,k,ns)+chi(i,j,k,ns)))*weiner(5+ns)
                     enddo
                     fluxx(i,j,k,5) = fluxx(i,j,k,5) +  soret

                     ! print*, "Hack 2, (stochflux) in x = ", fluxx(i,j,k,5), soret
                     ! stop

                  ! end if

               end do
            end do
         end do

         !!!!!!!!!!!!!!!!!!! y-flux !!!!!!!!!!!!!!!!!!!

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
                  
                  fluxy(i,j,k,2:5) = fluxy(i,j,k,2:5) + weiner(2:5)

                  ! Viscous heating:
                  phiflx =  weiner(2)*(prim(i,j-1,k,2)+prim(i,j,k,2)) +          &
                       weiner(3)*(prim(i,j-1,k,3)+prim(i,j,k,3)) + &
                       weiner(4)*(prim(i,j-1,k,4)+prim(i,j,k,4))
                  phiflx =  - 0.5d0*phiflx

                  fluxy(i,j,k,5) = fluxy(i,j,k,5) - phiflx

                  ! if(.not.single_component) then

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

                        ! if(lbc(2) == -1)then

                        !    if(i == lb(2)-1)then
                        !       weiner(5+ns)=0.d0
                        !    endif

                        ! elseif (lbc(2) == -2)then

                        !    if(i == lb(2)-1)then
                        !       weiner(5+ns)=weiner(5+ns)*sqrt(2.d0)
                        !    endif

                        ! endif

                        ! if(ubc(2) == -1)then

                        !    if(i == ub(2))then
                        !       weiner(5+ns)=0.d0
                        !    endif

                        ! elseif (ubc(2) == -2)then

                        !    if(i == ub(2))then
                        !       weiner(5+ns)=weiner(5+ns)*sqrt(2.d0)
                        !    endif

                        ! endif


                        fluxy(i,j,k,5+ns) = weiner(5+ns)

                     enddo

                     call get_enthalpies(meanT, hk)

                     soret = 0.d0

                     do ns = 1, nspecies
                        soret = soret + (hk(ns) + Runiv*meanT/molmass(ns) & 
                             *0.5d0*(chi(i,j-1,k,ns)+chi(i,j,k,ns)))*weiner(5+ns)
                     enddo
                     fluxy(i,j,k,5) = fluxy(i,j,k,5) +  soret

                  ! end if

               end do
            end do
         end do

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

                 ! if(.not.single_component) then

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

                        ! if(lbc(2) == -1)then

                        !    if(i == lb(2)-1)then
                        !       weiner(5+ns)=0.d0
                        !    endif

                        ! elseif (lbc(2) == -2)then

                        !    if(i == lb(2)-1)then
                        !       weiner(5+ns)=weiner(5+ns)*sqrt(2.d0)
                        !    endif

                        ! endif

                        ! if(ubc(2) == -1)then

                        !    if(i == ub(2))then
                        !       weiner(5+ns)=0.d0
                        !    endif

                        ! elseif (ubc(2) == -2)then

                        !    if(i == ub(2))then
                        !       weiner(5+ns)=weiner(5+ns)*sqrt(2.d0)
                        !    endif

                        ! endif


                        fluxz(i,j,k,5+ns) = weiner(5+ns)

                     enddo

                     call get_enthalpies(meanT, hk)

                     soret = 0.d0

                     do ns = 1, nspecies
                        soret = soret + (hk(ns) + Runiv*meanT/molmass(ns) & 
                             *0.5d0*(chi(i,j,k-1,ns)+chi(i,j,k,ns)))*weiner(5+ns)
                     enddo
                     fluxz(i,j,k,5) = fluxz(i,j,k,5) +  soret

                  ! end if

               end do
            end do
         end do

      else

         ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!! DL's tensor form !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! !!!!!!!!!!!!!!!!!!! x-flux !!!!!!!!!!!!!!!!!!!

         ! do k = lo(3),hi(3)
         !    do j = lo(2),hi(2)
         !       do i = lo(1)-1,hi(1)

         !          kappattF = (kappa(i,j,k)*prim(i,j,k,5)*prim(i,j,k,5)+kappa(i+1,j,k)*prim(i+1,j,k,5)*prim(i+1,j,k,5))
         !          etatF = (eta(i,j,k)*prim(i,j,k,5)+eta(i+1,j,k)*prim(i+1,j,k,5))
         !          velu = 0.5*(prim(i,j,k,2)+prim(i+1,j,k,2))
         !          velv = 0.5*(prim(i,j,k,3)+prim(i+1,j,k,3))
         !          velw = 0.5*(prim(i,j,k,4)+prim(i+1,j,k,4))

         !          do l = 2,4
         !             xflux(i+1,j,k,l) = xflux(i+1,j,k,l) + sqrt(sFac*etatF)*xsflux(i+1,j,k,l)
         !          end do
         !          xflux(i+1,j,k,5) = xflux(i+1,j,k,5) + sqrt(qFac*kappattF)*xsflux(i+1,j,k,5) + velu*sqrt(sFac*etatF)*xsflux(i+1,j,k,2) + velv*sqrt(sFac*etatF)*xsflux(i+1,j,k,3) + velw*sqrt(sFac*etatF)*xsflux(i+1,j,k,4)

         !       end do
         !    end do
         ! end do

         ! !!!!!!!!!!!!!!!!!!! y-flux !!!!!!!!!!!!!!!!!!!

         ! do k = lo(3),hi(3)
         !    do j = lo(2)-1,hi(2)
         !       do i = lo(1),hi(1)

         !          kappattF = (kappa(i,j,k)*prim(i,j,k,5)*prim(i,j,k,5)+kappa(i,j+1,k)*prim(i,j+1,k,5)*prim(i,j+1,k,5))
         !          etatF = (eta(i,j,k)*prim(i,j,k,5)+eta(i,j+1,k)*prim(i,j+1,k,5))
         !          velu = 0.5*(prim(i,j,k,2)+prim(i,j+1,k,2))            
         !          velv = 0.5*(prim(i,j,k,3)+prim(i,j+1,k,3))
         !          velw = 0.5*(prim(i,j,k,4)+prim(i,j+1,k,4))

         !          do l = 2,4
         !             yflux(i,j+1,k,l) = yflux(i,j+1,k,l) + sqrt(sFac*etatF)*ysflux(i,j+1,k,l)
         !          end do
         !          yflux(i,j+1,k,5) = yflux(i,j+1,k,5) + sqrt(qFac*kappattF)*ysflux(i,j+1,k,5) + velu*sqrt(sFac*etatF)*ysflux(i,j+1,k,2) + velv*sqrt(sFac*etatF)*ysflux(i,j+1,k,3) + velw*sqrt(sFac*etatF)*ysflux(i,j+1,k,4)

         !       end do
         !    end do
         ! end do

         ! !!!!!!!!!!!!!!!!!!! z-flux !!!!!!!!!!!!!!!!!!!

         ! do k = lo(3)-1,hi(3)
         !    do j = lo(2),hi(2)
         !       do i = lo(1),hi(1)

         !          kappattF = (kappa(i,j,k)*prim(i,j,k,5)*prim(i,j,k,5)+kappa(i,j,k+1)*prim(i,j,k+1,5)*prim(i,j,k+1,5))
         !          etatF = (eta(i,j,k)*prim(i,j,k,5)+eta(i,j,k+1)*prim(i,j,k+1,5))
         !          velu = 0.5*(prim(i,j,k,2)+prim(i,j,k+1,2))
         !          velv = 0.5*(prim(i,j,k,3)+prim(i,j,k+1,3))
         !          velw = 0.5*(prim(i,j,k,4)+prim(i,j,k+1,4))

         !          do l = 2,4
         !             zflux(i,j,k+1,l) = zflux(i,j,k+1,l) + sqrt(sFac*etatF)*zsflux(i,j,k+1,l)
         !          end do
         !          zflux(i,j,k+1,5) = zflux(i,j,k+1,5) + sqrt(qFac*kappattF)*zsflux(i,j,k+1,5) + velu*sqrt(sFac*etatF)*zsflux(i,j,k+1,2) + velv*sqrt(sFac*etatF)*zsflux(i,j,k+1,3) + velw*sqrt(sFac*etatF)*zsflux(i,j,k+1,4)

         !       end do
         !    end do
         ! end do

      endif

      !!!!!!!!!!!!! Enforce flux boundary conditions !!!!!!!!!!!!!
      
      call stoch_flux_bounds(lo,hi, fluxx, fluxy, &
#if (AMREX_SPACEDIM == 3)
                        fluxz &
#endif
                        )
      
  end subroutine stoch_flux

  subroutine stoch_flux_bounds(lo,hi, xflux, yflux, &
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
      
      !!!!!!!!!!!!!! x-flux BCs !!!!!!!!!!!!!!
      
      !!!!!!! GM: DL's work is insane in my membrane. Commenting for now... !!!!!!!
      
      !       !wall cell - hard wired for specular adiabatic for now
      !       if(lo(1) .eq. membrane_cell) then
      !         do k = lo(3),hi(3)
      !           do j = lo(2),hi(2)
      
      !               xflux(membrane_cell,j,k,2) = 0        
      !               xflux(membrane_cell,j,k,3) = 0 
      !               xflux(membrane_cell,j,k,4) = 0 
      !               xflux(membrane_cell,j,k,5) = 0
      
      ! !              xflux(membrane_cell,j,k,2) = 1.4142*xflux(membrane_cell,j,k,2)        
      ! !              xflux(membrane_cell,j,k,5) = 1.4142*xflux(membrane_cell,j,k,5)        
      
      !           end do
      !         end do
      !       endif
      
      !       !wall cell
      !       if(hi(1) .eq. membrane_cell-1) then
      !         do k = lo(3),hi(3)
      !           do j = lo(2),hi(2)
      
      !               xflux(membrane_cell,j,k,2) = 0  
      !               xflux(membrane_cell,j,k,3) = 0  
      !               xflux(membrane_cell,j,k,4) = 0  
      !               xflux(membrane_cell,j,k,5) = 0
      
      ! !              xflux(membrane_cell,j,k,2) = 1.4142*xflux(membrane_cell,j,k,2)        
      ! !              xflux(membrane_cell,j,k,5) = 1.4142*xflux(membrane_cell,j,k,5)        
      
      !           end do
      !         end do
      !       endif
      
      !if on lower bound and specular
      if((lo(1) .eq. 0) .and. (bc_lo(1) .eq. 1)) then
         do l = 2,5
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)

                  xflux(0,j,k,l) = 0        

               end do
            end do
         end do
      endif

      !if on upper bound and specular
      if((hi(1) .eq. n_cells(1)-1) .and. (bc_hi(1) .eq. 1)) then
         do l = 2,5
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)

                  xflux(hi(1)+1,j,k,l) = 0           

               end do
            end do
         end do
      endif

      !if on lower bound and diff
      if((lo(1) .eq. 0) .and. (bc_lo(1) .eq. 2)) then
         do l = 2,5
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)

                  xflux(0,j,k,l) = 1.4142*xflux(0,j,k,2)

               end do
            end do
         end do
      endif

      !if on upper bound and diff
      if((hi(1) .eq. n_cells(1)-1) .and. (bc_hi(1) .eq. 2)) then
         do l = 2,5
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)

                  xflux(hi(1)+1,j,k,l) = 1.4142*xflux(hi(1)+1,j,k,l)        

               end do
            end do
         end do
      endif

      !!!!!!!!!!!!!! y-flux BCs !!!!!!!!!!!!!!

      !if on lower bound and specular
      if((lo(2) .eq. 0) .and. (bc_lo(2) .eq. 1)) then
         do l = 2,5
            do k = lo(3),hi(3)
               do i = lo(1),hi(1)

                  yflux(i,0,k,l) = 0        

               end do
            end do
         end do
      endif

      !if on upper bound and specular
      if((hi(2) .eq. n_cells(2)-1) .and. (bc_hi(2) .eq. 1)) then
         do l = 2,5
            do k = lo(3),hi(3)
               do i = lo(1),hi(1)

                  yflux(i,hi(2)+1,k,l) = 0           

               end do
            end do
         end do
      endif

      !if on lower bound and diff
      if((lo(2) .eq. 0) .and. (bc_lo(2) .eq. 2)) then
         do l = 2,5
            do k = lo(3),hi(3)
               do i = lo(1),hi(1)

                  yflux(i,0,k,l) = 1.4142*yflux(i,0,k,4)

               end do
            end do
         end do
      endif

      !if on upper bound and diff
      if((hi(2) .eq. n_cells(2)-1) .and. (bc_hi(2) .eq. 2)) then
         do l = 2,5
            do k = lo(3),hi(3)
               do i = lo(1),hi(1)

                  yflux(i,hi(2)+1,k,l) = 1.4142*yflux(i,hi(2)+1,k,l)        

               end do
            end do
         end do
      endif

      !!!!!!!!!!!!!! z-flux BCs !!!!!!!!!!!!!!

      !if on lower bound and specular
      if((lo(3) .eq. 0) .and. (bc_lo(3) .eq. 1)) then
         do l = 2,5
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)

                  zflux(i,j,0,l) = 0

               end do
            end do
         end do
      endif

      !if on upper bound and specular
      if((hi(3) .eq. n_cells(3)-1) .and. (bc_hi(3) .eq. 1)) then
         do l = 2,5
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)

                  zflux(i,j,hi(3)+1,l) = 0           

               end do
            end do
         end do
      endif

      !if on lower bound and diff
      if((lo(3) .eq. 0) .and. (bc_lo(3) .eq. 2)) then
         do l = 2,5
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)

                  zflux(i,j,0,l) = 1.4142*zflux(i,j,0,5)

               end do
            end do
         end do
      endif

      !if on upper bound and diff
      if((hi(3) .eq. n_cells(3)-1) .and. (bc_hi(3) .eq. 2)) then
         do l = 2,5
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)

                  zflux(i,j,hi(3)+1,l) = 1.4142*zflux(i,j,hi(3)+1,l)        

               end do
            end do
         end do
      endif


    end subroutine stoch_flux_bounds

  subroutine diff_flux_sym(lo,hi, cons, prim, & 
                           eta, zeta, kappa, & 
                           chi, Dij, &
                           fluxx, fluxy, &
#if (AMREX_SPACEDIM == 3)
                           fluxz, &
#endif
                           cornux, cornvx, cornwx, & 
                           cornuy, cornvy, cornwy, & 
                           cornuz, cornvz, cornwz, & 
                           visccorn, dx) bind(C,name="diff_flux_sym")

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

               ! if(.not.single_component) then

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
                  
                  ! print*, "Hack: loc = ", i,j,k
                  ! print*, "Hack: chi = ", chi(i,j,k,:)
                  ! print*, "Hack: chim = ", chi(i-1,j,k,:)
                  ! print*, "Hack: Dij = ", Dij(i,j,k,:,:)
                  ! print*, "Hack: Dijm = ", Dij(i-1,j,k,:,:)
                  ! print*, "Hack: flux = ", fluxx(i,j,k,:)
                  ! stop

               ! end if

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

             ! if(.not.single_component) then
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

             ! end if

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

               ! if(.not.single_component) then
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

               ! end if

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

  end subroutine diff_flux_sym

end module flux_module

































