module flux_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, cell_depth, k_b, bc_lo, bc_hi, n_cells, membrane_cell, visc_type
  use conv_module, only : get_temperature, get_pressure_gas, get_temperature_gas, get_energy, get_density_gas, get_energy_gas
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

      real(amrex_real) :: conserved(nvars), primitive(nprimvars), wgt1, wgt2, vsqr, intenergy, massvec(nspecies), fracvec(nspecies), rho, temp, pt

      integer :: i,j,k,l

      wgt2 = 1.0/12.0 !fourth order interpolation
      wgt1 = 0.5 + wgt2 

      !wgt2 = 0
      !wgt1 = 0.5 + wgt2 !second order

      !wgt2 = (sqrt(7d0)-1d0)/4d0 !adjusted for correct variance fourth order interpolation - this apparently makes the overall spectrum worse
      !wgt1 = (sqrt(7d0)+1d0)/4d0


      !Interpolating conserved quantaties for conv term, apparently this has some advantge over interpolating primitives

      !x flux

    !print *, "flux: ", cons(1,0,0,1), cons(0,0,0,1), cons(-2,0,0,1), cons(-1,0,0,1)
     
      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1


            do l = 2,nprimvars             
               primitive(l) = wgt1*(prim(i,j,k,l)+prim(i-1,j,k,l)) -wgt2*(prim(i-2,j,k,l)+prim(i+1,j,k,l))  
            enddo

            temp = primitive(5)
            pt = primitive(6)
            
            call get_density_gas(pt,rho, temp)

            conserved(1) = rho

            call get_energy_gas(pt, intenergy)

            vsqr = primitive(2)**2 + primitive(3)**2 + primitive(4)**2

            conserved(5) = intenergy + 0.5*rho*vsqr

            xflux(i,j,k,1) = xflux(i,j,k,1) + conserved(1)*primitive(2)
            xflux(i,j,k,2) = xflux(i,j,k,2) + conserved(1)*(primitive(2)**2)+primitive(6)
            xflux(i,j,k,3) = xflux(i,j,k,3) + conserved(1)*primitive(2)*primitive(3)
            xflux(i,j,k,4) = xflux(i,j,k,4) + conserved(1)*primitive(2)*primitive(4)

            xflux(i,j,k,5) = xflux(i,j,k,5) + primitive(2)*conserved(5) + primitive(6)*primitive(2)
 
          end do
        end do
      end do


     !y flux
    
     do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1
         do i = lo(1),hi(1)

           do l = 2,nprimvars 
             primitive(l) = wgt1*(prim(i,j,k,l)+prim(i,j-1,k,l)) -wgt2*(prim(i,j-2,k,l)+prim(i,j+1,k,l))
           enddo

           temp = primitive(5)
           pt = primitive(6)

           call get_density_gas(pt,rho, temp)

           conserved(1) = rho

           call get_energy_gas(pt, intenergy)

           vsqr = primitive(2)**2 + primitive(3)**2 + primitive(4)**2

           conserved(5) = intenergy + 0.5*rho*vsqr

           yflux(i,j,k,1) = yflux(i,j,k,1) + conserved(1)*primitive(3)
           yflux(i,j,k,2) = yflux(i,j,k,2) + conserved(1)*primitive(2)*primitive(3)
           yflux(i,j,k,3) = yflux(i,j,k,3) + conserved(1)*primitive(3)**2+primitive(6)
           yflux(i,j,k,4) = yflux(i,j,k,4) + conserved(1)*primitive(4)*primitive(3)

           yflux(i,j,k,5) = yflux(i,j,k,5) + primitive(3)*conserved(5) + primitive(6)*primitive(3)

         end do
       end do
     end do

     !z flux
     
     do k = lo(3),hi(3)+1
       do j = lo(2),hi(2)
         do i = lo(1),hi(1)

           do l = 2,nprimvars 
             primitive(l) = wgt1*(prim(i,j,k,l)+prim(i,j,k-1,l)) -wgt2*(prim(i,j,k-2,l)+prim(i,j,k+1,l))
           enddo

           temp = primitive(5)
           pt = primitive(6)

           call get_density_gas(pt,rho, temp)

           conserved(1) = rho

           call get_energy_gas(pt, intenergy)

           vsqr = primitive(2)**2 + primitive(3)**2 + primitive(4)**2

           conserved(5) = intenergy + 0.5*rho*vsqr

           zflux(i,j,k,1) = zflux(i,j,k,1) + conserved(1)*primitive(4)
           zflux(i,j,k,2) = zflux(i,j,k,2) + conserved(1)*primitive(2)*primitive(4)
           zflux(i,j,k,3) = zflux(i,j,k,3) + conserved(1)*primitive(3)*primitive(4)
           zflux(i,j,k,4) = zflux(i,j,k,4) + conserved(1)*primitive(4)**2+primitive(6)
           zflux(i,j,k,5) = zflux(i,j,k,5) + primitive(4)*conserved(5) + primitive(6)*primitive(4)

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
                        rancorn, eta, zeta, kappa, dx, dt) bind(C,name="stoch_flux")

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

      real(amrex_real) ::etatF, kappattF, dtinv, volinv, sFac, qFac, velu, velv, velw, wgt1, wgt2, muxp, kxp, weiner(5), fweights(5), nweight, muzepp, muzemp, muzepm, muzemm, phiflx, muyp, muzp, kyp, kzp

      integer :: i,j,k,l

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
                  ! if (abs(visc_type) .eq. 3) then 
                  !    zetaxp = (zeta(i,j,k)*prim(i,j,k,5) + zeta(i-1,j,k)*prim(i-1,j,k,5))
                  ! else
                  !    zetaxp = 0.0
                  ! endif

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

               end do
            end do
         end do

         !!!!!!!!!!!!!!!!!!! y-flux !!!!!!!!!!!!!!!!!!!

         do k = lo(3),hi(3)
            do j = lo(2),hi(2)+1
               do i = lo(1),hi(1)
                  
                  muyp = eta(i,j,k)*prim(i,j,k,5) + eta(i,j-1,k)*prim(i,j-1,k,5)
                  kyp = kappa(i,j,k)*prim(i,j,k,5)**2 + kappa(i,j-1,k)*prim(i,j-1,k,5)**2
                  ! if (abs(visc_type) .eq. 3) then 
                  !    zetayp = (zeta(i,j,k)*prim(i,j,k,5) + zeta(i,j-1,k)*prim(i,j-1,k,5))
                  ! else
                  !    zetayp = 0.0
                  ! endif
                  
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

               end do
            end do
         end do

         !!!!!!!!!!!!!!!!!!! z-flux !!!!!!!!!!!!!!!!!!!

         do k = lo(3),hi(3)+1
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  
                  muzp = eta(i,j,k)*prim(i,j,k,5) + eta(i,j,k-1)*prim(i,j,k-1,5)
                  kzp = kappa(i,j,k)*prim(i,j,k,5)**2 + kappa(i,j,k-1)*prim(i,j,k-1,5)**2
                  ! if (abs(visc_type) .eq. 3) then 
                  !    zetazp = (zeta(i,j,k)*prim(i,j,k,5) + zeta(i,j,k-1)*prim(i,j,k-1,5))
                  ! else
                  !    zetazp = 0.0
                  ! endif

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

  subroutine diff_flux_sym(lo,hi, cons, prim, eta, zeta, kappa, fluxx, fluxy, &
#if (AMREX_SPACEDIM == 3)
                        fluxz, &
#endif
                        cornux, cornvx, cornwx, cornuy, cornvy, cornwy, cornuz, cornvz, cornwz, visccorn, dx) bind(C,name="diff_flux_sym")

      integer         , intent(in   ) :: lo(3),hi(3)
      real(amrex_real), intent(in   ) :: dx(3)
      real(amrex_real), intent(inout) :: fluxx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
      real(amrex_real), intent(inout) :: fluxy(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
      real(amrex_real), intent(inout) :: fluxz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
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

      real(amrex_real), intent(in   ) :: cons(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      real(amrex_real), intent(in   ) :: prim(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)

      real(amrex_real), intent(in   ) :: eta(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
      real(amrex_real), intent(in   ) :: zeta(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
      real(amrex_real), intent(in   ) :: kappa(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))

      integer :: i,j,k 
      real(amrex_real) :: half, dxinv(3), twothirds, muxp, kxp, tauxxp, tauyxp, tauzxp, divxp, phiflx, muyp, kyp, tauxyp, tauyyp, tauzyp, divyp, zetaxp, onetwelfth
#if (AMREX_SPACEDIM == 3)
      real(amrex_real) :: muzp, kzp, tauxzp, tauyzp, tauzzp, divzp
#endif

      dxinv = 1d0/dx
      twothirds = 2d0/3d0
      half = 0.5d0
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
             ! if (abs(visc_type) .eq. 3) then
             !    zetaxp = 0.125d0*(zeta(i,j-1,k-1) + zeta(i-1,j-1,k-1) + zeta(i,j,k-1) + zeta(i-1,j,k-1)+ &
             !         zeta(i,j-1,k) + zeta(i-1,j-1,k) + zeta(i,j,k) + zeta(i-1,j,k))
             ! else
             !    zetaxp = 0.0
             ! endif
             
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

































