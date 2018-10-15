module flux_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, cell_depth, k_b, bc_lo, bc_hi, n_cells
  use conv_module, only : get_temperature, get_pressure_gas
  implicit none

  private

  public :: diff_flux

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

      real(amrex_real), intent(in   ) :: cons(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)
      real(amrex_real), intent(in   ) :: prim(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nprimvars)

      real(amrex_real), intent(in   ) :: eta(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc)
      real(amrex_real), intent(in   ) :: zeta(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc)
      real(amrex_real), intent(in   ) :: kappa(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc)

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
            div = dux + dvy + dwz

            taux = zetaf*div + etaf*(2d0*dux - twothirds*div)
            tauy = etaf*(duy + dvx)
            tauz = etaf*(duz + dwx)

            !momentum fluxes
            xflux(i+1,j,k,2) = xflux(i+1,j,k,2) - taux
            xflux(i+1,j,k,3) = xflux(i+1,j,k,3) - tauy
            xflux(i+1,j,k,4) = xflux(i+1,j,k,4) - tauz
            
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

            !energy flux
            xflux(i+1,j,k,5) = xflux(i+1,j,k,5) - (u*taux + v*tauy + w*tauz) - kappaf*(dtx + dty + dtz)

!            if(i .eq. -1) then
! 
!              print *, i,j,k, "temp: ", prim(i,j,k,5)
!              print *, i+1,j,k, "temp: ", prim(i+1,j,k,5)

!              !print *, i,j,k, "dty: ", dty

!            endif

          end do
        end do
      end do

!      do k = lo(3),hi(3)
!        do j = lo(2),hi(2)
!          do i = lo(1)-1,hi(1)


!            totalE = totalE + cons(i,j,k,5)
!          enddo
!        enddo
!      enddo


!      print *, "fluxlo: ", fluxlo, " fluxhi: ", fluxhi
!      print *, "Total E: ", totalE

      duz = 0
      dvz = 0
      dwz = 0
      dtz = 0

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
            tauy = zetaf*div + etaf*(2d0*dvy - twothirds*div)                  
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
            tauz = zetaf*div + etaf*(2d0*dwz - twothirds*div)

            !momentum fluxes
            zflux(i,j,k+1,2) = zflux(i,j,k+1,2) - taux
            zflux(i,j,k+1,3) = zflux(i,j,k+1,3) - tauy
            zflux(i,j,k+1,4) = zflux(i,j,k+1,4) - tauz
            
            ! z derivative of T on z face
            dtz = (prim(i,j,k+1,5) - prim(i,j,k,5))*dxinv(1)

            ! y derivative of T on z face
            up    = 0.5*(prim(i,j+1,k,5)+prim(i,j+1,k+1,5))
            down  = 0.5*(prim(i,j-1,k,5)+prim(i,j-1,k+1,5))
            dty   = 0.5*(up - down)*dxinv(2)

            ! x derivative of T on z face
            up    = 0.5*(prim(i+1,j,k+1,5)+prim(i+1,j,k,5))
            down  = 0.5*(prim(i-1,j,k+1,5)+prim(i-1,j,k,5))
            dtx = 0.5*(up - down)*dxinv(3)

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

      real(amrex_real), intent(in   ) :: cons(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)
      real(amrex_real), intent(in   ) :: prim(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nprimvars)

      real(amrex_real) :: conserved(nvars), primitive(6), wgt1, wgt2, vsqr, intenergy, massvec(nspecies), fracvec(nspecies)

      integer :: i,j,k,l

      wgt2 = 1.0/12.0 !fourth order interpolation
      !wgt2 = 0  !second order interpolation
      wgt1 = 0.5 + wgt2 

      !Interpolating conserved quantaties for conv term, apparently this has some advantge over interpolating primitives

      !x flux
     
      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1)-1,hi(1)

            do l = 1,nvars 
              conserved(l) = wgt1*(cons(i+1,j,k,l)+cons(i,j,k,l)) -wgt2*(cons(i-1,j,k,l)+cons(i+2,j,k,l)) !Better to interpolate primitives here?
            enddo

            primitive(1) = conserved(1)  
            primitive(2) = conserved(2)/conserved(1)
            primitive(3) = conserved(3)/conserved(1)
            primitive(4) = conserved(4)/conserved(1)

            vsqr = primitive(2)**2 + primitive(3)**2 + primitive(4)**2

            intenergy = conserved(5) - 0.5*vsqr*conserved(1)
 
            fracvec = conserved(6:nvars)

            massvec = fracvec*conserved(1)

            call get_temperature(intenergy, massvec, primitive(5))

            call get_pressure_gas(primitive(6), fracvec, primitive(1), primitive(5))

            xflux(i+1,j,k,1) = xflux(i+1,j,k,1) + primitive(1)*primitive(2)
            xflux(i+1,j,k,2) = xflux(i+1,j,k,2) + primitive(1)*primitive(2)**2+primitive(6)
            xflux(i+1,j,k,3) = xflux(i+1,j,k,3) + primitive(1)*primitive(3)*primitive(2)
            xflux(i+1,j,k,4) = xflux(i+1,j,k,4) + primitive(1)*primitive(4)*primitive(2)
            xflux(i+1,j,k,5) = xflux(i+1,j,k,5) + primitive(2)*conserved(5) + primitive(6)*primitive(2)


          end do
        end do
      end do

      !y flux
     
      do k = lo(3),hi(3)
        do j = lo(2)-1,hi(2)
          do i = lo(1),hi(1)

            do l = 1,nvars 
              conserved(l) = wgt1*(cons(i,j+1,k,l)+cons(i,j,k,l)) -wgt2*(cons(i,j-1,k,l)+cons(i,j+2,k,l)) !Better to interpolate primitives here?
            enddo

            primitive(1) = conserved(1)  
            primitive(2) = conserved(2)/conserved(1)
            primitive(3) = conserved(3)/conserved(1)
            primitive(4) = conserved(4)/conserved(1)

            vsqr = primitive(2)**2 + primitive(3)**2 + primitive(4)**2

            intenergy = conserved(5) - 0.5*vsqr*conserved(1)
 
            fracvec = conserved(6:nvars)

            massvec = fracvec*conserved(1)

            call get_temperature(intenergy, massvec, primitive(5))

            call get_pressure_gas(primitive(6), fracvec, primitive(1), primitive(5))

            yflux(i,j+1,k,1) = yflux(i,j+1,k,1) + primitive(1)*primitive(3)
            yflux(i,j+1,k,2) = yflux(i,j+1,k,2) + primitive(1)*primitive(2)*primitive(3)
            yflux(i,j+1,k,3) = yflux(i,j+1,k,3) + primitive(1)*primitive(3)**2+primitive(6)
            yflux(i,j+1,k,4) = yflux(i,j+1,k,4) + primitive(1)*primitive(4)*primitive(3)
            yflux(i,j+1,k,5) = yflux(i,j+1,k,5) + primitive(3)*conserved(5) + primitive(6)*primitive(3)


          end do
        end do
      end do

      !z flux
     
      do k = lo(3)-1,hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)

            do l = 1,nvars 
              conserved(l) = wgt1*(cons(i,j,k+1,l)+cons(i,j,k,l)) -wgt2*(cons(i,j,k-1,l)+cons(i,j,k+2,l)) !Better to interpolate primitives here?
            enddo

            primitive(1) = conserved(1)  
            primitive(2) = conserved(2)/conserved(1)
            primitive(3) = conserved(3)/conserved(1)
            primitive(4) = conserved(4)/conserved(1)

            vsqr = primitive(2)**2 + primitive(3)**2 + primitive(4)**2

            intenergy = conserved(5) - 0.5*vsqr*conserved(1)
 
            fracvec = conserved(6:nvars)

            massvec = fracvec*conserved(1)

            call get_temperature(intenergy, massvec, primitive(5))

            call get_pressure_gas(primitive(6), fracvec, primitive(1), primitive(5))

            zflux(i,j,k+1,1) = zflux(i,j,k+1,1) + primitive(1)*primitive(4)
            zflux(i,j,k+1,2) = zflux(i,j,k+1,2) + primitive(1)*primitive(2)*primitive(4)
            zflux(i,j,k+1,3) = zflux(i,j,k+1,3) + primitive(1)*primitive(3)*primitive(4)
            zflux(i,j,k+1,4) = zflux(i,j,k+1,4) + primitive(1)*primitive(4)**2+primitive(6)
            zflux(i,j,k+1,5) = zflux(i,j,k+1,5) + primitive(4)*conserved(5) + primitive(6)*primitive(4)

          end do
        end do
      end do


  end subroutine hyp_flux

  subroutine stoch_flux(lo,hi, cons, prim, xflux, yflux, &
#if (AMREX_SPACEDIM == 3)
                        zflux, &
#endif
                        xsflux, ysflux, &
#if (AMREX_SPACEDIM == 3)
                        zsflux, &
#endif
                        eta, zeta, kappa, dx, dt) bind(C,name="stoch_flux")

      integer         , intent(in   ) :: lo(3),hi(3)
      real(amrex_real), intent(in   ) :: dx(3), dt
      real(amrex_real), intent(inout) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
      real(amrex_real), intent(inout) :: yflux(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
      real(amrex_real), intent(inout) :: zflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif

      real(amrex_real), intent(inout) :: xsflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
      real(amrex_real), intent(inout) :: ysflux(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
      real(amrex_real), intent(inout) :: zsflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif

      real(amrex_real), intent(in   ) :: cons(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)
      real(amrex_real), intent(in   ) :: prim(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nprimvars)


      real(amrex_real), intent(in   ) :: eta(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc)
      real(amrex_real), intent(in   ) :: zeta(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc)
      real(amrex_real), intent(in   ) :: kappa(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc)

      real(amrex_real) ::etatF, kappattF, dtinv, volinv, sFac, qFac, velF

      integer :: i,j,k,l

      dtinv = 1d0/dt
#if (AMREX_SPACEDIM == 3)
      volinv = 1d0/(dx(1)*dx(2)*dx(3))
#endif

#if (AMREX_SPACEDIM == 2)
      volinv = 1d0/(dx(1)*dx(2)*cell_depth)
#endif

      sFac = 4d0*k_b*volinv*dtinv/3d0
      qFac = k_b*volinv*dtinv

     
      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1)-1,hi(1)

            kappattF = (kappa(i,j,k)*prim(i,j,k,5)*prim(i,j,k,5)+kappa(i+1,j,k)*prim(i+1,j,k,5)*prim(i+1,j,k,5))
            etatF = (eta(i,j,k)*prim(i,j,k,5)+eta(i+1,j,k)*prim(i+1,j,k,5))
            velF = 0.5*(prim(i,j,k,2)+prim(i+1,j,k,2))

            xflux(i+1,j,k,2) = xflux(i+1,j,k,2) + sqrt(sFac*etatF)*xsflux(i+1,j,k,2)
            xflux(i+1,j,k,5) = xflux(i+1,j,k,5) + sqrt(qFac*kappattF)*xsflux(i+1,j,k,5) + velF*sqrt(sFac)*xsflux(i+1,j,k,2)

          end do
        end do
      end do

      !if on lower bound and not periodic
      if((lo(1) .eq. 0) .and. (bc_lo(1) .ne. -1)) then
        do k = lo(3),hi(3)
          do j = lo(2),hi(2)
            
              xflux(0,j,k,2) = 0        
              xflux(0,j,k,5) = 0        

          end do
        end do
      endif

      !if on upper bound and not periodic
      if((hi(1) .eq. n_cells(1)-1) .and. (bc_hi(1) .ne. -1)) then
        do k = lo(3),hi(3)
          do j = lo(2),hi(2)
            
              xflux(hi(1)+1,j,k,2) = 0        
              xflux(hi(1)+1,j,k,5) = 0        

          end do
        end do
      endif

      do k = lo(3),hi(3)
        do j = lo(2)-1,hi(2)
          do i = lo(1),hi(1)

            kappattF = (kappa(i,j,k)*prim(i,j,k,5)*prim(i,j,k,5)+kappa(i,j+1,k)*prim(i,j+1,k,5)*prim(i,j+1,k,5))
            etatF = (eta(i,j,k)*prim(i,j,k,5)+eta(i,j+1,k)*prim(i,j+1,k,5))
            velF = 0.5*(prim(i,j,k,3)+prim(i,j+1,k,3))


            yflux(i,j+1,k,3) = yflux(i,j+1,k,3) + sqrt(sFac*etatF)*ysflux(i,j+1,k,3)
            yflux(i,j+1,k,5) = yflux(i,j+1,k,5) + sqrt(qFac*kappattF)*ysflux(i,j+1,k,5) + velF*sqrt(sFac)*ysflux(i,j+1,k,3)

          end do
        end do
      end do

      !if on lower bound and not periodic
      if((lo(2) .eq. 0) .and. (bc_lo(2) .ne. -1)) then
        do k = lo(3),hi(3)
          do i = lo(1),hi(1)
            
              xflux(i,0,k,3) = 0        
              xflux(i,0,k,5) = 0        

          end do
        end do
      endif

      !if on upper bound and not periodic
      if((hi(2) .eq. n_cells(2)-1) .and. (bc_hi(2) .ne. -1)) then
        do k = lo(3),hi(3)
          do i = lo(1),hi(1)
            
              xflux(i,hi(2)+1,k,3) = 0        
              xflux(i,hi(2)+1,k,5) = 0        

          end do
        end do
      endif

      do k = lo(3)-1,hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)

            kappattF = (kappa(i,j,k)*prim(i,j,k,5)*prim(i,j,k,5)+kappa(i,j,k+1)*prim(i,j,k+1,5)*prim(i,j,k+1,5))
            etatF = (eta(i,j,k)*prim(i,j,k,5)+eta(i,j,k+1)*prim(i,j,k+1,5))
            velF = 0.5*(prim(i,j,k,4)+prim(i,j,k+1,4))

            zflux(i,j,k+1,4) = zflux(i,j,k+1,4) + sqrt(sFac*etatF)*zsflux(i,j,k+1,4)
            zflux(i,j,k+1,5) = zflux(i,j,k+1,5) + sqrt(qFac*kappattF)*zsflux(i,j,k+1,5) + velF*sqrt(sFac)*zsflux(i,j,k+1,4)

          end do
        end do
      end do

      !if on lower bound and not periodic
      if((lo(3) .eq. 0) .and. (bc_lo(3) .ne. -1)) then
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)
            
              zflux(i,j,0,4) = 0        
              zflux(i,j,0,5) = 0      

              print *, "Here!"

          end do
        end do
      endif

      !if on upper bound and not periodic
      if((hi(3) .eq. n_cells(3)-1) .and. (bc_hi(3) .ne. -1)) then
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)
            
              zflux(i,j,hi(3)+1,4) = 0        
              zflux(i,j,hi(3)+1,5) = 0        

              print *, "Here!"

          end do
        end do
      endif
      


  end subroutine stoch_flux

end module flux_module
