module time_step_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, cfl
  use GL_namelist_module
  implicit none

  private

  public :: rk2_stage1, rk2_stage2, initphi, integrate, setdt

contains

  subroutine rk2_stage1(lo,hi, phi, phin, rannums, integral, dx, dt) bind(C,name="rk2_stage1")

      integer         , intent(in   ) :: lo(2),hi(2)
      real(amrex_real), intent(in   ) :: dx(2), dt

      real(amrex_real), intent(inout) :: phi(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2))
      real(amrex_real), intent(inout) :: phin(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2))

      real(amrex_real), intent(in)    :: rannums(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2))
      real(amrex_real), intent(in)    :: integral

      integer :: i,j,k,l
      real(amrex_real) :: func,factor

       factor = sqrt(2.d0*noise_coef/(dt*dx(1)*dx(2)))

         do  j=lo(2),hi(2)
           do  i=lo(1),hi(1)

             func =  acoef+2.d0*bcoef*phi(i,j)+3.d0*ccoef*phi(i,j)**2+4.d0*dcoef*phi(i,j)**3

             phin(i,j) = phi(i,j) + dt*diff_coef*(phi(i+1,j)-2.d0*phi(i,j)+phi(i-1,j))/dx(1)**2  &
              + dt*diff_coef*(phi(i,j+1)-2.d0*phi(i,j)+phi(i,j-1))/dx(2)**2  &
              -dt*func-dt*umbrella*integral + dt*factor * rannums(i,j)


          enddo
        enddo

         do  j=lo(2),hi(2)
           do  i=lo(1),hi(1)

               phi(i,j)= phin(i,j)

          enddo
        enddo

  end subroutine rk2_stage1

  subroutine rk2_stage2(lo,hi, phi, phin, rannums, integral, dx, dt) bind(C,name="rk2_stage2")

      integer         , intent(in   ) :: lo(2),hi(2)
      real(amrex_real), intent(in   ) :: dx(2), dt

      real(amrex_real), intent(inout) :: phi(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2))
      real(amrex_real), intent(inout) :: phin(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2))

      real(amrex_real), intent(in)    :: rannums(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2))
      real(amrex_real), intent(in)    :: integral

      integer :: i,j,k,l
      real(amrex_real) :: func,factor

       factor = sqrt(2.d0*noise_coef/(dt*dx(1)*dx(2)))

         do  j=lo(2),hi(2)
           do  i=lo(1),hi(1)


             func =  acoef+2.d0*bcoef*phi(i,j)+3.d0*ccoef*phi(i,j)**2+4.d0*dcoef*phi(i,j)**3

             phin(i,j) = phi(i,j) + dt*diff_coef*(phi(i+1,j)-2.d0*phi(i,j)+phi(i-1,j))/dx(1)**2  &
              + dt*diff_coef*(phi(i,j+1)-2.d0*phi(i,j)+phi(i,j-1))/dx(2)**2  &
              -dt*func-dt*umbrella*integral + dt*factor * rannums(i,j)

          enddo
        enddo

  end subroutine rk2_stage2

  subroutine initphi(lo,hi, phi, dx) bind(C,name="initphi")

      integer         , intent(in   ) :: lo(2),hi(2)
      real(amrex_real), intent(in   ) :: dx(2)

      real(amrex_real), intent(inout) :: phi(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2))

      real(amrex_real) :: xloc,yloc
      integer :: i,j

         do  j=lo(2),hi(2)
           do  i=lo(1),hi(1)

              xloc = dx(1)*dfloat(i-1)
              yloc = dx(2)*dfloat(j-1)

           !  phi(i,j) =  (sin(2.d0*pi*xloc)*sin(2.d0*pi*yloc)) **2
              if( (xloc-.5d0)**2 + (yloc-.5d0)**2 .lt. rad**2)then
                 phi(i,j) = 1.d0
              else
                 phi(i,j) = 0.d0
              endif


          enddo
        enddo

  end subroutine initphi

  subroutine integrate(lo,hi, phi, dx, integral) bind(C,name="integrate")

      integer         , intent(in   ) :: lo(2),hi(2)
      real(amrex_real), intent(in   ) :: dx(2)
      real(amrex_real), intent(out  ) :: integral

      real(amrex_real), intent(inout) :: phi(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2))

      real(amrex_real) :: xloc,yloc
      integer :: i,j

         integral = 0.d0

         do  j=lo(2),hi(2)
           do  i=lo(1),hi(1)

              integral = integral + phi(i,j) - phi0

          enddo
        enddo

        integral = integral*dx(1)*dx(2)

  end subroutine integrate

  subroutine setdt (  dx, dt ) bind(C,name="setdt")

  double precision::  dt, dx(2)

  if(cfl.gt.0.d0)then

      dt = .25d0*cfl*min(dx(1)**2, dx(2)**2)/diff_coef
    
  endif

  end subroutine setdt

end module time_step_module
