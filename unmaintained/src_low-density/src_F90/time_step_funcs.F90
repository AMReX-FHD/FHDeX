module time_step_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars
  implicit none

  private

  public :: rk3_stage1, rk3_stage2, rk3_stage3, euler_step

contains

  subroutine euler_step(lo,hi, cu, xflux, yflux, &
#if (AMREX_SPACEDIM == 3)
                        zflux, &
#endif
                        dx, dt) bind(C,name="euler_step")

      integer         , intent(in   ) :: lo(3),hi(3)
      real(amrex_real), intent(in   ) :: dx(3), dt

      real(amrex_real), intent(inout) :: cu(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)


      real(amrex_real), intent(in   ) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
      real(amrex_real), intent(in   ) :: yflux(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
      real(amrex_real), intent(in   ) :: zflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif
      integer :: i,j,k,l
      real(amrex_real) :: dxinv(3)

      dxinv = 1d0/dx

      !print *, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3)

      do  l=1,nvars
        do  k=lo(3),hi(3)
         do  j=lo(2),hi(2)
           do  i=lo(1),hi(1)

                print *, "density before", cu(i,j,k,l)
                 cu(i,j,k,l) = cu(i,j,k,l)                                      &
                                - dt*(xflux(i+1,j,k,l)-xflux(i,j,k,l))*dxinv(1)  &
                                - dt*(yflux(i,j+1,k,l)-yflux(i,j,k,l))*dxinv(2)  &
#if (AMREX_SPACEDIM == 3)
                                - dt*(zflux(i,j,k+1,l)-zflux(i,j,k,l))*dxinv(3)
#endif
                                print *, "density after" ,cu(i,j,k,l)
           enddo
          enddo
        enddo
      enddo

  end subroutine euler_step

  subroutine rk3_stage1(lo,hi, cu, cup, xflux, yflux, &
#if (AMREX_SPACEDIM == 3)
                        zflux, &
#endif
                        dx, dt) bind(C,name="rk3_stage1")

      integer         , intent(in   ) :: lo(3),hi(3)
      real(amrex_real), intent(in   ) :: dx(3), dt

      real(amrex_real), intent(inout) :: cu(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      real(amrex_real), intent(inout) :: cup(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)


      real(amrex_real), intent(in   ) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
      real(amrex_real), intent(in   ) :: yflux(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
      real(amrex_real), intent(in   ) :: zflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif
      integer :: i,j,k,l
      real(amrex_real) :: dxinv(3)

      dxinv = 1d0/dx

      do  l=1,nvars
        do  k=lo(3),hi(3)
         do  j=lo(2),hi(2)
           do  i=lo(1),hi(1)


                 cup(i,j,k,l) = cu(i,j,k,l)                                      &
                                - dt*(xflux(i+1,j,k,l)-xflux(i,j,k,l))*dxinv(1)  &
                                - dt*(yflux(i,j+1,k,l)-yflux(i,j,k,l))*dxinv(2)  &
#if (AMREX_SPACEDIM == 3)
                                - dt*(zflux(i,j,k+1,l)-zflux(i,j,k,l))*dxinv(3)
#endif

           enddo
          enddo
        enddo
      enddo

  end subroutine rk3_stage1

  subroutine rk3_stage2(lo,hi, cu, cup, cup2, xflux, yflux, &
#if (AMREX_SPACEDIM == 3)
                        zflux, &
#endif
                        dx, dt) bind(C,name="rk3_stage2")

      integer         , intent(in   ) :: lo(3),hi(3)
      real(amrex_real), intent(in   ) :: dx(3), dt

      real(amrex_real), intent(in   ) :: cu(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      real(amrex_real), intent(in   ) :: cup(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      real(amrex_real), intent(inout) :: cup2(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)


      real(amrex_real), intent(in   ) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
      real(amrex_real), intent(in   ) :: yflux(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
      real(amrex_real), intent(in   ) :: zflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif
      integer :: i,j,k,l
      real(amrex_real) :: dxinv(3)

      dxinv = 1d0/dx

      do  l=1,nvars
        do  k=lo(3),hi(3)
         do  j=lo(2),hi(2)
           do  i=lo(1),hi(1)

                 cup2(i,j,k,l) =  0.25d0*(3.0d0*cu(i,j,k,l) + cup(i,j,k,l)              &
                                - dt*(xflux(i+1,j,k,l)-xflux(i,j,k,l))*dxinv(1)  &
                                - dt*(yflux(i,j+1,k,l)-yflux(i,j,k,l))*dxinv(2)  &
#if (AMREX_SPACEDIM == 3)
                                - dt*(zflux(i,j,k+1,l)-zflux(i,j,k,l))*dxinv(3)  &
#endif
                               )


           enddo
          enddo
        enddo
      enddo

  end subroutine rk3_stage2

  subroutine rk3_stage3(lo,hi, cu, cup, cup2, xflux, yflux, &
#if (AMREX_SPACEDIM == 3)
                        zflux, &
#endif
                        dx, dt) bind(C,name="rk3_stage3")

      integer         , intent(in   ) :: lo(3),hi(3)
      real(amrex_real), intent(in   ) :: dx(3), dt

      real(amrex_real), intent(inout) :: cu(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      real(amrex_real), intent(in   ) :: cup(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      real(amrex_real), intent(in   ) :: cup2(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)

      real(amrex_real), intent(in   ) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
      real(amrex_real), intent(in   ) :: yflux(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
      real(amrex_real), intent(in   ) :: zflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif
      integer :: i,j,k,l
      real(amrex_real) :: dxinv(3), twothirds

      dxinv = 1d0/dx
      twothirds = 2d0/3d0

      do  l=1,nvars
        do  k=lo(3),hi(3)
         do  j=lo(2),hi(2)
           do  i=lo(1),hi(1)

                 cu(i,j,k,l) =  twothirds*(0.5*cu(i,j,k,l) + cup2(i,j,k,l)              &
                                - dt*(xflux(i+1,j,k,l)-xflux(i,j,k,l))*dxinv(1)  &
                                - dt*(yflux(i,j+1,k,l)-yflux(i,j,k,l))*dxinv(2)  &
#if (AMREX_SPACEDIM == 3)
                                - dt*(zflux(i,j,k+1,l)-zflux(i,j,k,l))*dxinv(3)  &
#endif
                                )

           enddo
          enddo
        enddo
      enddo


  end subroutine rk3_stage3

end module time_step_module
