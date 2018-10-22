module time_step_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars
  implicit none

  private

  public :: rk3_stage1, rk3_stage2, rk3_stage3

contains

  subroutine rk3_stage1(lo,hi, cu, cup, source, xflux, yflux, &
#if (AMREX_SPACEDIM == 3)
                        zflux, &
#endif
                        dx, dt) bind(C,name="rk3_stage1")

      integer         , intent(in   ) :: lo(3),hi(3)
      real(amrex_real), intent(in   ) :: dx(3), dt

      real(amrex_real), intent(inout) :: cu(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)
      real(amrex_real), intent(inout) :: cup(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)

      real(amrex_real), intent(in)    :: source(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)

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

!                if(l .eq. 1) then
!                  print *, "Stage1"
!                  print *, i,j,k
!                  print *, "CUpre: ", cu(i,j,k,5)
!                endif

                 cup(i,j,k,l) = cu(i,j,k,l)                                      &
                                - dt*(xflux(i+1,j,k,l)-xflux(i,j,k,l))*dxinv(1)  & 
                                !- dt*(yflux(i,j+1,k,l)-yflux(i,j,k,l))*dxinv(2)  &
#if (AMREX_SPACEDIM == 3)
                               ! - dt*(zflux(i,j,k+1,l)-zflux(i,j,k,l))*dxinv(3)  &
#endif
                                + dt*source(i,j,k,l)


                        
!                if(l .eq. 1) then

!                  print *, "FluxL: ", zflux(i,j,k,5)*dt*dxinv(3)
!                  print *, "FluxR: ", -zflux(i,j,k+1,5)*dt*dxinv(3)
!                  print *, "CUpost: ", cu(i,j,k,5)

!                endif

           enddo
          enddo
        enddo
      enddo

     !call exit()

  end subroutine rk3_stage1

  subroutine rk3_stage2(lo,hi, cu, cup, cup2, source, xflux, yflux, &
#if (AMREX_SPACEDIM == 3)
                        zflux, &
#endif
                        dx, dt) bind(C,name="rk3_stage2")

      integer         , intent(in   ) :: lo(3),hi(3)
      real(amrex_real), intent(in   ) :: dx(3), dt

      real(amrex_real), intent(in   ) :: cu(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)
      real(amrex_real), intent(in   ) :: cup(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)
      real(amrex_real), intent(inout) :: cup2(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)

      real(amrex_real), intent(in)    :: source(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)

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
                               ! - dt*(yflux(i,j+1,k,l)-yflux(i,j,k,l))*dxinv(2)  &
#if (AMREX_SPACEDIM == 3)
                              !  - dt*(zflux(i,j,k+1,l)-zflux(i,j,k,l))*dxinv(3)  &
#endif
                                + dt*source(i,j,k,l))


           enddo
          enddo
        enddo
      enddo

  end subroutine rk3_stage2

  subroutine rk3_stage3(lo,hi, cu, cup, cup2, source, xflux, yflux, &
#if (AMREX_SPACEDIM == 3)
                        zflux, &
#endif
                        dx, dt) bind(C,name="rk3_stage3")

      integer         , intent(in   ) :: lo(3),hi(3)
      real(amrex_real), intent(in   ) :: dx(3), dt

      real(amrex_real), intent(inout) :: cu(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)
      real(amrex_real), intent(in   ) :: cup(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)
      real(amrex_real), intent(in   ) :: cup2(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)

      real(amrex_real), intent(in)    :: source(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nvars)

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
                              !  - dt*(yflux(i,j+1,k,l)-yflux(i,j,k,l))*dxinv(2)  &
#if (AMREX_SPACEDIM == 3)
                              !  - dt*(zflux(i,j,k+1,l)-zflux(i,j,k,l))*dxinv(3)  &
#endif
                                + dt*source(i,j,k,l))

           enddo
          enddo
        enddo
      enddo

  end subroutine rk3_stage3

end module time_step_module
