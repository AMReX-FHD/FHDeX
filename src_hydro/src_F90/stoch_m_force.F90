module stoch_m_force_module

  use amrex_error_module
  ! use common_namelist_module, only: visc_type

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 2)

      subroutine stoch_m_force(lo,hi, &
                               flux_cc, fluxcclo, fluxcchi, ncomp_cc, &
                               flux_nd, fluxndlo, fluxndhi, ncomp_nd, &
                               divx, divxlo, divxhi, &
                               divy, divylo, divyhi, &
                               dx, increment) bind(C,name="stoch_m_force")

      integer         , intent(in   ) :: lo(2),hi(2)
      integer         , intent(in   ) :: ncomp_cc, ncomp_nd
      integer         , intent(in   ) :: fluxcclo(2),fluxcchi(2), fluxndlo(2),fluxndhi(2)
      integer         , intent(in   ) :: divxlo(2),divxhi(2), divylo(2),divyhi(2)
      double precision, intent(in   ) :: flux_cc(fluxcclo(1):fluxcchi(1),fluxcclo(2):fluxcchi(2),ncomp_cc)
      double precision, intent(in   ) :: flux_nd(fluxndlo(1):fluxndhi(1),fluxndlo(2):fluxndhi(2),ncomp_nd)
      double precision, intent(inout) ::    divx(divxlo(1):divxhi(1),divxlo(2):divxhi(2))
      double precision, intent(inout) ::    divy(divylo(1):divyhi(1),divylo(2):divyhi(2))
      double precision, intent(in   ) :: dx(2)
      integer         , intent(in   ) :: increment

      integer :: i,j
      double precision :: dxinv

      dxinv = 1.d0/dx(1)

      if (increment == 1) then

         ! divergence on x-faces
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            divx(i,j) = divx(i,j) + (flux_cc(i,j,1) - flux_cc(i-1,j,1)) * dxinv + &
                                    (flux_nd(i,j+1,1) - flux_nd(i,j,1)) * dxinv
         end do
         end do

         ! divergence on y-faces
         do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            divy(i,j) = divy(i,j) + (flux_nd(i+1,j,2) - flux_nd(i,j,2)) * dxinv + &
                                    (flux_cc(i,j,2) - flux_cc(i,j-1,2)) * dxinv
         end do
         end do

      else

         ! divergence on x-faces
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            divx(i,j) = (flux_cc(i,j,1) - flux_cc(i-1,j,1)) * dxinv + &
                        (flux_nd(i,j+1,1) - flux_nd(i,j,1)) * dxinv
         end do
         end do

         ! divergence on y-faces
         do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            divy(i,j) = (flux_nd(i+1,j,2) - flux_nd(i,j,2)) * dxinv + &
                        (flux_cc(i,j,2) - flux_cc(i,j-1,2)) * dxinv
         end do
         end do

      end if

    end subroutine stoch_m_force

#endif

#if (AMREX_SPACEDIM == 3)

      subroutine stoch_m_force(lo,hi, &
                               flux_cc, fluxcclo, fluxcchi, ncomp_cc, &
                               flux_xy, fluxxylo, fluxxyhi, ncomp_xy, &
                               flux_xz, fluxxzlo, fluxxzhi, ncomp_xz, &
                               flux_yz, fluxyzlo, fluxyzhi, ncomp_yz, &
                               divx, divxlo, divxhi, &
                               divy, divylo, divyhi, &
                               divz, divzlo, divzhi, &
                               dx, increment) bind(C,name="stoch_m_force")

      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: ncomp_cc, ncomp_xy, ncomp_xz, ncomp_yz
      integer         , intent(in   ) :: fluxcclo(3),fluxcchi(3)
      integer         , intent(in   ) :: fluxxylo(3),fluxxyhi(3), fluxxzlo(3),fluxxzhi(3), fluxyzlo(3),fluxyzhi(3)
      integer         , intent(in   ) :: divxlo(3),divxhi(3), divylo(3),divyhi(3), divzlo(3),divzhi(3)
      double precision, intent(in   ) :: flux_cc(fluxcclo(1):fluxcchi(1),fluxcclo(2):fluxcchi(2),fluxcclo(3):fluxcchi(3),ncomp_cc)
      double precision, intent(in   ) :: flux_xy(fluxxylo(1):fluxxyhi(1),fluxxylo(2):fluxxyhi(2),fluxxylo(3):fluxxyhi(3),ncomp_xy)
      double precision, intent(in   ) :: flux_xz(fluxxzlo(1):fluxxzhi(1),fluxxzlo(2):fluxxzhi(2),fluxxzlo(3):fluxxzhi(3),ncomp_xz)
      double precision, intent(in   ) :: flux_yz(fluxyzlo(1):fluxyzhi(1),fluxyzlo(2):fluxyzhi(2),fluxyzlo(3):fluxyzhi(3),ncomp_yz)
      double precision, intent(inout) ::    divx(divxlo(1):divxhi(1),divxlo(2):divxhi(2),divxlo(3):divxhi(3))
      double precision, intent(inout) ::    divy(divylo(1):divyhi(1),divylo(2):divyhi(2),divylo(3):divyhi(3))
      double precision, intent(inout) ::    divz(divzlo(1):divzhi(1),divzlo(2):divzhi(2),divzlo(3):divzhi(3))
      double precision, intent(in   ) :: dx(3)
      integer         , intent(in   ) :: increment

      integer :: i,j,k
      double precision :: dxinv

      dxinv = 1.d0/dx(1)

      if (increment == 1) then

         ! divergence on x-faces
         do k=lo(3),hi(3)
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            divx(i,j,k) = divx(i,j,k) + (flux_cc(i,j,k,1) - flux_cc(i-1,j,k,1)) * dxinv + &
                                        (flux_xy(i,j+1,k,1) - flux_xy(i,j,k,1)) * dxinv + &
                                        (flux_xz(i,j,k+1,1) - flux_xz(i,j,k,1)) * dxinv
         end do
         end do
         end do

         ! divergence on y-faces
         do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            divy(i,j,k) = divy(i,j,k) + (flux_xy(i+1,j,k,2) - flux_xy(i,j,k,2)) * dxinv + &
                                        (flux_cc(i,j,k,2) - flux_cc(i,j-1,k,2)) * dxinv + &
                                        (flux_yz(i,j,k+1,1) - flux_yz(i,j,k,1)) * dxinv
         end do
         end do
         end do

         ! divergence on z-faces
         do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            divz(i,j,k) = divz(i,j,k) + (flux_xz(i+1,j,k,2) - flux_xz(i,j,k,2)) * dxinv + &
                                        (flux_yz(i,j+1,k,2) - flux_yz(i,j,k,2)) * dxinv + &
                                        (flux_cc(i,j,k,3) - flux_cc(i,j,k-1,3)) * dxinv
         end do
         end do
         end do

      else

         ! divergence on x-faces
         do k=lo(3),hi(3)
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            divx(i,j,k) = (flux_cc(i,j,k,1) - flux_cc(i-1,j,k,1)) * dxinv + &
                          (flux_xy(i,j+1,k,1) - flux_xy(i,j,k,1)) * dxinv + &
                          (flux_xz(i,j,k+1,1) - flux_xz(i,j,k,1)) * dxinv
         end do
         end do
         end do

         ! divergence on y-faces
         do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            divy(i,j,k) = (flux_xy(i+1,j,k,2) - flux_xy(i,j,k,2)) * dxinv + &
                          (flux_cc(i,j,k,2) - flux_cc(i,j-1,k,2)) * dxinv + &
                          (flux_yz(i,j,k+1,1) - flux_yz(i,j,k,1)) * dxinv
         end do
         end do
         end do

         ! divergence on z-faces
         do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            divz(i,j,k) = (flux_xz(i+1,j,k,2) - flux_xz(i,j,k,2)) * dxinv + &
                          (flux_yz(i,j+1,k,2) - flux_yz(i,j,k,2)) * dxinv + &
                          (flux_cc(i,j,k,3) - flux_cc(i,j,k-1,3)) * dxinv
         end do
         end do
         end do

      end if

    end subroutine stoch_m_force
    
#endif

end module stoch_m_force_module
