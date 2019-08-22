module mk_advective_m_fluxdiv_module

  use amrex_error_module
  ! use common_namelist_module, only: visc_type

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 2)

    subroutine mk_advective_m_fluxdiv(lo,hi, &
                                      umac, umaclo, umachi, &
                                      vmac, vmaclo, vmachi, &
                                      mx, mxlo, mxhi, &
                                      my, mylo, myhi, &
                                      m_updatex, mudxlo, mudxhi, &
                                      m_updatey, mudylo, mudyhi, &
                                      dx, increment) bind(C,name="mk_advective_m_fluxdiv")

      integer         , intent(in   ) :: lo(2),hi(2)
      integer         , intent(in   ) :: umaclo(2),umachi(2), vmaclo(2),vmachi(2)
      integer         , intent(in   ) :: mxlo(2),mxhi(2), mylo(2),myhi(2)
      integer         , intent(in   ) :: mudxlo(2),mudxhi(2), mudylo(2),mudyhi(2)
      double precision, intent(in   ) ::      umac(umaclo(1):umachi(1),umaclo(2):umachi(2))
      double precision, intent(in   ) ::      vmac(vmaclo(1):vmachi(1),vmaclo(2):vmachi(2))
      double precision, intent(in   ) ::        mx(mxlo(1):mxhi(1),mxlo(2):mxhi(2))
      double precision, intent(in   ) ::        my(mylo(1):myhi(1),mylo(2):myhi(2))
      double precision, intent(inout) :: m_updatex(mudxlo(1):mudxhi(1),mudxlo(2):mudxhi(2))
      double precision, intent(inout) :: m_updatey(mudylo(1):mudyhi(1),mudylo(2):mudyhi(2))
      double precision, intent(in   ) :: dx(2)
      integer         , intent(in   ) :: increment

      ! local
      integer :: i,j

      double precision :: dxinv
      double precision :: fluxx_hi, fluxx_lo, fluxy_hi, fluxy_lo

      dxinv = 1.d0/dx(1)

      !! Note: we are computing Eq (42) in 10.2140/camcos.2014.9.47
      !! Low mach number fluctuating hydrodynamics of diffusively mixing fluids.
      !! A. Donev, A. Nonaka, Y. Sun, T. G. Fai, A. Garcia and J. Bell

      ! computation is redundant; but this is prep for GPU to avoid using extra
      ! memory or moving/allocating temporary storage on the GPU

      !=============================
      ! mx fluxes and divergence
      !=============================
      if (increment==1) then
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            fluxx_hi = 0.25d0*(mx(i,j)+mx(i+1,j))*(umac(i,j)+umac(i+1,j))
            fluxx_lo = 0.25d0*(mx(i-1,j)+mx(i,j))*(umac(i-1,j)+umac(i,j))
            fluxy_hi = 0.25d0*(mx(i,j)+mx(i,j+1))*(vmac(i-1,j+1)+vmac(i,j+1))
            fluxy_lo = 0.25d0*(mx(i,j-1)+mx(i,j))*(vmac(i-1,j)+vmac(i,j))               
            m_updatex(i,j) = m_updatex(i,j) - ( fluxx_hi-fluxx_lo + fluxy_hi-fluxy_lo ) * dxinv
         end do
         end do
      else if (increment==0) then
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            fluxx_hi = 0.25d0*(mx(i,j)+mx(i+1,j))*(umac(i,j)+umac(i+1,j))
            fluxx_lo = 0.25d0*(mx(i-1,j)+mx(i,j))*(umac(i-1,j)+umac(i,j))
            fluxy_hi = 0.25d0*(mx(i,j)+mx(i,j+1))*(vmac(i-1,j+1)+vmac(i,j+1))
            fluxy_lo = 0.25d0*(mx(i,j-1)+mx(i,j))*(vmac(i-1,j)+vmac(i,j))               
            m_updatex(i,j) = -( fluxx_hi-fluxx_lo + fluxy_hi-fluxy_lo ) * dxinv
         end do
         end do
      end if

      !=============================
      ! my fluxes and divergence
      !=============================
      if (increment==1) then
         do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            fluxx_hi = 0.25d0*(my(i,j)+my(i+1,j))*(umac(i+1,j-1)+umac(i+1,j))
            fluxx_lo = 0.25d0*(my(i-1,j)+my(i,j))*(umac(i,j-1)+umac(i,j))
            fluxy_hi = 0.25d0*(my(i,j)+my(i,j+1))*(vmac(i,j)+vmac(i,j+1))
            fluxy_lo = 0.25d0*(my(i,j-1)+my(i,j))*(vmac(i,j-1)+vmac(i,j))
            m_updatey(i,j) = m_updatey(i,j) - ( fluxx_hi-fluxx_lo + fluxy_hi-fluxy_lo ) * dxinv
         end do
         end do
      else if (increment==0) then
         do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            fluxx_hi = 0.25d0*(my(i,j)+my(i+1,j))*(umac(i+1,j-1)+umac(i+1,j))
            fluxx_lo = 0.25d0*(my(i-1,j)+my(i,j))*(umac(i,j-1)+umac(i,j))
            fluxy_hi = 0.25d0*(my(i,j)+my(i,j+1))*(vmac(i,j)+vmac(i,j+1))
            fluxy_lo = 0.25d0*(my(i,j-1)+my(i,j))*(vmac(i,j-1)+vmac(i,j))
            m_updatey(i,j) = -( fluxx_hi-fluxx_lo + fluxy_hi-fluxy_lo ) * dxinv
         end do
         end do
      end if

    end subroutine mk_advective_m_fluxdiv

#endif

#if (AMREX_SPACEDIM == 3)

    subroutine mk_advective_m_fluxdiv(lo,hi, &
                                      umac, umaclo, umachi, &
                                      vmac, vmaclo, vmachi, &
                                      wmac, wmaclo, wmachi, &
                                      mx, mxlo, mxhi, &
                                      my, mylo, myhi, &
                                      mz, mzlo, mzhi, &
                                      m_updatex, mudxlo, mudxhi, &
                                      m_updatey, mudylo, mudyhi, &
                                      m_updatez, mudzlo, mudzhi, &
                                      dx, increment) bind(C,name="mk_advective_m_fluxdiv")

      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: umaclo(3),umachi(3), vmaclo(3),vmachi(3), wmaclo(3),wmachi(3)
      integer         , intent(in   ) :: mxlo(3),mxhi(3), mylo(3),myhi(3), mzlo(3),mzhi(3)
      integer         , intent(in   ) :: mudxlo(3),mudxhi(3), mudylo(3),mudyhi(3), mudzlo(3),mudzhi(3)
      double precision, intent(in   ) ::      umac(umaclo(1):umachi(1),umaclo(2):umachi(2),umaclo(3):umachi(3))
      double precision, intent(in   ) ::      vmac(vmaclo(1):vmachi(1),vmaclo(2):vmachi(2),vmaclo(3):vmachi(3))
      double precision, intent(in   ) ::      wmac(wmaclo(1):wmachi(1),wmaclo(2):wmachi(2),wmaclo(3):wmachi(3))
      double precision, intent(in   ) ::        mx(mxlo(1):mxhi(1),mxlo(2):mxhi(2),mxlo(3):mxhi(3))
      double precision, intent(in   ) ::        my(mylo(1):myhi(1),mylo(2):myhi(2),mylo(3):myhi(3))
      double precision, intent(in   ) ::        mz(mzlo(1):mzhi(1),mzlo(2):mzhi(2),mzlo(3):mzhi(3))
      double precision, intent(inout) :: m_updatex(mudxlo(1):mudxhi(1),mudxlo(2):mudxhi(2),mudxlo(3):mudxhi(3))
      double precision, intent(inout) :: m_updatey(mudylo(1):mudyhi(1),mudylo(2):mudyhi(2),mudylo(3):mudyhi(3))
      double precision, intent(inout) :: m_updatez(mudzlo(1):mudzhi(1),mudzlo(2):mudzhi(2),mudzlo(3):mudzhi(3))
      double precision, intent(in   ) :: dx(3)
      integer         , intent(in   ) :: increment

      ! local
      integer :: i,j,k

      double precision :: dxinv
      double precision :: fluxx_hi, fluxx_lo, fluxy_hi, fluxy_lo, fluxz_hi, fluxz_lo

      dxinv = 1.d0/dx(1)

      !=============================
      ! mx fluxes and divergence
      !=============================
      if (increment==1) then
         do k=lo(3),hi(3)
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            fluxx_hi = 0.25d0*(mx(i,j,k)+mx(i+1,j,k))*(umac(i,j,k)+umac(i+1,j,k))
            fluxx_lo = 0.25d0*(mx(i-1,j,k)+mx(i,j,k))*(umac(i-1,j,k)+umac(i,j,k))
            fluxy_hi = 0.25d0*(mx(i,j,k)+mx(i,j+1,k))*(vmac(i-1,j+1,k)+vmac(i,j+1,k))
            fluxy_lo = 0.25d0*(mx(i,j-1,k)+mx(i,j,k))*(vmac(i-1,j,k)+vmac(i,j,k))
            fluxz_hi = 0.25d0*(mx(i,j,k)+mx(i,j,k+1))*(wmac(i-1,j,k+1)+wmac(i,j,k+1))
            fluxz_lo = 0.25d0*(mx(i,j,k-1)+mx(i,j,k))*(wmac(i-1,j,k)+wmac(i,j,k))
            m_updatex(i,j,k) = m_updatex(i,j,k) - ( fluxx_hi-fluxx_lo + fluxy_hi-fluxy_lo + fluxz_hi-fluxz_lo ) * dxinv
         end do
         end do
         end do
      else if (increment==0) then
         do k=lo(3),hi(3)
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            fluxx_hi = 0.25d0*(mx(i,j,k)+mx(i+1,j,k))*(umac(i,j,k)+umac(i+1,j,k))
            fluxx_lo = 0.25d0*(mx(i-1,j,k)+mx(i,j,k))*(umac(i-1,j,k)+umac(i,j,k))
            fluxy_hi = 0.25d0*(mx(i,j,k)+mx(i,j+1,k))*(vmac(i-1,j+1,k)+vmac(i,j+1,k))
            fluxy_lo = 0.25d0*(mx(i,j-1,k)+mx(i,j,k))*(vmac(i-1,j,k)+vmac(i,j,k))
            fluxz_hi = 0.25d0*(mx(i,j,k)+mx(i,j,k+1))*(wmac(i-1,j,k+1)+wmac(i,j,k+1))
            fluxz_lo = 0.25d0*(mx(i,j,k-1)+mx(i,j,k))*(wmac(i-1,j,k)+wmac(i,j,k))
            m_updatex(i,j,k) = -( fluxx_hi-fluxx_lo + fluxy_hi-fluxy_lo + fluxz_hi-fluxz_lo ) * dxinv
         end do
         end do
         end do
      end if

      !=============================
      ! my fluxes and divergence
      !=============================
      if (increment==1) then
         do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            fluxx_hi = 0.25d0*(my(i,j,k)+my(i+1,j,k))*(umac(i+1,j-1,k)+umac(i+1,j,k))
            fluxx_lo = 0.25d0*(my(i-1,j,k)+my(i,j,k))*(umac(i,j-1,k)+umac(i,j,k))
            fluxy_hi = 0.25d0*(my(i,j,k)+my(i,j+1,k))*(vmac(i,j,k)+vmac(i,j+1,k))
            fluxy_lo = 0.25d0*(my(i,j-1,k)+my(i,j,k))*(vmac(i,j-1,k)+vmac(i,j,k))
            fluxz_hi = 0.25d0*(my(i,j,k)+my(i,j,k+1))*(wmac(i,j-1,k+1)+wmac(i,j,k+1))
            fluxz_lo = 0.25d0*(my(i,j,k-1)+my(i,j,k))*(wmac(i,j-1,k)+wmac(i,j,k))
            m_updatey(i,j,k) = m_updatey(i,j,k) - ( fluxx_hi-fluxx_lo + fluxy_hi-fluxy_lo + fluxz_hi-fluxz_lo ) * dxinv
         end do
         end do
         end do
      else if (increment==0) then
         do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            fluxx_hi = 0.25d0*(my(i,j,k)+my(i+1,j,k))*(umac(i+1,j-1,k)+umac(i+1,j,k))
            fluxx_lo = 0.25d0*(my(i-1,j,k)+my(i,j,k))*(umac(i,j-1,k)+umac(i,j,k))
            fluxy_hi = 0.25d0*(my(i,j,k)+my(i,j+1,k))*(vmac(i,j,k)+vmac(i,j+1,k))
            fluxy_lo = 0.25d0*(my(i,j-1,k)+my(i,j,k))*(vmac(i,j-1,k)+vmac(i,j,k))
            fluxz_hi = 0.25d0*(my(i,j,k)+my(i,j,k+1))*(wmac(i,j-1,k+1)+wmac(i,j,k+1))
            fluxz_lo = 0.25d0*(my(i,j,k-1)+my(i,j,k))*(wmac(i,j-1,k)+wmac(i,j,k))
            m_updatey(i,j,k) = -( fluxx_hi-fluxx_lo + fluxy_hi-fluxy_lo + fluxz_hi-fluxz_lo ) * dxinv
         end do
         end do
         end do
      end if

      !=============================
      ! mz fluxes and divergence
      !=============================
      if (increment==1) then
         do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            fluxx_hi = 0.25d0*(mz(i,j,k)+mz(i+1,j,k))*(umac(i+1,j,k-1)+umac(i+1,j,k))
            fluxx_lo = 0.25d0*(mz(i-1,j,k)+mz(i,j,k))*(umac(i,j,k-1)+umac(i,j,k))
            fluxy_hi = 0.25d0*(mz(i,j,k)+mz(i,j+1,k))*(vmac(i,j+1,k-1)+vmac(i,j+1,k))
            fluxy_lo = 0.25d0*(mz(i,j-1,k)+mz(i,j,k))*(vmac(i,j,k-1)+vmac(i,j,k))
            fluxz_hi = 0.25d0*(mz(i,j,k)+mz(i,j,k+1))*(wmac(i,j,k)+wmac(i,j,k+1))
            fluxz_lo = 0.25d0*(mz(i,j,k-1)+mz(i,j,k))*(wmac(i,j,k-1)+wmac(i,j,k))
            m_updatez(i,j,k) = m_updatez(i,j,k) - ( fluxx_hi-fluxx_lo + fluxy_hi-fluxy_lo + fluxz_hi-fluxz_lo ) * dxinv
         end do
         end do
         end do
      else if (increment==0) then
         do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            fluxx_hi = 0.25d0*(mz(i,j,k)+mz(i+1,j,k))*(umac(i+1,j,k-1)+umac(i+1,j,k))
            fluxx_lo = 0.25d0*(mz(i-1,j,k)+mz(i,j,k))*(umac(i,j,k-1)+umac(i,j,k))
            fluxy_hi = 0.25d0*(mz(i,j,k)+mz(i,j+1,k))*(vmac(i,j+1,k-1)+vmac(i,j+1,k))
            fluxy_lo = 0.25d0*(mz(i,j-1,k)+mz(i,j,k))*(vmac(i,j,k-1)+vmac(i,j,k))
            fluxz_hi = 0.25d0*(mz(i,j,k)+mz(i,j,k+1))*(wmac(i,j,k)+wmac(i,j,k+1))
            fluxz_lo = 0.25d0*(mz(i,j,k-1)+mz(i,j,k))*(wmac(i,j,k-1)+wmac(i,j,k))
            m_updatez(i,j,k) = -( fluxx_hi-fluxx_lo + fluxy_hi-fluxy_lo + fluxz_hi-fluxz_lo ) * dxinv
         end do
         end do
         end do
      end if

    end subroutine mk_advective_m_fluxdiv

#endif

end module mk_advective_m_fluxdiv_module
