module mk_advective_s_fluxdiv_module

  use amrex_error_module
  ! use common_namelist_module, only: visc_type

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 2)

    subroutine mk_advective_s_fluxdiv(lo,hi, &
                                      umac, umaclo, umachi, &
                                      vmac, vmaclo, vmachi, &
                                      mx, mxlo, mxhi, &
                                      my, mylo, myhi, &
                                      m_update, mudlo, mudhi, &
                                      dx, increment) bind(C,name="mk_advective_s_fluxdiv")

      integer         , intent(in   ) :: lo(2),hi(2)
      integer         , intent(in   ) :: umaclo(2),umachi(2), vmaclo(2),vmachi(2)
      integer         , intent(in   ) :: mxlo(2),mxhi(2), mylo(2),myhi(2)
      integer         , intent(in   ) :: mudlo(2),mudhi(2)
      double precision, intent(in   ) ::      umac(umaclo(1):umachi(1),umaclo(2):umachi(2))
      double precision, intent(in   ) ::      vmac(vmaclo(1):vmachi(1),vmaclo(2):vmachi(2))
      double precision, intent(in   ) ::        mx(mxlo(1):mxhi(1),mxlo(2):mxhi(2))
      double precision, intent(in   ) ::        my(mylo(1):myhi(1),mylo(2):myhi(2))
      double precision, intent(inout) :: m_update(mudlo(1):mudhi(1),mudlo(2):mudhi(2))
      double precision, intent(in   ) :: dx(2)
      integer         , intent(in   ) :: increment

      ! local
      integer :: i,j

      double precision :: m_fluxx
      double precision :: m_fluxy

      double precision :: dxinv

      dxinv = 1.d0/dx(1)

      !=============================
      ! fluxes and divergence
      !=============================
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
         m_fluxx = dxinv*(mx(i,j)*umac(i,j) - mx(i-1,j)*umac(i-1,j))
         m_fluxy = dxinv*(my(i,j)*vmac(i,j) - my(i,j-1)*vmac(i,j-1))
         m_update(i,j) = -( m_fluxx + m_fluxy )
      end do
      end do

    end subroutine mk_advective_s_fluxdiv

#endif

#if (AMREX_SPACEDIM == 3)

    subroutine mk_advective_s_fluxdiv(lo,hi, &
                                      umac, umaclo, umachi, &
                                      vmac, vmaclo, vmachi, &
                                      wmac, wmaclo, wmachi, &
                                      mx, mxlo, mxhi, &
                                      my, mylo, myhi, &
                                      mz, mzlo, mzhi, &
                                      m_updatex, mudxlo, mudxhi, &
                                      m_updatey, mudylo, mudyhi, &
                                      m_updatez, mudzlo, mudzhi, &
                                      dx, increment) bind(C,name="mk_advective_s_fluxdiv")

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

      double precision ::   mx_fluxx(lo(1):hi(1)+2,lo(2):hi(2)  ,lo(3):hi(3)  )
      double precision ::   mx_fluxy(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
      double precision ::   mx_fluxz(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)

      double precision ::   my_fluxx(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
      double precision ::   my_fluxy(lo(1):hi(1)  ,lo(2):hi(2)+2,lo(3):hi(3)  )
      double precision ::   my_fluxz(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)

      double precision ::   mz_fluxx(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
      double precision ::   mz_fluxy(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
      double precision ::   mz_fluxz(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+2)

      double precision :: dxinv

      dxinv = 1.d0/dx(1)

      !=============================
      ! mx fluxes and divergence
      !=============================
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)+2
         mx_fluxx(i,j,k) = 0.25d0*(mx(i-1,j,k)+mx(i,j,k))*(umac(i-1,j,k)+umac(i,j,k))
      end do
      end do
      end do

      do k=lo(3),hi(3)+1
      do j=lo(2),hi(2)+1
      do i=lo(1),hi(1)+1
         mx_fluxy(i,j,k) = 0.25d0*(mx(i,j-1,k)+mx(i,j,k))*(vmac(i-1,j,k)+vmac(i,j,k))
         mx_fluxz(i,j,k) = 0.25d0*(mx(i,j,k-1)+mx(i,j,k))*(wmac(i-1,j,k)+wmac(i,j,k))
      end do
      end do
      end do

      if (increment==1) then
         do k=lo(3),hi(3)
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            m_updatex(i,j,k) = m_updatex(i,j,k) - ( (mx_fluxx(i+1,j,k)-mx_fluxx(i,j,k)) * dxinv + &
                                                    (mx_fluxy(i,j+1,k)-mx_fluxy(i,j,k)) * dxinv + &
                                                    (mx_fluxz(i,j,k+1)-mx_fluxz(i,j,k)) * dxinv )
         end do
         end do
         end do
      else if (increment==0) then
         do k=lo(3),hi(3)
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            m_updatex(i,j,k) = -( (mx_fluxx(i+1,j,k)-mx_fluxx(i,j,k)) * dxinv + &
                                  (mx_fluxy(i,j+1,k)-mx_fluxy(i,j,k)) * dxinv + &
                                  (mx_fluxz(i,j,k+1)-mx_fluxz(i,j,k)) * dxinv )
         end do
         end do
         end do
      end if

      !=============================
      ! my fluxes and divergence
      !=============================
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)+2
      do i=lo(1),hi(1)
         my_fluxy(i,j,k) = 0.25d0*(my(i,j-1,k)+my(i,j,k))*(vmac(i,j-1,k)+vmac(i,j,k))
      end do
      end do
      end do

      do k=lo(3),hi(3)+1
      do j=lo(2),hi(2)+1
      do i=lo(1),hi(1)+1
         my_fluxx(i,j,k) = 0.25d0*(my(i-1,j,k)+my(i,j,k))*(umac(i,j-1,k)+umac(i,j,k))
         my_fluxz(i,j,k) = 0.25d0*(my(i,j,k-1)+my(i,j,k))*(wmac(i,j-1,k)+wmac(i,j,k))
      end do
      end do
      end do

      if (increment==1) then
         do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            m_updatey(i,j,k) = m_updatey(i,j,k) - ( (my_fluxx(i+1,j,k)-my_fluxx(i,j,k)) * dxinv + &
                                                    (my_fluxy(i,j+1,k)-my_fluxy(i,j,k)) * dxinv + &
                                                    (my_fluxz(i,j,k+1)-my_fluxz(i,j,k)) * dxinv )
         end do
         end do
         end do
      else if (increment==0) then
         do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            m_updatey(i,j,k) = -( (my_fluxx(i+1,j,k)-my_fluxx(i,j,k)) * dxinv + &
                                  (my_fluxy(i,j+1,k)-my_fluxy(i,j,k)) * dxinv + &
                                  (my_fluxz(i,j,k+1)-my_fluxz(i,j,k)) * dxinv )
         end do
         end do
         end do
      end if

      !=============================
      ! mz fluxes and divergence
      !=============================
      do k=lo(3),hi(3)+2
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
         mz_fluxz(i,j,k) = 0.25d0*(mz(i,j,k-1)+mz(i,j,k))*(wmac(i,j,k-1)+wmac(i,j,k))
      end do
      end do
      end do

      do k=lo(3),hi(3)+1
      do j=lo(2),hi(2)+1
      do i=lo(1),hi(1)+1
         mz_fluxx(i,j,k) = 0.25d0*(mz(i-1,j,k)+mz(i,j,k))*(umac(i,j,k-1)+umac(i,j,k))
         mz_fluxy(i,j,k) = 0.25d0*(mz(i,j-1,k)+mz(i,j,k))*(vmac(i,j,k-1)+vmac(i,j,k))
      end do
      end do
      end do

      if (increment==1) then
         do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            m_updatez(i,j,k) = m_updatez(i,j,k) - ( (mz_fluxx(i+1,j,k)-mz_fluxx(i,j,k)) * dxinv + &
                                                    (mz_fluxy(i,j+1,k)-mz_fluxy(i,j,k)) * dxinv + &
                                                    (mz_fluxz(i,j,k+1)-mz_fluxz(i,j,k)) * dxinv )
         end do
         end do
         end do
      else if (increment==0) then
         do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            m_updatez(i,j,k) = -( (mz_fluxx(i+1,j,k)-mz_fluxx(i,j,k)) * dxinv + &
                                  (mz_fluxy(i,j+1,k)-mz_fluxy(i,j,k)) * dxinv + &
                                  (mz_fluxz(i,j,k+1)-mz_fluxz(i,j,k)) * dxinv )
         end do
         end do
         end do
      end if

    end subroutine mk_advective_s_fluxdiv
    
#endif

end module mk_advective_s_fluxdiv_module
