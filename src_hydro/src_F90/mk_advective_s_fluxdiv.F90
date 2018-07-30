module mk_advective_s_fluxdiv_module

  use amrex_error_module
  use common_namelist_module, only: visc_coef, fixed_dt

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 2)

    subroutine mk_advective_s_fluxdiv(lo,hi, &
                                      umac, umaclo, umachi, &
                                      vmac, vmaclo, vmachi, &
                                      mx, mxlo, mxhi, &
                                      my, mylo, myhi, &
                                      m, mlo, mhi, &
                                      m_update, mudlo, mudhi, &
                                      dx, increment) bind(C,name="mk_advective_s_fluxdiv")

      integer         , intent(in   ) :: lo(2),hi(2)
      integer         , intent(in   ) :: umaclo(2),umachi(2), vmaclo(2),vmachi(2)
      integer         , intent(in   ) :: mlo(2),mhi(2)
      integer         , intent(in   ) :: mxlo(2),mxhi(2), mylo(2),myhi(2)
      integer         , intent(in   ) :: mudlo(2),mudhi(2)
      double precision, intent(in   ) ::      umac(umaclo(1):umachi(1),umaclo(2):umachi(2))
      double precision, intent(in   ) ::      vmac(vmaclo(1):vmachi(1),vmaclo(2):vmachi(2))
      double precision, intent(in   ) ::        mx(mxlo(1):mxhi(1),mxlo(2):mxhi(2))
      double precision, intent(in   ) ::        my(mylo(1):myhi(1),mylo(2):myhi(2))
      double precision, intent(in   ) ::         m(mlo(1):mhi(1),mlo(2):mhi(2))
      double precision, intent(inout) :: m_update(mudlo(1):mudhi(1),mudlo(2):mudhi(2))
      double precision, intent(in   ) :: dx(2)
      integer         , intent(in   ) :: increment

      ! local
      integer :: i,j

      double precision :: m_fluxx
      double precision :: m_fluxy

      double precision :: dxinv

      double precision :: tracer_visc_coef

      dxinv = 1.d0/dx(1)
      
      !! Coefficient of viscous term to reduce cell-reynolds number:
      ! tracer_visc_coef = 0.05d0*(0.5d0*dx(1)*dx(1)/(AMREX_SPACEDIM*fixed_dt))

      !=============================
      ! fluxes and divergence
      !=============================
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
         m_fluxx = dxinv*(mx(i,j)*umac(i,j) - mx(i-1,j)*umac(i-1,j))
         m_fluxy = dxinv*(my(i,j)*vmac(i,j) - my(i,j-1)*vmac(i,j-1))

         m_update(i,j) = -( m_fluxx + m_fluxy )
         
         !! Add viscous term
         ! m_update(i,j) = m_update(i,j) &
         !      + tracer_visc_coef*(dxinv**2)*(m(i-1,j)+m(i+1,j)+m(i,j-1)+m(i,j+1)-4.0d0*m(i,j))

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
                                      m, mlo, mhi, &
                                      m_update, mudlo, mudhi, &
                                      dx, increment) bind(C,name="mk_advective_s_fluxdiv")

      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: umaclo(3),umachi(3), vmaclo(3),vmachi(3), wmaclo(3),wmachi(3)
      integer         , intent(in   ) :: mxlo(3),mxhi(3), mylo(3),myhi(3), mzlo(3),mzhi(3)
      integer         , intent(in   ) :: mlo(3),mhi(3)
      integer         , intent(in   ) :: mudlo(3),mudhi(3)
      double precision, intent(in   ) ::      umac(umaclo(1):umachi(1),umaclo(2):umachi(2),umaclo(3):umachi(3))
      double precision, intent(in   ) ::      vmac(vmaclo(1):vmachi(1),vmaclo(2):vmachi(2),vmaclo(3):vmachi(3))
      double precision, intent(in   ) ::      wmac(wmaclo(1):wmachi(1),wmaclo(2):wmachi(2),wmaclo(3):wmachi(3))
      double precision, intent(in   ) ::        mx(mxlo(1):mxhi(1),mxlo(2):mxhi(2),mxlo(3):mxhi(3))
      double precision, intent(in   ) ::        my(mylo(1):myhi(1),mylo(2):myhi(2),mylo(3):myhi(3))
      double precision, intent(in   ) ::        mz(mzlo(1):mzhi(1),mzlo(2):mzhi(2),mzlo(3):mzhi(3))
      double precision, intent(in   ) ::        m(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
      double precision, intent(inout) :: m_update(mudlo(1):mudhi(1),mudlo(2):mudhi(2),mudlo(3):mudhi(3))
      double precision, intent(in   ) :: dx(3)
      integer         , intent(in   ) :: increment

      ! local
      integer :: i,j,k

      double precision :: m_fluxx
      double precision :: m_fluxy
      double precision :: m_fluxz

      double precision :: dxinv

      dxinv = 1.d0/dx(1)

      !=============================
      ! fluxes and divergence
      !=============================
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
         m_fluxx = dxinv*(mx(i  ,j  ,k  )*umac(i  ,j  ,k  ) - mx(i-1,j  ,k  )*umac(i-1,j  ,k  ))
         m_fluxy = dxinv*(my(i  ,j  ,k  )*vmac(i  ,j  ,k  ) - my(i  ,j-1,k  )*vmac(i  ,j-1,k  ))
         m_fluxz = dxinv*(mz(i  ,j  ,k  )*wmac(i  ,j  ,k  ) - mz(i  ,j  ,k-1)*wmac(i  ,j  ,k-1))

         m_update(i,j) = -( m_fluxx + m_fluxy + m_fluxz)
      end do
      end do
      end do

    end subroutine mk_advective_s_fluxdiv
    
#endif

end module mk_advective_s_fluxdiv_module
