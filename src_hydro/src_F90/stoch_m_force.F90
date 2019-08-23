module stoch_m_force_module

  use amrex_error_module
  ! use common_namelist_module, only: visc_type

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 2)

      subroutine mult_by_sqrt_eta_temp(lo,hi, &
                                  mflux_cc, mfluxcclo, mfluxcchi, ncomp_cc, &
                                  mflux_nd, mfluxndlo, mfluxndhi, ncomp_nd, &
                                  eta, etacclo, etacchi, &
                                  eta_nodal, etandlo, etandhi, &
                                  temperature, tempcclo, tempcchi, &
                                  temperature_nodal, tempndlo, tempndhi &
                                  ) bind(C,name="mult_by_sqrt_eta_temp")

      integer         , intent(in   ) :: lo(2),hi(2)
      integer         , intent(in   ) :: ncomp_cc, ncomp_nd
      integer         , intent(in   ) :: mfluxcclo(2),mfluxcchi(2), mfluxndlo(2),mfluxndhi(2)
      integer         , intent(in   ) :: etacclo(2),etacchi(2), etandlo(2),etandhi(2)
      integer         , intent(in   ) :: tempcclo(2),tempcchi(2), tempndlo(2),tempndhi(2)
      double precision, intent(inout) :: mflux_cc(mfluxcclo(1):mfluxcchi(1),mfluxcclo(2):mfluxcchi(2),ncomp_cc)
      double precision, intent(inout) :: mflux_nd(mfluxndlo(1):mfluxndhi(1),mfluxndlo(2):mfluxndhi(2),ncomp_nd)
      double precision, intent(in   ) :: eta(      etacclo(1):etacchi(1),etacclo(2):etacchi(2))
      double precision, intent(in   ) :: eta_nodal(etandlo(1):etandhi(1),etandlo(2):etandhi(2))
      double precision, intent(in   ) :: temperature(      tempcclo(1):tempcchi(1),tempcclo(2):tempcchi(2))
      double precision, intent(in   ) :: temperature_nodal(tempndlo(1):tempndhi(1),tempndlo(2):tempndhi(2))

      integer :: i,j

      do j=lo(2)-1,hi(2)+1
         do i=lo(1)-1,hi(1)+1
            mflux_cc(i,j,:) = mflux_cc(i,j,:) * sqrt(eta(i,j)*temperature(i,j))
         end do
      end do

      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)+1
            mflux_nd(i,j,:) = mflux_nd(i,j,:) * sqrt(eta_nodal(i,j)*temperature_nodal(i,j))
         end do
      end do

    end subroutine mult_by_sqrt_eta_temp

#endif

#if (AMREX_SPACEDIM == 3)

      subroutine mult_by_sqrt_eta_temp(lo,hi, &
                                  mflux_cc, mfluxcclo, mfluxcchi, ncomp_cc, &
                                  mflux_xy, mfluxxylo, mfluxxyhi, ncomp_xy, &
                                  mflux_xz, mfluxxzlo, mfluxxzhi, ncomp_xz, &
                                  mflux_yz, mfluxyzlo, mfluxyzhi, ncomp_yz, &
                                  eta,    etacclo, etacchi, &
                                  eta_xy, etaxylo, etaxyhi, &
                                  eta_xz, etaxzlo, etaxzhi, &
                                  eta_yz, etayzlo, etayzhi, &
                                  temperature,    tempcclo, tempcchi, &
                                  temperature_xy, tempxylo, tempxyhi, &
                                  temperature_xz, tempxzlo, tempxzhi, &
                                  temperature_yz, tempyzlo, tempyzhi &
                                  ) bind(C,name="mult_by_sqrt_eta_temp")

      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: ncomp_cc, ncomp_xy, ncomp_xz, ncomp_yz
      integer         , intent(in   ) :: mfluxcclo(3),mfluxcchi(3)
      integer         , intent(in   ) :: mfluxxylo(3),mfluxxyhi(3), mfluxxzlo(3),mfluxxzhi(3), mfluxyzlo(3),mfluxyzhi(3)
      integer         , intent(in   ) :: etacclo(3),etacchi(3)
      integer         , intent(in   ) :: etaxylo(3),etaxyhi(3), etaxzlo(3),etaxzhi(3), etayzlo(3),etayzhi(3)
      integer         , intent(in   ) :: tempcclo(3),tempcchi(3)
      integer         , intent(in   ) :: tempxylo(3),tempxyhi(3), tempxzlo(3),tempxzhi(3), tempyzlo(3),tempyzhi(3)
      double precision, intent(inout) :: mflux_cc(mfluxcclo(1):mfluxcchi(1),mfluxcclo(2):mfluxcchi(2),mfluxcclo(3):mfluxcchi(3),ncomp_cc)
      double precision, intent(inout) :: mflux_xy(mfluxxylo(1):mfluxxyhi(1),mfluxxylo(2):mfluxxyhi(2),mfluxxylo(3):mfluxxyhi(3),ncomp_xy)
      double precision, intent(inout) :: mflux_xz(mfluxxzlo(1):mfluxxzhi(1),mfluxxzlo(2):mfluxxzhi(2),mfluxxzlo(3):mfluxxzhi(3),ncomp_xz)
      double precision, intent(inout) :: mflux_yz(mfluxyzlo(1):mfluxyzhi(1),mfluxyzlo(2):mfluxyzhi(2),mfluxyzlo(3):mfluxyzhi(3),ncomp_yz)
      double precision, intent(in   ) :: eta(   etacclo(1):etacchi(1),etacclo(2):etacchi(2),etacclo(3):etacchi(3))
      double precision, intent(in   ) :: eta_xy(etaxylo(1):etaxyhi(1),etaxylo(2):etaxyhi(2),etaxylo(3):etaxyhi(3))
      double precision, intent(in   ) :: eta_xz(etaxzlo(1):etaxzhi(1),etaxzlo(2):etaxzhi(2),etaxzlo(3):etaxzhi(3))
      double precision, intent(in   ) :: eta_yz(etayzlo(1):etayzhi(1),etayzlo(2):etayzhi(2),etayzlo(3):etayzhi(3))
      double precision, intent(in   ) :: temperature(   tempcclo(1):tempcchi(1),tempcclo(2):tempcchi(2),tempcclo(3):tempcchi(3))
      double precision, intent(in   ) :: temperature_xy(tempxylo(1):tempxyhi(1),tempxylo(2):tempxyhi(2),tempxylo(3):tempxyhi(3))
      double precision, intent(in   ) :: temperature_xz(tempxzlo(1):tempxzhi(1),tempxzlo(2):tempxzhi(2),tempxzlo(3):tempxzhi(3))
      double precision, intent(in   ) :: temperature_yz(tempyzlo(1):tempyzhi(1),tempyzlo(2):tempyzhi(2),tempyzlo(3):tempyzhi(3))

      integer :: i,j,k

      do k=lo(3)-1,hi(3)+1
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               mflux_cc(i,j,k,:) = mflux_cc(i,j,k,:) * sqrt(eta(i,j,k)*temperature(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)+1
               mflux_xy(i,j,k,:) = mflux_xy(i,j,k,:) * sqrt(eta_xy(i,j,k)*temperature_xy(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
               mflux_xz(i,j,k,:) = mflux_xz(i,j,k,:) * sqrt(eta_xz(i,j,k)*temperature_xz(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
               mflux_yz(i,j,k,:) = mflux_yz(i,j,k,:) * sqrt(eta_yz(i,j,k)*temperature_yz(i,j,k))
            end do
         end do
      end do

    end subroutine mult_by_sqrt_eta_temp
    
#endif

end module stoch_m_force_module
