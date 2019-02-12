module correction_flux_module

  use amrex_error_module
  use common_namelist_module
  use multispec_namelist_module

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 2)

  subroutine correction_flux(lo,hi, &
       rho, rlo, rhi, nspecies, &
       rhotot, rtlo, rthi, &
       flux_x, fxlo, fxhi, &
       flux_y, fylo, fyhi) bind(C,name="correction_flux")

    integer         , intent(in   ) :: lo(2),hi(2), nspecies
    integer         , intent(in   ) :: rlo(2),rhi(2), rtlo(2),rthi(2)
    integer         , intent(in   ) :: fxlo(2),fxhi(2), fylo(2),fyhi(2)
    double precision, intent(in   ) :: rho(rlo(1):rhi(1),rlo(2):rhi(2),nspecies)
    double precision, intent(in   ) :: rhotot(rtlo(1):rthi(1),rtlo(2):rthi(2))
    double precision, intent(inout) :: flux_x(fxlo(1):fxhi(1),fxlo(2):fxhi(2),nspecies)
    double precision, intent(inout) :: flux_y(fylo(1):fyhi(1),fylo(2):fyhi(2),nspecies)

    ! local
    integer          :: i,j,n
    double precision :: sumx,sumy,corr,total_corr

    ! x-faces
    total_corr=0.0d0
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1

          ! free the data
          sumx = 0.d0
          corr = 0.d0

          ! sum the x-fluxes upto nspecies-1 
          do n=1, nspecies-1
             sumx = sumx + flux_x(i,j,n)
          end do

          ! caculate corr and print error if not zero 
          corr = flux_x(i,j,nspecies) + sumx

          ! correct x-flux for last species  
          flux_x(i,j,nspecies) = -sumx     
          total_corr = total_corr + abs(corr)

       end do
    end do

    ! y-faces
    total_corr=0.0d0
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)

          ! free the data
          sumy  = 0.d0
          corr = 0.d0

          ! sum the y-fluxes upto nspecies-1 
          do n=1, nspecies-1
             sumy = sumy + flux_y(i,j,n)
          end do

          ! caculate corr and print error if not zero 
          corr = flux_y(i,j,nspecies) + sumy

          ! correct y-flux for last species  
          flux_y(i,j,nspecies) = -sumy             
          total_corr = total_corr + abs(corr)

       end do
    end do

  end subroutine correction_flux

#endif

#if (AMREX_SPACEDIM == 3)

  subroutine correction_flux(lo,hi, &
       rho, rlo, rhi, nspecies, &
       rhotot, rtlo, rthi, &
       flux_x, fxlo, fxhi, &
       flux_y, fylo, fyhi, &
       flux_z, fzlo, fzhi) bind(C,name="correction_flux")

    integer         , intent(in   ) :: lo(3),hi(3), nspecies
    integer         , intent(in   ) :: rlo(3),rhi(3), rtlo(3),rthi(3)
    integer         , intent(in   ) :: fxlo(3),fxhi(3), fylo(3),fyhi(3), fzlo(3),fzhi(3)
    double precision, intent(in   ) :: rho(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3),nspecies)
    double precision, intent(in   ) :: rhotot(rtlo(1):rthi(1),rtlo(2):rthi(2))
    double precision, intent(inout) :: flux_x(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),nspecies)
    double precision, intent(inout) :: flux_y(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),nspecies)
    double precision, intent(inout) :: flux_z(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),nspecies)

    ! local
    integer         :: i,j,k,n
    double precision :: sumx,sumy,sumz,corr,total_corr

    ! x-faces
    total_corr=0.0d0
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1

             ! free the data
             sumx = 0.d0
             corr = 0.d0

             ! sum the x-fluxes upto nspecies-1 
             do n=1, nspecies-1
                sumx = sumx + flux_x(i,j,k,n)
             end do

             ! caculate corr and print error if not zero 
             corr = flux_x(i,j,k,nspecies) + sumx 

             ! correct x-flux for last species  
             flux_x(i,j,k,nspecies) = -sumx             
             total_corr = total_corr + abs(corr) 

          end do
       end do
    end do

    ! y-faces
    total_corr=0.0d0
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)

             ! free the data
             sumy  = 0.d0
             corr = 0.d0

             ! sum the y-fluxes upto nspecies-1 
             do n=1, nspecies-1
                sumy = sumy + flux_y(i,j,k,n)
             end do

             ! caculate corr and print error if not zero 
             corr = flux_y(i,j,k,nspecies) + sumy 

             ! correct y-flux for last species  
             flux_y(i,j,k,nspecies) = -sumy             
             total_corr = total_corr + abs(corr) 

          end do
       end do
    end do

    ! z-faces
    total_corr=0.0d0
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! free the data
             sumz  = 0.d0
             corr = 0.d0

             ! sum the z-fluxes upto nspecies-1 
             do n=1, nspecies-1
                sumz = sumz + flux_z(i,j,k,n)
             end do

             ! caculate corr and print error if not zero 
             corr = flux_z(i,j,k,nspecies) + sumz 

             ! correct z-flux for last species  
             flux_z(i,j,k,nspecies) = -sumz             
             total_corr = total_corr + abs(corr) 

          end do
       end do
    end do

  end subroutine correction_flux

#endif

end module correction_flux_module
