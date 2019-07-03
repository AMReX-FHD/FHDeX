module flux_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, molmass, cell_depth, k_b, runiv, bc_lo, bc_hi, n_cells, membrane_cell, visc_type, algorithm_type

  implicit none

  private

  public :: stoch_flux, stoch_flux_BC

contains


  subroutine stoch_flux(lo,hi, cons, fluxx, fluxy, &
#if (AMREX_SPACEDIM == 3)
       fluxz, &
#endif
       ranfluxx, ranfluxy, &
#if (AMREX_SPACEDIM == 3)
       ranfluxz, &
#endif
       dx, dt) bind(C,name="stoch_flux")

    integer         , intent(in   ) :: lo(3),hi(3)
    real(amrex_real), intent(in   ) :: dx(3), dt
    real(amrex_real), intent(inout) :: fluxx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
    real(amrex_real), intent(inout) :: fluxy(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
    real(amrex_real), intent(inout) :: fluxz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif

    real(amrex_real), intent(inout) :: ranfluxx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
    real(amrex_real), intent(inout) :: ranfluxy(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
    real(amrex_real), intent(inout) :: ranfluxz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif

    real(amrex_real), intent(in   ) :: cons(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)


    real(amrex_real) :: volinv, dtinv


    integer :: i,j,k,l
    integer :: ll, ns

    dtinv = 1d0/dt
#if (AMREX_SPACEDIM == 3)
    volinv = 1d0/(dx(1)*dx(2)*dx(3))
#endif

#if (AMREX_SPACEDIM == 2)
    volinv = 1d0/(dx(1)*dx(2)*cell_depth)
#endif

    if (abs(visc_type) .gt. 1) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! JB's tensor form !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!! x-flux !!!!!!!!!!!!!!!!!!!

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)+1




             end do
          end do
       end do

       ! print*, "Hack: got here (end) stochflux"

!!!!!!!!!!!!!!!!!!! y-flux !!!!!!!!!!!!!!!!!!!

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)

 
             end do
          end do
       end do

!!!!!!!!!!!!!!!!!!! z-flux !!!!!!!!!!!!!!!!!!!

       do k = lo(3),hi(3)+1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

 

             end do
          end do
       end do

    else

    endif

!!!!!!!!!!!!! Enforce flux boundary conditions !!!!!!!!!!!!!

    call stoch_flux_BC(lo,hi, fluxx, fluxy, &
#if (AMREX_SPACEDIM == 3)
         fluxz &
#endif
         )

  end subroutine stoch_flux

  subroutine stoch_flux_BC(lo,hi, xflux, yflux, &
#if (AMREX_SPACEDIM == 3)
       zflux &
#endif
       ) bind(C,name="stoch_flux_bounds")

    integer         , intent(in   ) :: lo(3),hi(3)
    real(amrex_real), intent(inout) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3), nvars)
    real(amrex_real), intent(inout) :: yflux(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3), nvars)
#if (AMREX_SPACEDIM == 3)
    real(amrex_real), intent(inout) :: zflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1, nvars)
#endif

    real(amrex_real) :: sqrtTwo

    integer :: i,j,k,l

    sqrtTwo = sqrt(2.0)

!!!!!!!!!!!!!! x-flux BCs !!!!!!!!!!!!!!

    !if on lower bound and specular
    if((lo(1) .eq. 0) .and. (bc_lo(1) .eq. 1)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)

                xflux(0,j,k,l) = 0        

             end do
          end do
       end do
    endif

    !if on upper bound and specular
    if((hi(1) .eq. n_cells(1)-1) .and. (bc_hi(1) .eq. 1)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)

                xflux(hi(1)+1,j,k,l) = 0           

             end do
          end do
       end do
    endif

    !if on lower bound and diff
    if((lo(1) .eq. 0) .and. (bc_lo(1) .eq. 2)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)

                xflux(0,j,k,l) = sqrtTwo*xflux(0,j,k,l)

             end do
          end do
       end do
    endif

    !if on upper bound and diff
    if((hi(1) .eq. n_cells(1)-1) .and. (bc_hi(1) .eq. 2)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)

                xflux(hi(1)+1,j,k,l) = sqrtTwo*xflux(hi(1)+1,j,k,l)        

             end do
          end do
       end do
    endif

!!!!!!!!!!!!!! y-flux BCs !!!!!!!!!!!!!!

    !if on lower bound and specular
    if((lo(2) .eq. 0) .and. (bc_lo(2) .eq. 1)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)

                yflux(i,0,k,l) = 0        

             end do
          end do
       end do
    endif

    !if on upper bound and specular
    if((hi(2) .eq. n_cells(2)-1) .and. (bc_hi(2) .eq. 1)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)

                yflux(i,hi(2)+1,k,l) = 0           

             end do
          end do
       end do
    endif

    !if on lower bound and diff
    if((lo(2) .eq. 0) .and. (bc_lo(2) .eq. 2)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)

                yflux(i,0,k,l) = sqrtTwo*yflux(i,0,k,l)

             end do
          end do
       end do
    endif

    !if on upper bound and diff
    if((hi(2) .eq. n_cells(2)-1) .and. (bc_hi(2) .eq. 2)) then
       do l = 2,nvars
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)

                yflux(i,hi(2)+1,k,l) = sqrtTwo*yflux(i,hi(2)+1,k,l)        

             end do
          end do
       end do
    endif

!!!!!!!!!!!!!! z-flux BCs !!!!!!!!!!!!!!

    !if on lower bound and specular
    if((lo(3) .eq. 0) .and. (bc_lo(3) .eq. 1)) then
       do l = 2,nvars
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                zflux(i,j,0,l) = 0

             end do
          end do
       end do
    endif

    !if on upper bound and specular
    if((hi(3) .eq. n_cells(3)-1) .and. (bc_hi(3) .eq. 1)) then
       do l = 2,nvars
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                zflux(i,j,hi(3)+1,l) = 0           

             end do
          end do
       end do
    endif

    !if on lower bound and diff
    if((lo(3) .eq. 0) .and. (bc_lo(3) .eq. 2)) then
       do l = 2,nvars
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                zflux(i,j,0,l) = sqrtTwo*zflux(i,j,0,l)

             end do
          end do
       end do
    endif

    !if on upper bound and diff
    if((hi(3) .eq. n_cells(3)-1) .and. (bc_hi(3) .eq. 2)) then
       do l = 2,nvars
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                zflux(i,j,hi(3)+1,l) = sqrtTwo*zflux(i,j,hi(3)+1,l)        

             end do
          end do
       end do
    endif


  end subroutine stoch_flux_BC

end module flux_module

































