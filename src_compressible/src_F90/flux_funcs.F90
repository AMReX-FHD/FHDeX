module flux_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, molmass, &
                                     cell_depth, k_b, runiv, n_cells, membrane_cell, &
                                     visc_type, algorithm_type, &
                                     bc_mass_lo, bc_mass_hi, bc_therm_lo, bc_therm_hi, &
                                     bc_vel_lo, bc_vel_hi
  implicit none

  private

  public :: stoch_flux_BC

contains

  subroutine stoch_flux(lo,hi, cons, prim, fluxx, fluxy, &
#if (AMREX_SPACEDIM == 3)
                        fluxz, &
#endif
                        ranfluxx, ranfluxy, &
#if (AMREX_SPACEDIM == 3)
                        ranfluxz, &
#endif
                        rancorn, & 
                        eta, zeta, kappa, &
                        chi, Dij, & 
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
    real(amrex_real), intent(in   ) :: cons(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nvars)
    real(amrex_real), intent(in   ) :: prim(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nprimvars)
    real(amrex_real), intent(in   ) :: rancorn(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
    real(amrex_real), intent(in   ) :: eta  (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
    real(amrex_real), intent(in   ) :: zeta (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
    real(amrex_real), intent(in   ) :: kappa(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
    real(amrex_real), intent(in   ) :: chi  (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nspecies)
    real(amrex_real), intent(in   ) :: Dij  (lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nspecies,nspecies)

    ! Enforce flux boundary conditions
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
    real(amrex_real), intent(inout) :: xflux(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),nvars)
    real(amrex_real), intent(inout) :: yflux(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),nvars)
#if (AMREX_SPACEDIM == 3)
    real(amrex_real), intent(inout) :: zflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,nvars)
#endif

    real(amrex_real) :: sqrtTwo

    integer :: bc_iter, bc_tmp, indx_lo, indx_hi
    integer :: mass_ind_lo,mass_ind_hi, therm_ind_lo,therm_ind_hi, vel_ind_lo,vel_ind_hi

    integer :: i,j,k,l

      !wall cell - hard wired for specular adiabatic for now
      if(lo(1) .eq. membrane_cell) then
        do k = lo(3),hi(3)
          do j = lo(2),hi(2)

              xflux(membrane_cell,j,k,2) = 0
              xflux(membrane_cell,j,k,3) = 0
              xflux(membrane_cell,j,k,4) = 0
              xflux(membrane_cell,j,k,5) = 0

!              xflux(membrane_cell,j,k,2) = 1.4142*xflux(membrane_cell,j,k,2)        
!              xflux(membrane_cell,j,k,5) = 1.4142*xflux(membrane_cell,j,k,5)        

          end do
        end do
      endif

      !wall cell
      if(hi(1) .eq. membrane_cell-1) then
        do k = lo(3),hi(3)
          do j = lo(2),hi(2)

              xflux(membrane_cell,j,k,2) = 0
              xflux(membrane_cell,j,k,3) = 0
              xflux(membrane_cell,j,k,4) = 0
              xflux(membrane_cell,j,k,5) = 0

!              xflux(membrane_cell,j,k,2) = 1.4142*xflux(membrane_cell,j,k,2)        
!              xflux(membrane_cell,j,k,5) = 1.4142*xflux(membrane_cell,j,k,5)        

          end do
        end do
      endif


!!!!!!!!!!!!!! y-flux BCs !!!!!!!!!!!!!!


  end subroutine stoch_flux_BC

end module flux_module

































