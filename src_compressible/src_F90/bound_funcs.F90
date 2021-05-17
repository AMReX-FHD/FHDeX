module bound_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, t_lo, t_hi, nprimvars, nvars, nspecies, n_cells, &
                                     algorithm_type, membrane_cell, MAX_SPECIES, &
                                     bc_mass_lo, bc_mass_hi, bc_therm_lo, &
                                     bc_therm_hi, bc_vel_lo, bc_vel_hi, &
                                     bc_Yk_x_lo, bc_Yk_x_hi, bc_Yk_y_lo, bc_Yk_y_hi, bc_Yk_z_lo, bc_Yk_z_hi, &
                                     bc_Xk_x_lo, bc_Xk_x_hi, bc_Xk_y_lo, bc_Xk_y_hi, bc_Xk_z_lo, bc_Xk_z_hi
  !use compressible_namelist_module, only : bc_Yk, bc_Xk
  use conv_module

  implicit none

  integer, parameter :: LOHI = 2

  private

contains

  subroutine setup_bc() bind(C,name="setup_bc")

    integer :: d

    ! for each direction, if bc_vel_lo/hi is periodic, then
    ! set the corresponding bc_mass_lo/hi and bc_therm_lo/hi to periodic
    do d=1,AMREX_SPACEDIM

       if (bc_vel_lo(d) .eq. -1) then
          bc_mass_lo(d) = -1
          bc_therm_lo(d) = -1
          bc_mass_hi(d) = -1
          bc_therm_hi(d) = -1
       end if
       
    enddo

  end subroutine setup_bc

  subroutine setup_cwall() bind(C,name="setup_cwall")
    
    integer :: ns

    real(amrex_real) :: sumx, sumy

    ! Compute Xk or Yk at the wall, depending on which is defined
    ! X walls
    if (bc_mass_lo(1).eq.2) then
       sumx = 0
       sumy = 0
       do ns=1,nspecies
          sumx = sumx + bc_Xk_x_lo(ns)
          sumy = sumy + bc_Yk_x_lo(ns)
       enddo
       if (abs(sumx-1).lt.1.d-10) then
          call get_massfrac(bc_Xk_x_lo(1:nspecies),bc_Yk_x_lo(1:nspecies))
       else if (abs(sumy-1).lt.1d-10) then
          call get_molfrac(bc_Yk_x_lo(1:nspecies),bc_Xk_x_lo(1:nspecies))
       endif
    endif

    if (bc_mass_hi(1).eq.2) then
       sumx = 0
       sumy = 0
       do ns=1,nspecies
          sumx = sumx + bc_Xk_x_hi(ns)
          sumy = sumy + bc_Yk_x_hi(ns)
       enddo
       if (abs(sumx-1).lt.1.d-10) then
          call get_massfrac(bc_Xk_x_hi(1:nspecies),bc_Yk_x_hi(1:nspecies))
       else if (abs(sumy-1).lt.1d-10) then
          call get_molfrac(bc_Yk_x_hi(1:nspecies),bc_Xk_x_hi(1:nspecies))
       endif
    endif

    ! Y walls
    if (bc_mass_lo(2).eq.2) then
       sumx = 0
       sumy = 0
       do ns=1,nspecies
          sumx = sumx + bc_Xk_y_lo(ns)
          sumy = sumy + bc_Yk_y_lo(ns)
       enddo
       if (abs(sumx-1).lt.1.d-10) then
          call get_massfrac(bc_Xk_y_lo(1:nspecies),bc_Yk_y_lo(1:nspecies))
       else if (abs(sumy-1).lt.1d-10) then
          call get_molfrac(bc_Yk_y_lo(1:nspecies),bc_Xk_y_lo(1:nspecies))
       endif
    endif

    if (bc_mass_hi(2).eq.2) then
       sumx = 0
       sumy = 0
       do ns=1,nspecies
          sumx = sumx + bc_Xk_y_hi(ns)
          sumy = sumy + bc_Yk_y_hi(ns)
       enddo
       if (abs(sumx-1).lt.1.d-10) then
          call get_massfrac(bc_Xk_y_hi(1:nspecies),bc_Yk_y_hi(1:nspecies))
       else if (abs(sumy-1).lt.1d-10) then
          call get_molfrac(bc_Yk_y_hi(1:nspecies),bc_Xk_y_hi(1:nspecies))
       endif
    endif

    ! Z walls
    if (bc_mass_lo(3).eq.2) then
       sumx = 0
       sumy = 0
       do ns=1,nspecies
          sumx = sumx + bc_Xk_z_lo(ns)
          sumy = sumy + bc_Yk_z_lo(ns)
       enddo
       if (abs(sumx-1).lt.1.d-10) then
          call get_massfrac(bc_Xk_z_lo(1:nspecies),bc_Yk_z_lo(1:nspecies))
       else if (abs(sumy-1).lt.1d-10) then
          call get_molfrac(bc_Yk_z_lo(1:nspecies),bc_Xk_z_lo(1:nspecies))
       endif
    endif

    if (bc_mass_hi(3).eq.2) then
       sumx = 0
       sumy = 0
       do ns=1,nspecies
          sumx = sumx + bc_Xk_z_hi(ns)
          sumy = sumy + bc_Yk_z_hi(ns)
       enddo
       if (abs(sumx-1).lt.1.d-10) then
          call get_massfrac(bc_Xk_z_hi(1:nspecies),bc_Yk_z_hi(1:nspecies))
       else if (abs(sumy-1).lt.1d-10) then
          call get_molfrac(bc_Yk_z_hi(1:nspecies),bc_Xk_z_hi(1:nspecies))
       endif
    endif

  end subroutine setup_cwall

end module bound_module

