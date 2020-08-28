module flux_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, molmass, &
                                     cell_depth, k_b, runiv, n_cells, membrane_cell, &
                                     visc_type, algorithm_type, &
                                     bc_mass_lo, bc_mass_hi, bc_therm_lo, bc_therm_hi, &
                                     bc_vel_lo, bc_vel_hi
  use conv_module, only : get_energy, &
                          get_density_gas, &
                          get_energy_gas, get_hc_gas
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

    sqrtTwo = sqrt(2.0)
    
    ! Set index ranges

    mass_ind_lo = 3 + AMREX_SPACEDIM
    mass_ind_hi = 2 + AMREX_SPACEDIM + nspecies

    therm_ind_lo = 2 + AMREX_SPACEDIM
    therm_ind_hi = 2 + AMREX_SPACEDIM

    vel_ind_lo = 2
    vel_ind_hi = 1 + AMREX_SPACEDIM

!!!!!!!!

    !! Template:

    ! do bc_iter = 1,3 ! iterate over 1) mass, 2) temperature, and 3) velocity BCs

    !    SELECT CASE (bc_iter)
    !    CASE (1) ! mass boundary conditions
    !       bc_tmp = bc_mass_()
    !       indx_lo = mass_ind_lo
    !       indx_hi = mass_ind_hi
    !    CASE (2) ! temperature boundary conditions
    !       bc_tmp = bc_therm_()
    !       indx_lo = therm_ind_lo
    !       indx_hi = therm_ind_hi
    !    CASE (3) ! velocity boundary conditions
    !       bc_tmp = bc_vel_()
    !       indx_lo = vel_ind_lo
    !       indx_hi = vel_ind_hi
    !    END SELECT

    !    if(bc_tmp .eq. 1) then ! neumann (0 scaling)
    !       do l = indx_lo,indx_hi
    !          if(l.ne.) then ! neumann if not normal velocity
    !          else           ! dirichlet if normal velocity
    !          endif
    !       enddo
    !    elseif(bc_tmp .eq. 2) then ! dirichlet (root 2 scaling) 
    !       do l = indx_lo,indx_hi
    !       enddo
    !    endif

    ! enddo

!!!!!!!!!!!!!! x-flux BCs !!!!!!!!!!!!!!

    !if on lower bound
    if(lo(1) .eq. 0) then !lower x bound
       
       do bc_iter = 1,3 ! iterate over 1) mass, 2) temperature, and 3) velocity BCs
          
          SELECT CASE (bc_iter)
          CASE (1) ! mass boundary conditions
             bc_tmp = bc_mass_lo(1)
             indx_lo = mass_ind_lo
             indx_hi = mass_ind_hi
          CASE (2) ! temperature boundary conditions
             bc_tmp = bc_therm_lo(1)
             indx_lo = therm_ind_lo
             indx_hi = therm_ind_hi
          CASE (3) ! velocity boundary conditions
             bc_tmp = bc_vel_lo(1)
             indx_lo = vel_ind_lo
             indx_hi = vel_ind_hi
          END SELECT

          if(bc_tmp .eq. 1) then ! neumann (0 scaling)

             do l = indx_lo,indx_hi

                if(l.ne.2) then ! neumann if not normal velocity
                   do k = lo(3),hi(3)
                      do j = lo(2),hi(2)
                         xflux(0,j,k,l) = 0
                      end do
                   end do
                else            ! dirichlet if normal velocity
                   do k = lo(3),hi(3)
                      do j = lo(2),hi(2)
                         xflux(0,j,k,l) = sqrtTwo*xflux(0,j,k,l)
                      end do
                   end do
                endif

             end do

          elseif(bc_tmp .eq. 2) then ! dirichlet (root 2 scaling) 

             do l = indx_lo,indx_hi

                do k = lo(3),hi(3)
                   do j = lo(2),hi(2)
                      xflux(0,j,k,l) = sqrtTwo*xflux(0,j,k,l)
                   end do
                end do

             end do

          endif

       enddo

    endif

    !if on upper bound
    if(hi(1) .eq. n_cells(1)-1) then

       do bc_iter = 1,3 ! iterate over 1) mass, 2) temperature, and 3) velocity BCs
          
          SELECT CASE (bc_iter)
          CASE (1) ! mass boundary conditions
             bc_tmp = bc_mass_hi(1)
             indx_lo = mass_ind_lo
             indx_hi = mass_ind_hi
          CASE (2) ! temperature boundary conditions
             bc_tmp = bc_therm_hi(1)
             indx_lo = therm_ind_lo
             indx_hi = therm_ind_hi
          CASE (3) ! velocity boundary conditions
             bc_tmp = bc_vel_hi(1)
             indx_lo = vel_ind_lo
             indx_hi = vel_ind_hi
          END SELECT

          if(bc_tmp .eq. 1) then ! neumann (0 scaling)

             do l = indx_lo,indx_hi

                if(l.ne.2) then ! neumann if not normal velocity
                   do k = lo(3),hi(3)
                      do j = lo(2),hi(2)
                         xflux(hi(1)+1,j,k,l) = 0
                      end do
                   end do
                else           ! dirichlet if normal velocity
                   do k = lo(3),hi(3)
                      do j = lo(2),hi(2)
                         xflux(hi(1)+1,j,k,l) = sqrtTwo*xflux(hi(1)+1,j,k,l)        
                      end do
                   end do
                endif

             enddo

          elseif(bc_tmp .eq. 2) then ! dirichlet (root 2 scaling) 

             do l = indx_lo,indx_hi

                do k = lo(3),hi(3)
                   do j = lo(2),hi(2)
                      xflux(hi(1)+1,j,k,l) = sqrtTwo*xflux(hi(1)+1,j,k,l)        
                   end do
                end do

             enddo

          endif

       enddo

    endif

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

    !if on lower bound
    if(lo(2) .eq. 0) then
       
       do bc_iter = 1,3 ! iterate over 1) mass, 2) temperature, and 3) velocity BCs

          SELECT CASE (bc_iter)
          CASE (1) ! mass boundary conditions
             bc_tmp = bc_mass_lo(2)
             indx_lo = mass_ind_lo
             indx_hi = mass_ind_hi
          CASE (2) ! temperature boundary conditions
             bc_tmp = bc_therm_lo(2)
             indx_lo = therm_ind_lo
             indx_hi = therm_ind_hi
          CASE (3) ! velocity boundary conditions
             bc_tmp = bc_vel_lo(2)
             indx_lo = vel_ind_lo
             indx_hi = vel_ind_hi
          END SELECT

          if(bc_tmp .eq. 1) then ! neumann (0 scaling)

             do l = indx_lo,indx_hi

                if(l.ne.3) then ! neumann if not normal velocity
                   do k = lo(3),hi(3)
                      do i = lo(1),hi(1)
                         yflux(i,0,k,l) = 0
                      end do
                   end do
                else           ! dirichlet if normal velocity
                   do k = lo(3),hi(3)
                      do i = lo(1),hi(1)
                         yflux(i,0,k,l) = sqrtTwo*yflux(i,0,k,l)
                      end do
                   end do
                endif

             enddo

          elseif(bc_tmp .eq. 2) then ! dirichlet (root 2 scaling) 

             do l = indx_lo,indx_hi

                do k = lo(3),hi(3)
                   do i = lo(1),hi(1)
                      yflux(i,0,k,l) = sqrtTwo*yflux(i,0,k,l)
                   end do
                end do

             enddo

          endif

       enddo

    endif

    !if on upper bound
    if(hi(2) .eq. n_cells(2)-1) then

       do bc_iter = 1,3 ! iterate over 1) mass, 2) temperature, and 3) velocity BCs

          SELECT CASE (bc_iter)
          CASE (1) ! mass boundary conditions
             bc_tmp = bc_mass_hi(2)
             indx_lo = mass_ind_lo
             indx_hi = mass_ind_hi
          CASE (2) ! temperature boundary conditions
             bc_tmp = bc_therm_hi(2)
             indx_lo = therm_ind_lo
             indx_hi = therm_ind_hi
          CASE (3) ! velocity boundary conditions
             bc_tmp = bc_vel_hi(2)
             indx_lo = vel_ind_lo
             indx_hi = vel_ind_hi
          END SELECT

          if(bc_tmp .eq. 1) then ! neumann (0 scaling)

             do l = indx_lo,indx_hi

                if(l.ne.3) then ! neumann if not normal velocity
                   do k = lo(3),hi(3)
                      do i = lo(1),hi(1)
                         yflux(i,hi(2)+1,k,l) = 0
                      end do
                   end do
                else           ! dirichlet if normal velocity
                   do k = lo(3),hi(3)
                      do i = lo(1),hi(1)
                         yflux(i,hi(2)+1,k,l) = sqrtTwo*yflux(i,hi(2)+1,k,l)        
                      end do
                   end do
                endif

             enddo

          elseif(bc_tmp .eq. 2) then ! dirichlet (root 2 scaling) 

             do l = indx_lo,indx_hi

                do k = lo(3),hi(3)
                   do i = lo(1),hi(1)
                      yflux(i,hi(2)+1,k,l) = sqrtTwo*yflux(i,hi(2)+1,k,l)        
                   end do
                end do

             enddo

          endif

       enddo

    endif

!!!!!!!!!!!!!! z-flux BCs !!!!!!!!!!!!!!

    !if on lower bound
    if(lo(3) .eq. 0) then

       do bc_iter = 1,3 ! iterate over 1) mass, 2) temperature, and 3) velocity BCs

          SELECT CASE (bc_iter)
          CASE (1) ! mass boundary conditions
             bc_tmp = bc_mass_lo(3)
             indx_lo = mass_ind_lo
             indx_hi = mass_ind_hi
          CASE (2) ! temperature boundary conditions
             bc_tmp = bc_therm_lo(3)
             indx_lo = therm_ind_lo
             indx_hi = therm_ind_hi
          CASE (3) ! velocity boundary conditions
             bc_tmp = bc_vel_lo(3)
             indx_lo = vel_ind_lo
             indx_hi = vel_ind_hi
          END SELECT

          if(bc_tmp .eq. 1) then ! neumann (0 scaling)

             do l = indx_lo,indx_hi

                if(l.ne.4) then ! neumann if not normal velocity
                   do j = lo(2),hi(2)
                      do i = lo(1),hi(1)
                         zflux(i,j,0,l) = 0
                      end do
                   end do
                else           ! dirichlet if normal velocity
                   do j = lo(2),hi(2)
                      do i = lo(1),hi(1)
                         zflux(i,j,0,l) = sqrtTwo*zflux(i,j,0,l)
                      end do
                   end do
                endif

             enddo

          elseif(bc_tmp .eq. 2) then ! dirichlet (root 2 scaling) 

             do l = indx_lo,indx_hi

                do j = lo(2),hi(2)
                   do i = lo(1),hi(1)
                      zflux(i,j,0,l) = sqrtTwo*zflux(i,j,0,l)
                   end do
                end do

             enddo

          endif

       enddo

    endif

    !if on upper bound
    if(hi(3) .eq. n_cells(3)-1) then

       do bc_iter = 1,3 ! iterate over 1) mass, 2) temperature, and 3) velocity BCs

          SELECT CASE (bc_iter)
          CASE (1) ! mass boundary conditions
             bc_tmp = bc_mass_hi(3)
             indx_lo = mass_ind_lo
             indx_hi = mass_ind_hi
          CASE (2) ! temperature boundary conditions
             bc_tmp = bc_therm_hi(3)
             indx_lo = therm_ind_lo
             indx_hi = therm_ind_hi
          CASE (3) ! velocity boundary conditions
             bc_tmp = bc_vel_hi(3)
             indx_lo = vel_ind_lo
             indx_hi = vel_ind_hi
          END SELECT

          if(bc_tmp .eq. 1) then ! neumann (0 scaling)

             do l = indx_lo,indx_hi

                if(l.ne.4) then ! neumann if not normal velocity
                   do j = lo(2),hi(2)
                      do i = lo(1),hi(1)
                         zflux(i,j,hi(3)+1,l) = 0
                      end do
                   end do
                else           ! dirichlet if normal velocity
                   do j = lo(2),hi(2)
                      do i = lo(1),hi(1)
                         zflux(i,j,hi(3)+1,l) = sqrtTwo*zflux(i,j,hi(3)+1,l)        
                      end do
                   end do
                endif

             enddo

          elseif(bc_tmp .eq. 2) then ! dirichlet (root 2 scaling) 

             do l = indx_lo,indx_hi

                do j = lo(2),hi(2)
                   do i = lo(1),hi(1)
                      zflux(i,j,hi(3)+1,l) = sqrtTwo*zflux(i,j,hi(3)+1,l)        
                   end do
                end do

             enddo

          endif

       enddo

    endif

  end subroutine stoch_flux_BC

end module flux_module

































