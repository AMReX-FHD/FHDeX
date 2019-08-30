module bound_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, t_lo, t_hi, nprimvars, nvars, nspecies, n_cells, &
       algorithm_type, membrane_cell, MAX_SPECIES, bc_mass_lo, bc_mass_hi, bc_therm_lo, &
       bc_therm_hi, bc_vel_lo, bc_vel_hi
  use compressible_namelist_module, only : bc_Yk, bc_Xk
  use conv_module
  use trans_module

  implicit none

  private

  public :: set_bc

contains

  subroutine set_bc(lo,hi, cons, prim) bind(C,name="set_bc")

    integer         , intent(in   ) :: lo(3),hi(3)

    real(amrex_real), intent(inout) :: prim(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)
    real(amrex_real), intent(inout) :: cons(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)

    integer :: i,j,k,l,idir,bcell

    real(amrex_real) :: massvec(nspecies), fracvec(nspecies), intenergy, temp, rho, pt
    real(amrex_real) :: Ywall(nspecies), Xwall(nspecies)

    Ywall = 0.0d0
    Xwall = 0.0d0

    if(lo(1) .eq. 0) then !lower x bound

       if(bc_mass_lo(1) .eq. 1) then ! wall

          if(algorithm_type.eq.2) then
             do k = lo(3)-ngc(3),hi(3)+ngc(3)
                do j = lo(2)-ngc(2),hi(2)+ngc(2)
                   do i = 1, ngc(1)          

                      prim(lo(1)-i,j,k,6:nprimvars) = prim(lo(1)-1+i,j,k,6:nprimvars)

                   enddo
                enddo
             enddo
          endif

       elseif(bc_mass_lo(1) .eq. 2) then ! reservoir

          if(algorithm_type.eq.2) then

             Ywall(1:nspecies) = bc_Yk(1,1,1:nspecies)
             Xwall(1:nspecies) = bc_Xk(1,1,1:nspecies)

             do k = lo(3)-ngc(3),hi(3)+ngc(3)
                do j = lo(2)-ngc(2),hi(2)+ngc(2)
                   do i = 1,ngc(1)

                      do l = 1, nspecies
                         prim(lo(1)-i,j,k,6+l)          = 2.d0*Ywall(l) - prim(lo(1)-1+i,j,k,6+l)
                         prim(lo(1)-i,j,k,6+nspecies+l) = 2.d0*Xwall(l) - prim(lo(1)-1+i,j,k,6+nspecies+l)
                      enddo

                   enddo
                enddo
             enddo
          endif

       endif

       if(bc_therm_lo(1) .eq. 1) then ! adiabatic

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = lo(2)-ngc(2),hi(2)+ngc(2)
                do i = 1, ngc(1)          

                   prim(lo(1)-i,j,k,5) = prim(lo(1)-1+i,j,k,5)
                   prim(lo(1)-i,j,k,6) = prim(lo(1)-1+i,j,k,6)

                enddo
             enddo
          enddo

       elseif(bc_therm_lo(1) .eq. 2) then ! isothermal

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = lo(2)-ngc(2),hi(2)+ngc(2)
                do i = 1,ngc(1)

                   prim(lo(1)-i,j,k,5) = -prim(lo(1)-1+i,j,k,5) + 2*t_lo(1)
                   prim(lo(1)-i,j,k,6) = prim(lo(1)-1+i,j,k,6)

                enddo
             enddo
          enddo

       endif

      if(bc_vel_lo(1) .eq. 1) then ! slip

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = lo(2)-ngc(2),hi(2)+ngc(2)
                do i = 1, ngc(1)    
      
                   cons(lo(1)-i,j,k,2) = -cons(lo(1)-1+i,j,k,2) 
                   cons(lo(1)-i,j,k,3) = cons(lo(1)-1+i,j,k,3) 
                   cons(lo(1)-i,j,k,4) = cons(lo(1)-1+i,j,k,4) 

                   prim(lo(1)-i,j,k,2) = -prim(lo(1)-1+i,j,k,2)
                   prim(lo(1)-i,j,k,3) = prim(lo(1)-1+i,j,k,3)
                   prim(lo(1)-i,j,k,4) = prim(lo(1)-1+i,j,k,4)

                   ! thermal & species (+pressure) BCs must be enforced first
                   fracvec = prim(lo(1)-i,j,k,6+1:6+nspecies)
                   temp = prim(lo(1)-i,j,k,5)
                   pt = prim(lo(1)-i,j,k,6)

                   call get_density(pt, rho, temp, fracvec)
                   call get_energy(intenergy, fracvec, temp)
                   ! call get_density_gas(pt, rho, temp)
                   ! call get_energy_gas(pt, intenergy)

                   ! total density depends on temperature
                   prim(lo(1)-i,j,k,1) = rho
                   cons(lo(1)-i,j,k,1) = rho
                   if(algorithm_type.eq.2) then
                      do l = 1, nspecies
                         cons(lo(1)-i,j,k,5+l) = rho*prim(lo(1)-i,j,k,6+l)
                      enddo
                   endif

                   ! must be last BC enforced: depends on rho, vel, & temp
                   cons(lo(1)-i,j,k,5) = rho*intenergy + 0.5*rho*(prim(lo(1)-i,j,k,2)**2 + & 
                        prim(lo(1)-i,j,k,3)**2 + prim(lo(1)-i,j,k,4)**2)

                enddo
             enddo
          enddo

       elseif(bc_vel_lo(1) .eq. 2) then ! no slip
 
          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = lo(2)-ngc(2),hi(2)+ngc(2)
                do i = 1,ngc(1)

                   prim(lo(1)-i,j,k,2) = -prim(lo(1)-1+i,j,k,2) ! + 2*vel_lo(1)
                   prim(lo(1)-i,j,k,3) = -prim(lo(1)-1+i,j,k,3) 
                   prim(lo(1)-i,j,k,4) = -prim(lo(1)-1+i,j,k,4)
                   
                   ! thermal & species (+pressure) BCs must be enforced first
                   fracvec = prim(lo(1)-i,j,k,6+1:6+nspecies)
                   temp = prim(lo(1)-i,j,k,5)
                   pt = prim(lo(1)-i,j,k,6)

                   call get_density(pt, rho, temp, fracvec)
                   call get_energy(intenergy, fracvec, temp)
                   ! call get_density_gas(pt, rho, temp)
                   ! call get_energy_gas(pt, intenergy)

                   prim(lo(1)-i,j,k,1) = rho
                   cons(lo(1)-i,j,k,1) = rho
                   if(algorithm_type.eq.2) then
                      do l = 1, nspecies
                         cons(lo(1)-i,j,k,5+l) = rho*prim(lo(1)-i,j,k,6+l)
                      enddo
                   endif

                   cons(lo(1)-i,j,k,2) = rho*prim(lo(1)-i,j,k,2)
                   cons(lo(1)-i,j,k,3) = rho*prim(lo(1)-i,j,k,3)
                   cons(lo(1)-i,j,k,4) = rho*prim(lo(1)-i,j,k,4)

                   ! must be last BC enforced: depends on rho, vel, & temp
                   cons(lo(1)-i,j,k,5) = rho*intenergy + 0.5*rho*(prim(lo(1)-i,j,k,2)**2 + & 
                        prim(lo(1)-i,j,k,3)**2 + prim(lo(1)-i,j,k,4)**2)

                enddo
             enddo
          enddo

       endif

    endif

    if(hi(1) .eq. (n_cells(1)-1)) then !upper x bound

       if(bc_mass_hi(1) .eq. 1) then ! wall

          if(algorithm_type.eq.2) then
             do k = lo(3)-ngc(3),hi(3)+ngc(3)
                do j = lo(2)-ngc(2),hi(2)+ngc(2)
                   do i = 1, ngc(1)

                      prim(hi(1)+i,j,k,6:nprimvars) = prim(hi(1)+1-i,j,k,6:nprimvars) 

                   enddo
                enddo
             enddo
          endif

       elseif(bc_mass_hi(1) .eq. 2) then ! reservoir

          if(algorithm_type.eq.2) then

             Ywall(1:nspecies) = bc_Yk(1,2,1:nspecies)
             Xwall(1:nspecies) = bc_Xk(1,2,1:nspecies)

             do k = lo(3)-ngc(3),hi(3)+ngc(3)
                do j = lo(2)-ngc(2),hi(2)+ngc(2)
                   do i = 1, ngc(1)

                      do l = 1, nspecies
                         prim(hi(1)+i,j,k,6+l)          = 2.d0*Ywall(l) - prim(hi(1)+1-i,j,k,6+l)
                         prim(hi(1)+i,j,k,6+nspecies+l) = 2.d0*Xwall(l) - prim(hi(1)+1-i,j,k,6+nspecies+l)
                      enddo

                   enddo
                enddo
             enddo

          endif

       endif

       if(bc_therm_hi(1) .eq. 1) then ! adiabatic

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = lo(2)-ngc(2),hi(2)+ngc(2)
                do i = 1, ngc(1)

                   prim(hi(1)+i,j,k,5) = prim(hi(1)+1-i,j,k,5)
                   prim(hi(1)+i,j,k,6) = prim(hi(1)+1-i,j,k,6)
                      
                enddo
             enddo
          enddo

       elseif(bc_therm_hi(1) .eq. 2) then ! isothermal

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = lo(2)-ngc(2),hi(2)+ngc(2)
                do i = 1, ngc(1)

                   prim(hi(1)+i,j,k,5) = -prim(hi(1)+1-i,j,k,5) + 2*t_hi(1)
                   prim(hi(1)+i,j,k,6) = prim(hi(1)+1-i,j,k,6)

                enddo
             enddo
          enddo

       endif

       if(bc_vel_hi(1) .eq. 1) then ! slip

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = lo(2)-ngc(2),hi(2)+ngc(2)
                do i = 1, ngc(1)

                   cons(hi(1)+i,j,k,2) = -cons(hi(1)+1-i,j,k,2) 
                   cons(hi(1)+i,j,k,3) = cons(hi(1)+1-i,j,k,3) 
                   cons(hi(1)+i,j,k,4) = cons(hi(1)+1-i,j,k,4)  

                   prim(hi(1)+i,j,k,2) = -prim(hi(1)+1-i,j,k,2)
                   prim(hi(1)+i,j,k,3) = prim(hi(1)+1-i,j,k,3)
                   prim(hi(1)+i,j,k,4) = prim(hi(1)+1-i,j,k,4)

                   ! thermal & species (+pressure) BCs must be enforced first
                   fracvec = prim(hi(1)+i,j,k,6+1:6+nspecies)
                   temp = prim(hi(1)+i,j,k,5)
                   pt = prim(hi(1)+i,j,k,6)

                   call get_density(pt, rho, temp, fracvec)
                   call get_energy(intenergy, fracvec, temp)
                   ! call get_density_gas(pt, rho, temp)
                   ! call get_energy_gas(pt, intenergy)

                   ! total density depends on temperature
                   prim(hi(1)+i,j,k,1) = rho
                   cons(hi(1)+i,j,k,1) = rho
                   if(algorithm_type.eq.2) then
                      do l = 1, nspecies
                         cons(hi(1)+i,j,k,5+l) = rho*prim(hi(1)+i,j,k,6+l)
                      enddo
                   endif

                   ! must be last BC enforced: depends on rho, vel, & temp
                   cons(hi(1)+i,j,k,5) = rho*intenergy + 0.5*rho*(prim(hi(1)+i,j,k,2)**2 + & 
                        prim(hi(1)+i,j,k,3)**2 + prim(hi(1)+i,j,k,4)**2)
                      
                enddo
             enddo
          enddo

       elseif(bc_vel_hi(1) .eq. 2) then ! no slip

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = lo(2)-ngc(2),hi(2)+ngc(2)
                do i = 1, ngc(1)

                   prim(hi(1)+i,j,k,2) = -prim(hi(1)+1-i,j,k,2) ! + 2*vel_hi(1)
                   prim(hi(1)+i,j,k,3) = -prim(hi(1)+1-i,j,k,3) 
                   prim(hi(1)+i,j,k,4) = -prim(hi(1)+1-i,j,k,4)
                   
                   ! thermal & species (+pressure) BCs must be enforced first
                   fracvec = prim(hi(1)+i,j,k,6+1:6+nspecies)
                   temp = prim(hi(1)+i,j,k,5)
                   pt = prim(hi(1)+i,j,k,6)

                   call get_density(pt, rho, temp, fracvec)
                   call get_energy(intenergy, fracvec, temp)
                   ! call get_density_gas(pt, rho, temp)
                   ! call get_energy_gas(pt, intenergy)

                   ! total density depends on temperature
                   prim(hi(1)+i,j,k,1) = rho
                   cons(hi(1)+i,j,k,1) = rho
                   if(algorithm_type.eq.2) then
                      do l = 1, nspecies
                         cons(hi(1)+i,j,k,5+l) = rho*prim(hi(1)+i,j,k,6+l)
                      enddo
                   endif

                   cons(hi(1)+i,j,k,2) = rho*prim(hi(1)+i,j,k,2)
                   cons(hi(1)+i,j,k,3) = rho*prim(hi(1)+i,j,k,3)
                   cons(hi(1)+i,j,k,4) = rho*prim(hi(1)+i,j,k,4)

                   ! must be last BC enforced: depends on rho, vel, & temp
                   cons(hi(1)+i,j,k,5) = rho*intenergy + 0.5*rho*(prim(hi(1)+i,j,k,2)**2 + & 
                        prim(hi(1)+i,j,k,3)**2 + prim(hi(1)+i,j,k,4)**2)
                   
                enddo
             enddo
          enddo

       endif

    endif


    if(lo(2) .eq. 0) then !lower y bound

       if(bc_mass_lo(2) .eq. 1) then ! wall

          if(algorithm_type.eq.2) then
             do k = lo(3)-ngc(3),hi(3)+ngc(3)
                do j = 1,ngc(2)
                   do i = lo(1)-ngc(1),hi(1)+ngc(1)

                      prim(i,lo(2)-j,k,6:nprimvars) = prim(i,lo(2)-1+j,k,6:nprimvars)

                   enddo
                enddo
             enddo
          endif

       elseif(bc_mass_lo(2) .eq. 2) then ! reservoir

          if(algorithm_type.eq.2) then

             Ywall(1:nspecies) = bc_Yk(2,1,1:nspecies)
             Xwall(1:nspecies) = bc_Xk(2,1,1:nspecies)

             do k = lo(3)-ngc(3),hi(3)+ngc(3)
                do j = 1,ngc(2)
                   do i = lo(1)-ngc(1),hi(1)+ngc(1)

                      do l = 1, nspecies
                         prim(i,lo(2)-j,k,6+l)          = 2.d0*Ywall(l) - prim(i,lo(2)-1+j,k,6+l)
                         prim(i,lo(2)-j,k,6+nspecies+l) = 2.d0*Xwall(l) - prim(i,lo(2)-1+j,k,6+nspecies+l)
                      enddo

                   enddo
                enddo
             enddo

          endif

       endif

       if(bc_therm_lo(2) .eq. 1) then ! adiabatic

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = 1,ngc(2)
                do i = lo(1)-ngc(1),hi(1)+ngc(1)

                   prim(i,lo(2)-j,k,5) = prim(i,lo(2)-1+j,k,5)
                   prim(i,lo(2)-j,k,6) = prim(i,lo(2)-1+j,k,6)

                enddo
             enddo
          enddo

       elseif(bc_therm_lo(2) .eq. 2) then ! isothermal

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = 1,ngc(2)
                do i = lo(1)-ngc(1),hi(1)+ngc(1)

                   prim(i,lo(2)-j,k,5) = -prim(i,lo(2)-1+j,k,5) + 2*t_lo(2)
                   prim(i,lo(2)-j,k,6) = prim(i,lo(2)-1+j,k,6)

                enddo
             enddo
          enddo

       endif

       if(bc_vel_lo(2) .eq. 1) then ! slip

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = 1,ngc(2)
                do i = lo(1)-ngc(1),hi(1)+ngc(1)          

                   cons(i,lo(2)-j,k,2) = cons(i,lo(2)-1+j,k,2) 
                   cons(i,lo(2)-j,k,3) = -cons(i,lo(2)-1+j,k,3) 
                   cons(i,lo(2)-j,k,4) = cons(i,lo(2)-1+j,k,4) 

                   prim(i,lo(2)-j,k,2) = prim(i,lo(2)-1+j,k,2)
                   prim(i,lo(2)-j,k,3) = -prim(i,lo(2)-1+j,k,3)
                   prim(i,lo(2)-j,k,4) = prim(i,lo(2)-1+j,k,4)

                   ! thermal & species (+pressure) BCs must be enforced first
                   fracvec = prim(i,lo(2)-j,k,6+1:6+nspecies)
                   temp = prim(i,lo(2)-j,k,5)
                   pt = prim(i,lo(2)-j,k,6)

                   call get_density(pt, rho, temp, fracvec)
                   call get_energy(intenergy, fracvec, temp)
                   ! call get_density_gas(pt, rho, temp)
                   ! call get_energy_gas(pt, intenergy)

                   prim(i,lo(2)-j,k,1) = rho
                   cons(i,lo(2)-j,k,1) = rho
                   if(algorithm_type.eq.2) then
                      do l = 1, nspecies
                         cons(i,lo(2)-j,k,5+l) = rho*prim(i,lo(2)-j,k,6+l)
                      enddo
                   endif

                   ! must be last BC enforced: depends on rho, vel, & temp
                   cons(i,lo(2)-j,k,5) = rho*intenergy + 0.5*rho*(prim(i,lo(2)-j,k,2)**2 + & 
                        prim(i,lo(2)-j,k,3)**2 + prim(i,lo(2)-j,k,4)**2)

                enddo
             enddo
          enddo

       elseif(bc_vel_lo(2) .eq. 2) then ! no slip

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = 1,ngc(2)
                do i = lo(1)-ngc(1),hi(1)+ngc(1)

                   prim(i,lo(2)-j,k,2) = -prim(i,lo(2)-1+j,k,2)
                   prim(i,lo(2)-j,k,3) = -prim(i,lo(2)-1+j,k,3) ! + 2*vel_lo(2) 
                   prim(i,lo(2)-j,k,4) = -prim(i,lo(2)-1+j,k,4)

                   ! thermal & species (+pressure) BCs must be enforced first
                   fracvec = prim(i,lo(2)-j,k,6+1:6+nspecies)
                   temp = prim(i,lo(2)-j,k,5)
                   pt = prim(i,lo(2)-j,k,6)

                   call get_density(pt, rho, temp, fracvec)
                   call get_energy(intenergy, fracvec, temp)
                   ! call get_density_gas(pt, rho, temp)
                   ! call get_energy_gas(pt, intenergy)

                   prim(i,lo(2)-j,k,1) = rho
                   cons(i,lo(2)-j,k,1) = rho
                   if(algorithm_type.eq.2) then
                      do l = 1, nspecies
                         cons(i,lo(2)-j,k,5+l) = rho*prim(i,lo(2)-j,k,6+l)
                      enddo
                   endif

                   cons(i,lo(2)-j,k,2) = rho*prim(i,lo(2)-j,k,2)
                   cons(i,lo(2)-j,k,3) = rho*prim(i,lo(2)-j,k,3)
                   cons(i,lo(2)-j,k,4) = rho*prim(i,lo(2)-j,k,4)

                   ! must be last BC enforced: depends on rho, vel, & temp
                   cons(i,lo(2)-j,k,5) = rho*intenergy + 0.5*rho*(prim(i,lo(2)-j,k,2)**2 + & 
                        prim(i,lo(2)-j,k,3)**2 + prim(i,lo(2)-j,k,4)**2)

                enddo
             enddo
          enddo

       endif

    endif

    if(hi(2) .eq. (n_cells(2)-1)) then !upper y bound

       if(bc_mass_hi(2) .eq. 1) then ! wall

          if(algorithm_type.eq.2) then
             do k = lo(3)-ngc(3),hi(3)+ngc(3)
                do j = 1,ngc(2)
                   do i = lo(1)-ngc(1),hi(1)+ngc(1)

                      prim(i,hi(2)+j,k,6:nprimvars) = prim(i,hi(2)+1-j,k,6:nprimvars)

                   enddo
                enddo
             enddo
          endif

       elseif(bc_mass_hi(2) .eq. 2) then ! reservoir

          if(algorithm_type.eq.2) then
             
             Ywall(1:nspecies) = bc_Yk(2,2,1:nspecies)
             Xwall(1:nspecies) = bc_Xk(2,2,1:nspecies)

             do k = lo(3)-ngc(3),hi(3)+ngc(3)
                do j = 1,ngc(2)
                   do i = lo(1)-ngc(1),hi(1)+ngc(1)

                      do l = 1, nspecies
                         prim(i,hi(2)+j,k,6+l)          = 2.d0*Ywall(l) - prim(i,hi(2)+1-j,k,6+l)
                         prim(i,hi(2)+j,k,6+nspecies+l) = 2.d0*Xwall(l) - prim(i,hi(2)+1-j,k,6+nspecies+l)
                      enddo

                   enddo
                enddo
             enddo
          endif

       endif

       if(bc_therm_hi(2) .eq. 1) then ! adiabatic

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = 1,ngc(2)
                do i = lo(1)-ngc(1),hi(1)+ngc(1)

                   prim(i,hi(2)+j,k,5) = prim(i,hi(2)+1-j,k,5)
                   prim(i,hi(2)+j,k,6) = prim(i,hi(2)+1-j,k,6)

                enddo
             enddo
          enddo

       elseif(bc_therm_hi(2) .eq. 2) then ! isothermal

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = 1,ngc(2)
                do i = lo(1)-ngc(1),hi(1)+ngc(1)

                   prim(i,hi(2)+j,k,5) = -prim(i,hi(2)+1-j,k,5) + 2*t_hi(2)
                   prim(i,hi(2)+j,k,6) = prim(i,hi(2)+1-j,k,6)

                enddo
             enddo
          enddo

       endif

       if(bc_vel_hi(2) .eq. 1) then ! slip

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = 1,ngc(2)
                do i = lo(1)-ngc(1),hi(1)+ngc(1)

                   cons(i,hi(2)+j,k,2) = cons(i,hi(2)+1-j,k,2) 
                   cons(i,hi(2)+j,k,3) = -cons(i,hi(2)+1-j,k,3) 
                   cons(i,hi(2)+j,k,4) = cons(i,hi(2)+1-j,k,4) 

                   prim(i,hi(2)+j,k,2) = prim(i,hi(2)+1-j,k,2)
                   prim(i,hi(2)+j,k,3) = -prim(i,hi(2)+1-j,k,3)
                   prim(i,hi(2)+j,k,4) = prim(i,hi(2)+1-j,k,4)
 
                   ! thermal & species (+pressure) BCs must be enforced first
                   fracvec = prim(i,hi(2)+j,k,6+1:6+nspecies)
                   temp = prim(i,hi(2)+j,k,5)
                   pt = prim(i,hi(2)+j,k,6)

                   call get_density(pt, rho, temp, fracvec)
                   call get_energy(intenergy, fracvec, temp)
                   ! call get_density_gas(pt, rho, temp)
                   ! call get_energy_gas(pt, intenergy)
                   
                   prim(i,hi(2)+j,k,1) = rho
                   cons(i,hi(2)+j,k,1) = rho
                   if(algorithm_type.eq.2) then
                      do l = 1, nspecies
                         cons(i,hi(2)+j,k,5+l) = rho*prim(i,hi(2)+j,k,6+l)
                      enddo
                   endif

                   ! must be last BC enforced: depends on rho, vel, & temp
                   cons(i,hi(2)+j,k,5) = rho*intenergy + 0.5*rho*(prim(i,hi(2)+j,k,2)**2 + & 
                        prim(i,hi(2)+j,k,3)**2 + prim(i,hi(2)+j,k,4)**2)

                enddo
             enddo
          enddo

       elseif(bc_vel_hi(2) .eq. 2) then ! no slip

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = 1,ngc(2)
                do i = lo(1)-ngc(1),hi(1)+ngc(1)

                   prim(i,hi(2)+j,k,2) = -prim(i,hi(2)+1-j,k,2)
                   prim(i,hi(2)+j,k,3) = -prim(i,hi(2)+1-j,k,3) ! + 2*vel_hi(2)
                   prim(i,hi(2)+j,k,4) = -prim(i,hi(2)+1-j,k,4)

                   ! thermal & species (+pressure) BCs must be enforced first
                   fracvec = prim(i,hi(2)+j,k,6+1:6+nspecies)
                   temp = prim(i,hi(2)+j,k,5)
                   pt = prim(i,hi(2)+j,k,6)

                   call get_density(pt, rho, temp, fracvec)
                   call get_energy(intenergy, fracvec, temp)
                   ! call get_density_gas(pt, rho, temp)
                   ! call get_energy_gas(pt, intenergy)
                   
                   prim(i,hi(2)+j,k,1) = rho
                   cons(i,hi(2)+j,k,1) = rho
                   if(algorithm_type.eq.2) then
                      do l = 1, nspecies
                         cons(i,hi(2)+j,k,5+l) = rho*prim(i,hi(2)+j,k,6+l)
                      enddo
                   endif

                   cons(i,hi(2)+j,k,2) = rho*prim(i,hi(2)+j,k,2)
                   cons(i,hi(2)+j,k,3) = rho*prim(i,hi(2)+j,k,3)
                   cons(i,hi(2)+j,k,4) = rho*prim(i,hi(2)+j,k,4)

                   ! must be last BC enforced: depends on rho, vel, & temp
                   cons(i,hi(2)+j,k,5) = rho*intenergy + 0.5*rho*(prim(i,hi(2)+j,k,2)**2 + & 
                        prim(i,hi(2)+j,k,3)**2 + prim(i,hi(2)+j,k,4)**2)
                   
                enddo
             enddo
          enddo

       endif

    endif

    if(n_cells(3).gt.1) then

       if(lo(3) .eq. 0) then !lower z bound

          if(bc_mass_lo(3) .eq. 1) then ! wall

             if(algorithm_type.eq.2) then
                do k = 1,ngc(3)
                   do j = lo(2)-ngc(2),hi(2)+ngc(2)
                      do i = lo(1)-ngc(1),hi(1)+ngc(1)

                         prim(i,j,lo(3)-k,6:nprimvars) = prim(i,j,lo(3)-1+k,6:nprimvars)

                      enddo
                   enddo
                enddo
             endif

          elseif(bc_mass_lo(3) .eq. 2) then ! reservoir

             if(algorithm_type.eq.2) then

                Ywall(1:nspecies) = bc_Yk(3,1,1:nspecies)
                Xwall(1:nspecies) = bc_Xk(3,1,1:nspecies)

                do k = 1,ngc(3)
                   do j = lo(2)-ngc(2),hi(2)+ngc(2)
                      do i = lo(1)-ngc(1),hi(1)+ngc(1)

                         do l = 1, nspecies
                            prim(i,j,lo(3)-k,6+l)          = 2.d0*Ywall(l) - prim(i,j,lo(3)-1+k,6+l)
                            prim(i,j,lo(3)-k,6+nspecies+l) = 2.d0*Xwall(l) - prim(i,j,lo(3)-1+k,6+nspecies+l)
                         enddo

                      enddo
                   enddo
                enddo
             endif


          endif

          if(bc_therm_lo(3) .eq. 1) then ! adiabatic

             do k = 1,ngc(3)
                do j = lo(2)-ngc(2),hi(2)+ngc(2)
                   do i = lo(1)-ngc(1),hi(1)+ngc(1)

                      prim(i,j,lo(3)-k,5) = prim(i,j,lo(3)-1+k,5)
                      prim(i,j,lo(3)-k,6) = prim(i,j,lo(3)-1+k,6)

                   enddo
                enddo
             enddo

          elseif(bc_therm_lo(3) .eq. 2) then ! isothermal

             do k = 1,ngc(3)
                do j = lo(2)-ngc(2),hi(2)+ngc(2)
                   do i = lo(1)-ngc(1),hi(1)+ngc(1)

                      prim(i,j,lo(3)-k,5) = -prim(i,j,lo(3)-1+k,5) + 2*t_lo(3)
                      prim(i,j,lo(3)-k,6) = prim(i,j,lo(3)-1+k,6)

                   enddo
                enddo
             enddo


          endif

          if(bc_vel_lo(3) .eq. 1) then ! slip

             do k = 1,ngc(3)
                do j = lo(2)-ngc(2),hi(2)+ngc(2)
                   do i = lo(1)-ngc(1),hi(1)+ngc(1)

                      cons(i,j,lo(3)-k,2) = cons(i,j,lo(3)-1+k,2) 
                      cons(i,j,lo(3)-k,3) = cons(i,j,lo(3)-1+k,3) 
                      cons(i,j,lo(3)-k,4) = -cons(i,j,lo(3)-1+k,4) 

                      prim(i,j,lo(3)-k,2) = prim(i,j,lo(3)-1+k,2)
                      prim(i,j,lo(3)-k,3) = prim(i,j,lo(3)-1+k,3)
                      prim(i,j,lo(3)-k,4) = -prim(i,j,lo(3)-1+k,4)

                      ! thermal & species (+pressure) BCs must be enforced first
                      fracvec = prim(i,j,lo(3)-k,6+1:6+nspecies)
                      temp = prim(i,j,lo(3)-k,5)
                      pt = prim(i,j,lo(3)-k,6)

                      call get_density(pt, rho, temp, fracvec)
                      call get_energy(intenergy, fracvec, temp)
                      ! call get_density_gas(pt, rho, temp)
                      ! call get_energy_gas(pt, intenergy)

                      prim(i,j,lo(3)-k,1) = rho
                      cons(i,j,lo(3)-k,1) = rho
                      if(algorithm_type.eq.2) then
                         do l = 1, nspecies
                            cons(i,j,lo(3)-k,5+l) = rho*prim(i,j,lo(3)-k,6+l)
                         enddo
                      endif

                      ! must be last BC enforced: depends on rho, vel, & temp
                      cons(i,j,lo(3)-k,5) = rho*intenergy + 0.5*rho*(prim(i,j,lo(3)-k,2)**2 + & 
                           prim(i,j,lo(3)-k,3)**2 + prim(i,j,lo(3)-k,4)**2)

                   enddo
                enddo
             enddo

          elseif(bc_vel_lo(3) .eq. 2) then ! no slip

             do k = 1,ngc(3)
                do j = lo(2)-ngc(2),hi(2)+ngc(2)
                   do i = lo(1)-ngc(1),hi(1)+ngc(1)

                      prim(i,j,lo(3)-k,2) = -prim(i,j,lo(3)-1+k,2) 
                      prim(i,j,lo(3)-k,3) = -prim(i,j,lo(3)-1+k,3) 
                      prim(i,j,lo(3)-k,4) = -prim(i,j,lo(3)-1+k,4) ! + 2*vel_lo(3)

                      ! thermal & species (+pressure) BCs must be enforced first
                      fracvec = prim(i,j,lo(3)-k,6+1:6+nspecies)
                      temp = prim(i,j,lo(3)-k,5)
                      pt = prim(i,j,lo(3)-k,6)

                      call get_density(pt, rho, temp, fracvec)
                      call get_energy(intenergy, fracvec, temp)
                      ! call get_density_gas(pt, rho, temp)
                      ! call get_energy_gas(pt, intenergy)

                      prim(i,j,lo(3)-k,1) = rho
                      cons(i,j,lo(3)-k,1) = rho
                      if(algorithm_type.eq.2) then
                         do l = 1, nspecies
                            cons(i,j,lo(3)-k,5+l) = rho*prim(i,j,lo(3)-k,6+l)
                         enddo
                      endif

                      cons(i,j,lo(3)-k,2) = rho*prim(i,j,lo(3)-k,2)
                      cons(i,j,lo(3)-k,3) = rho*prim(i,j,lo(3)-k,3)
                      cons(i,j,lo(3)-k,4) = rho*prim(i,j,lo(3)-k,4)

                      ! must be last BC enforced: depends on rho, vel, & temp
                      cons(i,j,lo(3)-k,5) = rho*intenergy + 0.5*rho*(prim(i,j,lo(3)-k,2)**2 + & 
                           prim(i,j,lo(3)-k,3)**2 + prim(i,j,lo(3)-k,4)**2)

                   enddo
                enddo
             enddo

          endif

       endif

       if(hi(3) .eq. (n_cells(3)-1)) then !upper z bound

          if(bc_mass_hi(3) .eq. 1) then ! wall

             if(algorithm_type.eq.2) then
                do k = 1,ngc(3)
                   do j = lo(2)-ngc(2),hi(2)+ngc(2)
                      do i = lo(1)-ngc(1),hi(1)+ngc(1)

                         prim(i,j,hi(3)+k,6:nprimvars) = prim(i,j,hi(3)+1-k,6:nprimvars)

                      enddo
                   enddo
                enddo
             endif

          elseif(bc_mass_hi(3) .eq. 2) then ! reservoir

             if(algorithm_type.eq.2) then

                Ywall(1:nspecies) = bc_Yk(3,2,1:nspecies)
                Xwall(1:nspecies) = bc_Xk(3,2,1:nspecies)

                do k = 1,ngc(3)
                   do j = lo(2)-ngc(2),hi(2)+ngc(2)
                      do i = lo(1)-ngc(1),hi(1)+ngc(1)

                         do l = 1, nspecies
                            prim(i,j,hi(3)+k,6+l)          = 2.d0*Ywall(l) - prim(i,j,hi(3)+1-k,6+l)
                            prim(i,j,hi(3)+k,6+nspecies+l) = 2.d0*Xwall(l) - prim(i,j,hi(3)+1-k,6+nspecies+l)
                         enddo

                      enddo
                   enddo
                enddo
             endif

          endif

          if(bc_therm_hi(3) .eq. 1) then ! adiabatic

             do k = 1,ngc(3)
                do j = lo(2)-ngc(2),hi(2)+ngc(2)
                   do i = lo(1)-ngc(1),hi(1)+ngc(1)

                      prim(i,j,hi(3)+k,5) = prim(i,j,hi(3)+1-k,5)
                      prim(i,j,hi(3)+k,6) = prim(i,j,hi(3)+1-k,6)

                   enddo
                enddo
             enddo

          elseif(bc_therm_hi(3) .eq. 2) then ! isothermal

             do k = 1,ngc(3)
                do j = lo(2)-ngc(2),hi(2)+ngc(2)
                   do i = lo(1)-ngc(1),hi(1)+ngc(1)

                      prim(i,j,hi(3)+k,5) = -prim(i,j,hi(3)+1-k,5) + 2*t_hi(3)
                      prim(i,j,hi(3)+k,6) = prim(i,j,hi(3)+1-k,6)

                   enddo
                enddo
             enddo

          endif

          if(bc_vel_hi(3) .eq. 1) then ! slip

             do k = 1,ngc(3)
                do j = lo(2)-ngc(2),hi(2)+ngc(2)
                   do i = lo(1)-ngc(1),hi(1)+ngc(1)

                      cons(i,j,hi(3)+k,2) = cons(i,j,hi(3)+1-k,2) 
                      cons(i,j,hi(3)+k,3) = cons(i,j,hi(3)+1-k,3) 
                      cons(i,j,hi(3)+k,4) = -cons(i,j,hi(3)+1-k,4) 

                      prim(i,j,hi(3)+k,2) = prim(i,j,hi(3)+1-k,2)
                      prim(i,j,hi(3)+k,3) = prim(i,j,hi(3)+1-k,3)
                      prim(i,j,hi(3)+k,4) = -prim(i,j,hi(3)+1-k,4)

                      ! thermal & species (+pressure) BCs must be enforced first
                      fracvec = prim(i,j,hi(3)+k,6+1:6+nspecies)
                      temp = prim(i,j,hi(3)+k,5)
                      pt = prim(i,j,hi(3)+k,6)

                      call get_density(pt, rho, temp, fracvec)
                      call get_energy(intenergy, fracvec, temp)
                      ! call get_density_gas(pt, rho, temp)
                      ! call get_energy_gas(pt, intenergy)

                      prim(i,j,hi(3)+k,1) = rho
                      cons(i,j,hi(3)+k,1) = rho
                      if(algorithm_type.eq.2) then
                         do l = 1, nspecies
                            cons(i,j,hi(3)+k,5+l) = rho*prim(i,j,hi(3)+k,6+l)
                         enddo
                      endif

                      ! must be last BC enforced: depends on rho, vel, & temp
                      cons(i,j,hi(3)+k,5) = rho*intenergy + 0.5*rho*(prim(i,j,hi(3)+k,2)**2 + & 
                           prim(i,j,hi(3)+k,3)**2 + prim(i,j,hi(3)+k,4)**2)

                   enddo
                enddo
             enddo

          elseif(bc_vel_hi(3) .eq. 2) then ! no slip

             do k = 1,ngc(3)
                do j = lo(2)-ngc(2),hi(2)+ngc(2)
                   do i = lo(1)-ngc(1),hi(1)+ngc(1)

                      prim(i,j,hi(3)+k,2) = -prim(i,j,hi(3)+1-k,2) 
                      prim(i,j,hi(3)+k,3) = -prim(i,j,hi(3)+1-k,3) 
                      prim(i,j,hi(3)+k,4) = -prim(i,j,hi(3)+1-k,4) ! + 2*vel_hi(3)

                      ! thermal & species (+pressure) BCs must be enforced first
                      fracvec = prim(i,j,hi(3)+k,6+1:6+nspecies)
                      temp = prim(i,j,hi(3)+k,5)
                      pt = prim(i,j,hi(3)+k,6)

                      call get_density(pt, rho, temp, fracvec)
                      call get_energy(intenergy, fracvec, temp)
                      ! call get_density_gas(pt, rho, temp)
                      ! call get_energy_gas(pt, intenergy)

                      prim(i,j,hi(3)+k,1) = rho
                      cons(i,j,hi(3)+k,1) = rho
                      if(algorithm_type.eq.2) then
                         do l = 1, nspecies
                            cons(i,j,hi(3)+k,5+l) = rho*prim(i,j,hi(3)+k,6+l)
                         enddo
                      endif

                      cons(i,j,hi(3)+k,2) = rho*prim(i,j,hi(3)+k,2)
                      cons(i,j,hi(3)+k,3) = rho*prim(i,j,hi(3)+k,3)
                      cons(i,j,hi(3)+k,4) = rho*prim(i,j,hi(3)+k,4)

                      ! must be last BC enforced: depends on rho, vel, & temp
                      cons(i,j,hi(3)+k,5) = rho*intenergy + 0.5*rho*(prim(i,j,hi(3)+k,2)**2 + & 
                           prim(i,j,hi(3)+k,3)**2 + prim(i,j,hi(3)+k,4)**2)

                   enddo
                enddo
             enddo

          endif

       endif

    endif

  end subroutine set_bc

  subroutine setup_bc() bind(C,name="setup_bc")

    integer :: d

    ! if bc_vel_lo/hi are periodic, set everything else to periodic
    do d=1,AMREX_SPACEDIM

       if (bc_vel_lo(d) .eq. -1) then
          bc_mass_lo(d) = -1
          bc_therm_lo(d) = -1
       end if

       if (bc_vel_hi(d) .eq. -1) then
          bc_mass_hi(d) = -1
          bc_therm_hi(d) = -1
       end if

    enddo

  end subroutine setup_bc

  subroutine setup_cwall() bind(C,name="setup_cwall")

    integer :: ns, d
    integer :: index, nsx, dx

    real(amrex_real) :: sumxt, sumyt, sumxb, sumyb

    ! Compute Xk or Yk at the wall, depending on which is defined
    do d=1,AMREX_SPACEDIM
       
       if(bc_mass_lo(d).eq.2)then

          sumxb = 0
          sumyb = 0

          do ns=1,nspecies

             sumxb = sumxb + bc_Xk(d,1,ns)
             sumyb = sumyb + bc_Yk(d,1,ns)

          enddo

          if(abs(sumxb-1).lt.1.d-10)then
             call get_massfrac(bc_Xk(d,1,1:nspecies),bc_Yk(d,1,1:nspecies))
          elseif(abs(sumyb-1).lt.1d-10)then
             call get_molfrac(bc_Yk(d,1,1:nspecies),bc_Xk(d,1,1:nspecies))
          endif

       endif
       
       if(bc_mass_hi(d).eq.2)then

          sumxt = 0
          sumyt = 0

          do ns=1,nspecies

             sumxt = sumxt + bc_Xk(d,2,ns)
             sumyt = sumyt + bc_Yk(d,2,ns)
 
          enddo

          if(abs(sumxt-1).lt.1.d-10)then
             call get_massfrac(bc_Xk(d,2,1:nspecies),bc_Yk(d,2,1:nspecies))
          elseif(abs(sumyt-1).lt.1d-10)then
             call get_molfrac(bc_Yk(d,2,1:nspecies),bc_Xk(d,2,1:nspecies))
          endif

       endif

    enddo

  end subroutine setup_cwall

end module bound_module

