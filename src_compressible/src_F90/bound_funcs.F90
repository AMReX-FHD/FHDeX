module bound_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, bc_lo, bc_hi, t_lo, t_hi, nprimvars, nvars, nspecies, n_cells, algorithm_type, membrane_cell
  use conv_module
  use trans_module

  implicit none

  private

  public :: set_bc

contains

  subroutine set_bc(lo,hi, cons, prim, eta, zeta, kappa, chi, Dij) bind(C,name="set_bc")

    integer         , intent(in   ) :: lo(3),hi(3)

    real(amrex_real), intent(inout) :: prim(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)
    real(amrex_real), intent(inout) :: cons(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
    real(amrex_real), intent(inout) :: eta(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
    real(amrex_real), intent(inout) :: zeta(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
    real(amrex_real), intent(inout) :: kappa(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
    real(amrex_real), intent(inout) :: chi(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nspecies)
    real(amrex_real), intent(inout) :: Dij(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3),nspecies,nspecies)

    integer :: i,j,k,l,idir,bcell

    real(amrex_real) :: massvec(nspecies), fracvec(nspecies), intenergy, temp, rho, pt

    if(lo(1) .eq. 0) then !lower x bound

       if(bc_lo(1) .eq. 1) then ! slip adiabatic

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = lo(2)-ngc(2),hi(2)+ngc(2)
                do i = 1, ngc(1)          

                   !print *, "Setting xLo slip: "

                   eta(lo(1)-i,j,k) = eta(lo(1)-1+i,j,k)
                   zeta(lo(1)-i,j,k) = zeta(lo(1)-1+i,j,k)
                   kappa(lo(1)-i,j,k) = kappa(lo(1)-1+i,j,k)

                   cons(lo(1)-i,j,k,1) = cons(lo(1)-1+i,j,k,1)
                   cons(lo(1)-i,j,k,2) = -cons(lo(1)-1+i,j,k,2) 
                   cons(lo(1)-i,j,k,3) = cons(lo(1)-1+i,j,k,3) 
                   cons(lo(1)-i,j,k,4) = cons(lo(1)-1+i,j,k,4) 
                   cons(lo(1)-i,j,k,5) = cons(lo(1)-1+i,j,k,5) 

                   prim(lo(1)-i,j,k,1) = prim(lo(1)-1+i,j,k,1)
                   prim(lo(1)-i,j,k,2) = -prim(lo(1)-1+i,j,k,2)
                   prim(lo(1)-i,j,k,3) = prim(lo(1)-1+i,j,k,3)
                   prim(lo(1)-i,j,k,4) = prim(lo(1)-1+i,j,k,4)
                   prim(lo(1)-i,j,k,5) = prim(lo(1)-1+i,j,k,5)
                   prim(lo(1)-i,j,k,6) = prim(lo(1)-1+i,j,k,6)
                   
                   if(algorithm_type.eq.2) then
                      chi(lo(1)-i,j,k,:) = chi(lo(1)-1+i,j,k,:)
                      Dij(lo(1)-i,j,k,:,:) = Dij(lo(1)-1+i,j,k,:,:)

                      cons(lo(1)-i,j,k,5:nvars) = cons(lo(1)-1+i,j,k,5:nvars)
                      prim(lo(1)-i,j,k,6:nprimvars) = prim(lo(1)-1+i,j,k,6:nprimvars)
                   endif

                enddo
             enddo
          enddo

       elseif(bc_lo(1) .eq. 2) then ! no slip thermal

          !print *, "Setting xLo thermal: "

          idir = 0 
          call setcwall(Xwall,Ywall,idir)

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = lo(2)-ngc(2),hi(2)+ngc(2)
                do i = 1,ngc(1)

                   eta(lo(1)-i,j,k) = eta(lo(1)-1+i,j,k)
                   zeta(lo(1)-i,j,k) = zeta(lo(1)-1+i,j,k)
                   kappa(lo(1)-i,j,k) = kappa(lo(1)-1+i,j,k)

                   prim(lo(1)-i,j,k,2) = -prim(lo(1)-1+i,j,k,2) 
                   prim(lo(1)-i,j,k,3) = -prim(lo(1)-1+i,j,k,3) 
                   prim(lo(1)-i,j,k,4) = -prim(lo(1)-1+i,j,k,4)
                   prim(lo(1)-i,j,k,5) = -prim(lo(1)-1+i,j,k,5) + 2*t_lo(1)
                   prim(lo(1)-i,j,k,6) = prim(lo(1)-1+i,j,k,6)

                   if(algorithm_type.eq.2) then
                      Dij(lo(1)-i,j,k,:,:) = Dij(lo(1)-1+i,j,k,:,:)
                      chi(lo(1)-i,j,k,:) = chi(lo(1)-1+i,j,k,:)

                      do l = 1, nspecies
                         prim(lo(1)-i,j,k,6+l)          = 2.d0*Ywall(l) - prim(lo(1)-1+i,j,k,6+l)
                         prim(lo(1)-i,j,k,6+nspecies+l) = 2.d0*Xwall(l) - prim(lo(1)-1+i,j,k,6+nspecies+l)
                      enddo
                   else
                      print *,'Error: need access to eos for low x bc ', bc_lo(1)
                      stop
                   end if
                   
                   fracvec = prim(lo(1)-i,j,k,6+1:6+nspecies)
                   temp = prim(lo(1)-i,j,k,5)
                   pt = prim(lo(1)-i,j,k,6)

                   call get_density(pt, rho, temp, fracvec)
                   call get_energy(intenergy, fracvec, temp)
                   ! call get_density_gas(pt, rho, temp)
                   ! call get_energy_gas(pt, intenergy)

                   prim(lo(1)-i,j,k,1) = rho
                   cons(lo(1)-i,j,k,1) = rho
                   cons(lo(1)-i,j,k,2) = rho*prim(lo(1)-i,j,k,2)
                   cons(lo(1)-i,j,k,3) = rho*prim(lo(1)-i,j,k,3)
                   cons(lo(1)-i,j,k,4) = rho*prim(lo(1)-i,j,k,4)
                   cons(lo(1)-i,j,k,5) = rho*intenergy + 0.5*rho*(prim(lo(1)-i,j,k,2)**2 + & 
                        prim(lo(1)-i,j,k,3)**2 + prim(lo(1)-i,j,k,4)**2)
                   
                   if(algorithm_type.eq.2) then
                      do l = 1, nspecies
                         cons(lo(1)-i,j,k,5+l) = rho*prim(lo(1)-i,j,k,6+l)
                      enddo
                   endif

                enddo
             enddo
          enddo


       endif
    endif

    if(hi(1) .eq. (n_cells(1)-1)) then !upper x bound

       if(bc_hi(1) .eq. 1) then ! slip adiabatic

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = lo(2)-ngc(2),hi(2)+ngc(2)
                do i = 1, ngc(1)

                   eta(hi(1)+i,j,k) = eta(hi(1)+1-i,j,k)
                   zeta(hi(1)+i,j,k) = zeta(hi(1)+1-i,j,k)
                   kappa(hi(1)+i,j,k) = kappa(hi(1)+1-i,j,k)    

                   cons(hi(1)+i,j,k,1) = cons(hi(1)+1-i,j,k,1)
                   cons(hi(1)+i,j,k,2) = -cons(hi(1)+1-i,j,k,2) 
                   cons(hi(1)+i,j,k,3) = cons(hi(1)+1-i,j,k,3) 
                   cons(hi(1)+i,j,k,4) = cons(hi(1)+1-i,j,k,4) 
                   cons(hi(1)+i,j,k,5) = cons(hi(1)+1-i,j,k,5) 

                   prim(hi(1)+i,j,k,1) = prim(hi(1)+1-i,j,k,1)
                   prim(hi(1)+i,j,k,2) = -prim(hi(1)+1-i,j,k,2)
                   prim(hi(1)+i,j,k,3) = prim(hi(1)+1-i,j,k,3)
                   prim(hi(1)+i,j,k,4) = prim(hi(1)+1-i,j,k,4)
                   prim(hi(1)+i,j,k,5) = prim(hi(1)+1-i,j,k,5)
                   prim(hi(1)+i,j,k,6) = prim(hi(1)+1-i,j,k,6)

                   if(algorithm_type.eq.2) then
                      chi(hi(1)+i,j,k,:) = chi(hi(1)+1-i,j,k,:)         
                      Dij(hi(1)+i,j,k,:,:) = Dij(hi(1)+1-i,j,k,:,:)  
                      
                      cons(hi(1)+i,j,k,5:nvars) = cons(hi(1)+1-i,j,k,5:nvars)
                      prim(hi(1)+i,j,k,6:nprimvars) = prim(hi(1)+1-i,j,k,6:nprimvars)
                   endif
                      
                enddo
             enddo
          enddo

       elseif(bc_hi(1) .eq. 2) then ! no slip thermal

          idir = 1
          call setcwall(Xwall,Ywall,idir)

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = lo(2)-ngc(2),hi(2)+ngc(2)
                do i = 1, ngc(1)

                   eta(hi(1)+i,j,k) = eta(hi(1)+1-i,j,k)
                   zeta(hi(1)+i,j,k) = zeta(hi(1)+1-i,j,k)
                   kappa(hi(1)+i,j,k) = kappa(hi(1)+1-i,j,k) 

                   prim(hi(1)+i,j,k,2) = -prim(hi(1)+1-i,j,k,2) 
                   prim(hi(1)+i,j,k,3) = -prim(hi(1)+1-i,j,k,3) 
                   prim(hi(1)+i,j,k,4) = -prim(hi(1)+1-i,j,k,4)
                   prim(hi(1)+i,j,k,5) = -prim(hi(1)+1-i,j,k,5) + 2*t_hi(1)
                   prim(hi(1)+i,j,k,6) = prim(hi(1)+1-i,j,k,6)

                   if(algorithm_type.eq.2) then
                      Dij(hi(1)+i,j,k,:,:) = Dij(hi(1)+1-i,j,k,:,:)
                      chi(hi(1)+i,j,k,:) = chi(hi(1)+1-i,j,k,:)

                      do l = 1, nspecies
                         prim(hi(1)+i,j,k,6+l)          = 2.d0*Ywall(l) - prim(hi(1)+1-i,j,k,6+l)
                         prim(hi(1)+i,j,k,6+nspecies+l) = 2.d0*Xwall(l) - prim(hi(1)+1-i,j,k,6+nspecies+l)
                      enddo
                   else
                      print *,'Error: need access to eos for low x bc ', bc_lo(1)
                      stop
                   end if
                   
                   fracvec = prim(hi(1)+i,j,k,6+1:6+nspecies)
                   temp = prim(hi(1)+i,j,k,5)
                   pt = prim(hi(1)+i,j,k,6)

                   call get_density(pt, rho, temp, fracvec)
                   call get_energy(intenergy, fracvec, temp)
                   ! call get_density_gas(pt, rho, temp)
                   ! call get_energy_gas(pt, intenergy)

                   prim(hi(1)+i,j,k,1) = rho
                   cons(hi(1)+i,j,k,1) = rho
                   cons(hi(1)+i,j,k,2) = rho*prim(hi(1)+i,j,k,2)
                   cons(hi(1)+i,j,k,3) = rho*prim(hi(1)+i,j,k,3)
                   cons(hi(1)+i,j,k,4) = rho*prim(hi(1)+i,j,k,4)
                   cons(hi(1)+i,j,k,5) = rho*intenergy + 0.5*rho*(prim(hi(1)+i,j,k,2)**2 + & 
                        prim(hi(1)+i,j,k,3)**2 + prim(hi(1)+i,j,k,4)**2)
                   
                   if(algorithm_type.eq.2) then
                      do l = 1, nspecies
                         cons(hi(1)+i,j,k,5+l) = rho*prim(hi(1)+i,j,k,6+l)
                      enddo
                   endif

                enddo
             enddo
          enddo


       endif
    endif

    !print *, "bc: ", prim(1,0,0,1), prim(0,0,0,1), prim(-2,0,0,1), prim(-1,0,0,1)


    if(lo(2) .eq. 0) then !lower y bound

       if(bc_lo(2) .eq. 1) then ! slip adiabatic

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = 1,ngc(2)
                do i = lo(1)-ngc(1),hi(1)+ngc(1)

                   eta(i,lo(2)-j,k) = eta(i,lo(2)-1+j,k)
                   zeta(i,lo(2)-j,k) = zeta(i,lo(2)-1+j,k)
                   kappa(i,lo(2)-j,k) = kappa(i,lo(2)-1+j,k)           

                   cons(i,lo(2)-j,k,1) = cons(i,lo(2)-1+j,k,1)
                   cons(i,lo(2)-j,k,2) = -cons(i,lo(2)-1+j,k,2) 
                   cons(i,lo(2)-j,k,3) = cons(i,lo(2)-1+j,k,3) 
                   cons(i,lo(2)-j,k,4) = cons(i,lo(2)-1+j,k,4) 
                   cons(i,lo(2)-j,k,5) = cons(i,lo(2)-1+j,k,5) 

                   prim(i,lo(2)-j,k,1) = prim(i,lo(2)-1+j,k,1)
                   prim(i,lo(2)-j,k,2) = -prim(i,lo(2)-1+j,k,2)
                   prim(i,lo(2)-j,k,3) = prim(i,lo(2)-1+j,k,3)
                   prim(i,lo(2)-j,k,4) = prim(i,lo(2)-1+j,k,4)
                   prim(i,lo(2)-j,k,5) = prim(i,lo(2)-1+j,k,5)
                   prim(i,lo(2)-j,k,6) = prim(i,lo(2)-1+j,k,6)

                   if(algorithm_type.eq.2) then
                      chi(i,lo(2)-j,k,:) = chi(i,lo(2)-1+j,k,:)
                      Dij(i,lo(2)-j,k,:,:) = Dij(i,lo(2)-1+j,k,:,:)  

                      cons(i,lo(2)-j,k,5:nvars) = cons(i,lo(2)-1+j,k,5:nvars)
                      prim(i,lo(2)-j,k,6:nprimvars) = prim(i,lo(2)-1+j,k,6:nprimvars)
                   endif

                enddo
             enddo
          enddo

       elseif(bc_lo(2) .eq. 2) then ! no slip thermal

          idir = 0
          call setcwall(Xwall,Ywall,idir)

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = 1,ngc(2)
                do i = lo(1)-ngc(1),hi(1)+ngc(1)

                   eta(i,lo(2)-j,k) = eta(i,lo(2)-1+j,k)
                   zeta(i,lo(2)-j,k) = zeta(i,lo(2)-1+j,k)
                   kappa(i,lo(2)-j,k) = kappa(i,lo(2)-1+j,k)

                   prim(i,lo(2)-j,k,2) = -prim(i,lo(2)-1+j,k,2) 
                   prim(i,lo(2)-j,k,3) = -prim(i,lo(2)-1+j,k,3) 
                   prim(i,lo(2)-j,k,4) = -prim(i,lo(2)-1+j,k,4)
                   prim(i,lo(2)-j,k,5) = -prim(i,lo(2)-1+j,k,5) + 2*t_lo(1)
                   prim(i,lo(2)-j,k,6) = prim(i,lo(2)-1+j,k,6)

                   if(algorithm_type.eq.2) then
                      Dij(i,lo(2)-j,k,:,:) = Dij(i,lo(2)-1+j,k,:,:)
                      chi(i,lo(2)-j,k,:) = chi(i,lo(2)-1+j,k,:)

                      do l = 1, nspecies
                         prim(i,lo(2)-j,k,6+l)          = 2.d0*Ywall(l) - prim(i,lo(2)-1+j,k,6+l)
                         prim(i,lo(2)-j,k,6+nspecies+l) = 2.d0*Xwall(l) - prim(i,lo(2)-1+j,k,6+nspecies+l)
                      enddo
                   else
                      print *,'Error: need access to eos for low x bc ', bc_lo(1)
                      stop
                   end if
                   
                   fracvec = prim(i,lo(2)-j,k,6+1:6+nspecies)
                   temp = prim(i,lo(2)-j,k,5)
                   pt = prim(i,lo(2)-j,k,6)

                   call get_density(pt, rho, temp, fracvec)
                   call get_energy(intenergy, fracvec, temp)
                   ! call get_density_gas(pt, rho, temp)
                   ! call get_energy_gas(pt, intenergy)

                   prim(i,lo(2)-j,k,1) = rho
                   cons(i,lo(2)-j,k,1) = rho
                   cons(i,lo(2)-j,k,2) = rho*prim(i,lo(2)-j,k,2)
                   cons(i,lo(2)-j,k,3) = rho*prim(i,lo(2)-j,k,3)
                   cons(i,lo(2)-j,k,4) = rho*prim(i,lo(2)-j,k,4)
                   cons(i,lo(2)-j,k,5) = rho*intenergy + 0.5*rho*(prim(i,lo(2)-j,k,2)**2 + & 
                        prim(i,lo(2)-j,k,3)**2 + prim(i,lo(2)-j,k,4)**2)
                   
                   if(algorithm_type.eq.2) then
                      do l = 1, nspecies
                         cons(i,lo(2)-j,k,5+l) = rho*prim(i,lo(2)-j,k,6+l)
                      enddo
                   endif

                enddo
             enddo
          enddo


       endif
    endif

    if(hi(2) .eq. (n_cells(2)-1)) then !upper y bound

       if(bc_hi(2) .eq. 1) then ! slip adiabatic

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = 1,ngc(2)
                do i = lo(1)-ngc(1),hi(1)+ngc(1)

                   eta(i,hi(2)+j,k) = eta(i,hi(2)+1-j,k)
                   zeta(i,hi(2)+j,k) = zeta(i,hi(2)+1-j,k)
                   kappa(i,hi(2)+j,k) = kappa(i,hi(2)+1-j,k)

                   cons(i,hi(2)+j,k,1) = cons(i,hi(2)+1-j,k,1)
                   cons(i,hi(2)+j,k,2) = -cons(i,hi(2)+1-j,k,2) 
                   cons(i,hi(2)+j,k,3) = cons(i,hi(2)+1-j,k,3) 
                   cons(i,hi(2)+j,k,4) = cons(i,hi(2)+1-j,k,4) 
                   cons(i,hi(2)+j,k,5) = cons(i,hi(2)+1-j,k,5) 

                   prim(i,hi(2)+j,k,1) = prim(i,hi(2)+1-j,k,1)
                   prim(i,hi(2)+j,k,2) = -prim(i,hi(2)+1-j,k,2)
                   prim(i,hi(2)+j,k,3) = prim(i,hi(2)+1-j,k,3)
                   prim(i,hi(2)+j,k,4) = prim(i,hi(2)+1-j,k,4)
                   prim(i,hi(2)+j,k,5) = prim(i,hi(2)+1-j,k,5)
                   prim(i,hi(2)+j,k,6) = prim(i,hi(2)+1-j,k,6)

                   if(algorithm_type.eq.2) then
                      chi(i,hi(2)+j,k,:) = chi(i,hi(2)+1-j,k,:)
                      Dij(i,hi(2)+j,k,:,:) = Dij(i,hi(2)+1-j,k,:,:)

                      cons(i,hi(2)+j,k,5:nvars) = cons(i,hi(2)+1-j,k,5:nvars)
                      prim(i,hi(2)+j,k,6:nprimvars) = prim(i,hi(2)+1-j,k,6:nprimvars)
                   endif

                enddo
             enddo
          enddo

       elseif(bc_hi(2) .eq. 2) then ! no slip thermal

          idir = 1
          call setcwall(Xwall,Ywall,idir)

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
             do j = 1,ngc(2)
                do i = lo(1)-ngc(1),hi(1)+ngc(1)

                   eta(i,hi(2)+j,k) = eta(i,hi(2)+1-j,k)              
                   zeta(i,hi(2)+j,k) = zeta(i,hi(2)+1-j,k)
                   kappa(i,hi(2)+j,k) = kappa(i,hi(2)+1-j,k)

                   prim(i,hi(2)+j,k,1) = prim(i,hi(2)+1-j,k,1)
                   prim(i,hi(2)+j,k,2) = -prim(i,hi(2)+1-j,k,2) 
                   prim(i,hi(2)+j,k,3) = -prim(i,hi(2)+1-j,k,3) 
                   prim(i,hi(2)+j,k,4) = -prim(i,hi(2)+1-j,k,4)
                   prim(i,hi(2)+j,k,5) = -prim(i,hi(2)+1-j,k,5) + 2*t_hi(2)

                   if(algorithm_type.eq.2) then
                      Dij(i,hi(2)+j,k,:,:) = Dij(i,hi(2)+1-j,k,:,:)
                      chi(i,hi(2)+j,k,:) = chi(i,hi(2)+1-j,k,:)

                      do l = 1, nspecies
                         prim(i,hi(2)+j,k,6+l)          = 2.d0*Ywall(l) - prim(i,hi(2)+1-j,k,6+l)
                         prim(i,hi(2)+j,k,6+nspecies+l) = 2.d0*Xwall(l) - prim(i,hi(2)+1-j,k,6+nspecies+l)
                      enddo
                   else
                      print *,'Error: need access to eos for low x bc ', bc_lo(1)
                      stop
                   end if
                   
                   fracvec = prim(i,hi(2)+j,k,6+1:6+nspecies)
                   temp = prim(i,hi(2)+j,k,5)
                   pt = prim(i,hi(2)+j,k,6)

                   call get_density(pt, rho, temp, fracvec)
                   call get_energy(intenergy, fracvec, temp)
                   ! call get_density_gas(pt, rho, temp)
                   ! call get_energy_gas(pt, intenergy)

                   prim(i,hi(2)+j,k,1) = rho
                   cons(i,hi(2)+j,k,1) = rho
                   cons(i,hi(2)+j,k,2) = rho*prim(i,hi(2)+j,k,2)
                   cons(i,hi(2)+j,k,3) = rho*prim(i,hi(2)+j,k,3)
                   cons(i,hi(2)+j,k,4) = rho*prim(i,hi(2)+j,k,4)
                   cons(i,hi(2)+j,k,5) = rho*intenergy + 0.5*rho*(prim(i,hi(2)+j,k,2)**2 + & 
                        prim(i,hi(2)+j,k,3)**2 + prim(i,hi(2)+j,k,4)**2)
                   
                   if(algorithm_type.eq.2) then
                      do l = 1, nspecies
                         cons(i,hi(2)+j,k,5+l) = rho*prim(i,hi(2)+j,k,6+l)
                      enddo
                   endif

                enddo
             enddo
          enddo


       endif
    endif

    if(lo(3) .eq. 0) then !lower z bound

       if(bc_lo(3) .eq. 1) then ! slip adiabatic

          do k = 1,ngc(3)
             do j = lo(2)-ngc(2),hi(2)+ngc(2)
                do i = lo(1)-ngc(1),hi(1)+ngc(1)

                   eta(i,j,lo(3)-k) = eta(i,j,lo(3)-1+k)
                   zeta(i,j,lo(3)-k) = zeta(i,j,lo(3)-1+k)
                   kappa(i,j,lo(3)-k) = kappa(i,j,lo(3)-1+k)  

                   cons(i,j,lo(3)-k,1) = cons(i,j,lo(3)-1+k,1)
                   cons(i,j,lo(3)-k,2) = -cons(i,j,lo(3)-1+k,2) 
                   cons(i,j,lo(3)-k,3) = cons(i,j,lo(3)-1+k,3) 
                   cons(i,j,lo(3)-k,4) = cons(i,j,lo(3)-1+k,4) 
                   cons(i,j,lo(3)-k,5) = cons(i,j,lo(3)-1+k,5) 

                   prim(i,j,lo(3)-k,1) = prim(i,j,lo(3)-1+k,1)
                   prim(i,j,lo(3)-k,2) = -prim(i,j,lo(3)-1+k,2)
                   prim(i,j,lo(3)-k,3) = prim(i,j,lo(3)-1+k,3)
                   prim(i,j,lo(3)-k,4) = prim(i,j,lo(3)-1+k,4)
                   prim(i,j,lo(3)-k,5) = prim(i,j,lo(3)-1+k,5)
                   prim(i,j,lo(3)-k,6) = prim(i,j,lo(3)-1+k,6)

                   if(algorithm_type.eq.2) then
                      chi(i,j,lo(3)-k,:) = chi(i,j,lo(3)-1+k,:)
                      Dij(i,j,lo(3)-k,:,:) = Dij(i,j,lo(3)-1+k,:,:)            

                      cons(i,j,lo(3)-k,5:nvars) = cons(i,j,lo(3)-1+k,5:nvars)
                      prim(i,j,lo(3)-k,6:nprimvars) = prim(i,j,lo(3)-1+k,6:nprimvars)
                   endif

                enddo
             enddo
          enddo

       elseif(bc_lo(3) .eq. 2) then ! no slip thermal

          idir = 0
          call setcwall(Xwall,Ywall,idir)

          do k = 1,ngc(3)
             do j = lo(2)-ngc(2),hi(2)+ngc(2)
                do i = lo(1)-ngc(1),hi(1)+ngc(1)

                   eta(i,j,lo(3)-k) = eta(i,j,lo(3)-1+k)
                   zeta(i,j,lo(3)-k) = zeta(i,j,lo(3)-1+k)
                   kappa(i,j,lo(3)-k) = kappa(i,j,lo(3)-1+k)

                   prim(i,j,lo(3)-k,2) = -prim(i,j,lo(3)-1+k,2) 
                   prim(i,j,lo(3)-k,3) = -prim(i,j,lo(3)-1+k,3) 
                   prim(i,j,lo(3)-k,4) = -prim(i,j,lo(3)-1+k,4)
                   prim(i,j,lo(3)-k,5) = -prim(i,j,lo(3)-1+k,5) + 2*t_lo(1)
                   prim(i,j,lo(3)-k,6) = prim(i,j,lo(3)-1+k,6)

                   if(algorithm_type.eq.2) then
                      Dij(i,j,lo(3)-k,:,:) = Dij(i,j,lo(3)-1+k,:,:)
                      chi(i,j,lo(3)-k) = chi(i,j,lo(3)-1+k,:)

                      do l = 1, nspecies
                         prim(i,j,lo(3)-k,6+l)          = 2.d0*Ywall(l) - prim(i,j,lo(3)-1+k,6+l)
                         prim(i,j,lo(3)-k,6+nspecies+l) = 2.d0*Xwall(l) - prim(i,j,lo(3)-1+k,6+nspecies+l)
                      enddo
                   else
                      print *,'Error: need access to eos for low x bc ', bc_lo(1)
                      stop
                   end if
                   
                   fracvec = prim(i,j,lo(3)-k,6+1:6+nspecies)
                   temp = prim(i,j,lo(3)-k,5)
                   pt = prim(i,j,lo(3)-k,6)

                   call get_density(pt, rho, temp, fracvec)
                   call get_energy(intenergy, fracvec, temp)
                   ! call get_density_gas(pt, rho, temp)
                   ! call get_energy_gas(pt, intenergy)

                   prim(i,j,lo(3)-k,1) = rho
                   cons(i,j,lo(3)-k,1) = rho
                   cons(i,j,lo(3)-k,2) = rho*prim(i,j,lo(3)-k,2)
                   cons(i,j,lo(3)-k,3) = rho*prim(i,j,lo(3)-k,3)
                   cons(i,j,lo(3)-k,4) = rho*prim(i,j,lo(3)-k,4)
                   cons(i,j,lo(3)-k,5) = rho*intenergy + 0.5*rho*(prim(i,j,lo(3)-k,2)**2 + & 
                        prim(i,j,lo(3)-k,3)**2 + prim(i,j,lo(3)-k,4)**2)
                   
                   if(algorithm_type.eq.2) then
                      do l = 1, nspecies
                         cons(i,j,lo(3)-k,5+l) = rho*prim(i,j,lo(3)-k,6+l)
                      enddo
                   endif

                enddo
             enddo
          enddo


       endif
    endif

    if(hi(3) .eq. (n_cells(3)-1)) then !upper z bound

       if(bc_hi(3) .eq. 1) then ! slip adiabatic

          do k = 1,ngc(3)
             do j = lo(2)-ngc(2),hi(2)+ngc(2)
                do i = lo(1)-ngc(1),hi(1)+ngc(1)

                   eta(i,j,hi(3)+k) = eta(i,j,hi(3)+1-k)           
                   zeta(i,j,hi(3)+k) = zeta(i,j,hi(3)+1-k)            
                   kappa(i,j,hi(3)+k) = kappa(i,j,hi(3)+1-k)           

                   cons(i,j,hi(3)+k,1) = cons(i,j,hi(3)+1-k,1)
                   cons(i,j,hi(3)+k,2) = -cons(i,j,hi(3)+1-k,2) 
                   cons(i,j,hi(3)+k,3) = cons(i,j,hi(3)+1-k,3) 
                   cons(i,j,hi(3)+k,4) = cons(i,j,hi(3)+1-k,4) 
                   cons(i,j,hi(3)+k,5) = cons(i,j,hi(3)+1-k,5) 

                   prim(i,j,hi(3)+k,1) = prim(i,j,hi(3)+1-k,1)
                   prim(i,j,hi(3)+k,2) = -prim(i,j,hi(3)+1-k,2)
                   prim(i,j,hi(3)+k,3) = prim(i,j,hi(3)+1-k,3)
                   prim(i,j,hi(3)+k,4) = prim(i,j,hi(3)+1-k,4)
                   prim(i,j,hi(3)+k,5) = prim(i,j,hi(3)+1-k,5)
                   prim(i,j,hi(3)+k,6) = prim(i,j,hi(3)+1-k,6)

                   if(algorithm_type.eq.2) then
                      chi(i,j,hi(3)+k,:) = chi(i,j,hi(3)+1-k,:)           
                      Dij(i,j,hi(3)+k,:,:) = Dij(i,j,hi(3)+1-k,:,:)           

                      cons(i,j,hi(3)+k,5:nvars) = cons(i,j,hi(3)+1-k,5:nvars)
                      prim(i,j,hi(3)+k,6:nprimvars) = prim(i,j,hi(3)+1-k,6:nprimvars)
                   endif

                enddo
             enddo
          enddo

       elseif(bc_hi(3) .eq. 2) then ! no slip thermal

          idir = 1
          call setcwall(Xwall,Ywall,idir)

          do k = 1,ngc(3)
             do j = lo(2)-ngc(2),hi(2)+ngc(2)
                do i = lo(1)-ngc(1),hi(1)+ngc(1)

                   eta(i,j,hi(3)+k) = eta(i,j,hi(3)+1-k)
                   zeta(i,j,hi(3)+k) = zeta(i,j,hi(3)+1-k)
                   kappa(i,j,hi(3)+k) = kappa(i,j,hi(3)+1-k)

                   prim(i,j,hi(3)+k,1) = prim(i,j,hi(3)+1-k,1)
                   prim(i,j,hi(3)+k,2) = -prim(i,j,hi(3)+1-k,2) 
                   prim(i,j,hi(3)+k,3) = -prim(i,j,hi(3)+1-k,3) 
                   prim(i,j,hi(3)+k,4) = -prim(i,j,hi(3)+1-k,4)
                   prim(i,j,hi(3)+k,5) = -prim(i,j,hi(3)+1-k,5) + 2*t_hi(3)

                   if(algorithm_type.eq.2) then
                      Dij(i,j,hi(3)+k,:,:) = Dij(i,j,hi(3)+1-k,:,:)
                      chi(i,j,hi(3)+k,:) = chi(i,j,hi(3)+1-k,:)

                      do l = 1, nspecies
                         prim(i,j,hi(3)+k,6+l)          = 2.d0*Ywall(l) - prim(i,j,hi(3)+1-k,6+l)
                         prim(i,j,hi(3)+k,6+nspecies+l) = 2.d0*Xwall(l) - prim(i,j,hi(3)+1-k,6+nspecies+l)
                      enddo
                   else
                      print *,'Error: need access to eos for low x bc ', bc_lo(1)
                      stop
                   end if
                   
                   fracvec = prim(i,j,hi(3)+k,6+1:6+nspecies)
                   temp = prim(i,j,hi(3)+k,5)
                   pt = prim(i,j,hi(3)+k,6)

                   call get_density(pt, rho, temp, fracvec)
                   call get_energy(intenergy, fracvec, temp)
                   ! call get_density_gas(pt, rho, temp)
                   ! call get_energy_gas(pt, intenergy)

                   prim(i,j,hi(3)+k,1) = rho
                   cons(i,j,hi(3)+k,1) = rho
                   cons(i,j,hi(3)+k,2) = rho*prim(i,j,hi(3)+k,2)
                   cons(i,j,hi(3)+k,3) = rho*prim(i,j,hi(3)+k,3)
                   cons(i,j,hi(3)+k,4) = rho*prim(i,j,hi(3)+k,4)
                   cons(i,j,hi(3)+k,5) = rho*intenergy + 0.5*rho*(prim(i,j,hi(3)+k,2)**2 + & 
                        prim(i,j,hi(3)+k,3)**2 + prim(i,j,hi(3)+k,4)**2)
                   
                   if(algorithm_type.eq.2) then
                      do l = 1, nspecies
                         cons(i,j,hi(3)+k,5+l) = rho*prim(i,j,hi(3)+k,6+l)
                      enddo
                   endif

                enddo
             enddo
          enddo


       endif
    endif

  end subroutine set_bc

end module bound_module

