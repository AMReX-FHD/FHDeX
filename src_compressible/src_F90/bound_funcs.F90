module bound_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, bc_lo, bc_hi, t_lo, t_hi, nprimvars, nvars, nspecies, n_cells, membrane_cell
  use conv_module

  implicit none

  private

  public :: set_bc

contains

  subroutine set_bc(lo,hi, cons, prim, eta, zeta, kappa) bind(C,name="set_bc")

      integer         , intent(in   ) :: lo(3),hi(3)

      real(amrex_real), intent(inout) :: prim(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)
      real(amrex_real), intent(inout) :: cons(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)
      real(amrex_real), intent(inout) :: eta(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
      real(amrex_real), intent(inout) :: zeta(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))
      real(amrex_real), intent(inout) :: kappa(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3))

      integer :: i,j,k,l,bcell

      real(amrex_real) :: massvec(nspecies), intenergy

      !Internal special case:

      bcell = membrane_cell

      if(lo(1) .eq. bcell) then !Interior rhs, apply slip adiabatic

        !print *, "Bcell!"

        do k = lo(3)-ngc(3),hi(3)+ngc(3)
          do j = lo(2)-ngc(2),hi(2)+ngc(2)
            do i = 1, ngc(1)

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
              
              do l = 1, nspecies
                cons(lo(1)-i,j,k,5+l) = cons(lo(1)-1+i,j,k,5+l)
                prim(lo(1)-i,j,k,6+l) = prim(lo(1)-1+i,j,k,6+l)
              enddo

            enddo
          enddo
        enddo

      endif


      if(hi(1) .eq. bcell-1) then !Interior lhs, apply slip adiabatic

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
              
              do l = 1, nspecies
                cons(hi(1)+i,j,k,5+l) = cons(hi(1)+1-i,j,k,5+l)
                prim(hi(1)+i,j,k,6+l) = prim(hi(1)+1-i,j,k,6+l)
              enddo

            enddo
          enddo
        enddo
      endif

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
                
                do l = 1, nspecies
                  cons(lo(1)-i,j,k,5+l) = cons(lo(1)-1+i,j,k,5+l)
                  prim(lo(1)-i,j,k,6+l) = prim(lo(1)-1+i,j,k,6+l)
                enddo

              enddo
            enddo
          enddo

        elseif(bc_lo(1) .eq. 2) then ! no slip thermal

          !print *, "Setting xLo thermal: "

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
            do j = lo(2)-ngc(2),hi(2)+ngc(2)
              do i = 1, ngc(1)

                eta(lo(1)-i,j,k) = eta(lo(1)-1+i,j,k)
                zeta(lo(1)-i,j,k) = zeta(lo(1)-1+i,j,k)
                kappa(lo(1)-i,j,k) = kappa(lo(1)-1+i,j,k)

                prim(lo(1)-i,j,k,1) = prim(lo(1)-1+i,j,k,1)
                prim(lo(1)-i,j,k,2) = -prim(lo(1)-1+i,j,k,2) 
                prim(lo(1)-i,j,k,3) = -prim(lo(1)-1+i,j,k,3) 
                prim(lo(1)-i,j,k,4) = -prim(lo(1)-1+i,j,k,4)
                prim(lo(1)-i,j,k,5) = -prim(lo(1)-1+i,j,k,5) + 2*t_lo(1)
                prim(lo(1)-i,j,k,6) = prim(lo(1)-1+i,j,k,6)


                do l = 1, nspecies
                  cons(lo(1)-i,j,k,5+l) = cons(lo(1)-1+i,j,k,5+l)
                  prim(lo(1)-i,j,k,6+l) = prim(lo(1)-1+i,j,k,6+l)
                enddo

                massvec = cons(lo(1)-i,j,k,6:nvars)*cons(lo(1)-i,j,k,1)

                call get_energy(intenergy, massvec, prim(lo(1)-i,j,k,5))

                cons(lo(1)-i,j,k,5) = intenergy + 0.5*cons(lo(1)-i,j,k,1)*(cons(lo(1)-i,j,k,2)**2 + cons(lo(1)-i,j,k,3)**2 + cons(lo(1)-i,j,k,4)**2)

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
                
                do l = 1, nspecies
                  cons(hi(1)+i,j,k,5+l) = cons(hi(1)+1-i,j,k,5+l)
                  prim(hi(1)+i,j,k,6+l) = prim(hi(1)+1-i,j,k,6+l)
                enddo

              enddo
            enddo
          enddo

        elseif(bc_hi(1) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
            do j = lo(2)-ngc(2),hi(2)+ngc(2)
              do i = 1, ngc(1)

                eta(hi(1)+i,j,k) = eta(hi(1)+1-i,j,k)
                zeta(hi(1)+i,j,k) = zeta(hi(1)+1-i,j,k)
                kappa(hi(1)+i,j,k) = kappa(hi(1)+1-i,j,k)

                prim(hi(1)+i,j,k,1) = prim(hi(1)+1-i,j,k,1)
                prim(hi(1)+i,j,k,2) = -prim(hi(1)+1-i,j,k,2) 
                prim(hi(1)+i,j,k,3) = -prim(hi(1)+1-i,j,k,3) 
                prim(hi(1)+i,j,k,4) = -prim(hi(1)+1-i,j,k,4)
                prim(hi(1)+i,j,k,5) = -prim(hi(1)+1-i,j,k,5) + 2*t_hi(1)
                prim(hi(1)+i,j,k,6) = prim(hi(1)+1-i,j,k,6)

                cons(hi(1)+i,j,k,1) = cons(hi(1)+1-i,j,k,1)
                cons(hi(1)+i,j,k,2) = -cons(hi(1)+1-i,j,k,2) 
                cons(hi(1)+i,j,k,3) = -cons(hi(1)+1-i,j,k,3) 
                cons(hi(1)+i,j,k,4) = -cons(hi(1)+1-i,j,k,4)

                do l = 1, nspecies
                  cons(hi(1)+i,j,k,5+l) = cons(hi(1)+1-i,j,k,5+l)
                  prim(hi(1)+i,j,k,6+l) = prim(hi(1)+1-i,j,k,6+l)
                enddo

                massvec = cons(hi(1)+i,j,k,6:nvars)*cons(hi(1)+i,j,k,1)

                call get_energy(intenergy, massvec, prim(hi(1)+i,j,k,5))

                cons(hi(1)+i,j,k,5) = intenergy + 0.5*cons(hi(1)+i,j,k,1)*(cons(hi(1)+i,j,k,2)**2 + cons(hi(1)+i,j,k,3)**2 + cons(hi(1)+i,j,k,4)**2)


              enddo
            enddo
          enddo


        endif
      endif


!      if(lo(2) .eq. 0) then !lower y bound

!        if(bc_lo(2) .eq. 1) then ! slip adiabatic

!          do k = lo(3)-ngc(3),hi(3)+ngc(3)
!            do j = 1,ngc(2)
!              do i = lo(1)-ngc(1),hi(1)+ngc(1)

!                eta(i,lo(2)-j,k) = eta(i,lo(2)-1+j,k)
!                zeta(i,lo(2)-j,k) = zeta(i,lo(2)-1+j,k)
!                kappa(i,lo(2)-j,k) = kappa(i,lo(2)-1+j,k)             

!                cons(i,lo(2)-j,k,1) = cons(i,lo(2)-1+j,k,1)
!                cons(i,lo(2)-j,k,2) = -cons(i,lo(2)-1+j,k,2) 
!                cons(i,lo(2)-j,k,3) = cons(i,lo(2)-1+j,k,3) 
!                cons(i,lo(2)-j,k,4) = cons(i,lo(2)-1+j,k,4) 
!                cons(i,lo(2)-j,k,5) = cons(i,lo(2)-1+j,k,5) 
!               
!                prim(i,lo(2)-j,k,1) = prim(i,lo(2)-1+j,k,1)
!                prim(i,lo(2)-j,k,2) = -prim(i,lo(2)-1+j,k,2)
!                prim(i,lo(2)-j,k,3) = prim(i,lo(2)-1+j,k,3)
!                prim(i,lo(2)-j,k,4) = prim(i,lo(2)-1+j,k,4)
!                prim(i,lo(2)-j,k,5) = prim(i,lo(2)-1+j,k,5)
!                prim(i,lo(2)-j,k,6) = prim(i,lo(2)-1+j,k,6)
!                
!                do l = 1, nspecies
!                  cons(i,lo(2)-j,k,5+l) = cons(i,lo(2)-1+j,k,5+l)
!                  prim(i,lo(2)-j,k,6+l) = prim(i,lo(2)-1+j,k,6+l)
!                enddo

!              enddo
!            enddo
!          enddo

!        elseif(bc_lo(2) .eq. 2) then ! no slip thermal

!          do k = lo(3)-ngc(3),hi(3)+ngc(3)
!            do j = 1,ngc(2)
!              do i = lo(1)-ngc(1),hi(1)+ngc(1)

!                eta(i,lo(2)-j,k) = eta(i,lo(2)-1+j,k)
!                zeta(i,lo(2)-j,k) = zeta(i,lo(2)-1+j,k)
!                kappa(i,lo(2)-j,k) = kappa(i,lo(2)-1+j,k)

!                prim(i,lo(2)-j,k,1) = prim(i,lo(2)-1+j,k,1)
!                prim(i,lo(2)-j,k,2) = -prim(i,lo(2)-1+j,k,2) 
!                prim(i,lo(2)-j,k,3) = -prim(i,lo(2)-1+j,k,3) 
!                prim(i,lo(2)-j,k,4) = -prim(i,lo(2)-1+j,k,4)
!                prim(i,lo(2)-j,k,5) = -prim(i,lo(2)-1+j,k,5) + 2*t_lo(2)

!                cons(i,lo(2)-j,k,1) = cons(i,lo(2)-1+j,k,1)
!                cons(i,lo(2)-j,k,2) = -cons(i,lo(2)-1+j,k,2) 
!                cons(i,lo(2)-j,k,3) = -cons(i,lo(2)-1+j,k,3) 
!                cons(i,lo(2)-j,k,4) = -cons(i,lo(2)-1+j,k,4)

!                do l = 1, nspecies
!                  cons(i,lo(2)-j,k,5+l) = cons(i,lo(2)-1+j,k,5+l)
!                  prim(i,lo(2)-j,k,6+l) = prim(i,lo(2)-1+j,k,6+l)
!                enddo

!                massvec = cons(i,lo(2)-j,k,6:nvars)*cons(i,lo(2)-j,k,1)

!                call get_energy(intenergy, massvec, prim(i,lo(2)-j,k,5))

!                cons(i,lo(2)-j,k,5) = intenergy + 0.5*cons(i,lo(2)-j,k,1)*(cons(i,lo(2)-j,k,2)**2 + cons(i,lo(2)-j,k,3)**2 + cons(i,lo(2)-j,k,4)**2)


!              enddo
!            enddo
!          enddo


!        endif
!      endif

!      if(hi(2) .eq. (n_cells(2)-1)) then !upper y bound

!        if(bc_hi(2) .eq. 1) then ! slip adiabatic

!          do k = lo(3)-ngc(3),hi(3)+ngc(3)
!            do j = 1,ngc(2)
!              do i = lo(1)-ngc(1),hi(1)+ngc(1)

!                eta(i,hi(2)+j,k) = eta(i,hi(2)+1-j,k)
!                zeta(i,hi(2)+j,k) = zeta(i,hi(2)+1-j,k)
!                kappa(i,hi(2)+j,k) = kappa(i,hi(2)+1-j,k)

!                cons(i,hi(2)+j,k,1) = cons(i,hi(2)+1-j,k,1)
!                cons(i,hi(2)+j,k,2) = -cons(i,hi(2)+1-j,k,2) 
!                cons(i,hi(2)+j,k,3) = cons(i,hi(2)+1-j,k,3) 
!                cons(i,hi(2)+j,k,4) = cons(i,hi(2)+1-j,k,4) 
!                cons(i,hi(2)+j,k,5) = cons(i,hi(2)+1-j,k,5) 
!               
!                prim(i,hi(2)+j,k,1) = prim(i,hi(2)+1-j,k,1)
!                prim(i,hi(2)+j,k,2) = -prim(i,hi(2)+1-j,k,2)
!                prim(i,hi(2)+j,k,3) = prim(i,hi(2)+1-j,k,3)
!                prim(i,hi(2)+j,k,4) = prim(i,hi(2)+1-j,k,4)
!                prim(i,hi(2)+j,k,5) = prim(i,hi(2)+1-j,k,5)
!                prim(i,hi(2)+j,k,6) = prim(i,hi(2)+1-j,k,6)
!                
!                do l = 1, nspecies
!                  cons(i,hi(2)+j,k,5+l) = cons(i,hi(2)+1-j,k,5+l)
!                  prim(i,hi(2)+j,k,6+l) = prim(i,hi(2)+1-j,k,6+l)
!                enddo

!              enddo
!            enddo
!          enddo

!        elseif(bc_hi(2) .eq. 2) then ! no slip thermal

!          do k = lo(3)-ngc(3),hi(3)+ngc(3)
!            do j = 1,ngc(2)
!              do i = lo(1)-ngc(1),hi(1)+ngc(1)

!                eta(i,hi(2)+j,k) = eta(i,hi(2)+1-j,k)              
!                zeta(i,hi(2)+j,k) = zeta(i,hi(2)+1-j,k)
!                kappa(i,hi(2)+j,k) = kappa(i,hi(2)+1-j,k)

!                prim(i,hi(2)+j,k,1) = prim(i,hi(2)+1-j,k,1)
!                prim(i,hi(2)+j,k,2) = -prim(i,hi(2)+1-j,k,2) 
!                prim(i,hi(2)+j,k,3) = -prim(i,hi(2)+1-j,k,3) 
!                prim(i,hi(2)+j,k,4) = -prim(i,hi(2)+1-j,k,4)
!                prim(i,hi(2)+j,k,5) = -prim(i,hi(2)+1-j,k,5) + 2*t_hi(2)

!                cons(i,hi(2)+j,k,1) = cons(i,hi(2)+1-j,k,1)
!                cons(i,hi(2)+j,k,2) = -cons(i,hi(2)+1-j,k,2) 
!                cons(i,hi(2)+j,k,3) = -cons(i,hi(2)+1-j,k,3) 
!                cons(i,hi(2)+j,k,4) = -cons(i,hi(2)+1-j,k,4)

!                do l = 1, nspecies
!                  cons(i,hi(2)+j,k,5+l) = cons(i,hi(2)+1-j,k,5+l)
!                  prim(i,hi(2)+j,k,6+l) = prim(i,hi(2)+1-j,k,6+l)
!                enddo

!                massvec = cons(i,hi(2)+j,k,6:nvars)*cons(i,hi(2)+j,k,1)

!                call get_energy(intenergy, massvec, prim(i,hi(2)+j,k,5))

!                cons(i,hi(2)+j,k,5) = intenergy + 0.5*cons(i,hi(2)+j,k,1)*(cons(i,hi(2)+j,k,2)**2 + cons(i,hi(2)+j,k,3)**2 + cons(i,hi(2)+j,k,4)**2)


!              enddo
!            enddo
!          enddo


!        endif
!      endif

!      if(lo(3) .eq. 0) then !lower z bound

!        if(bc_lo(3) .eq. 1) then ! slip adiabatic

!          do k = 1,ngc(3)
!            do j = lo(2)-ngc(2),hi(2)+ngc(2)
!              do i = lo(1)-ngc(1),hi(1)+ngc(1)

!                eta(i,j,lo(3)-k) = eta(i,j,lo(3)-1+k)
!                zeta(i,j,lo(3)-k) = zeta(i,j,lo(3)-1+k)
!                kappa(i,j,lo(3)-k) = kappa(i,j,lo(3)-1+k)              

!                cons(i,j,lo(3)-k,1) = cons(i,j,lo(3)-1+k,1)
!                cons(i,j,lo(3)-k,2) = -cons(i,j,lo(3)-1+k,2) 
!                cons(i,j,lo(3)-k,3) = cons(i,j,lo(3)-1+k,3) 
!                cons(i,j,lo(3)-k,4) = cons(i,j,lo(3)-1+k,4) 
!                cons(i,j,lo(3)-k,5) = cons(i,j,lo(3)-1+k,5) 
!               
!                prim(i,j,lo(3)-k,1) = prim(i,j,lo(3)-1+k,1)
!                prim(i,j,lo(3)-k,2) = -prim(i,j,lo(3)-1+k,2)
!                prim(i,j,lo(3)-k,3) = prim(i,j,lo(3)-1+k,3)
!                prim(i,j,lo(3)-k,4) = prim(i,j,lo(3)-1+k,4)
!                prim(i,j,lo(3)-k,5) = prim(i,j,lo(3)-1+k,5)
!                prim(i,j,lo(3)-k,6) = prim(i,j,lo(3)-1+k,6)
!                
!                do l = 1, nspecies
!                  cons(i,j,lo(3)-k,5+l) = cons(i,j,lo(3)-1+k,5+l)
!                  prim(i,j,lo(3)-k,6+l) = prim(i,j,lo(3)-1+k,6+l)
!                enddo

!              enddo
!            enddo
!          enddo

!        elseif(bc_lo(3) .eq. 2) then ! no slip thermal

!          do k = 1,ngc(3)
!            do j = lo(2)-ngc(2),hi(2)+ngc(2)
!              do i = lo(1)-ngc(1),hi(1)+ngc(1)

!                eta(i,j,lo(3)-k) = eta(i,j,lo(3)-1+k)             
!                zeta(i,j,lo(3)-k) = zeta(i,j,lo(3)-1+k)               
!                kappa(i,j,lo(3)-k) = kappa(i,j,lo(3)-1+k)               

!                prim(i,j,lo(3)-k,1) = prim(i,j,lo(3)-1+k,1)
!                prim(i,j,lo(3)-k,2) = -prim(i,j,lo(3)-1+k,2) 
!                prim(i,j,lo(3)-k,3) = -prim(i,j,lo(3)-1+k,3) 
!                prim(i,j,lo(3)-k,4) = -prim(i,j,lo(3)-1+k,4)
!                prim(i,j,lo(3)-k,5) = -prim(i,j,lo(3)-1+k,5) + 2*t_lo(3)

!                cons(i,j,lo(3)-k,1) = cons(i,j,lo(3)-1+k,1)
!                cons(i,j,lo(3)-k,2) = -cons(i,j,lo(3)-1+k,2) 
!                cons(i,j,lo(3)-k,3) = -cons(i,j,lo(3)-1+k,3) 
!                cons(i,j,lo(3)-k,4) = -cons(i,j,lo(3)-1+k,4)

!                do l = 1, nspecies
!                  cons(i,j,lo(3)-k,5+l) = cons(i,j,lo(3)-1+k,5+l)
!                  prim(i,j,lo(3)-k,6+l) = prim(i,j,lo(3)-1+k,6+l)
!                enddo

!                massvec = cons(i,j,lo(3)-k,6:nvars)*cons(i,j,lo(3)-k,1)

!                call get_energy(intenergy, massvec, prim(i,j,lo(3)-k,5))

!                cons(i,j,lo(3)-k,5) = intenergy + 0.5*cons(i,j,lo(3)-k,1)*(cons(i,j,lo(3)-k,2)**2 + cons(i,j,lo(3)-k,3)**2 + cons(i,j,lo(3)-k,4)**2)


!              enddo
!            enddo
!          enddo


!        endif
!      endif

!      if(hi(3) .eq. (n_cells(3)-1)) then !upper z bound

!        if(bc_hi(3) .eq. 1) then ! slip adiabatic

!          do k = 1,ngc(3)
!            do j = lo(2)-ngc(2),hi(2)+ngc(2)
!              do i = lo(1)-ngc(1),hi(1)+ngc(1)

!                eta(i,j,hi(3)+k) = eta(i,j,hi(3)+1-k)           
!                zeta(i,j,hi(3)+k) = zeta(i,j,hi(3)+1-k)            
!                kappa(i,j,hi(3)+k) = kappa(i,j,hi(3)+1-k)           

!                cons(i,j,hi(3)+k,1) = cons(i,j,hi(3)+1-k,1)
!                cons(i,j,hi(3)+k,2) = -cons(i,j,hi(3)+1-k,2) 
!                cons(i,j,hi(3)+k,3) = cons(i,j,hi(3)+1-k,3) 
!                cons(i,j,hi(3)+k,4) = cons(i,j,hi(3)+1-k,4) 
!                cons(i,j,hi(3)+k,5) = cons(i,j,hi(3)+1-k,5) 
!               
!                prim(i,j,hi(3)+k,1) = prim(i,j,hi(3)+1-k,1)
!                prim(i,j,hi(3)+k,2) = -prim(i,j,hi(3)+1-k,2)
!                prim(i,j,hi(3)+k,3) = prim(i,j,hi(3)+1-k,3)
!                prim(i,j,hi(3)+k,4) = prim(i,j,hi(3)+1-k,4)
!                prim(i,j,hi(3)+k,5) = prim(i,j,hi(3)+1-k,5)
!                prim(i,j,hi(3)+k,6) = prim(i,j,hi(3)+1-k,6)
!                
!                do l = 1, nspecies
!                  cons(i,j,hi(3)+k,5+l) = cons(i,j,hi(3)+1-k,5+l)
!                  prim(i,j,hi(3)+k,6+l) = prim(i,j,hi(3)+1-k,6+l)
!                enddo

!              enddo
!            enddo
!          enddo

!        elseif(bc_hi(3) .eq. 2) then ! no slip thermal

!          do k = 1,ngc(3)
!            do j = lo(2)-ngc(2),hi(2)+ngc(2)
!              do i = lo(1)-ngc(1),hi(1)+ngc(1)

!                eta(i,j,hi(3)+k) = eta(i,j,hi(3)+1-k)
!                zeta(i,j,hi(3)+k) = zeta(i,j,hi(3)+1-k)
!                kappa(i,j,hi(3)+k) = kappa(i,j,hi(3)+1-k)

!                prim(i,j,hi(3)+k,1) = prim(i,j,hi(3)+1-k,1)
!                prim(i,j,hi(3)+k,2) = -prim(i,j,hi(3)+1-k,2) 
!                prim(i,j,hi(3)+k,3) = -prim(i,j,hi(3)+1-k,3) 
!                prim(i,j,hi(3)+k,4) = -prim(i,j,hi(3)+1-k,4)
!                prim(i,j,hi(3)+k,5) = -prim(i,j,hi(3)+1-k,5) + 2*t_hi(3)

!                cons(i,j,hi(3)+k,1) = cons(i,j,hi(3)+1-k,1)
!                cons(i,j,hi(3)+k,2) = -cons(i,j,hi(3)+1-k,2) 
!                cons(i,j,hi(3)+k,3) = -cons(i,j,hi(3)+1-k,3) 
!                cons(i,j,hi(3)+k,4) = -cons(i,j,hi(3)+1-k,4)

!                do l = 1, nspecies
!                  cons(i,j,hi(3)+k,5+l) = cons(i,j,hi(3)+1-k,5+l)
!                  prim(i,j,hi(3)+k,6+l) = prim(i,j,hi(3)+1-k,6+l)
!                enddo

!                massvec = cons(i,j,hi(3)+k,6:nvars)*cons(i,j,hi(3)+k,1)

!                call get_energy(intenergy, massvec, prim(i,j,hi(3)+k,5))

!                cons(i,j,hi(3)+k,5) = intenergy + 0.5*cons(i,j,hi(3)+k,1)*(cons(i,j,hi(3)+k,2)**2 + cons(i,j,hi(3)+k,3)**2 + cons(i,j,hi(3)+k,4)**2)


!              enddo
!            enddo
!          enddo


!        endif
!      endif

  end subroutine set_bc

end module bound_module

