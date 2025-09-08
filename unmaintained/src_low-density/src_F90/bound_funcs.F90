module bound_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, bc_vel_lo, bc_vel_hi, t_lo, t_hi, nprimvars, nvars, nspecies, n_cells, membrane_cell

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

      real(amrex_real) :: massvec(nspecies), fracvec(nspecies), intenergy, temp, rho, pt

      !Internal special case:

      bcell = membrane_cell

      if(bcell .ge. 0) then !Apply if membrane

         if(lo(1) .eq. bcell) then !Interior rhs, apply slip adiabatic

            !print *, "Bcell!"

            do k = lo(3)-ngc(3),hi(3)+ngc(3)
               do j = lo(2)-ngc(2),hi(2)+ngc(2)
                  do i = 1, ngc(1)

                     cons(lo(1)-i,j,k,1) = cons(lo(1)-1+i,j,k,1)

                  enddo
               enddo
            enddo

         endif


         if(hi(1) .eq. bcell-1) then !Interior lhs, apply slip adiabatic

            do k = lo(3)-ngc(3),hi(3)+ngc(3)
               do j = lo(2)-ngc(2),hi(2)+ngc(2)
                  do i = 1, ngc(1)

                     cons(hi(1)+i,j,k,1) = cons(hi(1)+1-i,j,k,1)

                  enddo
               enddo
            enddo
         endif

      endif

      if(lo(1) .eq. 0) then !lower x bound

        if(bc_vel_lo(1) .eq. 1) then ! slip adiabatic

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
            do j = lo(2)-ngc(2),hi(2)+ngc(2)
              do i = 1, ngc(1)

                cons(lo(1)-i,j,k,1) = cons(lo(1)-1+i,j,k,1)

              enddo
            enddo
          enddo

        elseif(bc_vel_lo(1) .eq. 2) then ! no slip thermal

          !print *, "Setting xLo thermal: "

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
            do j = lo(2)-ngc(2),hi(2)+ngc(2)
              do i = 1, ngc(1)

                cons(lo(1)-i,j,k,1) = cons(lo(1)-1+i,j,k,1)

              enddo
            enddo
          enddo


        endif
      endif

      if(hi(1) .eq. (n_cells(1)-1)) then !upper x bound

        if(bc_vel_hi(1) .eq. 1) then ! slip adiabatic

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
            do j = lo(2)-ngc(2),hi(2)+ngc(2)
              do i = 1, ngc(1)

                cons(hi(1)+i,j,k,1) = cons(hi(1)+1-i,j,k,1)

              enddo
            enddo
          enddo

        elseif(bc_vel_hi(1) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
            do j = lo(2)-ngc(2),hi(2)+ngc(2)
              do i = 1, ngc(1)

                 cons(hi(1)+i,j,k,1) = cons(hi(1)+1-i,j,k,1)

              enddo
            enddo
          enddo


        endif
      endif

       !print *, "bc: ", prim(1,0,0,1), prim(0,0,0,1), prim(-2,0,0,1), prim(-1,0,0,1)


      if(lo(2) .eq. 0) then !lower y bound

        if(bc_vel_lo(2) .eq. 1) then ! slip adiabatic

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
            do j = 1,ngc(2)
              do i = lo(1)-ngc(1),hi(1)+ngc(1)

                cons(i,lo(2)-j,k,1) = cons(i,lo(2)-1+j,k,1)

              enddo
            enddo
          enddo

        elseif(bc_vel_lo(2) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
            do j = 1,ngc(2)
              do i = lo(1)-ngc(1),hi(1)+ngc(1)

                cons(i,lo(2)-j,k,1) = cons(i,lo(2)-1+j,k,1)

              enddo
            enddo
          enddo


        endif
      endif

      if(hi(2) .eq. (n_cells(2)-1)) then !upper y bound

        if(bc_vel_hi(2) .eq. 1) then ! slip adiabatic

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
            do j = 1,ngc(2)
              do i = lo(1)-ngc(1),hi(1)+ngc(1)

                cons(i,hi(2)+j,k,1) = cons(i,hi(2)+1-j,k,1)

              enddo
            enddo
          enddo

        elseif(bc_vel_hi(2) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
            do j = 1,ngc(2)
              do i = lo(1)-ngc(1),hi(1)+ngc(1)

                cons(i,hi(2)+j,k,1) = cons(i,hi(2)+1-j,k,1)

              enddo
            enddo
          enddo


        endif
      endif

      if(lo(3) .eq. 0) then !lower z bound

        if(bc_vel_lo(3) .eq. 1) then ! slip adiabatic

          do k = 1,ngc(3)
            do j = lo(2)-ngc(2),hi(2)+ngc(2)
              do i = lo(1)-ngc(1),hi(1)+ngc(1)

                cons(i,j,lo(3)-k,1) = cons(i,j,lo(3)-1+k,1)

              enddo
            enddo
          enddo

        elseif(bc_vel_lo(3) .eq. 2) then ! no slip thermal

          do k = 1,ngc(3)
            do j = lo(2)-ngc(2),hi(2)+ngc(2)
              do i = lo(1)-ngc(1),hi(1)+ngc(1)

                cons(i,j,lo(3)-k,1) = cons(i,j,lo(3)-1+k,1)

              enddo
            enddo
          enddo


        endif
      endif

      if(hi(3) .eq. (n_cells(3)-1)) then !upper z bound

        if(bc_vel_hi(3) .eq. 1) then ! slip adiabatic

          do k = 1,ngc(3)
            do j = lo(2)-ngc(2),hi(2)+ngc(2)
              do i = lo(1)-ngc(1),hi(1)+ngc(1)

                cons(i,j,hi(3)+k,1) = cons(i,j,hi(3)+1-k,1)

              enddo
            enddo
          enddo

        elseif(bc_vel_hi(3) .eq. 2) then ! no slip thermal

          do k = 1,ngc(3)
            do j = lo(2)-ngc(2),hi(2)+ngc(2)
              do i = lo(1)-ngc(1),hi(1)+ngc(1)

                cons(i,j,hi(3)+k,1) = cons(i,j,hi(3)+1-k,1)

              enddo
            enddo
          enddo


        endif
      endif

  end subroutine set_bc

end module bound_module

