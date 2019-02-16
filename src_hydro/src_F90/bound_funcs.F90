module bound_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, bc_lo, bc_hi, n_cells

  implicit none

  private

  public :: set_bc
  public :: fab_physbc

contains


  subroutine fab_physbc(lo, hi, pressure) bind(C, name="fab_physbc")

    integer,          intent(in   ) :: lo(3), hi(3)
    real(amrex_real), intent(inout) :: pressure(lo(1)-ngc(1):hi(1)+ngc(1), &
         &                                      lo(2)-ngc(2):hi(2)+ngc(2), &
         &                                      lo(3)-ngc(3):hi(3)+ngc(3))


    ! ** loop indices
    integer :: i,j,k


    !____________________________________________________________________________
    ! Apply BC to X faces

    if (lo(1) .eq. 0) then ! lower bound
       if(bc_lo(1) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc(3), hi(3)+ngc(3)
             do j = lo(2)-ngc(2), hi(2)+ngc(2)
                do i = 1, ngc(1)

                   pressure(lo(1)-i, j, k) = pressure(lo(1)-1+i, j, k)

                end do
             end do
          end do

       end if
    end if

    if(hi(1) .eq. (n_cells(1)-1)) then ! upper bound
       if(bc_hi(1) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc(3), hi(3)+ngc(3)
             do j = lo(2)-ngc(2), hi(2)+ngc(2)
                do i = 1, ngc(1)

                   pressure(hi(1)+i, j, k) = pressure(hi(1)+1-i, j, k)

                end do
             end do
          end do

       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Y faces

    if(lo(2) .eq. 0) then ! lower bound
       if(bc_lo(2) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc(3), hi(3)+ngc(3)
             do j = 1, ngc(2)
                do i = lo(1)-ngc(1), hi(1)+ngc(1)

                   pressure(i, lo(2)-j, k) = pressure(i, lo(2)-1+j, k)

                end do
             end do
          end do

       end if
    end if

    if(hi(2) .eq. (n_cells(2)-1)) then ! upper bound
       if(bc_hi(2) .eq. 2) then ! no slip thermal

          do k = lo(3)-ngc(3), hi(3)+ngc(3)
             do j = 1, ngc(2)
                do i = lo(1)-ngc(1), hi(1)+ngc(1)

                   pressure(i, hi(2)+j, k) = pressure(i, hi(2)+1-j, k)

                end do
             end do
          end do

       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Z faces

    if(lo(3) .eq. 0) then ! lower bound
       if(bc_lo(3) .eq. 2) then ! no slip thermal

          do k = 1, ngc(3)
             do j = lo(2)-ngc(2), hi(2)+ngc(2)
                do i = lo(1)-ngc(1), hi(1)+ngc(1)

                   pressure(i, j, lo(3)-k) = pressure(i, j, lo(3)-1+k)

                end do
             end do
          end do

       end if
    end if

    if(hi(3) .eq. (n_cells(3)-1)) then ! upper bound
       if(bc_hi(3) .eq. 2) then ! no slip thermal

          do k = 1, ngc(3)
             do j = lo(2)-ngc(2), hi(2)+ngc(2)
                do i = lo(1)-ngc(1), hi(1)+ngc(1)

                   pressure(i, j, hi(3)+k) = pressure(i, j, hi(3)+1-k)

                end do
             end do
          end do

       end if
    end if

  end subroutine fab_physbc



  subroutine set_bc(lo, hi, u_mac, v_mac,  w_mac ) bind(C, name="set_bc")

    integer,          intent(in   ) :: lo(3), hi(3)

    real(amrex_real), intent(inout) :: u_mac(lo(1)-ngc(1):hi(1)+ngc(1), &
         &                                   lo(2)-ngc(2):hi(2)+ngc(2), &
         &                                   lo(3)-ngc(3):hi(3)+ngc(3))
    real(amrex_real), intent(inout) :: v_mac(lo(1)-ngc(1):hi(1)+ngc(1), &
         &                                   lo(2)-ngc(2):hi(2)+ngc(2), &
         &                                   lo(3)-ngc(3):hi(3)+ngc(3))
    real(amrex_real), intent(inout) :: w_mac(lo(1)-ngc(1):hi(1)+ngc(1), &
         &                                   lo(2)-ngc(2):hi(2)+ngc(2), &
         &                                   lo(3)-ngc(3):hi(3)+ngc(3))


    ! ** loop indices
    integer :: i,j,k


    !____________________________________________________________________________
    ! Apply BC to X faces

    if(lo(1) .eq. 0) then ! lower bound
       if(bc_lo(1) .eq. 2) then ! no slip thermal

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)

                u_mac(lo(1), j, k) = 0

             end do
          end do

          do k = lo(3)-ngc(3), hi(3)+ngc(3)
             do j = lo(2)-ngc(2), hi(2)+ngc(2)
                do i = 1, ngc(1)

                   ! Indices about x-face are symmetric
                   u_mac(lo(1)-i, j, k) = -u_mac(lo(1)+i, j, k)

                   v_mac(lo(1)-i, j, k) = -v_mac(lo(1)-1+i, j, k)
                   w_mac(lo(1)-i, j, k) = -w_mac(lo(1)-1+i, j, k)

                end do
             end do
          end do

       end if
    end if

    if(hi(1) .eq. (n_cells(1)-1)) then ! upper bound
       if(bc_hi(1) .eq. 2) then ! no slip thermal

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)

                u_mac(hi(1), j, k) = 0

             end do
          end do

          do k = lo(3)-ngc(3), hi(3)+ngc(3)
             do j = lo(2)-ngc(2), hi(2)+ngc(2)
                do i = 1, ngc(1)

                   ! Indices about x-face are symmetric
                   u_mac(hi(1)+i, j, k) = -u_mac(hi(1)-i, j, k)

                   v_mac(hi(1)+i, j, k) = -v_mac(hi(1)+1-i, j, k)
                   w_mac(hi(1)+i, j, k) = -w_mac(hi(1)+1-i, j, k)

                end do
             end do
          end do

       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Y faces

    if(lo(2) .eq. 0) then ! lower bound
       if(bc_lo(2) .eq. 2) then ! no slip thermal

          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

                v_mac(i, lo(2), k) = 0

             end do
          end do

          do k = lo(3)-ngc(3), hi(3)+ngc(3)
             do j = 1, ngc(2)
                do i = lo(1)-ngc(1), hi(1)+ngc(1)

                   ! Indices about y-face are symmetric
                   v_mac(i, lo(2)-j, k) = -v_mac(i, lo(2)+j, k)

                   u_mac(i, lo(2)-j, k) = -u_mac(i, lo(2)-1+j, k)
                   w_mac(i, lo(2)-j, k) = -w_mac(i, lo(2)-1+j, k)

                end do
             end do
          end do

       end if
    end if

    if(hi(2) .eq. (n_cells(2)-1)) then ! upper bound
       if(bc_hi(2) .eq. 2) then ! no slip thermal

          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

                v_mac(i, hi(2), k) = 0

             end do
          end do

          do k = lo(3)-ngc(3), hi(3)+ngc(3)
             do j = 1, ngc(2)
                do i = lo(1)-ngc(1), hi(1)+ngc(1)

                   ! Indices about y-face are symmetric
                   v_mac(i, hi(2)+j, k) = -v_mac(i, hi(2)-j, k)

                   u_mac(i, hi(2)+j, k) = -u_mac(i, hi(2)+1-j, k)
                   w_mac(i, hi(2)+j, k) = -w_mac(i, hi(2)+1-j, k)

                end do
             end do
          end do

       end if
    end if


    !____________________________________________________________________________
    ! Apply BC to Z faces

    if(lo(3) .eq. 0) then ! lower bound
       if(bc_lo(3) .eq. 2) then ! no slip thermal

          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                w_mac(i, j, lo(3)) = 0

             end do
          end do

          do k = 1, ngc(3)
             do j = lo(2)-ngc(2), hi(2)+ngc(2)
                do i = lo(1)-ngc(1), hi(1)+ngc(1)

                   ! Indices about z-face are symmetric
                   w_mac(i, j, lo(3)-k) = -w_mac(i, j, lo(3)+k)

                   u_mac(i, j, lo(3)-k) = -u_mac(i, j, lo(3)-1+k)
                   v_mac(i, j, lo(3)-k) = -v_mac(i, j, lo(3)-1+k)

                end do
             end do
          end do

       end if
    end if

    if(hi(3) .eq. (n_cells(3)-1)) then ! upper bound
       if(bc_hi(3) .eq. 2) then ! no slip thermal

          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                w_mac(i, j, hi(3)) = 0

             end do
          end do

          do k = 1, ngc(3)
             do j = lo(2)-ngc(2), hi(2)+ngc(2)
                do i = lo(1)-ngc(1), hi(1)+ngc(1)

                   ! Indices about z-face are symmetric
                   w_mac(i, j, hi(3)+k) = -w_mac(i, j, hi(3)-k)

                   u_mac(i, j, hi(3)+k) = -u_mac(i, j, hi(3)+1-k)
                   v_mac(i, j, hi(3)+k) = -v_mac(i, j, hi(3)+1-k)

                end do
             end do
          end do

       end if
    end if

  end subroutine set_bc

end module bound_module
