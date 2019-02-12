module bound_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, bc_lo, bc_hi, n_cells

  implicit none

  private

  public :: set_bc

contains

  subroutine set_bc(lo, hi, u_mac, v_mac,  w_mac ) bind(C,name="set_bc")

  integer, intent(in) :: lo(3), hi(3)

  real(amrex_real), intent(inout) :: u_mac(lo(1)-ngc(1):hi(1)+ngc(1), &
                                           lo(2)-ngc(2):hi(2)+ngc(2), &
                                           lo(3)-ngc(3):hi(3)+ngc(3))

  real(amrex_real), intent(inout) :: v_mac(lo(1)-ngc(1):hi(1)+ngc(1), &
                                           lo(2)-ngc(2):hi(2)+ngc(2), &
                                           lo(3)-ngc(3):hi(3)+ngc(3))

  real(amrex_real), intent(inout) :: w_mac(lo(1)-ngc(1):hi(1)+ngc(1), &
                                           lo(2)-ngc(2):hi(2)+ngc(2), &
                                           lo(3)-ngc(3):hi(3)+ngc(3))

  integer :: i,j,k

  if(hi(1) .eq. (n_cells(1)-1)) then !upper bound x

    if(bc_hi(1) .eq. 2) then ! no slip thermal

          do k = lo(3),hi(3)
            do j = lo(2),hi(2)

             u_mac(hi(1),j,k) = 0        

            end do
          end do    

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
            do j = lo(2)-ngc(2),hi(2)+ngc(2)
              do i = 1, ngc(1)

                u_mac(hi(1)+i,j,k)=-u_mac(hi(1)-i,j,k)                
                v_mac(hi(1)+i,j,k)=-v_mac(hi(1)+1-i,j,k)                
                w_mac(hi(1)+i,j,k)=-w_mac(hi(1)+1-i,j,k)                

              end do
            end do
          end do

    endif

  end if

  if(hi(2) .eq. (n_cells(2)-1)) then !upper bound y

    if(bc_hi(2) .eq. 2) then ! no slip thermal

          do k = lo(3),hi(3)
            do i = lo(1),hi(1)

             v_mac(i,hi(2),k) = 0        

            end do
          end do    

          do k = lo(3)-ngc(3),hi(3)+ngc(3)
            do j = lo(2)-ngc(2),hi(2)+ngc(2)
              do i = 1, ngc(1)

                v_mac(hi(1)+i,j,k)=-v_mac(hi(1)-i,j,k)                
                u_mac(hi(1)+i,j,k)=-u_mac(hi(1)+1-i,j,k)                
                w_mac(hi(1)+i,j,k)=-w_mac(hi(1)+1-i,j,k)                

              end do
            end do
          end do

    endif

  end if

  end subroutine set_bc

end module bound_module
