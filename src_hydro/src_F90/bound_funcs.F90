module bound_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, bc_lo, bc_hi, n_cells

  implicit none

  private

  public :: set_bc

contains

  subroutine set_bc(lo, hi, u_mac, v_mac, w_mac) bind(C,name="set_bc")

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


  end subroutine set_bc

end module bound_module
