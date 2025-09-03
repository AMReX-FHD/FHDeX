module bc_fill_module

! since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
! for e.g., #if (BL_SPACEDIM == 1) statements.

  implicit none

  public

contains

  subroutine confill(con,con_lo,con_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="confill")

    use amrex_fort_module, only : bl_spacedim, amrex_real
    use amrex_filcc_module, only : amrex_filccn

    implicit none

    integer      :: con_lo(3),con_hi(3)
    integer      :: bc(bl_spacedim,2)
    integer      :: domlo(3), domhi(3)
    real(amrex_real) :: delta(3), xlo(3), time
    real(amrex_real) :: con(con_lo(1):con_hi(1),con_lo(2):con_hi(2),con_lo(3):con_hi(3))

    call amrex_filccn(con_lo, con_hi, con, con_lo, con_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine confill

end module bc_fill_module
