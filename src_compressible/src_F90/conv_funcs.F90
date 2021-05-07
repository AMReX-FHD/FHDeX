module conv_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : molmass, nspecies, hcv, hcp, runiv, dof
  implicit none

  private

  public :: get_molfrac, get_massfrac, get_hc_gas

contains

  subroutine get_energy(energy, massvec, temp)  bind(C,name="get_energy")    

    !This function originaly had a reference to e0 - check this.

    real(amrex_real), intent(in   ) :: temp
    real(amrex_real), intent(inout) :: energy, massvec(nspecies)

    integer :: i
    real(amrex_real) :: cvmix, e0

    cvmix = 0.0d0; e0 = 0.0d0

    do i = 1, nspecies
       cvmix = cvmix + massvec(i)*hcv(i)
       ! e0 = e0 + massvec(i)*e0ref(i)
    enddo

    energy = e0 + temp*cvmix 

  end subroutine get_energy
  
  subroutine get_molfrac(Yk, Xk)

    real(amrex_real), intent(inout) :: Xk(nspecies)
    real(amrex_real), intent(in   ) :: Yk(nspecies)

    integer :: ns
    real(amrex_real) :: molmix

    molmix = 0.0d0
    do ns = 1, nspecies
       molmix = molmix + Yk(ns)/molmass(ns)
    enddo
    molmix = 1.0d0/molmix
    do ns = 1, nspecies
       Xk(ns) = Yk(ns)*(molmix/molmass(ns))
    enddo

  end subroutine get_molfrac

  subroutine get_massfrac(Xk, Yk)

    real(amrex_real), intent(inout) :: Yk(nspecies)
    real(amrex_real), intent(in   ) :: Xk(nspecies)

    integer :: ns
    real(amrex_real) :: molmix

    molmix = 0.0d0
    do ns = 1, nspecies
       molmix = molmix + Xk(ns)*molmass(ns)
    enddo
    do ns = 1, nspecies
       Yk(ns) = Xk(ns)*(molmass(ns)/molmix)
    enddo

  end subroutine get_massfrac

  subroutine get_hc_gas() bind(C,name="get_hc_gas")

    integer :: i

    do i=1, nspecies
       if(hcv(i) .lt. 0) then   
          hcv(i) = 0.5d0*dof(i)*Runiv/molmass(i)
       endif
       if (hcp(i) .lt. 0) then
          hcp(i) = 0.5d0*(2+dof(i))*Runiv/molmass(i)
       end if
    enddo

  end subroutine get_hc_gas

end module conv_module

