module conv_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, diameter, max_species, &
                                     molmass, k_b, nspecies, hcv, hcp, runiv, dof
  implicit none

  private

  public :: get_temperature, get_density, get_energy, get_molfrac, &
       get_massfrac, get_hc_gas, get_pressure_gas, get_density_gas, &
       get_energy_gas

contains

  subroutine get_temperature(energy, massfrac, temp)     

    !This function originaly had a reference to e0 - check this.

    real(amrex_real), intent(inout) :: temp
    real(amrex_real), intent(in   ) :: energy, massfrac(nspecies)

    integer :: i
    real(amrex_real) :: cvmix, e0

    cvmix = 0.0d0; e0 = 0.0d0

    do i = 1, nspecies
       cvmix = cvmix + massfrac(i)*hcv(i)
       ! e0 = e0 + massfrac(i)*e0ref(i)
    enddo

    temp = (energy-e0)/cvmix 

  end subroutine get_temperature

  subroutine get_density(pressure, density, temp, massfrac)  bind(C,name="get_density")    

    real(amrex_real), intent(in   ) :: temp, pressure, massfrac(nspecies)
    real(amrex_real), intent(inout) :: density

    integer :: i
    real(amrex_real) :: molmix

    molmix = 0.0d0
    do i = 1, nspecies
       molmix = molmix + massfrac(i)/molmass(i)
    enddo
    molmix = 1.0d0/molmix

    density = pressure/(runiv/molmix)/temp

  end subroutine get_density

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

  subroutine get_pressure_gas(pressure, fracvec, density, temp)  bind(C,name="get_pressure_gas")    

    real(amrex_real), intent(in   ) :: temp, fracvec(nspecies), density
    real(amrex_real), intent(inout) :: pressure

    integer :: i
    real(amrex_real) :: molmix

    molmix = 0.0d0
    do i = 1, nspecies
       molmix = molmix + fracvec(i)/molmass(i)
    enddo
    molmix = 1.0d0/molmix

    pressure = density*(runiv/molmix)*temp 

  end subroutine get_pressure_gas

  subroutine get_energy_gas(pressure, intenergy)  bind(C,name="get_energy_gas")    

    real(amrex_real), intent(in   ) :: pressure
    real(amrex_real), intent(inout) :: intenergy

    ! FIXME: this energy not scaled by density
    intenergy = pressure*3d0/2d0
    write(6,*) "called get_energy_gas, which is wrong"
    stop

  end subroutine get_energy_gas

  subroutine get_density_gas(pressure, density, temp)  bind(C,name="get_density_gas")    

    real(amrex_real), intent(in   ) :: temp, pressure
    real(amrex_real), intent(inout) :: density

    integer :: i
    real(amrex_real) :: avm

    avm = 0.0d0

    do i = 1, nspecies
       avm = avm + (1d0/nspecies)*molmass(i)

    enddo

    density = avm*pressure/(temp*runiv)

  end subroutine get_density_gas

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

