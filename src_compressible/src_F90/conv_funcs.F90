module conv_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, diameter, max_species, &
       molmass, k_b, nspecies, hcv, hcp, runiv, dof
  implicit none

  private

  public :: cons_to_prim, get_temperature, get_density, get_energy, get_molfrac, &
       get_massfrac, get_enthalpies, get_hc_gas, get_pressure_gas, get_density_gas, &
       get_temperature_gas, get_energy_gas

contains

  subroutine cons_to_prim(lo,hi, cons, prim) bind(C,name="cons_to_prim")

    integer         , intent(in   ) :: lo(3),hi(3)

    real(amrex_real), intent(inout) :: prim(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nprimvars)
    real(amrex_real), intent(in   ) :: cons(lo(1)-ngc(1):hi(1)+ngc(1),lo(2)-ngc(2):hi(2)+ngc(2),lo(3)-ngc(3):hi(3)+ngc(3), nvars)

    integer :: i,j,k,ns

    real(amrex_real) :: vsqr, sumYk, Yk(nspecies), Yk_fixed(nspecies), Xk(nspecies), intenergy

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             prim(i,j,k,1) = cons(i,j,k,1)
             prim(i,j,k,2) = cons(i,j,k,2)/cons(i,j,k,1)
             prim(i,j,k,3) = cons(i,j,k,3)/cons(i,j,k,1)
             prim(i,j,k,4) = cons(i,j,k,4)/cons(i,j,k,1)

             vsqr = prim(i,j,k,2)**2 + prim(i,j,k,3)**2 + prim(i,j,k,4)**2
             intenergy = cons(i,j,k,5)/cons(i,j,k,1) - 0.5*vsqr

             sumYk = 0.d0
             do ns = 1, nspecies
                Yk(ns) = cons(i,j,k,5+ns)/cons(i,j,k,1)
                Yk_fixed(ns) = max(0.d0,min(1.d0,Yk(ns)))
                sumYk = sumYk + Yk_fixed(ns)

                ! if((i.eq.0).and.(j.eq.0).and.(k.eq.0)) then
                !    print *, "massfrac: Hack = ", i, j, k, ns, Yk(ns), Yk_fixed(ns)
                ! endif
             enddo

             Yk_fixed(:) = Yk_fixed(:)/sumYk

             ! Yk_fixed(:) = Yk(:)

             ! update temperature in-place using internal energy
             call get_temperature(intenergy, Yk_fixed, prim(i,j,k,5))

             ! HACK: OVERRIDE
             ! prim(i,j,k,6) = 2.0*cons(i,j,k,1)*intenergy/3.0

             ! compute mole fractions from mass fractions
             call get_molfrac(Yk, Xk)

             ! mass fractions
             do ns = 1, nspecies
                prim(i,j,k,6+ns) = Yk(ns)
                prim(i,j,k,6+nspecies+ns) = Xk(ns)
             enddo

             call get_pressure_gas(prim(i,j,k,6), Yk, prim(i,j,k,1), prim(i,j,k,5))

          enddo
       enddo
    enddo

  end subroutine cons_to_prim

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

  subroutine get_enthalpies(T, hk) bind(C,name="get_enthalpies")

    real(amrex_real), intent(inout) :: hk(1:nspecies)
    real(amrex_real), intent(in   ) :: T

    integer :: ns

    !! FIXME: should have enthalpy reference eventually

    do ns = 1, nspecies 
       ! hk(ns) = h0ref(ns) + hcp(ns)*T
       hk(ns) = 0.0 + hcp(ns)*T
    enddo

  end subroutine get_enthalpies

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

  end subroutine get_energy_gas

  subroutine get_temperature_gas(pressure, fracvec, density, temp)  bind(C,name="get_temperature_gas")    

    real(amrex_real), intent(in   ) :: pressure, fracvec(nspecies), density
    real(amrex_real), intent(inout) :: temp

    integer :: i
    real(amrex_real) :: avm

    avm = 0.0d0

    do i = 1, nspecies
       avm = avm + fracvec(i)*molmass(i)

    enddo

    temp = avm*pressure/(runiv*density)

  end subroutine get_temperature_gas

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
          hcp(i) = 0.5d0*(2+dof(i))*Runiv/molmass(i)

          !print *, hcv(i)

       endif

    enddo

    !print *, hcv

  end subroutine get_hc_gas

end module conv_module

