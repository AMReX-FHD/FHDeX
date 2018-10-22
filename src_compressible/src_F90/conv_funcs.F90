module conv_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, diameter, max_species, molmass, k_b, nspecies, hcv, hcp, runiv, dof
  implicit none

  private

  public :: cons_to_prim, get_temperature, get_energy, get_hc_gas, get_pressure_gas

contains

  subroutine cons_to_prim(lo,hi, cons, prim) bind(C,name="cons_to_prim")

      integer         , intent(in   ) :: lo(3),hi(3)

      real(amrex_real), intent(inout) :: prim(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nprimvars)
      real(amrex_real), intent(in   ) :: cons(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nprimvars)

      integer :: i,j,k

      real(amrex_real) :: vsqr, massvec(nspecies), intenergy

      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)

            prim(i,j,k,1) = cons(i,j,k,1)
            prim(i,j,k,2) = cons(i,j,k,2)/cons(i,j,k,1)
            prim(i,j,k,3) = cons(i,j,k,3)/cons(i,j,k,1)
            prim(i,j,k,4) = cons(i,j,k,4)/cons(i,j,k,1)

            vsqr = prim(i,j,k,2)**2 + prim(i,j,k,3)**2 + prim(i,j,k,4)**2

            intenergy = cons(i,j,k,5) - 0.5*vsqr*cons(i,j,k,1)

            massvec = cons(i,j,k,6:nvars)*cons(i,j,k,1)

            call get_temperature(intenergy, massvec, prim(i,j,k,5))

!            if((i .eq. 36) .and. (j .eq. 0) .and. (k .eq. 0)) then

              !print *, "primcalc: ", i,j,k, " energy: ", cons(i,j,k,5), " temp: ", prim(i,j,k,5)

!            endif

          enddo
        enddo
      enddo

  end subroutine cons_to_prim

  subroutine get_temperature(energy, massvec, temp)     

    !This function originaly had a reference to e0 - check this.

    real(amrex_real), intent(inout) :: temp
    real(amrex_real), intent(in   ) :: energy, massvec(nspecies)

    integer :: i
    real(amrex_real) :: cvmix

    cvmix = 0.0d0

    do i = 1, nspecies
      cvmix = cvmix + massvec(i)*hcv(i)

    enddo

    temp = (energy)/cvmix 

  end subroutine

  subroutine get_energy(energy, massvec, temp)  bind(C,name="get_energy")    

    !This function originaly had a reference to e0 - check this.

    real(amrex_real), intent(in   ) :: temp
    real(amrex_real), intent(inout) :: energy, massvec(nspecies)

    integer :: i
    real(amrex_real) :: cvmix

    cvmix = 0.0d0

    do i = 1, nspecies
      cvmix = cvmix + massvec(i)*hcv(i)

    enddo

    energy = temp*cvmix 

  end subroutine

  subroutine get_pressure_gas(pressure, fracvec, density, temp)  bind(C,name="get_pressure_gas")    

    real(amrex_real), intent(in   ) :: temp, fracvec(nspecies), density
    real(amrex_real), intent(inout) :: pressure

    integer :: i
    real(amrex_real) :: avm

    avm = 0.0d0

    do i = 1, nspecies
      avm = avm + fracvec(i)*molmass(i)

    enddo

    pressure = temp*runiv*density/avm 

  end subroutine

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

  end subroutine

end module conv_module

