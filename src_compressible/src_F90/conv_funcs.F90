module conv_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, diameter, max_species, molmass, k_b
  implicit none

  private

  public :: cons_to_prim

contains

  subroutine cons_to_prim(lo,hi, cons, prim) bind(C,name="cons_to_prim")

      integer         , intent(in   ) :: lo(3),hi(3)

      real(amrex_real), intent(inout) :: prim(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nprimvars)
      real(amrex_real), intent(in   ) :: cons(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nprimvars)

      integer :: i,j,k

      real(amrex_real) :: vsqr

      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)

            prim(i,j,k,1) = cons(i,j,k,1)
            prim(i,j,k,2) = cons(i,j,k,2)/cons(i,j,k,1)
            prim(i,j,k,3) = cons(i,j,k,3)/cons(i,j,k,1)
            prim(i,j,k,4) = cons(i,j,k,4)/cons(i,j,k,1)

            vsqr = prim(i,j,k,2)**2 + prim(i,j,k,3)**2 + prim(i,j,k,4)**2


          enddo
        enddo
      enddo


!          rho = con(i,j,k,1)
!          uvel = con(i,j,k,2) / rho
!          vvel = con(i,j,k,3) / rho
!          wvel = con(i,j,k,4) / rho
!          vmag = (uvel**2 + vvel**2 + wvel**2 )
!          eint = con(i,j,k,5)/rho -0.5d0*vmag


  end subroutine cons_to_prim

!  subroutine get_temperature(energy, Yk, IWRK, RWRK, temp,ierr)     

!    real(amrex_real), intent(inout) :: temp
!    real(amrex_real), intent(in   ) :: energy

!    real(kind=8) :: temp,Yk(1:nspecies),RWRK,eintmix
!    integer :: IWRK, i , ierr
!    real(amrex_real), :: cvmix, e0                                  

!    cvmix = 0.0d0; e0 = 0.0d0

!    do i = 1, nspecies
!      cvmix = cvmix + Yk(ns)*cvgas(ns)
!      e0 = e0 + Yk(ns)*e0ref(ns)
!    enddo

!    temp = (energy-e0)/cvmix 

!  end subroutine

end module conv_module



!    do ns =  1,nspecies
!       e0ref(ns) = 0.d0
!!      R_g(ns) = Runiv / molecular_weight(ns)
!       molecular_mass(ns) = molecular_weight(ns) / AVOGADRO
!       cvgas(ns) = 0.5d0*(3+int_deg_free(ns))*Runiv/molecular_weight(ns)
!       cpgas(ns) = 0.5d0*(5+int_deg_free(ns))*Runiv/molecular_weight(ns)

!    enddo
