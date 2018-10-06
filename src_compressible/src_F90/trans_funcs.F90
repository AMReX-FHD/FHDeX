module trans_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, diameter, max_species, molmass, k_b
  implicit none

  private

  public :: trans_coeffs

contains

  subroutine trans_coeffs(lo,hi, prim, eta, zeta, kappa) bind(C,name="trans_coeffs")

      integer         , intent(in   ) :: lo(3),hi(3)

      real(amrex_real), intent(in   ) :: prim(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc, nprimvars)

      real(amrex_real), intent(inout) :: eta(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc)
      real(amrex_real), intent(inout) :: zeta(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc)
      real(amrex_real), intent(inout) :: kappa(lo(1)-ngc:hi(1)+ngc,lo(2)-ngc:hi(2)+ngc,lo(3)-ngc:hi(3)+ngc)

      integer :: i,j,k 
      real(amrex_real) :: mgrams(MAX_SPECIES)
      real(amrex_real) :: vconst, tconst, R, gamma1, gamma2, rootT

      mgrams = molmass/(6.022140857e23)

      R = k_B/mgrams(1)

      gamma1 = 1.2700
      gamma2 = 1.9223

      !single species, hard sphere for now. Check these defs
      vconst = (gamma1/4)*sqrt(R/3.14159265359)*(mgrams(1)/(diameter(1)**2))  ! Per Sone (Sone Grad-Hilbert, Bird Chapman Enskog? Check this)
      tconst = (5*gamma2/8)*R*sqrt(R/3.14159265359)*(mgrams(1)/(diameter(1)**2))

      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1)

            rootT = sqrt(prim(i,j,k,5))

            zeta(i,j,k) = 0 !no bulk viscosity for now

            eta(i,j,k) = rootT*vconst
            kappa(i,j,k) = rootT*tconst

          enddo
        enddo
      enddo


  end subroutine trans_coeffs

end module trans_module
