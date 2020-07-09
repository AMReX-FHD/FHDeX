module multispec_module

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : n_cells, ngc, k_b, Runiv, nprimvars, nspecies, molmass, diameter, hcp, hcv
  use conv_module, only : get_molfrac

  implicit none

  private

  public :: cholesky_decomp

contains

  subroutine cholesky_decomp(a,np,sqda) bind(C,name="cholesky_decomp")  

    integer :: np

    real(amrex_real), intent(inout) :: a(np,np)
    real(amrex_real), intent(inout) :: sqda(np,np)

    real(amrex_real) :: p(np)

    integer :: i, j, k
    real(amrex_real) :: sum1

    integer :: ising

    sum1 = 0.d0

    do i = 1, np
       ising = 0

       do j = i, np
          sum1 = a(i,j)

          do k = i-1, 1, -1
             sum1 = sum1 - a(i,k)*a(j,k)
          enddo

          if(i.eq.j) then
             if(sum1.le.0.d0) then
                p(i) = 0.d0
                ising = 1
             else
                p(i) = sqrt(sum1)
             endif
          else
             if(ising.eq.0)then
                a(j,i) = sum1/p(i)
             else
                a(j,i) = 0.d0
             endif
          endif

       enddo
    enddo

    sqda = 0.0d0

    do i = 1, np
       do j = i-1, 1, -1
          sqdA(i,j) = a(i,j)
       enddo
       sqdA(i,i) = p(i)
    enddo

  end subroutine cholesky_decomp

end module multispec_module
