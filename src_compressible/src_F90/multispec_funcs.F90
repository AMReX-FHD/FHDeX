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

    real(amrex_real) :: p(np), sqda2(np,np), dij(np,np)
    real(amrex_real) :: dd(np,np)
    real(amrex_real) :: yy(np), mwmix 

    integer :: i, j, k, ii, jj
    real(amrex_real) :: sum1
    real(amrex_real) :: small_number
    real(amrex_real) :: sum

    integer :: idiag,ising

    small_number = 0.d0
    sum1 = 0.d0

    ! NOTE: For idiag=1, please refer to original LLNS code
    idiag = 0

    do i = 1, np
       ising = 0

       do j = i, np
          sum1 = a(i,j)

          do k = i-1, 1, -1
             sum1 = sum1 - a(i,k)*a(j,k)
          enddo

          if(i.eq.j) then
             if(sum1.le.small_number) then
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
