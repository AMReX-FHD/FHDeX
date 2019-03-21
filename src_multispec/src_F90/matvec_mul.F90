module matvec_mul_module

  use amrex_error_module
  use common_namelist_module
  use multispec_namelist_module

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 2)

  subroutine matvec_mul(tlo, thi, &
                        xp, xlo, xhi, nc, &
                        ap, alo, ahi, nc2) bind(C,name="matvec_mul")

    integer,          intent(in   ) :: tlo(2),thi(2), nc, nc2
    integer,          intent(in   ) :: xlo(2),xhi(2), alo(2),ahi(2)
    double precision, intent(inout) :: xp(xlo(1):xhi(1),xlo(2):xhi(2),nc) ! last dimension for nc
    double precision, intent(in   ) :: ap(alo(1):ahi(1),alo(2):ahi(2),nc2)

    ! local
    integer :: i,j

    do j=tlo(2),thi(2)
       do i=tlo(1),thi(1)
          call matvec_mul_comp(xp(i,j,:), ap(i,j,:))
       end do
    end do

    contains 
    
    ! Use contained subroutine to do rank conversion and mat-vec mult
    subroutine matvec_mul_comp(xp_ij, ap_ij)

        double precision, dimension(nc),    intent(inout) :: xp_ij
        double precision, dimension(nc,nc), intent(in)    :: ap_ij
        
        xp_ij = matmul(ap_ij, xp_ij)
        
    end subroutine matvec_mul_comp 

  end subroutine matvec_mul

#endif

#if (AMREX_SPACEDIM == 3)
 
  subroutine matvec_mul(tlo, thi, &
                        xp, xlo, xhi, &
                        ap, alo, ahi, nc) bind(C,name="matvec_mul")

    integer,          intent(in   ) :: tlo(3),thi(3), nc
    integer,          intent(in   ) :: xlo(3),xhi(3), alo(3),ahi(3)
    double precision, intent(inout) :: xp(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3),nc) ! last dimension for nc
    double precision, intent(in   ) :: ap(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3),nc)

    ! local
    integer :: i,j,k
    
    do k=tlo(3),thi(3)
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)
             call matvec_mul_comp(xp(i,j,k,:), ap(i,j,k,:))
          end do
       end do
    end do
    contains 
    
    ! Use contained subroutine to do rank conversion and mat-vec mult
    subroutine matvec_mul_comp(xp_ij, ap_ij)

        double precision, dimension(nc),    intent(inout) :: xp_ij
        double precision, dimension(nc,nc), intent(in)    :: ap_ij  
        
        xp_ij = matmul(ap_ij, xp_ij)
 
    end subroutine matvec_mul_comp 

  end subroutine matvec_mul

#endif

end module matvec_mul_module
