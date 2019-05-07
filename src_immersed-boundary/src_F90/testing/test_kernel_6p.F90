#include "AMReX_fort_mod.F90"
#include "ib_fort_utils.F90"

program test_kernel_6p

    use amrex_fort_module, only: amrex_real
    use ib_fort_utils,     only: kernel_6p


    ! ** loop iteration
    integer                     :: i
    integer, parameter          :: ilo=-100, ihi=100


    ! ** distance from kernel center
    real(amrex_real)            :: x

    !** 8 because I want to go one cell out form the edge (6) in both directions
    real(amrex_real), parameter :: dx=8./(ihi-ilo)



    !*****************************************************************************!
    !                                                                             !
    !    TEST PROGRAM                                                             !
    !                                                                             !
    !*****************************************************************************!


    ! ** output information on test case
    write(*,*) "Testing 6-point kernel"
    write(*,*) "dx=", dx


    ! ** loop over test points and ouput (distance, function-value)
    do i = ilo, ihi
        x = i*dx
        write(*,*) x, kernel_6p(x)
    end do


end program test_kernel_6p
