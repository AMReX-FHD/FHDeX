module debug
        implicit none

contains

subroutine print_potential(lo, hi, pot, potlo, pothi, xstart, ystart, zstart) bind(C, name="print_potential")
        integer, intent(in) :: lo(3), hi(3), potlo(3), pothi(3), xstart, ystart, zstart 
        double precision, intent(in) :: pot(potlo(1):pothi(1),potlo(2):pothi(2), potlo(3):pothi(3))
        integer :: i,j,k
        print*, "POTENTIAL MFAB: "
        do k=zstart,zstart+5
        do j=ystart,ystart+5
        do i=xstart,xstart+5
                print *, "(", i, j, k, ")", pot(i,ystart,k)
        enddo
        enddo
        enddo
end subroutine

end module debug
