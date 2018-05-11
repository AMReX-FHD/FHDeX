
!Computes divergence at cell centres from velcocities at cell faces
#if AMREX_SPACEDIM == 2
subroutine compute_div(lo, hi, velx, velxlo, velxhi, vely, velylo, velyhi, div, divlo, divhi, dx, increment) bind(C, name="compute_div")

  use amrex_fort_module, only : amrex_real
  implicit none

	integer, intent(in) :: lo(2), hi(2), velxlo(2), velxhi(2), velylo(2), velyhi(2), divlo(2), divhi(2), increment
  real(amrex_real), intent(inout) :: velx(velxlo(1):velxhi(1),velxlo(2):velxhi(2))
  real(amrex_real), intent(inout) :: vely(velylo(1):velyhi(1),velylo(2):velyhi(2))
  real(amrex_real), intent(inout) :: div(divlo(1):divhi(1),divlo(2):divhi(2))
	real(amrex_real), intent(in   ) :: dx(2)

  ! local variables
  integer i,j
  double precision dxInv, dyInv

  dxInv = 1.d0/dx(1)
  dyInv = 1.d0/dx(2)

  if (increment .eq. 0) then
  	!compute divergence
    do j = lo(2), hi(2)
    	do i = lo(1), hi(1)
       div(i,j) = (-velx(i,j) + velx(i+1,j))*dxInv + (-vely(i,j) + vely(i,j+1))*dyInv
    	end do
    end do
  else
  	!increment divergence
    do j = lo(2), hi(2)
    	do i = lo(1), hi(1)
       div(i,j) = div(i,j) + (-velx(i,j) + velx(i+1,j))*dxInv + (-vely(i,j) + vely(i,j+1))*dyInv
    	end do
    end do
  endif

end subroutine compute_div
#endif

#if AMREX_SPACEDIM == 3
subroutine compute_div(lo, hi, velx, velxlo, velxhi, vely, velylo, velyhi, velz, velzlo, velzhi, div, divlo, divhi, dx, increment) bind(C, name="compute_div")

  use amrex_fort_module, only : amrex_real
  implicit none

	integer, intent(in) :: lo(3), hi(3), velxlo(3), velxhi(3), velylo(3), velyhi(3), velzlo(3), velzhi(3), divlo(3), divhi(3), increment
  real(amrex_real), intent(inout) :: velx(velxlo(1):velxhi(1),velxlo(2):velxhi(2),velxlo(3):velxhi(3))
  real(amrex_real), intent(inout) :: vely(velylo(1):velyhi(1),velylo(2):velyhi(2),velylo(3):velyhi(3))
  real(amrex_real), intent(inout) :: velz(velzlo(1):velzhi(1),velzlo(2):velzhi(2),velzlo(3):velzhi(3))
  real(amrex_real), intent(inout) :: div(divlo(1):divhi(1),divlo(2):divhi(2),divlo(3):divhi(3))
	real(amrex_real), intent(in   ) :: dx(3) 

  ! local variables
  integer i,j,k
  double precision dxInv, dyInv, dzInv

  dxInv = 1.d0/dx(1)
  dyInv = 1.d0/dx(2)
  dzInv = 1.d0/dx(3)

  if (increment .eq. 0) then
	  !compute divergence
	  do k = lo(3), hi(3)
		  do j = lo(2), hi(2)
			  do i = lo(1), hi(1)
		     div(i,j,k) = (-velx(i,j,k) + velx(i+1,j,k))*dxInv + (-vely(i,j,k) + vely(i,j+1,k))*dyInv + (-velz(i,j,k) + vely(i,j,k+1))*dzInv
			  end do
		  end do
	  end do
  else
	  !compute divergence
	  do k = lo(3), hi(3)
		  do j = lo(2), hi(2)
			  do i = lo(1), hi(1)
		     div(i,j,k) = (-velx(i,j,k) + velx(i+1,j,k))*dxInv + (-vely(i,j,k) + vely(i,j+1,k))*dyInv + (-velz(i,j,k) + vely(i,j,k+1))*dzInv
			  end do
		  end do
	  end do
  endif

end subroutine compute_div
#endif

!Computes gradient of cell centred scalar field at cell faces
#if AMREX_SPACEDIM == 2
subroutine compute_grad(lo, hi, velx, velxlo, velxhi, vely, velylo, velyhi, phi, philo, phihi, dx) bind(C, name="compute_grad")

  use amrex_fort_module, only : amrex_real
  implicit none

	integer, intent(in) :: lo(2), hi(2), velxlo(2), velxhi(2), velylo(2), velyhi(2), philo(2), phihi(2)
  real(amrex_real), intent(inout) :: velx(velxlo(1):velxhi(1),velxlo(2):velxhi(2))
  real(amrex_real), intent(inout) :: vely(velylo(1):velyhi(1),velylo(2):velyhi(2))
  real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2))
	real(amrex_real), intent(in   ) :: dx(2) 

  ! local variables
  integer i,j
  double precision dxInv, dyInv;

  dxInv = 1.d0/dx(1)
  dyInv = 1.d0/dx(2)

  !Currently assuming periodic

	!x faces
  do j = lo(2), hi(2)
  	do i = lo(1), hi(1)+1
     velx(i,j) = (phi(i,j) - phi(i-1,j))*dxInv
  	end do
  end do

	!y faces
  do j = lo(2), hi(2)+1
  	do i = lo(1), hi(1)
     vely(i,j) = (phi(i,j) - phi(i,j-1))*dyInv
  	end do
  end do

end subroutine compute_grad
#endif

#if AMREX_SPACEDIM == 3
subroutine compute_grad(lo, hi, velx, velxlo, velxhi, vely, velylo, velyhi, velz, velzlo, velzhi, phi, philo, phihi, dx) bind(C, name="compute_grad")

  use amrex_fort_module, only : amrex_real
  implicit none

	integer, intent(in) :: lo(3), hi(3), velxlo(3), velxhi(3), velylo(3), velyhi(3), velzlo(3), velzhi(3), philo(3), phihi(3)
  real(amrex_real), intent(inout) :: velx(velxlo(1):velxhi(1),velxlo(2):velxhi(2),velxlo(3):velxhi(3))
  real(amrex_real), intent(inout) :: vely(velylo(1):velyhi(1),velylo(2):velyhi(2),velylo(3):velyhi(3))
  real(amrex_real), intent(inout) :: velz(velzlo(1):velzhi(1),velzlo(2):velzhi(2),velzlo(3):velzhi(3))
  real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
	real(amrex_real), intent(in   ) :: dx(3) 

  ! local variables
  integer i,j,k
  double precision dxInv, dyInv, dzInv;

  dxInv = 1.d0/dx(1)
  dyInv = 1.d0/dx(2)
  dzInv = 1.d0/dx(3)

  !Currently assuming periodic

	!x faces
  do k = lo(3), hi(3)
    do j = lo(2), hi(2)
    	do i = lo(1), hi(1)+1
       velx(i,j,k) = (phi(i,j,k) - phi(i-1,j,k))*dxInv
    	end do
    end do
  end do

	!y faces
  do k = lo(3), hi(3)
    do j = lo(2), hi(2)+1
    	do i = lo(1), hi(1)
       vely(i,j,k) = (phi(i,j,k) - phi(i,j-1,k))*dyInv
    	end do
    end do
  end do

	!z faces
  do k = lo(3), hi(3)+1
    do j = lo(2), hi(2)
    	do i = lo(1), hi(1)
       vely(i,j,k) = (phi(i,j,k) - phi(i,j,k-1))*dzInv
    	end do
    end do
  end do


end subroutine compute_grad
#endif



