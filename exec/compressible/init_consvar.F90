subroutine init_consvar(lo, hi, cu, culo, cuhi, ncomp, dx &
     , reallo, realhi) bind(C, name="init_consvar")

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, bc_lo, bc_hi, n_cells, membrane_cell

  implicit none

  integer         , intent(in   ) :: lo(3), hi(3), culo(3), cuhi(3), ncomp
  real(amrex_real), intent(inout) :: cu(culo(1):cuhi(1),culo(2):cuhi(2),culo(3):cuhi(3), nvars)
  real(amrex_real), intent(in   ) :: dx(3) 
  real(amrex_real), intent(in   ) :: reallo(3), realhi(3)

  integer          :: i,j,k,l
  double precision :: pos(3),center(3),partdom,itVec(3),relpos(3),rad,rad2

  double precision :: L_hlf, k1, k1_inv, k2, k2_inv, r_a, r_b
  double precision :: pi, freq, amp, width1, width2, perturb, slope, fun, fun_ptrb

  center = (realhi - reallo)/2d0
  L_hlf = (realhi(1) - reallo(1))/2d0

  !! IC parameters
  pi = acos(-1.d0)

  ! k1 & k2 determine steepness of profile:
  k1 = 1d-2*L_hlf
  k2 = k1
  k1_inv = 1/k1
  k2_inv = 1/k2

  ! Vortex:
  ! [r_a r_b] defines radial bounds of bump:
  ! r_a = 0.5d0*L_hlf
  ! r_b = L_hlf - r_a

  ! Stream:
  width1 = L_hlf/2.0d0

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           itVec(1) = dble(i)*dx(1)
           itVec(2) = dble(j)*dx(2)
           itVec(3) = dble(k)*dx(3)

           pos = reallo + itVec
           relpos = pos - center

           ! rad2 = DOT_PRODUCT(relpos,relpos)
           ! rad = SQRT(rad2)

           ! Circle
           ! fun = 0.5d0*(1d0+tanh(k2_inv*(r_a-rad)))

           ! Stream:
           ! perturb = 0d0
           ! fun = 0.25d0*(1d0+tanh(k1_inv*(rad - (-width1/1.5d0+perturb)))) &
           !      *(1d0+tanh(k2_inv*((width1/1.5d0+perturb) - rad)))

           ! Total density must be pre-set
           if (relpos(3) .gt. 0) then
              cu(i,j,k,5+1) = cu(i,j,k,1)*0.4
              cu(i,j,k,5+2) = cu(i,j,k,1)*0.4
              cu(i,j,k,5+3) = cu(i,j,k,1)*0.1
              cu(i,j,k,5+4) = cu(i,j,k,1)*0.1
           else
              cu(i,j,k,5+1) = cu(i,j,k,1)*0.1
              cu(i,j,k,5+2) = cu(i,j,k,1)*0.1
              cu(i,j,k,5+3) = cu(i,j,k,1)*0.4
              cu(i,j,k,5+4) = cu(i,j,k,1)*0.4
           endif

        end do
     end do
  end do

end subroutine init_consvar
