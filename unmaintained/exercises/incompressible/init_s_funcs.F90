subroutine init_s_vel(lo, hi, phic, phiclo, phichi, dx &
                      , reallo, realhi) bind(C, name="init_s_vel")

  use amrex_fort_module, only : amrex_real

  implicit none

  integer         , intent(in   ) :: lo(3), hi(3), phiclo(3), phichi(3)
  real(amrex_real), intent(inout) :: phic(phiclo(1):phichi(1),phiclo(2):phichi(2),phiclo(3):phichi(3))
  real(amrex_real), intent(in   ) :: dx(3)
  real(amrex_real), intent(in   ) :: reallo(3), realhi(3)

  integer          :: i,j,k
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
  r_a = 0.5d0*L_hlf
  ! r_b = L_hlf - r_a

  ! Stream:
  ! freq = 3.d0*pi/L_hlf
  ! amp = 2.0d-1*L_hlf
  width1 = L_hlf/2.0d0

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           itVec(1) = dble(i)*dx(1)
           itVec(2) = dble(j)*dx(2)
           itVec(3) = dble(k)*dx(3)

           pos = reallo + itVec
           relpos = pos - center
           rad2 = DOT_PRODUCT(relpos,relpos)
           rad = SQRT(rad2)

           ! Circle
           ! phic(i,j,k) = 0.5d0*(1d0+tanh(k2_inv*(r_a-rad)))

           ! Stream:
           ! perturb = amp*sin(freq*relpos(1))
           ! slope = amp*freq*cos(freq*relpos(1))
           perturb = 0d0
           fun_ptrb = 0.25d0*(1d0+tanh(k1_inv*(relpos(2) - (-width1/2.d0+perturb)))) &
                            *(1d0+tanh(k2_inv*((width1/2.d0+perturb) - relpos(2))))
           phic(i,j,k) = fun_ptrb

        end do
     end do
  end do

end subroutine init_s_vel
