subroutine init_consvar(lo, hi, cu, culo, cuhi, pu, pulo, puhi, dx &
     , reallo, realhi) bind(C, name="init_consvar")

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, bc_lo, bc_hi, n_cells, prob_type, membrane_cell
  use conv_module, only : get_energy

  implicit none

  integer         , intent(in   ) :: lo(3), hi(3), culo(3), cuhi(3), pulo(3), puhi(3)
  real(amrex_real), intent(inout) :: cu(culo(1):cuhi(1),culo(2):cuhi(2),culo(3):cuhi(3), nvars)
  real(amrex_real), intent(in   ) :: pu(pulo(1):puhi(1),pulo(2):puhi(2),pulo(3):puhi(3), nprimvars)
  real(amrex_real), intent(in   ) :: dx(3)
  real(amrex_real), intent(in   ) :: reallo(3), realhi(3)

  integer          :: i,j,k,l
  double precision :: pos(3),center(3),itVec(3),relpos(3)
  double precision :: L_hlf, pi
  double precision :: massvec(nspecies), intEnergy

  center = (realhi - reallo)/2d0
  L_hlf = (realhi(1) - reallo(1))/2d0

  pi = acos(-1.d0)

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           itVec(1) = dble(i)*dx(1)
           itVec(2) = dble(j)*dx(2)
           itVec(3) = dble(k)*dx(3)

           pos = reallo + itVec
           relpos = pos - center

           ! Total density must be pre-set
           
           if (prob_type.eq.2) then ! Rayleigh-Taylor
              if (relpos(3) .ge. 0) then
                 massvec = (/0.4, 0.4, 0.1, 0.1/)
                 do l = 1,nspecies
                    cu(i,j,k,5+l) = cu(i,j,k,1)*massvec(l)
                 enddo
              else
                 massvec = (/0.1, 0.1, 0.4, 0.4/)
                 do l = 1,nspecies
                    cu(i,j,k,5+l) = cu(i,j,k,1)*massvec(l)
                 enddo
              endif

              call get_energy(intEnergy, massvec, pu(i,j,k,5))
              cu(i,j,k,5) = cu(i,j,k,1)*intEnergy + 0.5*cu(i,j,k,1)*(pu(i,j,k,2)**2 + &
                   pu(i,j,k,3)**2 + pu(i,j,k,4)**2)
           endif


        end do
     end do
  end do

end subroutine init_consvar
