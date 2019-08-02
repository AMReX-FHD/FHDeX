subroutine init_consvar(lo, hi, cu, culo, cuhi, pu, pulo, puhi, dx &
     , reallo, realhi) bind(C, name="init_consvar")

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, bc_lo, bc_hi, n_cells, prob_type, molmass, Runiv, grav, membrane_cell
  use conv_module, only : get_energy, get_pressure_gas

  implicit none

  integer         , intent(in   ) :: lo(3), hi(3), culo(3), cuhi(3), pulo(3), puhi(3)
  real(amrex_real), intent(inout) :: cu(culo(1):cuhi(1),culo(2):cuhi(2),culo(3):cuhi(3), nvars)
  real(amrex_real), intent(in   ) :: pu(pulo(1):puhi(1),pulo(2):puhi(2),pulo(3):puhi(3), nprimvars)
  real(amrex_real), intent(in   ) :: dx(3)
  real(amrex_real), intent(in   ) :: reallo(3), realhi(3)

  integer          :: i,j,k,l
  double precision :: pos(3),center(3),itVec(3),relpos(3)
  double precision :: L_hlf, pi
  double precision :: massvec(nspecies), intEnergy, pamb, molmix, rgasmix, alpha

  center = (realhi - reallo)/2d0
  L_hlf = (realhi(1) - reallo(1))/2d0

  pi = acos(-1.d0)

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           itVec(1) = (dble(i)+0.5)*dx(1)
           itVec(2) = (dble(j)+0.5)*dx(2)
           itVec(3) = (dble(k)+0.5)*dx(3)

           pos = reallo + itVec
           relpos = pos - center

           ! Total density must be pre-set
           
           if (prob_type.eq.2) then ! Rayleigh-Taylor
              if (relpos(3) .ge. 0) then
                 massvec = (/0.4, 0.4, 0.1, 0.1/)
              else
                 massvec = (/0.1, 0.1, 0.4, 0.4/)
              endif
              
              call get_pressure_gas(pamb, massvec, cu(i,j,k,1), pu(i,j,k,5))
              molmix = 0.0d0
              do l = 1, nspecies
                 molmix = molmix + massvec(l)/molmass(l)
              enddo
              molmix = 1.0d0/molmix
              rgasmix = Runiv/molmix
              alpha = grav(3)/(rgasmix*pu(i,j,k,5))
              
              ! rho = exponential in z-dir to init @ hydrostatic eqm.
              ! must satisfy system: dP/dz = -rho*g & P = rhogasmix*rho*T
              ! Assumes temp=const
              cu(i,j,k,1) = pamb*exp(alpha*pos(3))/(rgasmix*pu(i,j,k,5)) 
              
              do l = 1,nspecies
                 cu(i,j,k,5+l) = cu(i,j,k,1)*massvec(l)
              enddo

              call get_energy(intEnergy, massvec, pu(i,j,k,5))
              cu(i,j,k,5) = cu(i,j,k,1)*intEnergy + 0.5*cu(i,j,k,1)*(pu(i,j,k,2)**2 + &
                   pu(i,j,k,3)**2 + pu(i,j,k,4)**2)
           endif


        end do
     end do
  end do

end subroutine init_consvar
