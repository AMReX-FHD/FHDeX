subroutine init_consvar(lo, hi, cu, culo, cuhi, pu, pulo, puhi, dx, &
                        reallo, realhi) bind(C, name="init_consvar")

  use amrex_fort_module, only : amrex_real
  use common_namelist_module, only : ngc, nvars, nprimvars, nspecies, &
                                     n_cells, prob_type, molmass, Runiv, grav, membrane_cell, rho0, prob_lo, prob_hi, t_lo, t_hi
  use compressible_namelist_module, only : bc_Yk
  use conv_module, only : get_energy, get_pressure_gas

  implicit none

  integer         , intent(in   ) :: lo(3), hi(3), culo(3), cuhi(3), pulo(3), puhi(3)
  real(amrex_real), intent(inout) :: cu(culo(1):cuhi(1),culo(2):cuhi(2),culo(3):cuhi(3),nvars)
  real(amrex_real), intent(in   ) :: pu(pulo(1):puhi(1),pulo(2):puhi(2),pulo(3):puhi(3),nprimvars)
  real(amrex_real), intent(in   ) :: dx(3)
  real(amrex_real), intent(in   ) :: reallo(3), realhi(3)

  integer          :: i,j,k,l
  double precision :: pos(3),center(3),itVec(3),relpos(3)
  double precision :: L_hlf, pi, Lf, pres, velscale, mach,x,y,z, hy
  double precision :: massvec(nspecies), intEnergy, pamb, molmix, rgasmix, alpha
  double precision :: Ygrad

  center = (realhi - reallo)/2.d0
  L_hlf = (realhi(1) - reallo(1))/2.d0
  Lf = realhi(1) - reallo(1)
  mach = 0.3
  velscale = 30565.2d0*mach

  hy = (prob_hi(2) - prob_lo(2))/3d0

  pi = acos(-1.d0)
  !write(6,*)"into init"

  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     
     itVec(1) = (dble(i)+0.5)*dx(1)
     itVec(2) = (dble(j)+0.5)*dx(2)
     itVec(3) = (dble(k)+0.5)*dx(3)

     pos = reallo + itVec
     relpos = pos - center

     ! Total density must be pre-set

     if (prob_type .eq. 2) then ! Rayleigh-Taylor

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

     elseif (prob_type .eq. 3) then ! diffusion barrier

        do l = 1,nspecies
           Ygrad = (bc_Yk(2,2,l) - bc_Yk(2,1,l))/(realhi(2) - reallo(2))
           massvec(l) = Ygrad*pos(2) + bc_Yk(2,1,l)
           cu(i,j,k,5+l) = cu(i,j,k,1)*massvec(l)
        enddo

        call get_energy(intEnergy, massvec, pu(i,j,k,5))
        cu(i,j,k,5) = cu(i,j,k,1)*intEnergy + 0.5*cu(i,j,k,1)*(pu(i,j,k,2)**2 + &
             pu(i,j,k,3)**2 + pu(i,j,k,4)**2)

     elseif (prob_type .eq. 4) then ! Taylor Green Vortex

        x=itVec(1)
        y=itVec(2)
        z=itVec(3)

        cu(i,j,k,1) = 1.784d-3
        cu(i,j,k,2) =  velscale*cu(i,j,k,1)*sin(2.d0*pi*x/Lf)*cos(2.d0*pi*y/Lf)*cos(2.d0*pi*z/Lf)
        cu(i,j,k,3) = -velscale*cu(i,j,k,1)*cos(2.d0*pi*x/Lf)*sin(2.d0*pi*y/Lf)*cos(2.d0*pi*z/Lf)
        cu(i,j,k,4) = 0.d0
        pres = 1.01325d6+cu(i,j,k,1)*velscale**2*cos(2.d0*pi*x/Lf)*cos(4.d0*pi*y/Lf)*(cos(4.d0*pi*z/Lf)+2.d0)
        cu(i,j,k,5) = pres/(5.d0/3.d0-1.d0) + 0.5*(cu(i,j,k,2)**2 + &
             cu(i,j,k,3)**2 + cu(i,j,k,4)**2) / cu(i,j,k,1)
        cu(i,j,k,6) = cu(i,j,k,1)
        cu(i,j,k,7) = 0.d0
        

     elseif (prob_type .eq. 5) then ! Taylor Green Vortex

       

        cu(i,j,k,1) = rho0
        cu(i,j,k,2) = 0
        cu(i,j,k,3) = 0
        cu(i,j,k,4) = 0
        if((prob_lo(2) + itVec(2)) < hy) then
          massvec = (/1.0,0.0/) 
          call get_energy(intEnergy, massvec, t_lo(2));
          cu(i,j,k,5) = cu(i,j,k,1)*intEnergy
          cu(i,j,k,6) = cu(i,j,k,1)
          cu(i,j,k,7) = 0
        elseif((prob_lo(2) + itVec(2)) < 2*hy) then
          massvec = (/0.0,1.0/)
          call get_energy(intEnergy, massvec, t_hi(2));
          cu(i,j,k,5) = cu(i,j,k,1)*intEnergy
          cu(i,j,k,6) = 0
          cu(i,j,k,7) = cu(i,j,k,1)
        else
          massvec = (/1.0,0.0/) 
          call get_energy(intEnergy, massvec, t_lo(2));
          cu(i,j,k,5) = cu(i,j,k,1)*intEnergy
          cu(i,j,k,6) = cu(i,j,k,1)
          cu(i,j,k,7) = 0
        endif
      endif

  end do
  end do
  end do

end subroutine init_consvar
