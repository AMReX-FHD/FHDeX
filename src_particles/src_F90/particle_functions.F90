module cell_sorted_particle_module
  use amrex_fort_module, only: amrex_real, amrex_particle_real
  use iso_c_binding ,    only: c_int
  
  implicit none
  private
  
  public particle_t, remove_particle_from_cell
  
  type, bind(C) :: particle_t
#if (BL_SPACEDIM == 3)
     real(amrex_particle_real) :: pos(3)
#endif
#if (BL_SPACEDIM == 2)
     real(amrex_particle_real) :: pos(2)
#endif


     real(amrex_particle_real) :: vel(3)     
     real(amrex_particle_real) :: mass
     real(amrex_particle_real) :: R
     real(amrex_particle_real) :: radius
     real(amrex_particle_real) :: accel_factor
     real(amrex_particle_real) :: drag_factor
     real(amrex_particle_real) :: angular_vel1
     real(amrex_particle_real) :: angular_vel2
     real(amrex_particle_real) :: angular_vel3

     integer(c_int)            :: id         
     integer(c_int)            :: cpu        
     integer(c_int)            :: sorted     
     integer(c_int)            :: species
  end type particle_t

contains
  
  subroutine remove_particle_from_cell(cell_parts, cell_np, new_np, i)
    
    use iso_c_binding, only: c_int
    
    implicit none
    
    integer(c_int), intent(inout) :: cell_parts(cell_np)
    integer(c_int), intent(in   ) :: cell_np
    integer(c_int), intent(inout) :: new_np
    integer(c_int), intent(in   ) :: i 

    cell_parts(i) = cell_parts(new_np)
    new_np = new_np - 1
        
  end subroutine remove_particle_from_cell
  
end module cell_sorted_particle_module

subroutine move_particles_dsmc(particles, np, lo, hi, &
     cell_part_ids, cell_part_cnt, clo, chi, plo, dx, dt, surfaces, ns) &
     bind(c,name="move_particles_dsmc")
  
  use amrex_fort_module, only: amrex_real
  use iso_c_binding, only: c_ptr, c_int, c_f_pointer
  use cell_sorted_particle_module, only: particle_t, remove_particle_from_cell
  use surfaces_module
  
  implicit none

  type(particle_t), intent(inout), target :: particles(np)
  type(surface_t), intent(in), target :: surfaces(ns)
  integer(c_int), intent(in) :: np, ns
  integer(c_int), intent(in) :: lo(3), hi(3)
  integer(c_int), intent(in) :: clo(3), chi(3)
  type(c_ptr), intent(inout) :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
  integer(c_int), intent(inout) :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
  real(amrex_real), intent(in) :: plo(3)
  real(amrex_real), intent(in) :: dx(3)
  real(amrex_real), intent(in) :: dt
  
  integer :: i, j, k, p, cell_np, new_np, intsurf, intside
  integer :: cell(3)
  integer(c_int), pointer :: cell_parts(:)
  type(particle_t), pointer :: part
  type(surface_t), pointer :: surf
  real(amrex_real) inv_dx(3), proj(3), runtime, inttime, adj

  adj = 0.999999

  inv_dx = 1.d0/dx

  

  !part => particles(1)

  !print *, surf%ux, surf%uy, surf%uz
  !call find_intersect_3d(part%vel(1),part%vel(2),part%vel(3),part%pos(1),part%pos(2),part%pos(3), dt,surfaces,ns)
  
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           cell_np = cell_part_cnt(i,j,k)
           call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])

           new_np = cell_np
           p = 1
           do while (p <= new_np)
              
              part => particles(cell_parts(p))

              runtime = dt

              do while (runtime .gt. 0)

                !print *, "Particle ", p, " prevel ", part%vel, " prepos ", part%pos
#if (BL_SPACEDIM == 3)                
                call find_intersect(part%vel(1),part%vel(2),part%vel(3),part%pos(1),part%pos(2),part%pos(3),runtime, surfaces, ns, intsurf, inttime, intside)
#endif
                !inttime = runtime
                !intsurf = -1

                ! move the particle in a straight line, adj factor prevents double detection of boundary intersection
                part%pos(1) = part%pos(1) + inttime*part%vel(1)*adj
                part%pos(2) = part%pos(2) + inttime*part%vel(2)*adj
#if (BL_SPACEDIM == 3)
                part%pos(3) = part%pos(3) + inttime*part%vel(3)*adj
#endif
                runtime = runtime - inttime

                if(intsurf .gt. 0) then

                  surf => surfaces(intsurf)

                  if(intside .eq. 0) then !lhs
#if (BL_SPACEDIM == 3)
                    !print *, "Applying bc. prevel: ", part%vel
                    call apply_bc(part%vel(1),part%vel(2),part%vel(3), surf%lnx, surf%lny,surf%lnz, surf%costhetaleft, surf%sinthetaleft, surf%cosphileft, surf%sinphileft, surf%porosityleft, surf%specularityleft, surf%temperatureleft, part%R)
                    !print *, "Bc applied. postvel: ", part%vel
#endif
#if (BL_SPACEDIM == 2)
                    call apply_bc(part%vel(1),part%vel(2), surf%lnx, surf%lny, surf%costhetaleft, surf%sinthetaleft, surf%porosityleft, surf%specularityleft, surf%temperatureleft, part%R)
#endif

                  else !rhs

#if (BL_SPACEDIM == 3)
                   call apply_bc(part%vel(1),part%vel(2),part%vel(3), surf%lnx, surf%lny,surf%lnz, surf%costhetaright, surf%sinthetaright, surf%cosphiright, surf%sinphiright, surf%porosityright, surf%specularityright, surf%temperatureright, part%R)
#endif
#if (BL_SPACEDIM == 2)
                   call apply_bc(part%vel(1),part%vel(2), surf%lnx, surf%lny, surf%costhetaright, surf%sinthetaright, surf%porosityright, surf%specularityright, surf%temperatureright, part%R)
#endif
                  endif

                endif
                !print *, "Particle ", p, " postvel ", part%vel, " postpos ", part%pos

              end do

              ! if it has changed cells, remove from vector.
              ! otherwise continue

              cell(1) = floor((part%pos(1) - plo(1))*inv_dx(1))              
              cell(2) = floor((part%pos(2) - plo(2))*inv_dx(2))
#if (BL_SPACEDIM == 3)
              cell(3) = floor((part%pos(3) - plo(3))*inv_dx(3))
#else
              cell(3) = 0
#endif
              if ((cell(1) /= i) .or. (cell(2) /= j) .or. (cell(3) /= k)) then
                 part%sorted = 0
                 call remove_particle_from_cell(cell_parts, cell_np, new_np, p)  
              else
                 p = p + 1
              end if
           end do

           cell_part_cnt(i,j,k) = new_np
           
        end do
     end do
  end do
  
end subroutine move_particles_dsmc

subroutine move_particles_fhd(particles, np, lo, hi, &
     cell_part_ids, cell_part_cnt, clo, chi, plo, dx, dt, plof, dxf, &
                                     velx, velxlo, velxhi, &
                                     vely, velylo, velyhi, &
#if (BL_SPACEDIM == 3)
                                     velz, velzlo, velzhi, &
#endif
                                     coordsx, coordsxlo, coordsxhi, &
                                     coordsy, coordsylo, coordsyhi, &
#if (BL_SPACEDIM == 3)
                                     coordsz, coordszlo, coordszhi, &
#endif
                                     beta, betalo, betahi, &

                                     rho, rholo, rhohi, &

                                     sourcex, sourcexlo, sourcexhi, &
                                     sourcey, sourceylo, sourceyhi, &
#if (BL_SPACEDIM == 3)
                                     sourcez, sourcezlo, sourcezhi, &
#endif
                                     surfaces, ns)bind(c,name="move_particles_fhd")
  use amrex_fort_module, only: amrex_real
  use iso_c_binding, only: c_ptr, c_int, c_f_pointer
  use cell_sorted_particle_module, only: particle_t, remove_particle_from_cell
  use common_namelist_module, only: visc_type, k_B
  use rng_functions_module
  use surfaces_module
  
  implicit none

  integer,          intent(in   )         :: np, ns, lo(3), hi(3), clo(3), chi(3), velxlo(3), velxhi(3), velylo(3), velyhi(3)
  integer,          intent(in   )         :: sourcexlo(3), sourcexhi(3), sourceylo(3), sourceyhi(3), rholo(3), rhohi(3)
  integer,          intent(in   )         :: coordsxlo(3), coordsxhi(3), coordsylo(3), coordsyhi(3)
  integer,          intent(in   )         :: betalo(3), betahi(3)
#if (AMREX_SPACEDIM == 3)
  integer,          intent(in   )         :: velzlo(3), velzhi(3), sourcezlo(3), sourcezhi(3), coordszlo(3), coordszhi(3)
#endif
  type(particle_t), intent(inout), target :: particles(np)
  type(surface_t),  intent(in),    target :: surfaces(ns)

  double precision, intent(in   )         :: dx(3), dxf(3), dt, plo(3), plof(3)

  double precision, intent(in   ) :: velx(velxlo(1):velxhi(1),velxlo(2):velxhi(2),velxlo(3):velxhi(3))
  double precision, intent(in   ) :: vely(velylo(1):velyhi(1),velylo(2):velyhi(2),velylo(3):velyhi(3))
#if (AMREX_SPACEDIM == 3)
  double precision, intent(in   ) :: velz(velzlo(1):velzhi(1),velzlo(2):velzhi(2),velzlo(3):velzhi(3))
#endif

  double precision, intent(in   ) :: coordsx(coordsxlo(1):coordsxhi(1),coordsxlo(2):coordsxhi(2),coordsxlo(3):coordsxhi(3),1:AMREX_SPACEDIM)
  double precision, intent(in   ) :: coordsy(coordsylo(1):coordsyhi(1),coordsylo(2):coordsyhi(2),coordsylo(3):coordsyhi(3),1:AMREX_SPACEDIM)
#if (AMREX_SPACEDIM == 3)
  double precision, intent(in   ) :: coordsz(coordszlo(1):coordszhi(1),coordszlo(2):coordszhi(2),coordszlo(3):coordszhi(3),1:AMREX_SPACEDIM)
#endif

  double precision, intent(in   ) :: beta(betalo(1):betahi(1),betalo(2):betahi(2),betalo(3):betahi(3))

  double precision, intent(in   ) :: rho(rholo(1):rhohi(1),rholo(2):rhohi(2),rholo(3):rhohi(3))

  double precision, intent(inout) :: sourcex(sourcexlo(1):sourcexhi(1),sourcexlo(2):sourcexhi(2),sourcexlo(3):sourcexhi(3))
  double precision, intent(inout) :: sourcey(sourceylo(1):sourceyhi(1),sourceylo(2):sourceyhi(2),sourceylo(3):sourceyhi(3))
#if (AMREX_SPACEDIM == 3)
  double precision, intent(inout) :: sourcez(sourcezlo(1):sourcezhi(1),sourcezlo(2):sourcezhi(2),sourcezlo(3):sourcezhi(3))
#endif

  type(c_ptr),      intent(inout) :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
  integer(c_int),   intent(inout) :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
  
  integer :: i, j, k, p, cell_np, new_np, intside, intsurf
  integer :: ni, nj, nk, fi, fj, fk
  integer(c_int), pointer :: cell_parts(:)
  type(particle_t), pointer :: part
  real(amrex_real) dxinv(3), dxfinv(3), onemdxf(3), ixf(3), localvel(3), localbeta, bfac(3), deltap(3), std, normalrand(3), nodalp, tempvel(3), intold, inttime, runerr, runtime

#if (BL_SPACEDIM == 3)
  double precision c000,c001,c010,c011,c100,c101,c110,c111, ctotal
  double precision r000,r001,r010,r011,r100,r101,r110,r111
#endif

#if (BL_SPACEDIM == 2)
  double precision c00,c01,c10,c11, ctotal
  double precision r00,r01,r10,r11
#endif

  !Zero source terms
  do k = sourcexlo(3), sourcexhi(3)
    do j = sourcexlo(2), sourcexhi(2)
      do i = sourcexlo(1), sourcexhi(1)

        sourcex(i,j,k) = 0;

      enddo
    enddo
  enddo

  do k = sourceylo(3), sourceyhi(3)
    do j = sourceylo(2), sourceyhi(2)
      do i = sourceylo(1), sourceyhi(1)

        sourcey(i,j,k) = 0;

      enddo
    enddo
  enddo

#if (BL_SPACEDIM == 3)
  do k = sourcezlo(3), sourcezhi(3)
    do j = sourcezlo(2), sourcezhi(2)
      do i = sourcezlo(1), sourcezhi(1)

        sourcez(i,j,k) = 0;

      enddo
    enddo
  enddo
#endif

  dxinv = 1.d0/dx

  dxfinv = 1.d0/dxf
  onemdxf = 1.d0 - dxf
  


  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           cell_np = cell_part_cnt(i,j,k)
           call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])

           new_np = cell_np
           p = 1

           do while (p <= new_np)

              runtime = dt
              part => particles(cell_parts(p))

              do while (runtime .gt. 0)

                !find fluid cell
                fi = floor((part%pos(1) - plof(1))*dxfinv(1))
                fj = floor((part%pos(2) - plof(2))*dxfinv(2))
#if (BL_SPACEDIM == 3)
                fk = floor((part%pos(3) - plof(3))*dxfinv(3))
#else
                fk = 0
#endif

                !Interpolate fluid fields
                ixf(1) = (part%pos(1) - coordsx(fi,fj,fk,1))*dxfInv(1)
                ixf(2) = (part%pos(2) - coordsy(fi,fj,fk,2))*dxfInv(2)
#if (BL_SPACEDIM == 3)
                ixf(3) = (part%pos(3) - coordsz(fi,fj,fk,3))*dxfInv(3)
#endif

#if (BL_SPACEDIM == 3)
                c000 = onemdxf(1)*onemdxf(2)*onemdxf(3)
                c001 = onemdxf(1)*onemdxf(2)*ixf(3)
                c010 = onemdxf(1)*onemdxf(3)*ixf(2)
                c011 = onemdxf(1)*ixf(2)*ixf(3)
                c101 = onemdxf(2)*onemdxf(3)*ixf(1)
                c100 = onemdxf(2)*ixf(1)*ixf(3)
                c110 = onemdxf(3)*ixf(1)*ixf(2)
                c111 = ixf(1)*ixf(2)*ixf(3)

                ctotal = (c000 + c001 + c010 + c011 + c100 + c101 + c110 + c111)

                r000 = c000/ctotal
                r001 = c001/ctotal
                r010 = c010/ctotal
                r011 = c011/ctotal
                r100 = c100/ctotal
                r101 = c101/ctotal
                r110 = c110/ctotal
                r111 = c111/ctotal

                if (visc_type .gt. 0) then
                  localbeta = beta(fi,fj,fk)
                else
                  !3d visc
                  localbeta = r000*beta(fi,fj,fk) + r001*beta(fi,fj,fk+1) + r010*beta(fi,fj+1,fk) + r011*beta(fi,fj+1,fk+1) + r100*beta(fi+1,fj,fk) + r101*beta(fi+1,fj,fk+1) + r110*beta(fi+1,fj+1,fk) + r111*beta(fi+1,fj+1,fk+1)
                endif

                localvel(1) = r000*velx(fi,fj,fk) + r001*velx(fi,fj,fk+1) + r010*velx(fi,fj+1,fk) + r011*velx(fi,fj+1,fk+1) + r100*velx(fi+1,fj,fk) + r101*velx(fi+1,fj,fk+1) + r110*velx(fi+1,fj+1,fk) + r111*velx(fi+1,fj+1,fk+1)
                localvel(2) = r000*vely(fi,fj,fk) + r001*vely(fi,fj,fk+1) + r010*vely(fi,fj+1,fk) + r011*vely(fi,fj+1,fk+1) + r100*vely(fi+1,fj,fk) + r101*vely(fi+1,fj,fk+1) + r110*vely(fi+1,fj+1,fk) + r111*vely(fi+1,fj+1,fk+1)
                localvel(3) = r000*velz(fi,fj,fk) + r001*velz(fi,fj,fk+1) + r010*velz(fi,fj+1,fk) + r011*velz(fi,fj+1,fk+1) + r100*velz(fi+1,fj,fk) + r101*velz(fi+1,fj,fk+1) + r110*velz(fi+1,fj+1,fk) + r111*velz(fi+1,fj+1,fk+1)

#endif

#if (BL_SPACEDIM == 2)

                c00 = onemdxf(1)*onemdxf(2)
                c01 = onemdxf(1)*ixf(2)
                c10 = ixf(1)*onemdxf(2)
                c11 = ixf(1)*ixf(2)

                ctotal = (c00 + c01 + c10 + c11)*2

                r00 = c00/ctotal
                r01 = c01/ctotal
                r10 = c10/ctotal
                r11 = c11/ctotal

                if (visc_type .gt. 0) then
                  localbeta = beta(fi,fj,fk)
                else
                  localbeta = beta(fi,fj,fk)*c00 + beta(fi,fj+1,fk)*c01 + beta(fi+1,fj,fk)*c10 + beta(fi+1,fj+1,fk)*c11
                endif
                !2d xvel
                localvel(1) = velx(fi,fj,fk)*c00 + velx(fi,fj+1,fk)*c01 + velx(fi+1,fj,fk)*c10 + velx(fi+1,fj+1,fk)*c11
                localvel(2) = vely(fi,fj,fk)*c00 + vely(fi,fj+1,fk)*c01 + vely(fi+1,fj,fk)*c10 + vely(fi+1,fj+1,fk)*c11
#endif

                !Brownian forcing

                call get_particle_normal(normalrand(1))
                call get_particle_normal(normalrand(2))
                call get_particle_normal(normalrand(3))

                runerr = 1d0;

                tempvel(1) = part%vel(1)
                tempvel(2) = part%vel(2)
#if (BL_SPACEDIM == 3)
                tempvel(3) = part%vel(3)
#endif

                inttime = runtime
                intold = inttime

                do while (runerr .gt. 0.1)                
                  !Assuming T=1 for now
                 std = sqrt(-part%drag_factor*localbeta*k_B*2d0*intold*1d0)/part%mass
                 std = 0

                 bfac(1) = std*normalrand(1)
                 bfac(2) = std*normalrand(2)
                 bfac(3) = std*normalrand(3)

                !Semi-implicit Euler velocity and position update

                 deltap(1) = tempvel(1)
                 deltap(2) = tempvel(2)
#if (BL_SPACEDIM == 3)
                 deltap(3) = tempvel(3)
#endif

                tempvel(1) = part%accel_factor*localbeta*(part%vel(1)-localvel(1))*intold + bfac(1) + part%vel(1)
                tempvel(2) = part%accel_factor*localbeta*(part%vel(2)-localvel(2))*intold + bfac(2) + part%vel(2)
#if (BL_SPACEDIM == 3)
                tempvel(3) = part%accel_factor*localbeta*(part%vel(3)-localvel(3))*intold + bfac(3) + part%vel(3)
#endif
#if (BL_SPACEDIM == 3)                
                call find_intersect(tempvel(1),tempvel(2),tempvel(3),part%pos(1),part%pos(2),part%pos(3),intold, surfaces, ns, intsurf, inttime, intside)
#endif
                runerr = abs((inttime - intold)/intold)

                intold = inttime

               end do

               part%vel(1) = tempvel(1)
               part%vel(2) = tempvel(2)
#if (BL_SPACEDIM == 3)
               part%vel(3) = tempvel(3)
#endif

               deltap(1) = part%mass*(part%vel(1) - deltap(1))
               deltap(2) = part%mass*(part%vel(2) - deltap(2))
#if (BL_SPACEDIM == 3)
               deltap(3) = part%mass*(part%vel(3) - deltap(3))
#endif

#if (BL_SPACEDIM == 3)
              !distribute x momentum change 
               nodalp = r000*deltap(1)
               sourcex(fi,fj,fk) = sourcex(fi,fj,fk) + nodalp
               sourcex(fi,fj-1,fk) = sourcex(fi,fj-1,fk) + nodalp
               sourcex(fi,fj,fk-1) = sourcex(fi,fj,fk-1) + nodalp
               sourcex(fi,fj-1,fk-1) = sourcex(fi,fj-1,fk-1) + nodalp

               nodalp = r001*deltap(1)
               sourcex(fi,fj,fk+1) = sourcex(fi,fj,fk+1) + nodalp
               sourcex(fi,fj-1,fk+1) = sourcex(fi,fj-1,fk+1) + nodalp
               sourcex(fi,fj,fk) = sourcex(fi,fj,fk) + nodalp
               sourcex(fi,fj-1,fk) = sourcex(fi,fj-1,fk) + nodalp

               nodalp = r010*deltap(1)
               sourcex(fi,fj+1,fk) = sourcex(fi,fj+1,fk) + nodalp
               sourcex(fi,fj,fk) = sourcex(fi,fj,fk) + nodalp
               sourcex(fi,fj+1,fk-1) = sourcex(fi,fj+1,fk-1) + nodalp
               sourcex(fi,fj,fk-1) = sourcex(fi,fj,fk-1) + nodalp

               nodalp = r011*deltap(1)
               sourcex(fi,fj+1,fk+1) = sourcex(fi,fj+1,fk+1) + nodalp
               sourcex(fi,fj,fk+1) = sourcex(fi,fj,fk+1) + nodalp
               sourcex(fi,fj+1,fk) = sourcex(fi,fj+1,fk) + nodalp
               sourcex(fi,fj,fk) = sourcex(fi,fj,fk) + nodalp

               nodalp = r100*deltap(1)
               sourcex(fi+1,fj,fk) = sourcex(fi+1,fj,fk) + nodalp
               sourcex(fi+1,fj-1,fk) = sourcex(fi+1,fj-1,fk) + nodalp
               sourcex(fi+1,fj,fk-1) = sourcex(fi+1,fj,fk-1) + nodalp
               sourcex(fi+1,fj-1,fk-1) = sourcex(fi+1,fj-1,fk-1) + nodalp

               nodalp = r101*deltap(1)
               sourcex(fi+1,fj,fk+1) = sourcex(fi+1,fj,fk+1) + nodalp
               sourcex(fi+1,fj-1,fk+1) = sourcex(fi+1,fj-1,fk+1) + nodalp
               sourcex(fi+1,fj,fk) = sourcex(fi+1,fj,fk) + nodalp
               sourcex(fi+1,fj-1,fk) = sourcex(fi+1,fj-1,fk) + nodalp

               nodalp = r110*deltap(1)
               sourcex(fi+1,fj+1,fk) = sourcex(fi+1,fj+1,fk) + nodalp
               sourcex(fi+1,fj,fk) = sourcex(fi+1,fj,fk) + nodalp
               sourcex(fi+1,fj+1,fk-1) = sourcex(fi+1,fj+1,fk-1) + nodalp
               sourcex(fi+1,fj,fk-1) = sourcex(fi+1,fj,fk-1) + nodalp

               nodalp = r111*deltap(1)
               sourcex(fi+1,fj+1,fk+1) = sourcex(fi+1,fj+1,fk+1) + nodalp
               sourcex(fi+1,fj,fk+1) = sourcex(fi+1,fj,fk+1) + nodalp
               sourcex(fi+1,fj+1,fk) = sourcex(fi+1,fj+1,fk) + nodalp
               sourcex(fi+1,fj,fk) = sourcex(fi+1,fj,fk) + nodalp

              !distribute y momentum change 
               nodalp = r000*deltap(2)
               sourcey(fi,fj,fk) = sourcey(fi,fj,fk) + nodalp
               sourcey(fi-1,fj,fk) = sourcey(fi-1,fj,fk) + nodalp
               sourcey(fi,fj,fk-1) = sourcey(fi,fj,fk-1) + nodalp
               sourcey(fi-1,fj,fk-1) = sourcey(fi-1,fj,fk-1) + nodalp

               nodalp = r001*deltap(2)
               sourcey(fi,fj,fk+1) = sourcey(fi,fj,fk+1) + nodalp
               sourcey(fi-1,fj,fk+1) = sourcey(fi-1,fj,fk+1) + nodalp
               sourcey(fi,fj,fk) = sourcey(fi,fj,fk) + nodalp
               sourcey(fi-1,fj,fk) = sourcey(fi-1,fj,fk) + nodalp

               nodalp = r010*deltap(2)
               sourcey(fi,fj+1,fk) = sourcey(fi,fj+1,fk) + nodalp
               sourcey(fi-1,fj+1,fk) = sourcey(fi-1,fj+1,fk) + nodalp
               sourcey(fi,fj+1,fk-1) = sourcey(fi,fj+1,fk-1) + nodalp
               sourcey(fi-1,fj+1,fk-1) = sourcey(fi-1,fj+1,fk-1) + nodalp

               nodalp = r011*deltap(2)
               sourcey(fi,fj+1,fk+1) = sourcey(fi,fj+1,fk+1) + nodalp
               sourcey(fi-1,fj+1,fk+1) = sourcey(fi-1,fj+1,fk+1) + nodalp
               sourcey(fi,fj+1,fk) = sourcey(fi,fj+1,fk) + nodalp
               sourcey(fi-1,fj+1,fk) = sourcey(fi-1,fj+1,fk) + nodalp

               nodalp = r100*deltap(2)
               sourcey(fi+1,fj,fk) = sourcey(fi+1,fj,fk) + nodalp
               sourcey(fi,fj,fk) = sourcey(fi,fj,fk) + nodalp
               sourcey(fi+1,fj,fk-1) = sourcey(fi+1,fj,fk-1) + nodalp
               sourcey(fi,fj,fk-1) = sourcey(fi,fj,fk-1) + nodalp

               nodalp = r101*deltap(2)
               sourcey(fi+1,fj,fk+1) = sourcey(fi+1,fj,fk+1) + nodalp
               sourcey(fi,fj,fk+1) = sourcey(fi,fj,fk+1) + nodalp
               sourcey(fi+1,fj,fk) = sourcey(fi+1,fj,fk) + nodalp
               sourcey(fi,fj,fk) = sourcey(fi,fj,fk) + nodalp

               nodalp = r110*deltap(2)
               sourcey(fi+1,fj+1,fk) = sourcey(fi+1,fj+1,fk) + nodalp
               sourcey(fi,fj+1,fk) = sourcey(fi,fj+1,fk) + nodalp
               sourcey(fi+1,fj+1,fk-1) = sourcey(fi+1,fj+1,fk-1) + nodalp
               sourcey(fi,fj+1,fk-1) = sourcey(fi,fj+1,fk-1) + nodalp

               nodalp = r111*deltap(2)
               sourcey(fi+1,fj+1,fk+1) = sourcey(fi+1,fj+1,fk+1) + nodalp
               sourcey(fi,fj+1,fk+1) = sourcey(fi,fj+1,fk+1) + nodalp
               sourcey(fi+1,fj+1,fk) = sourcey(fi+1,fj+1,fk) + nodalp
               sourcey(fi,fj+1,fk) = sourcey(fi,fj+1,fk) + nodalp

              !distribute z momentum change 
               nodalp = r000*deltap(3)
               sourcez(fi,fj,fk) = sourcez(fi,fj,fk) + nodalp
               sourcez(fi-1,fj,fk) = sourcez(fi-1,fj,fk) + nodalp
               sourcez(fi,fj-1,fk) = sourcez(fi,fj-1,fk) + nodalp
               sourcez(fi-1,fj-1,fk) = sourcez(fi-1,fj-1,fk) + nodalp

               nodalp = r001*deltap(3)
               sourcez(fi,fj,fk+1) = sourcez(fi,fj,fk+1) + nodalp
               sourcez(fi-1,fj,fk+1) = sourcez(fi-1,fj,fk+1) + nodalp
               sourcez(fi,fj-1,fk+1) = sourcez(fi,fj-1,fk+1) + nodalp
               sourcez(fi-1,fj-1,fk+1) = sourcez(fi-1,fj-1,fk+1) + nodalp

               nodalp = r010*deltap(3)
               sourcez(fi,fj+1,fk) = sourcez(fi,fj+1,fk) + nodalp
               sourcez(fi-1,fj+1,fk) = sourcez(fi-1,fj+1,fk) + nodalp
               sourcez(fi,fj,fk) = sourcez(fi,fj,fk) + nodalp
               sourcez(fi-1,fj,fk) = sourceZ(fi-1,fj,fk) + nodalp

               nodalp = r011*deltap(3)
               sourcez(fi,fj+1,fk+1) = sourcez(fi,fj+1,fk+1) + nodalp
               sourcez(fi-1,fj+1,fk+1) = sourcez(fi-1,fj+1,fk+1) + nodalp
               sourcez(fi,fj,fk+1) = sourcez(fi,fj,fk+1) + nodalp
               sourcez(fi-1,fj,fk+1) = sourcez(fi-1,fj,fk+1) + nodalp

               nodalp = r100*deltap(3)
               sourcez(fi+1,fj,fk) = sourcez(fi+1,fj,fk) + nodalp
               sourcez(fi,fj,fk) = sourcez(fi,fj,fk) + nodalp
               sourcez(fi+1,fj-1,fk) = sourcez(fi+1,fj-1,fk) + nodalp
               sourcez(fi,fj-1,fk) = sourcez(fi,fj-1,fk) + nodalp

               nodalp = r101*deltap(3)
               sourcez(fi+1,fj,fk+1) = sourcez(fi+1,fj,fk+1) + nodalp
               sourcez(fi,fj,fk+1) = sourcez(fi,fj,fk+1) + nodalp
               sourcez(fi+1,fj-1,fk+1) = sourcez(fi+1,fj-1,fk+1) + nodalp
               sourcez(fi,fj-1,fk+1) = sourcez(fi,fj-1,fk+1) + nodalp

               nodalp = r110*deltap(3)
               sourcez(fi+1,fj+1,fk) = sourcez(fi+1,fj+1,fk) + nodalp
               sourcez(fi,fj+1,fk) = sourcez(fi,fj+1,fk) + nodalp
               sourcez(fi+1,fj,fk) = sourcez(fi+1,fj,fk) + nodalp
               sourcez(fi,fj,fk) = sourcez(fi,fj,fk) + nodalp

               nodalp = r111*deltap(3)
               sourcez(fi+1,fj+1,fk+1) = sourcez(fi+1,fj+1,fk+1) + nodalp
               sourcez(fi,fj+1,fk+1) = sourcez(fi,fj+1,fk+1) + nodalp
               sourcez(fi+1,fj,fk+1) = sourcez(fi+1,fj,fk+1) + nodalp
               sourcez(fi,fj,fk+1) = sourcez(fi,fj,fk+1) + nodalp
#endif

#if (BL_SPACEDIM == 2)
              !distribute x momentum change 
               nodalp = r00*deltap(1)
               sourcex(fi,fj,fk) = sourcex(fi,fj,fk) + nodalp
               sourcex(fi,fj-1,fk) = sourcex(fi,fj-1,fk) + nodalp

               nodalp = r01*deltap(1)
               sourcex(fi,fj+1,fk) = sourcex(fi,fj+1,fk) + nodalp
               sourcex(fi,fj,fk) = sourcex(fi,fj,fk) + nodalp

               nodalp = r10*deltap(1)
               sourcex(fi+1,fj,fk) = sourcex(fi+1,fj,fk) + nodalp
               sourcex(fi+1,fj-1,fk) = sourcex(fi+1,fj-1,fk) + nodalp

               nodalp = r11*deltap(1)
               sourcex(fi+1,fj+1,fk) = sourcex(fi+1,fj+1,fk) + nodalp
               sourcex(fi+1,fj,fk) = sourcex(fi+1,fj,fk) + nodalp

              !distribute y momentum change
               nodalp = r00*deltap(2)
               sourcey(fi,fj,fk) = sourcey(fi,fj,fk) + nodalp
               sourcey(fi-1,fj,fk) = sourcey(fi-1,fj,fk) + nodalp

               nodalp = r01*deltap(2)
               sourcey(fi,fj+1,fk) = sourcey(fi,fj+1,fk) + nodalp
               sourcey(fi-1,fj+1,fk) = sourcey(fi-1,fj+1,fk) + nodalp

               nodalp = r10*deltap(2)
               sourcey(fi+1,fj,fk) = sourcey(fi+1,fj,fk) + nodalp
               sourcey(fi,fj,fk) = sourcey(fi,fj,fk) + nodalp

               nodalp = r11*deltap(2)
               sourcey(fi+1,fj+1,fk) = sourcey(fi+1,fj+1,fk) + nodalp
               sourcey(fi,fj+1,fk) = sourcey(fi,fj+1,fk) + nodalp
#endif
                ! move the particle in a straight line
                part%pos(1) = part%pos(1) + dt*part%vel(1)
                part%pos(2) = part%pos(2) + dt*part%vel(2)
#if (BL_SPACEDIM == 3)
                part%pos(3) = part%pos(3) + dt*part%vel(3)
#endif

                runtime = runtime - inttime

              enddo

              ! if it has changed cells, remove from vector.
              ! otherwise continue
              ni = floor((part%pos(1) - plo(1))*dxinv(1))           
              nj = floor((part%pos(2) - plo(2))*dxinv(2))
#if (BL_SPACEDIM == 3)
              nk = floor((part%pos(3) - plo(3))*dxinv(3))
#else
              nk = 0
#endif

              if ((ni /= i) .or. (nj /= j) .or. (nk /= k)) then
                 part%sorted = 0
                 call remove_particle_from_cell(cell_parts, cell_np, new_np, p)  
              else
                 p = p + 1
              end if
           end do

           cell_part_cnt(i,j,k) = new_np
    
        end do
     end do
  end do
  
end subroutine move_particles_fhd

