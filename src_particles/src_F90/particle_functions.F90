subroutine move_particles_dsmc(particles, np, lo, hi, &
     cell_part_ids, cell_part_cnt, clo, chi, plo, phi, dx, dt, surfaces, ns) &
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
  real(amrex_real), intent(in) :: plo(3), phi(3), dx(3)
  real(amrex_real), intent(in) :: dt
  
  integer :: i, j, k, p, cell_np, new_np, intsurf, intside, push, intcount
  integer :: cell(3)
  integer(c_int), pointer :: cell_parts(:)
  type(particle_t), pointer :: part
  type(surface_t), pointer :: surf
  real(amrex_real) inv_dx(3), runtime, inttime, adjalt, adj, inv_dt, domsize(3), posalt(3), prex, postx

  adj = 0.9999999
  adjalt = 2d0*(1d0 - adj)

  inv_dx = 1.d0/dx
  inv_dt = 1.d0/dt
  
  domsize = phi - plo

  do p = 1, ns

    surf => surfaces(p)  

    surf%fxleft = 0
    surf%fyleft = 0
    surf%fzleft = 0

    surf%fxright = 0
    surf%fyright = 0
    surf%fzright = 0

  enddo        
  
  intcount = 0

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           cell_np = cell_part_cnt(i,j,k)
           call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])

           new_np = cell_np
           p = 1
           do while (p <= new_np)

              !if(cell_parts(p) .eq. 1320) then
                
                !print *, "Accessing ", p, " which is ", cell_parts(p), " and has id ", i,j,k

            !  endif
              part => particles(cell_parts(p))

              runtime = dt

              do while (runtime .gt. 0)

                prex = part%vel(1)

                call find_intersect(part,runtime, surfaces, ns, intsurf, inttime, intside, phi, plo)

                if(intsurf .eq. 6) then

                  intcount = intcount + 1

                endif
                
                 
                !print *, "Parrticle ", p, " intersect ", inttime, intsurf 

                !print *, "Prepos ", part%pos, " prevel ", part%vel

                !call sleep(1)

                posalt(1) = inttime*part%vel(1)*adjalt
                posalt(2) = inttime*part%vel(2)*adjalt
#if (BL_SPACEDIM == 3)
                posalt(3) = inttime*part%vel(3)*adjalt
#endif

                ! move the particle in a straight line, adj factor prevents double detection of boundary intersection
                part%pos(1) = part%pos(1) + inttime*part%vel(1)*adj
                part%pos(2) = part%pos(2) + inttime*part%vel(2)*adj
#if (BL_SPACEDIM == 3)
                part%pos(3) = part%pos(3) + inttime*part%vel(3)*adj
#endif
                runtime = runtime - inttime

                if(intsurf .gt. 0) then

                  surf => surfaces(intsurf)

                  call apply_bc(surf, part, intside, domsize, push)

                    if(push .eq. 1) then
                      
                      part%pos(1) = part%pos(1) + posalt(1)
                      part%pos(2) = part%pos(2) + posalt(2)
#if (BL_SPACEDIM == 3)
                      part%pos(3) = part%pos(3) + posalt(3)
#endif
                    endif
                    
                endif

                !print *, "Postpos ", part%pos, " postvel ", part%vel
                postx = part%vel(1)

              end do

              if(prex .ne. postx) then
                print *, "Oops!"
              endif

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

  !print *, "intcount: ", intcount

  do p = 1, ns

    surf => surfaces(p)  

    surf%fxleft = surf%fxleft*inv_dt
    surf%fyleft = surf%fyleft*inv_dt
    surf%fzleft = surf%fzleft*inv_dt

    surf%fxright = surf%fxright*inv_dt
    surf%fyright = surf%fyright*inv_dt
    surf%fzright = surf%fzright*inv_dt

  enddo  
  
end subroutine move_particles_dsmc

subroutine redirect(part)

    use cell_sorted_particle_module, only: particle_t

    implicit none

    type(particle_t), intent(inout) :: part

    double precision speed

    speed = sqrt(part%vel(1)**2 + part%vel(2)**2 + part%vel(3)**2)

    if(speed .ne. 0) then
      part%dir = part%vel/speed
    endif

end subroutine redirect

subroutine get_interpolation_weights(cc, rr, ixf, onemdxf)

  use amrex_fort_module, only: amrex_real

  double precision, intent(in   )  :: ixf(3), onemdxf(3)

  double precision, intent(inout)  :: cc(0:7)
  double precision, intent(inout)  :: rr(0:7)

#if (BL_SPACEDIM == 3)
  cc(0) = onemdxf(1)*onemdxf(2)*onemdxf(3)
  cc(4) = onemdxf(1)*onemdxf(2)*ixf(3)
  cc(2) = onemdxf(1)*onemdxf(3)*ixf(2)
  cc(6) = onemdxf(1)*ixf(2)*ixf(3)
  cc(5) = onemdxf(2)*onemdxf(3)*ixf(1)
  cc(1) = onemdxf(2)*ixf(1)*ixf(3)
  cc(3) = onemdxf(3)*ixf(1)*ixf(2)
  cc(7) = ixf(1)*ixf(2)*ixf(3)

  rr = cc/sum(cc)
#endif

#if (BL_SPACEDIM == 2)
  cc(0) = onemdxf(1)*onemdxf(2)
  cc(2) = onemdxf(1)*ixf(2)
  cc(1) = ixf(1)*onemdxf(2)
  cc(3) = ixf(1)*ixf(2)

  rr = cc/(sum(cc)*2)
#endif

end subroutine get_interpolation_weights

#if (BL_SPACEDIM == 3)
subroutine get_local_properties(cc, rr, fi, velx, velxlo, velxhi, vely, velylo, velyhi, velz, velzlo, velzhi, beta, localvel, localbeta, betalo, betahi)
#endif
#if (BL_SPACEDIM == 2)
subroutine get_local_properties(cc, rr, fi, velx, velxlo, velxhi, vely, velylo, velyhi, beta, localvel, localbeta, betalo, betahi)
#endif

  use amrex_fort_module, only: amrex_real
  use common_namelist_module, only: visc_type, k_B
  
  implicit none

  integer,          intent(in   )         :: fi(3), velxlo(3), velxhi(3), velylo(3), velyhi(3), betalo(3), betahi(3)
#if (AMREX_SPACEDIM == 3)
  integer,          intent(in   )         :: velzlo(3), velzhi(3)
#endif
  double precision, intent(inout)         :: rr(0:7), localvel(3), localbeta, cc(0:7)

  double precision, intent(in   ) :: velx(velxlo(1):velxhi(1),velxlo(2):velxhi(2),velxlo(3):velxhi(3))
  double precision, intent(in   ) :: vely(velylo(1):velyhi(1),velylo(2):velyhi(2),velylo(3):velyhi(3))
#if (AMREX_SPACEDIM == 3)
  double precision, intent(in   ) :: velz(velzlo(1):velzhi(1),velzlo(2):velzhi(2),velzlo(3):velzhi(3))
#endif

  double precision, intent(in   ) :: beta(betalo(1):betahi(1),betalo(2):betahi(2),betalo(3):betahi(3))

#if (BL_SPACEDIM == 3)

  if (visc_type .gt. 0) then
    localbeta = beta(fi(1),fi(2),fi(3))
  else
    !3d visc
    localbeta = rr(0)*beta(fi(1),fi(2),fi(3)) + rr(4)*beta(fi(1),fi(2),fi(3)+1) + rr(2)*beta(fi(1),fi(2)+1,fi(3)) + rr(6)*beta(fi(1),fi(2)+1,fi(3)+1) + rr(1)*beta(fi(1)+1,fi(2),fi(3)) + rr(5)*beta(fi(1)+1,fi(2),fi(3)+1) + rr(3)*beta(fi(1)+1,fi(2)+1,fi(3)) + rr(7)*beta(fi(1)+1,fi(2)+1,fi(3)+1)
  endif

  localvel(1) = rr(0)*velx(fi(1),fi(2),fi(3)) + rr(4)*velx(fi(1),fi(2),fi(3)+1) + rr(2)*velx(fi(1),fi(2)+1,fi(3)) + rr(6)*velx(fi(1),fi(2)+1,fi(3)+1) + rr(1)*velx(fi(1)+1,fi(2),fi(3)) + rr(5)*velx(fi(1)+1,fi(2),fi(3)+1) + rr(3)*velx(fi(1)+1,fi(2)+1,fi(3)) + rr(7)*velx(fi(1)+1,fi(2)+1,fi(3)+1)
  localvel(2) = rr(0)*vely(fi(1),fi(2),fi(3)) + rr(4)*vely(fi(1),fi(2),fi(3)+1) + rr(2)*vely(fi(1),fi(2)+1,fi(3)) + rr(6)*vely(fi(1),fi(2)+1,fi(3)+1) + rr(1)*vely(fi(1)+1,fi(2),fi(3)) + rr(5)*vely(fi(1)+1,fi(2),fi(3)+1) + rr(3)*vely(fi(1)+1,fi(2)+1,fi(3)) + rr(7)*vely(fi(1)+1,fi(2)+1,fi(3)+1)
  localvel(3) = rr(0)*velz(fi(1),fi(2),fi(3)) + rr(4)*velz(fi(1),fi(2),fi(3)+1) + rr(2)*velz(fi(1),fi(2)+1,fi(3)) + rr(6)*velz(fi(1),fi(2)+1,fi(3)+1) + rr(1)*velz(fi(1)+1,fi(2),fi(3)) + rr(5)*velz(fi(1)+1,fi(2),fi(3)+1) + rr(3)*velz(fi(1)+1,fi(2)+1,fi(3)) + rr(7)*velz(fi(1)+1,fi(2)+1,fi(3)+1)
#endif

#if (BL_SPACEDIM == 2)

!CHECK THIS!

  if (visc_type .gt. 0) then
    localbeta = beta(fi(1),fi(2),fi(3))
  else
    localbeta = beta(fi(1),fi(2),fi(3))*rr(0) + beta(fi(1),fi(2)+1,fi(3))*rr(2) + beta(fi(1)+1,fi(2),fi(3))*rr(1) + beta(fi(1)+1,fi(2)+1,fi(3))*rr(3)
  endif
    !2d xvel
  localvel(1) = velx(fi(1),fi(2),fi(3))*rr(0) + velx(fi(1),fi(2)+1,fi(3))*rr(2) + velx(fi(1)+1,fi(2),fi(3))*rr(1) + velx(fi(1)+1,fi(2)+1,fi(3))*rr(3)
  localvel(2) = vely(fi(1),fi(2),fi(3))*rr(0) + vely(fi(1),fi(2)+1,fi(3))*rr(2) + vely(fi(1)+1,fi(2),fi(3))*rr(1) + vely(fi(1)+1,fi(2)+1,fi(3))*rr(3)
#endif

end subroutine get_local_properties


#if (BL_SPACEDIM == 3)
subroutine distribute_momentum(deltap, rr, fi ,sourcex, sourcexlo, sourcexhi, sourcey, sourceylo, sourceyhi, sourcez, sourcezlo, sourcezhi)
#endif
#if (BL_SPACEDIM == 2)
subroutine distribute_momentum(deltap, rr, fi ,sourcex, sourcexlo, sourcexhi, sourcey, sourceylo, sourceyhi)
#endif

  integer,          intent(in   ) :: fi(3), sourcexlo(3), sourcexhi(3), sourceylo(3), sourceyhi(3)
#if (AMREX_SPACEDIM == 3)
  integer,          intent(in   ) :: sourcezlo(3), sourcezhi(3)
#endif

  double precision, intent(in   ) :: rr(0:7)

  double precision, intent(inout) :: deltap(3)

  double precision, intent(inout) :: sourcex(sourcexlo(1):sourcexhi(1),sourcexlo(2):sourcexhi(2),sourcexlo(3):sourcexhi(3))
  double precision, intent(inout) :: sourcey(sourceylo(1):sourceyhi(1),sourceylo(2):sourceyhi(2),sourceylo(3):sourceyhi(3))
#if (AMREX_SPACEDIM == 3)
  double precision, intent(inout) :: sourcez(sourcezlo(1):sourcezhi(1),sourcezlo(2):sourcezhi(2),sourcezlo(3):sourcezhi(3))
#endif

  double precision nodalp

#if (BL_SPACEDIM == 3)
  !distribute x momentum change 
   nodalp = rr(0)*deltap(1)
   sourcex(fi(1),fi(2),fi(3)) = sourcex(fi(1),fi(2),fi(3)) + nodalp
   sourcex(fi(1),fi(2)-1,fi(3)) = sourcex(fi(1),fi(2)-1,fi(3)) + nodalp
   sourcex(fi(1),fi(2),fi(3)-1) = sourcex(fi(1),fi(2),fi(3)-1) + nodalp
   sourcex(fi(1),fi(2)-1,fi(3)-1) = sourcex(fi(1),fi(2)-1,fi(3)-1) + nodalp

   nodalp = rr(4)*deltap(1)
   sourcex(fi(1),fi(2),fi(3)+1) = sourcex(fi(1),fi(2),fi(3)+1) + nodalp
   sourcex(fi(1),fi(2)-1,fi(3)+1) = sourcex(fi(1),fi(2)-1,fi(3)+1) + nodalp
   sourcex(fi(1),fi(2),fi(3)) = sourcex(fi(1),fi(2),fi(3)) + nodalp
   sourcex(fi(1),fi(2)-1,fi(3)) = sourcex(fi(1),fi(2)-1,fi(3)) + nodalp

   nodalp = rr(2)*deltap(1)
   sourcex(fi(1),fi(2)+1,fi(3)) = sourcex(fi(1),fi(2)+1,fi(3)) + nodalp
   sourcex(fi(1),fi(2),fi(3)) = sourcex(fi(1),fi(2),fi(3)) + nodalp
   sourcex(fi(1),fi(2)+1,fi(3)-1) = sourcex(fi(1),fi(2)+1,fi(3)-1) + nodalp
   sourcex(fi(1),fi(2),fi(3)-1) = sourcex(fi(1),fi(2),fi(3)-1) + nodalp

   nodalp = rr(6)*deltap(1)
   sourcex(fi(1),fi(2)+1,fi(3)+1) = sourcex(fi(1),fi(2)+1,fi(3)+1) + nodalp
   sourcex(fi(1),fi(2),fi(3)+1) = sourcex(fi(1),fi(2),fi(3)+1) + nodalp
   sourcex(fi(1),fi(2)+1,fi(3)) = sourcex(fi(1),fi(2)+1,fi(3)) + nodalp
   sourcex(fi(1),fi(2),fi(3)) = sourcex(fi(1),fi(2),fi(3)) + nodalp

   nodalp = rr(1)*deltap(1)
   sourcex(fi(1)+1,fi(2),fi(3)) = sourcex(fi(1)+1,fi(2),fi(3)) + nodalp
   sourcex(fi(1)+1,fi(2)-1,fi(3)) = sourcex(fi(1)+1,fi(2)-1,fi(3)) + nodalp
   sourcex(fi(1)+1,fi(2),fi(3)-1) = sourcex(fi(1)+1,fi(2),fi(3)-1) + nodalp
   sourcex(fi(1)+1,fi(2)-1,fi(3)-1) = sourcex(fi(1)+1,fi(2)-1,fi(3)-1) + nodalp

   nodalp = rr(5)*deltap(1)
   sourcex(fi(1)+1,fi(2),fi(3)+1) = sourcex(fi(1)+1,fi(2),fi(3)+1) + nodalp
   sourcex(fi(1)+1,fi(2)-1,fi(3)+1) = sourcex(fi(1)+1,fi(2)-1,fi(3)+1) + nodalp
   sourcex(fi(1)+1,fi(2),fi(3)) = sourcex(fi(1)+1,fi(2),fi(3)) + nodalp
   sourcex(fi(1)+1,fi(2)-1,fi(3)) = sourcex(fi(1)+1,fi(2)-1,fi(3)) + nodalp

   nodalp = rr(3)*deltap(1)
   sourcex(fi(1)+1,fi(2)+1,fi(3)) = sourcex(fi(1)+1,fi(2)+1,fi(3)) + nodalp
   sourcex(fi(1)+1,fi(2),fi(3)) = sourcex(fi(1)+1,fi(2),fi(3)) + nodalp
   sourcex(fi(1)+1,fi(2)+1,fi(3)-1) = sourcex(fi(1)+1,fi(2)+1,fi(3)-1) + nodalp
   sourcex(fi(1)+1,fi(2),fi(3)-1) = sourcex(fi(1)+1,fi(2),fi(3)-1) + nodalp

   nodalp = rr(7)*deltap(1)
   sourcex(fi(1)+1,fi(2)+1,fi(3)+1) = sourcex(fi(1)+1,fi(2)+1,fi(3)+1) + nodalp
   sourcex(fi(1)+1,fi(2),fi(3)+1) = sourcex(fi(1)+1,fi(2),fi(3)+1) + nodalp
   sourcex(fi(1)+1,fi(2)+1,fi(3)) = sourcex(fi(1)+1,fi(2)+1,fi(3)) + nodalp
   sourcex(fi(1)+1,fi(2),fi(3)) = sourcex(fi(1)+1,fi(2),fi(3)) + nodalp

  !distribute y momentum change 
   nodalp = rr(0)*deltap(2)
   sourcey(fi(1),fi(2),fi(3)) = sourcey(fi(1),fi(2),fi(3)) + nodalp
   sourcey(fi(1)-1,fi(2),fi(3)) = sourcey(fi(1)-1,fi(2),fi(3)) + nodalp
   sourcey(fi(1),fi(2),fi(3)-1) = sourcey(fi(1),fi(2),fi(3)-1) + nodalp
   sourcey(fi(1)-1,fi(2),fi(3)-1) = sourcey(fi(1)-1,fi(2),fi(3)-1) + nodalp

   nodalp = rr(4)*deltap(2)
   sourcey(fi(1),fi(2),fi(3)+1) = sourcey(fi(1),fi(2),fi(3)+1) + nodalp
   sourcey(fi(1)-1,fi(2),fi(3)+1) = sourcey(fi(1)-1,fi(2),fi(3)+1) + nodalp
   sourcey(fi(1),fi(2),fi(3)) = sourcey(fi(1),fi(2),fi(3)) + nodalp
   sourcey(fi(1)-1,fi(2),fi(3)) = sourcey(fi(1)-1,fi(2),fi(3)) + nodalp

   nodalp = rr(2)*deltap(2)
   sourcey(fi(1),fi(2)+1,fi(3)) = sourcey(fi(1),fi(2)+1,fi(3)) + nodalp
   sourcey(fi(1)-1,fi(2)+1,fi(3)) = sourcey(fi(1)-1,fi(2)+1,fi(3)) + nodalp
   sourcey(fi(1),fi(2)+1,fi(3)-1) = sourcey(fi(1),fi(2)+1,fi(3)-1) + nodalp
   sourcey(fi(1)-1,fi(2)+1,fi(3)-1) = sourcey(fi(1)-1,fi(2)+1,fi(3)-1) + nodalp

   nodalp = rr(6)*deltap(2)
   sourcey(fi(1),fi(2)+1,fi(3)+1) = sourcey(fi(1),fi(2)+1,fi(3)+1) + nodalp
   sourcey(fi(1)-1,fi(2)+1,fi(3)+1) = sourcey(fi(1)-1,fi(2)+1,fi(3)+1) + nodalp
   sourcey(fi(1),fi(2)+1,fi(3)) = sourcey(fi(1),fi(2)+1,fi(3)) + nodalp
   sourcey(fi(1)-1,fi(2)+1,fi(3)) = sourcey(fi(1)-1,fi(2)+1,fi(3)) + nodalp

   nodalp = rr(1)*deltap(2)
   sourcey(fi(1)+1,fi(2),fi(3)) = sourcey(fi(1)+1,fi(2),fi(3)) + nodalp
   sourcey(fi(1),fi(2),fi(3)) = sourcey(fi(1),fi(2),fi(3)) + nodalp
   sourcey(fi(1)+1,fi(2),fi(3)-1) = sourcey(fi(1)+1,fi(2),fi(3)-1) + nodalp
   sourcey(fi(1),fi(2),fi(3)-1) = sourcey(fi(1),fi(2),fi(3)-1) + nodalp

   nodalp = rr(5)*deltap(2)
   sourcey(fi(1)+1,fi(2),fi(3)+1) = sourcey(fi(1)+1,fi(2),fi(3)+1) + nodalp
   sourcey(fi(1),fi(2),fi(3)+1) = sourcey(fi(1),fi(2),fi(3)+1) + nodalp
   sourcey(fi(1)+1,fi(2),fi(3)) = sourcey(fi(1)+1,fi(2),fi(3)) + nodalp
   sourcey(fi(1),fi(2),fi(3)) = sourcey(fi(1),fi(2),fi(3)) + nodalp

   nodalp = rr(3)*deltap(2)
   sourcey(fi(1)+1,fi(2)+1,fi(3)) = sourcey(fi(1)+1,fi(2)+1,fi(3)) + nodalp
   sourcey(fi(1),fi(2)+1,fi(3)) = sourcey(fi(1),fi(2)+1,fi(3)) + nodalp
   sourcey(fi(1)+1,fi(2)+1,fi(3)-1) = sourcey(fi(1)+1,fi(2)+1,fi(3)-1) + nodalp
   sourcey(fi(1),fi(2)+1,fi(3)-1) = sourcey(fi(1),fi(2)+1,fi(3)-1) + nodalp

   nodalp = rr(7)*deltap(2)
   sourcey(fi(1)+1,fi(2)+1,fi(3)+1) = sourcey(fi(1)+1,fi(2)+1,fi(3)+1) + nodalp
   sourcey(fi(1),fi(2)+1,fi(3)+1) = sourcey(fi(1),fi(2)+1,fi(3)+1) + nodalp
   sourcey(fi(1)+1,fi(2)+1,fi(3)) = sourcey(fi(1)+1,fi(2)+1,fi(3)) + nodalp
   sourcey(fi(1),fi(2)+1,fi(3)) = sourcey(fi(1),fi(2)+1,fi(3)) + nodalp

  !distribute z momentum change 
   nodalp = rr(0)*deltap(3)
   sourcez(fi(1),fi(2),fi(3)) = sourcez(fi(1),fi(2),fi(3)) + nodalp
   sourcez(fi(1)-1,fi(2),fi(3)) = sourcez(fi(1)-1,fi(2),fi(3)) + nodalp
   sourcez(fi(1),fi(2)-1,fi(3)) = sourcez(fi(1),fi(2)-1,fi(3)) + nodalp
   sourcez(fi(1)-1,fi(2)-1,fi(3)) = sourcez(fi(1)-1,fi(2)-1,fi(3)) + nodalp

   nodalp = rr(4)*deltap(3)
   sourcez(fi(1),fi(2),fi(3)+1) = sourcez(fi(1),fi(2),fi(3)+1) + nodalp
   sourcez(fi(1)-1,fi(2),fi(3)+1) = sourcez(fi(1)-1,fi(2),fi(3)+1) + nodalp
   sourcez(fi(1),fi(2)-1,fi(3)+1) = sourcez(fi(1),fi(2)-1,fi(3)+1) + nodalp
   sourcez(fi(1)-1,fi(2)-1,fi(3)+1) = sourcez(fi(1)-1,fi(2)-1,fi(3)+1) + nodalp

   nodalp = rr(2)*deltap(3)
   sourcez(fi(1),fi(2)+1,fi(3)) = sourcez(fi(1),fi(2)+1,fi(3)) + nodalp
   sourcez(fi(1)-1,fi(2)+1,fi(3)) = sourcez(fi(1)-1,fi(2)+1,fi(3)) + nodalp
   sourcez(fi(1),fi(2),fi(3)) = sourcez(fi(1),fi(2),fi(3)) + nodalp
   sourcez(fi(1)-1,fi(2),fi(3)) = sourceZ(fi(1)-1,fi(2),fi(3)) + nodalp

   nodalp = rr(6)*deltap(3)
   sourcez(fi(1),fi(2)+1,fi(3)+1) = sourcez(fi(1),fi(2)+1,fi(3)+1) + nodalp
   sourcez(fi(1)-1,fi(2)+1,fi(3)+1) = sourcez(fi(1)-1,fi(2)+1,fi(3)+1) + nodalp
   sourcez(fi(1),fi(2),fi(3)+1) = sourcez(fi(1),fi(2),fi(3)+1) + nodalp
   sourcez(fi(1)-1,fi(2),fi(3)+1) = sourcez(fi(1)-1,fi(2),fi(3)+1) + nodalp

   nodalp = rr(1)*deltap(3)
   sourcez(fi(1)+1,fi(2),fi(3)) = sourcez(fi(1)+1,fi(2),fi(3)) + nodalp
   sourcez(fi(1),fi(2),fi(3)) = sourcez(fi(1),fi(2),fi(3)) + nodalp
   sourcez(fi(1)+1,fi(2)-1,fi(3)) = sourcez(fi(1)+1,fi(2)-1,fi(3)) + nodalp
   sourcez(fi(1),fi(2)-1,fi(3)) = sourcez(fi(1),fi(2)-1,fi(3)) + nodalp

   nodalp = rr(5)*deltap(3)
   sourcez(fi(1)+1,fi(2),fi(3)+1) = sourcez(fi(1)+1,fi(2),fi(3)+1) + nodalp
   sourcez(fi(1),fi(2),fi(3)+1) = sourcez(fi(1),fi(2),fi(3)+1) + nodalp
   sourcez(fi(1)+1,fi(2)-1,fi(3)+1) = sourcez(fi(1)+1,fi(2)-1,fi(3)+1) + nodalp
   sourcez(fi(1),fi(2)-1,fi(3)+1) = sourcez(fi(1),fi(2)-1,fi(3)+1) + nodalp

   nodalp = rr(3)*deltap(3)
   sourcez(fi(1)+1,fi(2)+1,fi(3)) = sourcez(fi(1)+1,fi(2)+1,fi(3)) + nodalp
   sourcez(fi(1),fi(2)+1,fi(3)) = sourcez(fi(1),fi(2)+1,fi(3)) + nodalp
   sourcez(fi(1)+1,fi(2),fi(3)) = sourcez(fi(1)+1,fi(2),fi(3)) + nodalp
   sourcez(fi(1),fi(2),fi(3)) = sourcez(fi(1),fi(2),fi(3)) + nodalp

   nodalp = rr(7)*deltap(3)
   sourcez(fi(1)+1,fi(2)+1,fi(3)+1) = sourcez(fi(1)+1,fi(2)+1,fi(3)+1) + nodalp
   sourcez(fi(1),fi(2)+1,fi(3)+1) = sourcez(fi(1),fi(2)+1,fi(3)+1) + nodalp
   sourcez(fi(1)+1,fi(2),fi(3)+1) = sourcez(fi(1)+1,fi(2),fi(3)+1) + nodalp
   sourcez(fi(1),fi(2),fi(3)+1) = sourcez(fi(1),fi(2),fi(3)+1) + nodalp
#endif

#if (BL_SPACEDIM == 2)
  !distribute x momentum change 
   nodalp = rr(0)*deltap(1)
   sourcex(fi(1),fi(2),fi(3)) = sourcex(fi(1),fi(2),fi(3)) + nodalp
   sourcex(fi(1),fi(2)-1,fi(3)) = sourcex(fi(1),fi(2)-1,fi(3)) + nodalp

   nodalp = rr(2)*deltap(1)
   sourcex(fi(1),fi(2)+1,fi(3)) = sourcex(fi(1),fi(2)+1,fi(3)) + nodalp
   sourcex(fi(1),fi(2),fi(3)) = sourcex(fi(1),fi(2),fi(3)) + nodalp

   nodalp = rr(1)*deltap(1)
   sourcex(fi(1)+1,fi(2),fi(3)) = sourcex(fi(1)+1,fi(2),fi(3)) + nodalp
   sourcex(fi(1)+1,fi(2)-1,fi(3)) = sourcex(fi(1)+1,fi(2)-1,fi(3)) + nodalp

   nodalp = rr(3)*deltap(1)
   sourcex(fi(1)+1,fi(2)+1,fi(3)) = sourcex(fi(1)+1,fi(2)+1,fi(3)) + nodalp
   sourcex(fi(1)+1,fi(2),fi(3)) = sourcex(fi(1)+1,fi(2),fi(3)) + nodalp

  !distribute y momentum change
   nodalp = rr(0)*deltap(2)
   sourcey(fi(1),fi(2),fi(3)) = sourcey(fi(1),fi(2),fi(3)) + nodalp
   sourcey(fi(1)-1,fi(2),fi(3)) = sourcey(fi(1)-1,fi(2),fi(3)) + nodalp

   nodalp = rr(2)*deltap(2)
   sourcey(fi(1),fi(2)+1,fi(3)) = sourcey(fi(1),fi(2)+1,fi(3)) + nodalp
   sourcey(fi(1)-1,fi(2)+1,fi(3)) = sourcey(fi(1)-1,fi(2)+1,fi(3)) + nodalp

   nodalp = rr(1)*deltap(2)
   sourcey(fi(1)+1,fi(2),fi(3)) = sourcey(fi(1)+1,fi(2),fi(3)) + nodalp
   sourcey(fi(1),fi(2),fi(3)) = sourcey(fi(1),fi(2),fi(3)) + nodalp

   nodalp = rr(3)*deltap(2)
   sourcey(fi(1)+1,fi(2)+1,fi(3)) = sourcey(fi(1)+1,fi(2)+1,fi(3)) + nodalp
   sourcey(fi(1),fi(2)+1,fi(3)) = sourcey(fi(1),fi(2)+1,fi(3)) + nodalp
#endif

end subroutine distribute_momentum

subroutine move_particles_fhd(particles, np, lo, hi, &
     cell_part_ids, cell_part_cnt, clo, chi, plo, phi, dx, dt, plof, dxf, &
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

  double precision, intent(in   )         :: dx(3), dxf(3), dt, plo(3), phi(3), plof(3)

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
  
  integer :: i, j, k, p, cell_np, new_np, intside, intsurf, push, loopcount
  integer :: ni(3), fi(3)
  integer(c_int), pointer :: cell_parts(:)
  type(particle_t), pointer :: part
  type(surface_t), pointer :: surf
  real(amrex_real) dxinv(3), dxfinv(3), onemdxf(3), ixf(3), localvel(3), localbeta, bfac(3), deltap(3), std, normalrand(3), nodalp, tempvel(3), intold, inttime, runerr, runtime, adj, adjalt, domsize(3), posalt(3), propvec(3)

  double precision  :: cc(0:7)
  double precision  :: rr(0:7)

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

   
  domsize = phi - plo

  adj = 0.999999
  adjalt = 2d0*(1d0 - adj)

  dxinv = 1.d0/dx

  dxfinv = 1.d0/dxf
  onemdxf = 1.d0 - dxf
  
  do p = 1, ns

    surf => surfaces(p)  

    surf%fxleft = 0
    surf%fyleft = 0
    surf%fzleft = 0

    surf%fxright = 0
    surf%fyright = 0
    surf%fzright = 0

  enddo

  !print *, "Starting"        

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

              !if(part%id .eq. 5) then
              !  print *, "Outer loop start vel: ", part%vel, " Outer loop start pos: ", part%pos/phi
              !endif


                !find fluid cell
              fi(1) = floor((part%pos(1) - plof(1))*dxfinv(1))
              fi(2) = floor((part%pos(2) - plof(2))*dxfinv(2))
#if (BL_SPACEDIM == 3)
              fi(3) = floor((part%pos(3) - plof(3))*dxfinv(3))
#else
              fi(3) = 0
#endif
              !Interpolate fluid fields. ixf is the particle position in local cell coordinates. fi is the fluid cell
              ixf(1) = (part%pos(1) - coordsx(fi(1),fi(2),fi(3),1))*dxfInv(1)
              ixf(2) = (part%pos(2) - coordsy(fi(1),fi(2),fi(3),2))*dxfInv(2)
#if (BL_SPACEDIM == 3)
              ixf(3) = (part%pos(3) - coordsz(fi(1),fi(2),fi(3),3))*dxfInv(3)
#endif

              !dxf is the size of the fluid cell


              call get_interpolation_weights(cc, rr, ixf, onemdxf)

#if (BL_SPACEDIM == 3)
              call get_local_properties(cc, rr, fi, velx, velxlo, velxhi, vely, velylo, velyhi, velz, velzlo, velzhi, beta, localvel, localbeta, betalo, betahi)
#endif
#if (BL_SPACEDIM == 2)
              call get_local_properties(cc, rr, fi, velx, velxlo, velxhi, vely, velylo, velyhi, beta, localvel, localbeta, betalo, betahi)
#endif

                !Brownian forcing

              call get_particle_normal(normalrand(1))
              call get_particle_normal(normalrand(2))
              call get_particle_normal(normalrand(3))

              runerr = 1d0;

              deltap(1) = part%vel(1)
              deltap(2) = part%vel(1)
#if (BL_SPACEDIM == 3)
              deltap(3) = part%vel(1)
#endif

!this velocity update introduces a timetep error when a boundary intersection occurs. Need to use some kind of iteration to get a consistent intersection time and velocity update.

              std = sqrt(-part%drag_factor*localbeta*k_B*2d0*runtime*293d0)/part%mass

              bfac(1) = std*normalrand(1)
              bfac(2) = std*normalrand(2)
              bfac(3) = std*normalrand(3)

              !print *, "brownian: ", bfac, " propusive: ", part%dir*part%propulsion*runtime

                !print *, "Position 1: ", part%pos, " Vel 1: ", part%vel, "localvel: ", localbeta

              part%vel(1) = part%accel_factor*localbeta*(part%vel(1)-localvel(1))*runtime + bfac(1) + part%vel(1)
              part%vel(2) = part%accel_factor*localbeta*(part%vel(2)-localvel(2))*runtime + bfac(2) + part%vel(2)
#if (BL_SPACEDIM == 3)
              part%vel(3) = part%accel_factor*localbeta*(part%vel(3)-localvel(3))*runtime + bfac(3) + part%vel(3)
#endif
              call redirect(part)
  
              part%vel = part%dir*part%propulsion*runtime + part%vel

              deltap(1) = part%mass*(part%vel(1) - deltap(1))
              deltap(2) = part%mass*(part%vel(2) - deltap(2))
#if (BL_SPACEDIM == 3)
              deltap(3) = part%mass*(part%vel(3) - deltap(3))
#endif

#if (BL_SPACEDIM == 3)
              call distribute_momentum(deltap, rr, fi ,sourcex, sourcexlo, sourcexhi, sourcey, sourceylo, sourceyhi, sourcez, sourcezlo, sourcezhi)
#endif
#if (BL_SPACEDIM == 2)
              call distribute_momentum(deltap, rr, fi ,sourcex, sourcexlo, sourcexhi, sourcey, sourceylo, sourceyhi)
#endif

              do while (runtime .gt. 0)

                !if(part%id .eq. 5) then
                !  print *, "Starting vel: ", part%vel, " Starting pos: ", part%pos/phi
                !endif


                call find_intersect(part,runtime, surfaces, ns, intsurf, inttime, intside, phi, plo)

                posalt(1) = inttime*part%vel(1)*adjalt
                posalt(2) = inttime*part%vel(2)*adjalt
#if (BL_SPACEDIM == 3)
                posalt(3) = inttime*part%vel(3)*adjalt
#endif

                ! move the particle in a straight line, adj factor prevents double detection of boundary intersection
                part%pos(1) = part%pos(1) + inttime*part%vel(1)*adj
                part%pos(2) = part%pos(2) + inttime*part%vel(2)*adj
#if (BL_SPACEDIM == 3)
                part%pos(3) = part%pos(3) + inttime*part%vel(3)*adj
#endif

                runtime = runtime - inttime

                if(intsurf .gt. 0) then

                  surf => surfaces(intsurf)
                  !if(intsurf .eq. 7) then
                 !   print *, part%id, " prevel: ", part%vel
                  !endif
                  call apply_bc(surf, part, intside, domsize, push)
                  !if(intsurf .eq. 7) then
                 !   print *, part%id, " postvel: ", part%vel
                 ! endif


                    if(push .eq. 1) then
                      
                      part%pos(1) = part%pos(1) + posalt(1)
                      part%pos(2) = part%pos(2) + posalt(2)
#if (BL_SPACEDIM == 3)
                      part%pos(3) = part%pos(3) + posalt(3)
#endif
                    endif
                    
                endif

               ! if(part%id .eq. 5) then
               !   print *, "Final vel: ", part%vel, " Final pos: ", part%pos/phi
              !  endif

              enddo

              ! if it has changed cells, remove from vector.
              ! otherwise continue
              ni(1) = floor((part%pos(1) - plo(1))*dxinv(1))
              ni(2) = floor((part%pos(2) - plo(2))*dxinv(2))
#if (BL_SPACEDIM == 3)
              ni(3) = floor((part%pos(3) - plo(3))*dxinv(3))
#else
              ni(3) = 0
#endif

              if ((ni(1) /= i) .or. (ni(2) /= j) .or. (ni(3) /= k)) then
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

  !print *, "Ending"        
  
end subroutine move_particles_fhd

