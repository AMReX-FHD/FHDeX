subroutine repulsive_force(part1,part2,dx, dr2) &
    bind(c,name="repulsive_force")

  use amrex_fort_module, only: amrex_real
  use iso_c_binding, only: c_ptr, c_int, c_f_pointer
  use cell_sorted_particle_module, only: particle_t

  implicit none
  type(particle_t), intent(inout) :: part1 
  type(particle_t), intent(inout) :: part2
  real(amrex_real), intent(in) :: dx, dr2

  !here using a (1/r)^4 potential interation--make sure sign is correct
  
  part1%force = part1%force + dx*part1%q*part2%q/(dr2*dr2*dr2)
  part2%force = part2%force - dx*part1%q*part2%q/(dr2*dr2*dr2)

end subroutine

subroutine force_function(part1,part2,domsize) &
    bind(c,name="force_function")

  use amrex_fort_module, only: amrex_real
  use iso_c_binding, only: c_ptr, c_int, c_f_pointer
  use cell_sorted_particle_module, only: particle_t
  use common_namelist_module, only: diameter

  implicit none
  type(particle_t), intent(inout) :: part1 !is this defined correctly?
  type(particle_t), intent(inout) :: part2
  real(amrex_real), intent(in) :: domsize(3)

  integer :: i,j,k,images
  real(amrex_real) :: dx(3), dr, dr2, permittivity, cutoff


  !here calculate forces as a function of distance

  dx = part1%pos-part2%pos

  !can we do smarter vector maniputions below? also--probably better ways to do this other than logic with ghost cells

  do while (i <= 3)

      if(dx(i) .gt. domsize(i)*.5) then !correct for boxsize; no particles farther than L/2

          dx(i) = dx(i) - domsize(i)

      end if

      if(dx(i) .lt. -1*domsize(i)*.5) then !correct for boxsize; no particles farther than L/2

          dx(i) = dx(i) + domsize(i)

      end if

      i = i + 1

   end do

  permittivity = 1 !for now we are keeping this at one (think this is true for cgs units)

  dr2 = dot_product(dx,dx)
  dr = sqrt(dr2)

  cutoff = 4*part1%radius

  !repulsive interaction
  if (dr .lt. cutoff) then
 
    call repulsive_force(part1,part2,dx,dr2) 

  end if

  !electrostatic -- need to determine how many images we should be adding
  images = 5 !change this to an input
  do while (i <= images)
  
    !change dx, dy, dz, dr2 for each image
     dx = dx + i*domsize     
     dr2 = dot_product(dx,dx)

     part1%force = part1%force + permittivity*(dx/abs(dx))*part1%q*part2%q/dr2
     part2%force = part2%force - permittivity*(dx/abs(dx))*part1%q*part2%q/dr2

     i = i + 1

     !make sure the above is doing the same thing as below!
!     part1%force(1) = part1%force(1) + (dx/abs(dx))*part1%q*part2%q/dr2
!     part2%force(1) = part2%force(1) - (dx/abs(dx))*part1%q*part2%q/dr2
   
!     part1%force(2) = part1%force(2) + (dy/abs(dy))*part1%q*part2%q/dr2
!     part2%force(2) = part2%force(2) - (dy/abs(dy))*part1%q*part2%q/dr2
!   
!     part1%force(3) = part1%force(3) + (dz/abs(dz))*part1%q*part2%q/dr2
!     part2%force(3) = part2%force(3) - (dz/abs(dz))*part1%q*part2%q/dr2

  end do
    

end subroutine force_function

subroutine calculate_force(particles, np, lo, hi, &
     cell_part_ids, cell_part_cnt, clo, chi, plo, phi, dx) &
     bind(c,name="calculate_force")
  
  use amrex_fort_module, only: amrex_real
  use iso_c_binding, only: c_ptr, c_int, c_f_pointer
  use cell_sorted_particle_module, only: particle_t, remove_particle_from_cell
  
  implicit none
  !this is everything we pass in
  type(particle_t), intent(inout), target :: particles(np)
  integer(c_int), intent(in) :: np 
  integer(c_int), intent(in) :: lo(3), hi(3)
  integer(c_int), intent(in) :: clo(3), chi(3) 
  type(c_ptr), intent(inout) :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
  integer(c_int), intent(inout) :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
  real(amrex_real), intent(in) :: plo(3), phi(3), dx(3) 
  
  integer :: i, j, k, p, n
  type(particle_t), pointer :: part
  type(particle_t), pointer :: part2 !added this
  real(amrex_real) inv_dx(3), domsize(3)

  inv_dx = 1.d0/dx
  
  domsize = phi - plo
  
  !zero all the forces
  p = 1
  do while (p <= np)

     part => particles(p) !this defines one particle--we can access all the data by doing part%something

     part%force=0

     p = p + 1

  end do

  !calculate N^2 interactions
  p = 1

  do while (p < np)

    part => particles(p) !this defines one particle--we can access all the data by doing part%something

    n = p + 1

    do while (n <= np) 

       part2 => particles(n) !this defines one particle--we can access all the data by doing part%something

       call force_function(part,part2,domsize)

       n = n + 1

    end do

    p = p + 1

  end do 
  
end subroutine calculate_force

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

         !below is surface specific stuff
         !       call find_intersect(part,runtime, surfaces, ns, intsurf, inttime, intside, phi, plo)

          !      if(intsurf .eq. 6) then

           !       intcount = intcount + 1

            !    endif
                
                 
                !print *, "Parrticle ", p, " intersect ", inttime, intsurf 

                !print *, "Prepos ", part%pos, " prevel ", part%vel

                !call sleep(1)
                 !inttime is the timestep for making a move here
    !            posalt(1) = inttime*part%vel(1)*adjalt
    !            posalt(2) = inttime*part%vel(2)*adjalt
!#if (BL_SPACEDIM == 3)
!                posalt(3) = inttime*part%vel(3)*adjalt
!#endif

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

!subroutine redirect(part)

!    use cell_sorted_particle_module, only: particle_t

!    implicit none

!    type(particle_t), intent(inout) :: part

!    double precision speed

!    speed = sqrt(part%vel(1)**2 + part%vel(2)**2 + part%vel(3)**2)

!    if(speed .ne. 0) then
!      part%dir = part%vel/speed
!    endif

!end subroutine redirect

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
              !call redirect(part)
  
              !part%vel = part%dir*part%propulsion*runtime + part%vel

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


subroutine peskin_3pt(r,w)

  !This isn't three point! Fill in correct values later

  double precision, intent(in   ) :: r
  double precision, intent(inout) :: w

  double precision rr
  rr = r*r

  if(r .le. -2) then

    w = 0

  elseif(r .le. -1) then

    w = 0.125*(5 + 2*r - sqrt(-7 - 12*r - 4*rr))

  elseif(r .le. 0) then

    w = 0.125*(3 + 2*r + sqrt(1 - 4*r - 4*rr))

  elseif(r .le. 1) then

    w = 0.125*(3 - 2*r + sqrt(1 + 4*r - 4*rr))

  elseif(r .le. 2) then

    w = 0.125*(5 - 2*r - sqrt(-7 + 12*r - 4*rr))

  else

    w = 0

  endif


end subroutine peskin_3pt


subroutine peskin_4pt(r,w)

  double precision, intent(in   ) :: r
  double precision, intent(inout) :: w

  double precision rr 
  rr = r*r

  if(r .le. -2) then

    w = 0

  elseif(r .le. -1) then

    w = 0.125*(5 + 2*r - sqrt(-7 - 12*r - 4*rr))

  elseif(r .le. 0) then

    w = 0.125*(3 + 2*r + sqrt(1 - 4*r - 4*rr))

  elseif(r .le. 1) then

    w = 0.125*(3 - 2*r + sqrt(1 + 4*r - 4*rr))

  elseif(r .le. 2) then

    w = 0.125*(5 - 2*r - sqrt(-7 + 12*r - 4*rr))

  else

    w = 0

  endif


end subroutine peskin_4pt

subroutine bspline_6pt(r,w)

  !Doesn't work?

  double precision, intent(in   ) :: r
  double precision, intent(inout) :: w

  double precision r1, r2, r3, r4, r5
 
  r1 = abs(r)
  r2 = r1*r1
  r3 = r2*r1
  r4 = r3*r1
  r5 = r4*r1

  if(r1 .le. 1) then

    w = 0.55 - 0.5*r2 + 0.25*r4 - (1d0/12d0)*r5

  elseif(r1 .le. 2) then

    w = 0.425 + 0.625*r - 1.75*r2 + 1.25*r3 - 0.375*r4 + (1d0/24d0)*r5

  elseif(r1 .le. 3) then

    w = 2.025 - 3.375*r + 2.25*r2 - 0.75*r3 + 0.125*r4 - (1d0/120d0)*r5;

  else

    w = 0

  endif


end subroutine bspline_6pt

subroutine peskin_6pt(r,w)

  double precision, intent(in   ) :: r
  double precision, intent(inout) :: w

  double precision alpha, beta, gamm, K, R1, R2, R3, discr
  integer sgn
 
  K = 59d0/60 - sqrt(29d0)/20

  R1 = r - ceiling(r) + 1
  R2 = R1*R1
  R3 = R2*R1
  alpha = 28
  beta  = 9d0/4 - 1.5 * (K + R2) + (22./3-7*K)*R1 - 7./3*R3;
  gamm = 0.25 * ( 0.5*(161./36 - 59./6*K + 5*K*K)*R2 + 1./3*(-109./24 + 5*K)*R2*R2 + 5./18*R3*R3  );
  discr = beta*beta - 4 * alpha * gamm;

  if((3/2 - K) .gt. 0) then

    sgn = 1

  else

    sgn = -1

  endif

  if(r .le. -3) then

    w = 0

  elseif(r .le. -2) then

    w = 1./(2*alpha) * ( -beta + sgn * sqrt(discr) );

  elseif(r .le. -1) then

    w = -3./(2*alpha) * ( -beta + sgn * sqrt(discr) ) - 1./16 + 1./8*( K+(r+2)*(r+2) ) + 1./12*(3*K-1)*(r+2) + 1./12*(r+2)*(r+2)*(r+2); 

  elseif(r .le. 0) then

    w = 2./(2*alpha) * ( -beta + sgn * sqrt(discr) ) + 1./4 + 1./6*(4-3*K)*(r+1) - 1./6*(r+1)*(r+1)*(r+1);

  elseif(r .le. 1) then

    w = 2./(2*alpha) * ( -beta + sgn * sqrt(discr) ) + 5./8 - 1./4 * ( K+r*r );

  elseif(r .le. 2) then

    w = -3./(2*alpha) * ( -beta + sgn * sqrt(discr) ) + 1./4 - 1./6*(4-3*K)*(r-1) + 1./6*(r-1)*(r-1)*(r-1);

  elseif(r .le. 3) then

    w = 1./(2*alpha) * ( -beta + sgn * sqrt(discr) ) - 1./16 + 1./8*(K+(r-2)*(r-2)) - 1./12*(3*K-1)*(r-2) - 1./12*(r-2)*(r-2)*(r-2); 	

  else

    w = 0

  endif


end subroutine peskin_6pt


subroutine get_weights(dxf, dxfinv, weights, indicies, &
                              coordsu, coordsulo, coordsuhi, &
                              coordsv, coordsvlo, coordsvhi, &
#if (BL_SPACEDIM == 3)
                              coordsw, coordswlo, coordswhi, &
#endif
                              part, ks, lo, hi, plof)

  use amrex_fort_module, only: amrex_real
  use cell_sorted_particle_module, only: particle_t
  use common_namelist_module

  implicit none

  double precision, intent(in   ) :: dxf(3), dxfinv(3), plof(3)
  integer,          intent(in   ) :: ks, coordsulo(3), coordsvlo(3), coordswlo(3), coordsuhi(3), coordsvhi(3), coordswhi(3), lo(3), hi(3) 
  type(particle_t), intent(in   ) :: part
  double precision, intent(inout) :: weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3)
  integer         , intent(inout) :: indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3)

  double precision, intent(in   ) :: coordsu(coordsulo(1):coordsuhi(1),coordsulo(2):coordsuhi(2),coordsulo(3):coordsuhi(3),1:AMREX_SPACEDIM)
  double precision, intent(in   ) :: coordsv(coordsvlo(1):coordsvhi(1),coordsvlo(2):coordsvhi(2),coordsvlo(3):coordsvhi(3),1:AMREX_SPACEDIM)
#if (AMREX_SPACEDIM == 3)
  double precision, intent(in   ) :: coordsw(coordswlo(1):coordswhi(1),coordswlo(2):coordswhi(2),coordswlo(3):coordswhi(3),1:AMREX_SPACEDIM)
#endif

  integer :: fi(3), fn(3),i, j, k
  double precision :: xx,yy,zz, w1, w2, w3, fr(3), fd(3), wcheck(3)

  !find fluid cell

  fr(1) = (part%pos(1) - plof(1))*dxfinv(1)
  fr(2) = (part%pos(2) - plof(2))*dxfinv(2)
  fr(3) = (part%pos(3) - plof(3))*dxfinv(3)

  fi(1) = floor(fr(1))
  fi(2) = floor(fr(2))
  fi(3) = floor(fr(3))

  fd(1) = fr(1) - fi(1)
  fd(2) = fr(2) - fi(2)
  fd(3) = fr(3) - fi(3)

  if(fd(1) .lt. 0.5) then
    fn(1) = -1
  else
    fn(1) = 0
  endif

  if(fd(2) .lt. 0.5) then
    fn(2) = -1
  else
    fn(2) = 0
  endif  

  if(fd(3) .lt. 0.5) then
    fn(3) = -1
  else
    fn(3) = 0
  endif    

  wcheck = 0

  do k = -(ks-1), ks
    do j = -(ks-1), ks
      do i = -(ks-1), ks

        xx = part%pos(1) - coordsu(fi(1)+i,fi(2)+j+fn(2),fi(3)+k+fn(3),1)
        yy = part%pos(2) - coordsu(fi(1)+i,fi(2)+j+fn(2),fi(3)+k+fn(3),2)
        zz = part%pos(3) - coordsu(fi(1)+i,fi(2)+j+fn(2),fi(3)+k+fn(3),3)

        if(pkernel_fluid .eq. 4) then
          call peskin_4pt(xx*dxfinv(1),w1)
          call peskin_4pt(yy*dxfinv(2),w2)
          call peskin_4pt(zz*dxfinv(3),w3)
        elseif(pkernel_fluid .eq. 6) then
          call peskin_6pt(xx*dxfinv(1),w1)
          call peskin_6pt(yy*dxfinv(2),w2)
          call peskin_6pt(zz*dxfinv(3),w3)
        endif 

        weights(i,j,k,1) = w1*w2*w3

        indicies(i,j,k,1,1) = fi(1)+i
        indicies(i,j,k,1,2) = fi(2)+j+fn(2)
        indicies(i,j,k,1,3) = fi(3)+k+fn(3)

        wcheck(1) = wcheck(1) + weights(i,j,k,1)

        !print*, weights(i,j,k,1)

        xx = part%pos(1) - coordsv(fi(1)+i+fn(1),fi(2)+j,fi(3)+k+fn(3),1)
        yy = part%pos(2) - coordsv(fi(1)+i+fn(1),fi(2)+j,fi(3)+k+fn(3),2)
        zz = part%pos(3) - coordsv(fi(1)+i+fn(1),fi(2)+j,fi(3)+k+fn(3),3)

        if(pkernel_fluid .eq. 4) then
          call peskin_4pt(xx*dxfinv(1),w1)
          call peskin_4pt(yy*dxfinv(2),w2)
          call peskin_4pt(zz*dxfinv(3),w3)
        elseif(pkernel_fluid .eq. 6) then
          call peskin_6pt(xx*dxfinv(1),w1)
          call peskin_6pt(yy*dxfinv(2),w2)
          call peskin_6pt(zz*dxfinv(3),w3)
        endif 

        weights(i,j,k,2) = w1*w2*w3

        indicies(i,j,k,2,1) = fi(1)+i+fn(1)
        indicies(i,j,k,2,2) = fi(2)+j
        indicies(i,j,k,2,3) = fi(3)+k+fn(3)

        wcheck(2) = wcheck(2) + weights(i,j,k,2)


        xx = part%pos(1) - coordsw(fi(1)+i+fn(1),fi(2)+j+fn(2),fi(3)+k,1)
        yy = part%pos(2) - coordsw(fi(1)+i+fn(1),fi(2)+j+fn(2),fi(3)+k,2)
        zz = part%pos(3) - coordsw(fi(1)+i+fn(1),fi(2)+j+fn(2),fi(3)+k,3)

        if(pkernel_fluid .eq. 4) then
          call peskin_4pt(xx*dxfinv(1),w1)
          call peskin_4pt(yy*dxfinv(2),w2)
          call peskin_4pt(zz*dxfinv(3),w3)
        elseif(pkernel_fluid .eq. 6) then
          call peskin_6pt(xx*dxfinv(1),w1)
          call peskin_6pt(yy*dxfinv(2),w2)
          call peskin_6pt(zz*dxfinv(3),w3)
        endif 

        weights(i,j,k,3) = w1*w2*w3

        indicies(i,j,k,3,1) = fi(1)+i+fn(1)
        indicies(i,j,k,3,2) = fi(2)+j+fn(2)
        indicies(i,j,k,3,3) = fi(3)+k

        wcheck(3) = wcheck(3) + weights(i,j,k,3)

      enddo
    enddo
  enddo


 ! print*, "Total: ", wcheck

!              !Interpolate fluid fields. ixf is the particle position in local cell coordinates. fi is the fluid cell
!              ixf(1) = (part%pos(1) - coordsx(fi(1),fi(2),fi(3),1))*dxfInv(1)
!              ixf(2) = (part%pos(2) - coordsy(fi(1),fi(2),fi(3),2))*dxfInv(2)
!#if (BL_SPACEDIM == 3)
!              ixf(3) = (part%pos(3) - coordsz(fi(1),fi(2),fi(3),3))*dxfInv(3)
!#endif



end subroutine get_weights

subroutine get_weights_scalar_cc(dx, dxinv, weights, indicies, &
                              coords, coordslo, coordshi, &
                              part, ks, lo, hi, plof, store)

  use amrex_fort_module, only: amrex_real
  use cell_sorted_particle_module, only: particle_t
  use common_namelist_module

  implicit none

  double precision, intent(in   ) :: dx(3), dxinv(3), plof(3)
  integer,          intent(in   ) :: ks, coordslo(3), coordshi(3), lo(3), hi(3), store
  type(particle_t), intent(in   ) :: part
  double precision, intent(inout) :: weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3)
  integer         , intent(inout) :: indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3)

  double precision, intent(in   ) :: coords(coordslo(1):coordshi(1),coordslo(2):coordshi(2),coordslo(3):coordshi(3),1:AMREX_SPACEDIM)

  integer :: fi(3), fn(3),i, j, k
  double precision :: xx,yy,zz, w1, w2, w3, fr(3), fd(3), wcheck

  !find scalar cell

  fr(1) = (part%pos(1) - plof(1))*dxinv(1)
  fr(2) = (part%pos(2) - plof(2))*dxinv(2)
  fr(3) = (part%pos(3) - plof(3))*dxinv(3)

  fi(1) = floor(fr(1))
  fi(2) = floor(fr(2))
  fi(3) = floor(fr(3))

  fd = fr - fi

  if(fd(1) .lt. 0.5) then
    fn(1) = -1
  else
    fn(1) = 0
  endif

  if(fd(2) .lt. 0.5) then
    fn(2) = -1
  else
    fn(2) = 0
  endif  

  if(fd(3) .lt. 0.5) then
    fn(3) = -1
  else
    fn(3) = 0
  endif    

  wcheck = 0

  do k = -(ks-1), ks
    do j = -(ks-1), ks
      do i = -(ks-1), ks

        xx = part%pos(1) - coords(fi(1)+i+fn(1),fi(2)+j+fn(2),fi(3)+k+fn(3),store)
        yy = part%pos(2) - coords(fi(1)+i+fn(1),fi(2)+j+fn(2),fi(3)+k+fn(3),store)
        zz = part%pos(3) - coords(fi(1)+i+fn(1),fi(2)+j+fn(2),fi(3)+k+fn(3),store)

        if(pkernel_fluid .eq. 4) then
          call peskin_4pt(xx*dxinv(1),w1)
          call peskin_4pt(yy*dxinv(2),w2)
          call peskin_4pt(zz*dxinv(3),w3)
        elseif(pkernel_fluid .eq. 6) then
          call peskin_6pt(xx*dxinv(1),w1)
          call peskin_6pt(yy*dxinv(2),w2)
          call peskin_6pt(zz*dxinv(3),w3)
        endif 

        weights(i,j,k,1) = w1*w2*w3

        indicies(i,j,k,1,1) = fi(1)+i+fn(1)
        indicies(i,j,k,1,2) = fi(2)+j+fn(2)
        indicies(i,j,k,1,3) = fi(3)+k+fn(3)

        wcheck = wcheck + weights(i,j,k,store)

      enddo
    enddo
  enddo


  print*, "Total: ", wcheck

!              !Interpolate fluid fields. ixf is the particle position in local cell coordinates. fi is the fluid cell
!              ixf(1) = (part%pos(1) - coordsx(fi(1),fi(2),fi(3),1))*dxfInv(1)
!              ixf(2) = (part%pos(2) - coordsy(fi(1),fi(2),fi(3),2))*dxfInv(2)
!#if (BL_SPACEDIM == 3)
!              ixf(3) = (part%pos(3) - coordsz(fi(1),fi(2),fi(3),3))*dxfInv(3)
!#endif



end subroutine get_weights_scalar_cc

subroutine spread_op(weights, indicies, &
                              sourceu, sourceulo, sourceuhi, &
                              sourcev, sourcevlo, sourcevhi, &
#if (BL_SPACEDIM == 3)
                              sourcew, sourcewlo, sourcewhi, &
#endif
                              velu, velulo, veluhi, &
                              velv, velvlo, velvhi, &
#if (BL_SPACEDIM == 3)
                              velw, velwlo, velwhi, &
#endif
                              part, ks, dxf)

  use amrex_fort_module, only: amrex_real
  use cell_sorted_particle_module, only: particle_t

  implicit none

  integer,          intent(in   ) :: ks, sourceulo(3), sourcevlo(3), sourcewlo(3), sourceuhi(3), sourcevhi(3), sourcewhi(3), velulo(3), velvlo(3), velwlo(3), veluhi(3), velvhi(3), velwhi(3)
  double precision, intent(in   ) :: dxf(3)
  type(particle_t), intent(in   ) :: part
  double precision, intent(inout) :: weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3)
  integer         , intent(inout) :: indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3)

  double precision, intent(inout) :: sourceu(sourceulo(1):sourceuhi(1),sourceulo(2):sourceuhi(2),sourceulo(3):sourceuhi(3))
  double precision, intent(inout) :: sourcev(sourcevlo(1):sourcevhi(1),sourcevlo(2):sourcevhi(2),sourcevlo(3):sourcevhi(3))
#if (AMREX_SPACEDIM == 3)
  double precision, intent(inout) :: sourcew(sourcewlo(1):sourcewhi(1),sourcewlo(2):sourcewhi(2),sourcewlo(3):sourcewhi(3))
#endif

  double precision, intent(in   ) :: velu(velulo(1):veluhi(1),velulo(2):veluhi(2),velulo(3):veluhi(3))
  double precision, intent(in   ) :: velv(velvlo(1):velvhi(1),velvlo(2):velvhi(2),velvlo(3):velvhi(3))
#if (AMREX_SPACEDIM == 3)
  double precision, intent(in   ) :: velw(velwlo(1):velwhi(1),velwlo(2):velwhi(2),velwlo(3):velwhi(3))
#endif

  integer :: i, j, k, ii, jj, kk
  double precision :: uloc, vloc, wloc, volinv

  volinv = 1/(dxf(1)*dxf(2)*dxf(3))

  uloc = 0
  vloc = 0
  wloc = 0

  do k = -(ks-1), ks
    do j = -(ks-1), ks
      do i = -(ks-1), ks

        ii = indicies(i,j,k,1,1)
        jj = indicies(i,j,k,1,2)
        kk = indicies(i,j,k,1,3)

        uloc = uloc + velu(ii,jj,kk)*weights(i,j,k,1)


        ii = indicies(i,j,k,2,1)
        jj = indicies(i,j,k,2,2)
        kk = indicies(i,j,k,2,3)

        vloc = vloc + velv(ii,jj,kk)*weights(i,j,k,2)

        ii = indicies(i,j,k,3,1)
        jj = indicies(i,j,k,3,2)
        kk = indicies(i,j,k,3,3)

        wloc = wloc + velw(ii,jj,kk)*weights(i,j,k,3)

      enddo
    enddo
  enddo

  !print*, "Fluid vel: ", uloc, wloc, vloc

  do k = -(ks-1), ks
    do j = -(ks-1), ks
      do i = -(ks-1), ks

        ii = indicies(i,j,k,1,1)
        jj = indicies(i,j,k,1,2)
        kk = indicies(i,j,k,1,3)

        sourceu(ii,jj,kk) = (part%vel(1)-uloc)*(1d-2)*weights(i,j,k,1)*part%drag_factor*volinv+part%force(1)*weights(i,j,k,1)*volinv

        ii = indicies(i,j,k,2,1)
        jj = indicies(i,j,k,2,2)
        kk = indicies(i,j,k,2,3)

        sourcev(ii,jj,kk) = (part%vel(2)-vloc)*(1d-2)*weights(i,j,k,2)*part%drag_factor*volinv+part%force(2)*weights(i,j,k,2)*volinv

        ii = indicies(i,j,k,3,1)
        jj = indicies(i,j,k,3,2)
        kk = indicies(i,j,k,3,3)

        sourcew(ii,jj,kk) = (part%vel(3)-wloc)*(1d-2)*weights(i,j,k,3)*part%drag_factor*volinv+part%force(3)*weights(i,j,k,3)*volinv

      enddo
    enddo
  enddo


end subroutine spread_op

subroutine spread_op_scalar_cc(weights, indicies, &
                              source, sourcelo, sourcehi, &
                              part, ks, dx, mq, store)

  use amrex_fort_module, only: amrex_real
  use cell_sorted_particle_module, only: particle_t

  implicit none

  integer,          intent(in   ) :: ks, sourcelo(3), sourcehi(3), store
  double precision, intent(in   ) :: dx(3), mq
  type(particle_t), intent(in   ) :: part
  double precision, intent(inout) :: weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3)
  integer         , intent(inout) :: indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3)

  double precision, intent(inout) :: source(sourcelo(1):sourcehi(1),sourcelo(2):sourcehi(2),sourcelo(3):sourcehi(3))

  integer :: i, j, k, ii, jj, kk
  double precision :: volinv

  volinv = 1/(dx(1)*dx(2)*dx(3))


  !print*, "Fluid vel: ", uloc, wloc, vloc

  do k = -(ks-1), ks
    do j = -(ks-1), ks
      do i = -(ks-1), ks

        ii = indicies(i,j,k,1,1)
        jj = indicies(i,j,k,1,2)
        kk = indicies(i,j,k,1,3)

        source(ii,jj,kk) = mq*weights(i,j,k,store)

      enddo
    enddo
  enddo


end subroutine spread_op_scalar_cc

subroutine inter_op(weights, indicies, &
                              velu, velulo, veluhi, &
                              velv, velvlo, velvhi, &
#if (BL_SPACEDIM == 3)
                              velw, velwlo, velwhi, &
#endif
                              part, ks, dxf)

  use amrex_fort_module, only: amrex_real
  use cell_sorted_particle_module, only: particle_t

  implicit none

  integer,          intent(in   ) :: ks, velulo(3), velvlo(3), velwlo(3), veluhi(3), velvhi(3), velwhi(3)
  double precision, intent(in   ) :: dxf(3)
  type(particle_t), intent(inout) :: part
  double precision, intent(inout) :: weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3)
  integer         , intent(inout) :: indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3)

  double precision, intent(in   ) :: velu(velulo(1):veluhi(1),velulo(2):veluhi(2),velulo(3):veluhi(3))
  double precision, intent(in   ) :: velv(velvlo(1):velvhi(1),velvlo(2):velvhi(2),velvlo(3):velvhi(3))
#if (AMREX_SPACEDIM == 3)
  double precision, intent(in   ) :: velw(velwlo(1):velwhi(1),velwlo(2):velwhi(2),velwlo(3):velwhi(3))
#endif

  integer :: i, j, k, ii, jj, kk


  part%vel = 0

  do k = -(ks-1), ks
    do j = -(ks-1), ks
      do i = -(ks-1), ks

        ii = indicies(i,j,k,1,1)
        jj = indicies(i,j,k,1,2)
        kk = indicies(i,j,k,1,3)

        part%vel(1) = part%vel(1) + weights(i,j,k,1)*velu(ii,jj,kk)

        ii = indicies(i,j,k,2,1)
        jj = indicies(i,j,k,2,2)
        kk = indicies(i,j,k,2,3)

        part%vel(2) = part%vel(2) + weights(i,j,k,2)*velv(ii,jj,kk)

        ii = indicies(i,j,k,3,1)
        jj = indicies(i,j,k,3,2)
        kk = indicies(i,j,k,3,3)

        part%vel(3) = part%vel(3) + weights(i,j,k,3)*velw(ii,jj,kk)

      enddo
    enddo
  enddo

  part%vel(1) = part%vel(1)


  !print*, "Intervel: ", part%vel

  part%multi = part%vel(1)

end subroutine inter_op

subroutine move_ions_fhd(particles, np, lo, hi, &
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

                                     sourcex, sourcexlo, sourcexhi, &
                                     sourcey, sourceylo, sourceyhi, &
#if (BL_SPACEDIM == 3)
                                     sourcez, sourcezlo, sourcezhi, &
#endif
                                     surfaces, ns, sw)bind(c,name="move_ions_fhd")
  use amrex_fort_module, only: amrex_real
  use iso_c_binding, only: c_ptr, c_int, c_f_pointer
  use cell_sorted_particle_module, only: particle_t, remove_particle_from_cell
  use common_namelist_module, only: visc_type, k_B, pkernel_fluid
  use rng_functions_module
  use surfaces_module
  
  implicit none

  integer,          intent(in   )         :: np, ns, lo(3), hi(3), clo(3), chi(3), velxlo(3), velxhi(3), velylo(3), velyhi(3), sw
  integer,          intent(in   )         :: sourcexlo(3), sourcexhi(3), sourceylo(3), sourceyhi(3)
  integer,          intent(in   )         :: coordsxlo(3), coordsxhi(3), coordsylo(3), coordsyhi(3)
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

  double precision, intent(inout) :: sourcex(sourcexlo(1):sourcexhi(1),sourcexlo(2):sourcexhi(2),sourcexlo(3):sourcexhi(3))
  double precision, intent(inout) :: sourcey(sourceylo(1):sourceyhi(1),sourceylo(2):sourceyhi(2),sourceylo(3):sourceyhi(3))
#if (AMREX_SPACEDIM == 3)
  double precision, intent(inout) :: sourcez(sourcezlo(1):sourcezhi(1),sourcezlo(2):sourcezhi(2),sourcezlo(3):sourcezhi(3))
#endif

  type(c_ptr),      intent(inout) :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
  integer(c_int),   intent(inout) :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
  
  integer :: i, j, k, p, cell_np, new_np, intside, intsurf, push, loopcount, pointcount, ks
  integer :: ni(3), fi(3)
  integer(c_int), pointer :: cell_parts(:)
  type(particle_t), pointer :: part
  type(surface_t), pointer :: surf
  real(amrex_real) dxinv(3), dxfinv(3), onemdxf(3), ixf(3), localvel(3), localbeta, bfac(3), deltap(3), std, normalrand(3), nodalp, tempvel(3), intold, inttime, runerr, runtime, adj, adjalt, domsize(3), posalt(3), propvec(3), norm(3), diffest, diffav, distav, diffinst, veltest

  double precision, allocatable :: weights(:,:,:,:)
  integer, allocatable :: indicies(:,:,:,:,:)

  if(pkernel_fluid .eq. 4) then
    ks = 2
  else
    ks = 3
  endif
  
  allocate(weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3))
  allocate(indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3))

  domsize = phi - plo

  adj = 0.999999
  adjalt = 2d0*(1d0 - adj)

  dxinv = 1.d0/dx

  dxfinv = 1.d0/dxf
  onemdxf = 1.d0 - dxf
  

  !print *, "Starting"


  diffav = 0
  distav = 0
  diffinst = 0
  veltest = 0

  
  call calculate_force(particles, np, lo, hi, cell_part_ids, cell_part_cnt, clo, chi, plo, phi, dx)


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

              call get_weights(dxf, dxfinv, weights, indicies, &
                              coordsx, coordsxlo, coordsxhi, &
                              coordsy, coordsylo, coordsyhi, &
#if (BL_SPACEDIM == 3)
                              coordsz, coordszlo, coordszhi, &
#endif
                              part, ks, lo, hi, plof)
              if(sw .ne. 2) then
        
               !print*, "INTERPOLATE"
      
               call inter_op(weights, indicies, &
                                velx, velxlo, velxhi, &
                                vely, velylo, velyhi, &
#if (BL_SPACEDIM == 3)
                                velz, velzlo, velzhi, &
#endif
                                part, ks, dxf)

              endif

                !print*, "MOVE"
              !print*, "OldPos: ", part%pos
              !print*, "MoveVel: ", part%vel

              runtime = dt

              do while (runtime .gt. 0)

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

                  call apply_bc(surf, part, intside, domsize, push)

                    if(push .eq. 1) then
                      
                      part%pos(1) = part%pos(1) + posalt(1)
                      part%pos(2) = part%pos(2) + posalt(2)
#if (BL_SPACEDIM == 3)
                      part%pos(3) = part%pos(3) + posalt(3)
#endif
                    endif
                    
                endif

              end do

              part%abspos = part%abspos + dt*part%vel

              distav = distav + dt*sqrt(part%vel(1)**2+part%vel(2)**2+part%vel(3)**2)

              part%travel_time = part%travel_time + dt

              norm = part%abspos - part%origin

              diffest = (norm(1)**2 + norm(2)**2 + norm(3)**2)/(6*part%travel_time)

              diffinst = diffinst + diffest

              if(part%step_count .ge. 50) then
                part%diff_av = (part%diff_av*(part%step_count-50) + diffest)/((part%step_count-50) + 1)

                diffav = diffav + part%diff_av
              endif

              part%step_count = part%step_count + 1

              veltest = veltest + part%multi

              !print*, "Diff est: ", diffest , ", av: ", part%diff_av

              !print *, "AbsPos: ", part%abspos
              !print *, "RelPos: ", part%pos

              !print*, "NewPos: ", part%pos


              if(sw .ne. 1) then

              !  print*, "SPREAD"
                call spread_op(weights, indicies, &
                                sourcex, sourcexlo, sourcexhi, &
                                sourcey, sourceylo, sourceyhi, &
#if (BL_SPACEDIM == 3)
                                sourcez, sourcezlo, sourcezhi, &
#endif
                                velx, velxlo, velxhi, &
                                vely, velylo, velyhi, &
#if (BL_SPACEDIM == 3)
                                velz, velzlo, velzhi, &
#endif
                                part, ks, dxf)

              endif

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

  !print *, "Diffav: ", diffav/np, " Diffinst: ", diffinst/np, " Distav: ", distav/np
  !print *, "veltest: ", veltest/np

  deallocate(weights)
  deallocate(indicies)
  
end subroutine move_ions_fhd

subroutine collect_charge(particles, np, lo, hi, &
     cell_part_ids, cell_part_cnt, clo, chi, plo, phi, dx, dt, ploes, dxes, &
                                     cellcenters, cellcenterslo, cellcentershi, &
                                     charge, chargelo, chargehi)bind(c,name="collect_charge")

  use amrex_fort_module, only: amrex_real
  use iso_c_binding, only: c_ptr, c_int, c_f_pointer
  use cell_sorted_particle_module, only: particle_t, remove_particle_from_cell
  use common_namelist_module, only: visc_type, k_B, pkernel_es
  use rng_functions_module
  use surfaces_module
  
  implicit none

  integer,          intent(in   )         :: np, lo(3), hi(3), clo(3), chi(3)
  integer,          intent(in   )         :: chargelo(3), chargehi(3)
  integer,          intent(in   )         :: cellcenterslo(3), cellcentershi(3)

  type(particle_t), intent(inout), target :: particles(np)

  double precision, intent(in   )         :: dx(3), dxes(3), dt, plo(3), phi(3), ploes(3)

  double precision, intent(in   ) :: cellcenters(cellcenterslo(1):cellcentershi(1),cellcenterslo(2):cellcentershi(2),cellcenterslo(3):cellcentershi(3),1:AMREX_SPACEDIM)
  double precision, intent(inout) :: charge(chargelo(1):chargehi(1),chargelo(2):chargehi(2),chargelo(3):chargehi(3))

  type(c_ptr),      intent(inout) :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
  integer(c_int),   intent(inout) :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
  
  integer :: i, j, k, p, cell_np, new_np, intside, intsurf, push, loopcount, pointcount, ks, store
  integer :: ni(3), fi(3)
  integer(c_int), pointer :: cell_parts(:)
  type(particle_t), pointer :: part
  type(surface_t), pointer :: surf
  real(amrex_real) dxinv(3), dxesinv(3), onemdxf(3), ixf(3), localvel(3), localbeta, bfac(3), deltap(3), std, normalrand(3), nodalp, tempvel(3), intold, inttime, runerr, runtime, domsize(3), posalt(3), propvec(3), norm(3), diffest, diffav, distav, diffinst, qm

  double precision, allocatable :: weights(:,:,:,:)
  integer, allocatable :: indicies(:,:,:,:,:)

  if(pkernel_es .eq. 4) then
    ks = 2
  else
    ks = 3
  endif
  
  allocate(weights(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3))
  allocate(indicies(-(ks-1):ks,-(ks-1):ks,-(ks-1):ks,3,3))

  domsize = phi - plo

  dxinv = 1.d0/dx

  dxesinv = 1.d0/dxes


  diffav = 0
  distav = 0
  diffinst = 0

  store = 1
  qm = 1

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           cell_np = cell_part_cnt(i,j,k)
           call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])

           new_np = cell_np
           p = 1

           do while (p <= new_np)

!              runtime = dt
              part => particles(cell_parts(p))


              call get_weights_scalar_cc(dxes, dxesinv, weights, indicies, &
                              cellcenters, cellcenterslo, cellcentershi, &
                              part, ks, lo, hi, ploes,store)
              !if(sw .ne. 2) then
        
               !print*, "INTERPOLATE"
      
!               call inter_op(weights, indicies, &
!                                velx, velxlo, velxhi, &
!                                vely, velylo, velyhi, &
!#if (BL_SPACEDIM == 3)
!                                velz, velzlo, velzhi, &
!#endif
!                                part, ks, dxf)

              !endif

                !print*, "MOVE"
              !print*, "OldPos: ", part%pos
              !print*, "MoveVel: ", part%vel

!              runtime = dt

!              do while (runtime .gt. 0)

!                call find_intersect(part,runtime, surfaces, ns, intsurf, inttime, intside, phi, plo)

!                posalt(1) = inttime*part%vel(1)*adjalt
!                posalt(2) = inttime*part%vel(2)*adjalt
!#if (BL_SPACEDIM == 3)
!                posalt(3) = inttime*part%vel(3)*adjalt
!#endif

!                ! move the particle in a straight line, adj factor prevents double detection of boundary intersection
!                part%pos(1) = part%pos(1) + inttime*part%vel(1)*adj
!                part%pos(2) = part%pos(2) + inttime*part%vel(2)*adj
!#if (BL_SPACEDIM == 3)
!                part%pos(3) = part%pos(3) + inttime*part%vel(3)*adj
!#endif
!                runtime = runtime - inttime

!                if(intsurf .gt. 0) then

!                  surf => surfaces(intsurf)

!                  call apply_bc(surf, part, intside, domsize, push)

!                    if(push .eq. 1) then
!                      
!                      part%pos(1) = part%pos(1) + posalt(1)
!                      part%pos(2) = part%pos(2) + posalt(2)
!#if (BL_SPACEDIM == 3)
!                      part%pos(3) = part%pos(3) + posalt(3)
!#endif
!                    endif
!                    
!                endif

!              end do

!              part%abspos = part%abspos + dt*part%vel

!              distav = distav + dt*sqrt(part%vel(1)**2+part%vel(2)**2+part%vel(3)**2)

!              part%travel_time = part%travel_time + dt

!              norm = part%abspos - part%origin

!              diffest = (norm(1)**2 + norm(2)**2 + norm(3)**2)/(6*part%travel_time)

!              diffinst = diffinst + diffest

!              if(part%step_count .ge. 50) then
!                part%diff_av = (part%diff_av*(part%step_count-50) + diffest)/((part%step_count-50) + 1)

!                diffav = diffav + part%diff_av
!              endif

!              part%step_count = part%step_count + 1

              !print*, "Diff est: ", diffest , ", av: ", part%diff_av

              !print *, "AbsPos: ", part%abspos
              !print *, "RelPos: ", part%pos

              !print*, "NewPos: ", part%pos

              !  print*, "SPREAD"
                call spread_op_scalar_cc(weights, indicies, &
                                charge, chargelo, chargehi, &
                                part, ks, dxes, qm, store)



              ! if it has changed cells, remove from vector.
              ! otherwise continue
!              ni(1) = floor((part%pos(1) - plo(1))*dxinv(1))
!              ni(2) = floor((part%pos(2) - plo(2))*dxinv(2))
!#if (BL_SPACEDIM == 3)
!              ni(3) = floor((part%pos(3) - plo(3))*dxinv(3))
!#else
!              ni(3) = 0
!#endif

!              if ((ni(1) /= i) .or. (ni(2) /= j) .or. (ni(3) /= k)) then
!                 part%sorted = 0
!                 call remove_particle_from_cell(cell_parts, cell_np, new_np, p)  
!              else
!                 p = p + 1
!              end if
!           end do

!           cell_part_cnt(i,j,k) = new_np

            p = p + 1

          enddo
    
        end do
     end do
  end do

  !print *, "Diffav: ", diffav/np, " Diffinst: ", diffinst/np, " Distav: ", distav/np

  deallocate(weights)
  deallocate(indicies)
  
end subroutine collect_charge

