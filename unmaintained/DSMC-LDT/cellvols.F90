module cellvols_module
  use amrex_fort_module, only: amrex_real, amrex_particle_real
  use iso_c_binding ,    only: c_int

  implicit none

  public get_cell_vols

contains

  subroutine get_cell_vols(vols, vlo, vhi, dx, samples, plo)bind(c,name="get_cell_vols")

    use amrex_fort_module, only: amrex_real
    use rng_functions_module

    implicit none

    integer,          intent(in   )         :: vlo(3), vhi(3), samples
    double precision, intent(in   )         :: dx(3), plo(3)

    double precision, intent(inout)         :: vols(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

    integer :: i, j, k, l
    double precision :: pos(3), vcount, scount, rad

    do k = vlo(3), vhi(3)
      do j = vlo(2), vhi(2)
        do i = vlo(1), vhi(1)

          vcount = 0
          scount = 0
          do l = 1, samples

            pos(1) = plo(1) + i*dx(1) + get_uniform_func()*dx(1)
            pos(2) = plo(2) + j*dx(2) + get_uniform_func()*dx(2)
            !pos(3) = plo(3) + k*dx(3) + get_uniform_func()*dx(3)

            rad = sqrt(pos(1)**2 + pos(2)**2)

            if(rad .le. 2.5e-5) then
              vcount = vcount + 1
            endif
              scount = scount + 1

          enddo

          if((vcount/scount) .ne. 0) then
            vols(i,j,k) = (vcount/scount)*vols(i,j,k)

          endif

          !print *, i,j,k,vcount, samples, vols(i,j,k)
        enddo
      enddo
    enddo

  end subroutine get_cell_vols

end module cellvols_module
