module convert_stag_module

  use amrex_error_module

  implicit none

  private

contains

  subroutine average_cc_to_face(lo,hi, &
                                cc, c_lo, c_hi, nc_c, &
                                facex, x_lo, x_hi, nc_x, &
                                facey, y_lo, y_hi, nc_y, &
#if (AMREX_SPACEDIM == 3)
                                facez, z_lo, z_hi, nc_z, &
#endif
                                cc_comp, face_comp, ncomp) bind (C,name="average_cc_to_face")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: c_lo(3), c_hi(3), nc_c
    integer         , intent(in   ) :: x_lo(3), x_hi(3), nc_x
    integer         , intent(in   ) :: y_lo(3), y_hi(3), nc_y
#if (AMREX_SPACEDIM == 3)
    integer         , intent(in   ) :: z_lo(3), z_hi(3), nc_z
#endif
    double precision, intent(in   ) ::    cc(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3),0:nc_c-1)
    double precision, intent(inout) :: facex(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3),0:nc_x-1)
    double precision, intent(inout) :: facey(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),0:nc_y-1)
#if (AMREX_SPACEDIM == 3)
    double precision, intent(inout) :: facez(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3),0:nc_z-1)
#endif
    integer         , intent(in   ) :: cc_comp, face_comp, ncomp

    ! local
    integer :: i,j,k
    integer :: f_comp,c_comp

    do f_comp = face_comp, face_comp+ncomp-1
       c_comp = cc_comp + f_comp - face_comp

       ! x-faces
       do k=lo(3),hi(3)
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          facex(i,j,k,f_comp) = 0.5d0*(cc(i,j,k,c_comp)+cc(i-1,j,k,c_comp))
       end do
       end do
       end do


       ! y -faces
       do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          facey(i,j,k,f_comp) = 0.5d0*(cc(i,j,k,c_comp)+cc(i,j-1,k,c_comp))
       end do
       end do
       end do


#if (AMREX_SPACEDIM == 3)
       ! z-faces
       do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          facez(i,j,k,f_comp) = 0.5d0*(cc(i,j,k,c_comp)+cc(i,j,k-1,c_comp))
       end do
       end do
       end do
#endif

    end do

  end subroutine average_cc_to_face

  subroutine shift_face_to_cc(lo,hi, &
                                face, f_lo, f_hi, nc_f, &
                                cc, c_lo, c_hi, nc_c, &
                                face_comp, cc_comp, ncomp, av_dim) bind (C,name="shift_face_to_cc")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: f_lo(3), f_hi(3), nc_f
    integer         , intent(in   ) :: c_lo(3), c_hi(3), nc_c
    double precision, intent(in   ) :: face(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),0:nc_f-1)
    double precision, intent(inout) ::   cc(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3),0:nc_c-1)
    integer         , intent(in   ) :: face_comp, cc_comp, ncomp, av_dim

    ! local
    integer :: i,j,k
    integer :: f_comp,c_comp

    do f_comp = face_comp, face_comp+ncomp-1
       c_comp = cc_comp + f_comp - face_comp

       do k=lo(3),hi(3)
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          select case (av_dim)
          case (0)
             cc(i,j,k,c_comp) = face(i,j,k,f_comp)
          case (1)
             cc(i,j,k,c_comp) = face(i,j,k,f_comp)
          case (2)
             cc(i,j,k,c_comp) = face(i,j,k,f_comp)
          case default
             call amrex_error("invalid av_dim in shift_face_to_cc")
          end select
       end do
       end do
       end do

    end do

  end subroutine shift_face_to_cc

end module convert_stag_module
