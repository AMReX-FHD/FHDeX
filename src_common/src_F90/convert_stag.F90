module convert_stag_module

  use amrex_error_module

  implicit none

  private

contains

  subroutine average_face_to_cc(lo,hi, &
                                face, f_lo, f_hi, nc_f, &
                                cc, c_lo, c_hi, nc_c, &
                                face_comp, cc_comp, ncomp, av_dim) bind (C,name="average_face_to_cc")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: f_lo(3), f_hi(3), nc_f
    integer         , intent(in   ) :: c_lo(3), c_hi(3), nc_c
    double precision, intent(in   ) :: face(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),nc_f)
    double precision, intent(inout) ::   cc(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3),nc_c)
    integer         , intent(in   ) :: face_comp, cc_comp, ncomp, av_dim

    integer :: i,j,k
    integer :: f_comp,c_comp

    do f_comp=face_comp,face_comp+ncomp-1
       c_comp = cc_comp + f_comp - face_comp

       do k=lo(3),hi(3)
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          select case (av_dim)
          case (1)
             cc(i,j,k,c_comp) = 0.5d0*(face(i+1,j,k,f_comp)+face(i,j,k,f_comp))
          case (2)
             cc(i,j,k,c_comp) = 0.5d0*(face(i,j+1,k,f_comp)+face(i,j,k,f_comp))
          case (3)
             cc(i,j,k,c_comp) = 0.5d0*(face(i,j,k+1,f_comp)+face(i,j,k,f_comp))
          case default
             call amrex_error("invalid av_dim in average_face_to_cc")
          end select
       end do
       end do
       end do
   
    end do

  end subroutine average_face_to_cc

end module convert_stag_module
