module stag_apply_op_module

  use amrex_error_module

  implicit none

  private

contains

#if (AMREX_SPACEDIM == 2)

  subroutine stag_apply_op(lo, hi, &
                           betacc, betacclo, betacchi, &
                           gammacc, gammacclo, gammacchi, &
                           betanodal, betanodallo, betanodalhi, &  
                           gammanodal, gammanodallo, gammanodalhi, &  
                           velxin, velxinlo, velxinhi, &
                           velyin, velyinlo, velyinhi, &
                           velxout, velxoutlo, velxouthi, &
                           velyout, velyoutlo, velyouthi, &
                           alpha, visctype) &
                           bind (C,name="stag_apply_op")

    integer         , intent(in   ) :: lo(2), hi(2), betacclo(2), betacchi(2), gammacclo(2), gammacchi(2), betanodallo(2), betanodalhi(2), gammanodallo(2), gammanodalhi(2)
    integer         , intent(in   ) :: velxinlo(2), velxinhi(2), velyinlo(2), velyinhi(2), velxoutlo(2), velxouthi(2), velyoutlo(2), velyouthi(2)
    integer         , intent(in   ) :: visctype
    double precision, intent(in   ) :: alpha
    double precision, intent(in   ) :: betacc(betacclo(1):betacchi(1),betacclo(2):betacchi(2))
    double precision, intent(in   ) :: gammacc(gammacclo(1):gammacchi(1),gammacclo(2):gammacchi(2))
    double precision, intent(in   ) :: betanodal(betanodallo(1):betanodalhi(1),betanodallo(2):betanodalhi(2))
    double precision, intent(in   ) :: gammanodal(gammanodallo(1):gammanodalhi(1),gammanodallo(2):gammanodalhi(2))

    double precision, intent(in   ) :: velxin(velxinlo(1):velxinhi(1),velxinlo(2):velxinhi(2))
    double precision, intent(in   ) :: velyin(velyinlo(1):velyinhi(1),velyinlo(2):velyinhi(2))
    double precision, intent(inout) :: velxout(velxoutlo(1):velxouthi(1),velxoutlo(2):velxouthi(2))
    double precision, intent(inout) :: velyout(velyoutlo(1):velyouthi(1),velyoutlo(2):velyouthi(2))

    ! local
    integer :: i,j,k
  
  end subroutine stag_apply_op


#endif

#if (AMREX_SPACEDIM == 3)

  subroutine stag_apply_op(lo, hi, &
                           betacc, betacclo, betacchi, &
                           gammacc, gammacclo, gammacchi, &
                           betanodal, betanodallo, betanodalhi, &  
                           gammanodal, gammanodallo, gammanodalhi, &  
                           velxin, velxinlo, velxinhi, &
                           velyin, velyinlo, velyinhi, &
                           velzin, velzinlo, velzinhi, &
                           velxout, velxoutlo, velxouthi, &
                           velyout, velyoutlo, velyouthi, &
                           velzout, velzoutlo, velzouthi, &
                           alpha, visctype) &
                           bind (C,name="stag_apply_op")

    integer         , intent(in   ) :: lo(3), hi(3), betaCClo(3), betaCChi(3), gammaCClo(3), gammaCChi(3), betaNodallo(3), betaNodalhi(3), gammaNodallo(3), gammaNodalhi(3)
    integer         , intent(in   ) :: velxin_lo(3), velxin_hi(3), velyin_lo(3), velyin_hi(3), velxout_lo(3), velxout_hi(3), velyout_lo(3), velyout_hi(3)
    integer         , intent(in   ) :: visctype
    double precision, intent(in   ) :: alpha

    double precision, intent(in   ) :: betacc(betacclo(1):betacchi(1),betacclo(2):betacchi(2),betacclo(3):betacchi(3))
    double precision, intent(in   ) :: gammacc(gammacclo(1):gammacchi(1),gammacclo(2):gammacchi(2),gammacclo(3):gammacchi(3))
    double precision, intent(in   ) :: betanodal(betanodallo(1):betanodalhi(1),betanodallo(2):betanodalhi(2),betanodallo(3):betanodalhi(3))
    double precision, intent(in   ) :: gammanodal(gammanodallo(1):gammanodalhi(1),gammanodallo(2):gammanodalhi(2),gammanodallo(3):gammanodalhi(3))

    double precision, intent(in   ) :: velxin(velxinlo(1):velxinhi(1),velxinlo(2):velxinhi(2),velxinlo(3):velxinhi(3))
    double precision, intent(in   ) :: velyin(velyinlo(1):velyinhi(1),velyinlo(2):velyinhi(2),velyinlo(3):velyinhi(3))
    double precision, intent(in   ) :: velzin(velzinlo(1):velzinhi(1),velzinlo(2):velzinhi(2),velzinlo(3):velzinhi(3))
    double precision, intent(inout) :: velxout(velxoutlo(1):velxouthi(1),velxoutlo(2):velxouthi(2),velxoutlo(3):velxouthi(3))
    double precision, intent(inout) :: velyout(velyoutlo(1):velyouthi(1),velyoutlo(2):velyouthi(2),velyoutlo(3):velyouthi(3))
    double precision, intent(inout) :: velzout(velzoutlo(1):velzouthi(1),velzoutlo(2):velzouthi(2),velzoutlo(3):velzouthi(3))

    ! local
    integer :: i,j,k
  
  end subroutine stag_apply_op


#endif

end module stag_apply_op_module
