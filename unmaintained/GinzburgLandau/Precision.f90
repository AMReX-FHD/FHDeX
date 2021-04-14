MODULE Precision
IMPLICIT NONE
PRIVATE
INTEGER * 1 :: INTEGER_var_1
INTEGER * 2 :: INTEGER_var_2
INTEGER * 4 :: INTEGER_var_4
INTEGER * 8 :: INTEGER_var_8
REAL * 4 :: REAL_var_4
REAL * 8 :: REAL_var_8
!REAL * 16 :: REAL_var_16
LOGICAL * 1 :: LOGICAL_var_1
LOGICAL * 2 :: LOGICAL_var_2
LOGICAL * 4 :: LOGICAL_var_4
LOGICAL * 8 :: LOGICAL_var_8
CHARACTER * 1 :: CHARACTER_var_1
INTEGER, PARAMETER, PUBLIC :: i_byte = KIND (INTEGER_var_1), i_short = KIND (INTEGER_var_2), i_sp = KIND (INTEGER_var_4), i_dp = &
& KIND (INTEGER_var_8), i_word = KIND (0)
INTEGER, PARAMETER, PUBLIC :: r_sp = KIND (0.0E0), r_dp = KIND (0.0D0), r_word = KIND (0.0)
INTEGER, PARAMETER, PUBLIC :: r_qp = KIND (REAL_var_8)
INTEGER, PARAMETER, PUBLIC :: l_byte = KIND (LOGICAL_var_1), l_short = KIND (LOGICAL_var_2), l_word = KIND (.TRUE.)
INTEGER, PARAMETER, PUBLIC :: c_byte = KIND (CHARACTER_var_1), c_ascii = KIND (' ')
INTEGER, PARAMETER, PUBLIC :: i_32 = KIND (INTEGER_var_4), i_64 = KIND (INTEGER_var_8)
INTEGER, PARAMETER, PUBLIC :: r_32 = KIND (REAL_var_4), r_64 = KIND (REAL_var_8)
INTEGER, PARAMETER, PUBLIC :: r_wp = r_dp, i_wp = i_word, l_wp = l_word
END MODULE Precision
