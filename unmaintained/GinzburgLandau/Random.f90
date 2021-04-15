MODULE Random_Numbers
USE Precision
!USE MarsenneTwister
IMPLICIT NONE
PRIVATE
PUBLIC :: RandomUniform, RandomNormal, RandomSeeds, UnpredictableSeeds, RandomBits, SaveRandomSeeds, ResetRandomSeeds, &
& RestoreRandomSeeds
INTEGER, SAVE, PUBLIC :: RNG_algorithm = 0
INTEGER (KIND=i_32), DIMENSION (4), SAVE :: seeds32 = (/ 153587801_i_32, - 759022222_i_32, 1288503317_i_32, 1718083407_i_32 /)
INTEGER (KIND=i_64), DIMENSION (5), SAVE :: seeds64 = (/ 153587801_i_64, - 759022222_i_64, 1288503317_i_64, - 1718083407_i_64, - &
& 123456789_i_64 /)
INTEGER (KIND=i_32), PARAMETER :: MAX_i_32 = HUGE (1_i_32), MIN_i_32 = - HUGE (1_i_32) - 1
INTEGER (KIND=i_64), PARAMETER :: MAX_i_64 = HUGE (1_i_64), MIN_i_64 = - HUGE (1_i_64) - 1
REAL (KIND=r_sp), PARAMETER :: normalization_r_sp = - 1.0_r_sp / MIN_i_32
REAL (KIND=r_dp), PARAMETER :: normalization_r_dp = - 1.0_r_dp / MIN_i_32
TYPE, PUBLIC :: Random_Real_Seeds
INTEGER (KIND=i_32), DIMENSION (4) :: seeds = 0
END TYPE
TYPE (Random_Real_Seeds), ALLOCATABLE, DIMENSION (:), PUBLIC :: rng_seeds
INTERFACE RandomBits
MODULE PROCEDURE RandomBitScalar_i_sp
MODULE PROCEDURE RandomBitScalar_i_dp
MODULE PROCEDURE RandomBitArray_i_sp
MODULE PROCEDURE RandomBitArray_i_dp
END INTERFACE
INTERFACE RandomUniform
MODULE PROCEDURE UniformScalar_i_sp
MODULE PROCEDURE UniformScalar_i_dp
MODULE PROCEDURE UniformArrayContiguous_1_i_sp
MODULE PROCEDURE UniformArrayContiguous_1_i_dp
MODULE PROCEDURE UniformArrayLoops_1_i_sp
MODULE PROCEDURE UniformArrayLoops_1_i_dp
MODULE PROCEDURE UniformArrayContiguous_2_i_sp
MODULE PROCEDURE UniformArrayContiguous_2_i_dp
MODULE PROCEDURE UniformArrayLoops_2_i_sp
MODULE PROCEDURE UniformArrayLoops_2_i_dp
MODULE PROCEDURE UniformArrayContiguous_3_i_sp
MODULE PROCEDURE UniformArrayContiguous_3_i_dp
MODULE PROCEDURE UniformArrayLoops_3_i_sp
MODULE PROCEDURE UniformArrayLoops_3_i_dp
MODULE PROCEDURE UniformScalar_r_sp
MODULE PROCEDURE UniformScalar_r_dp
MODULE PROCEDURE UniformArrayContiguous_1_r_sp
MODULE PROCEDURE UniformArrayContiguous_1_r_dp
MODULE PROCEDURE UniformArrayLoops_1_r_sp
MODULE PROCEDURE UniformArrayLoops_1_r_dp
MODULE PROCEDURE UniformArrayContiguous_2_r_sp
MODULE PROCEDURE UniformArrayContiguous_2_r_dp
MODULE PROCEDURE UniformArrayLoops_2_r_sp
MODULE PROCEDURE UniformArrayLoops_2_r_dp
MODULE PROCEDURE UniformArrayContiguous_3_r_sp
MODULE PROCEDURE UniformArrayContiguous_3_r_dp
MODULE PROCEDURE UniformArrayLoops_3_r_sp
MODULE PROCEDURE UniformArrayLoops_3_r_dp
END INTERFACE
INTERFACE RandomNormal
MODULE PROCEDURE NormalScalar_r_sp
MODULE PROCEDURE NormalScalar_r_dp
MODULE PROCEDURE NormalArrayLoops_1_r_sp
MODULE PROCEDURE NormalArrayLoops_1_r_dp
MODULE PROCEDURE NormalArrayLoops_2_r_sp
MODULE PROCEDURE NormalArrayLoops_2_r_dp
MODULE PROCEDURE NormalArrayLoops_3_r_sp
MODULE PROCEDURE NormalArrayLoops_3_r_dp
END INTERFACE
CONTAINS
SUBROUTINE InitializeSeeds32 (i_seeds)
IMPLICIT NONE
INTEGER (KIND=i_32), DIMENSION (4), INTENT (IN) :: i_seeds
seeds32 = i_seeds
IF (IAND(seeds32(1),-2_i_32) == 0) seeds32 (1) = i_seeds (1) - 1023_i_32
IF (IAND(seeds32(2),-8_i_32) == 0) seeds32 (2) = i_seeds (2) - 1023_i_32
IF (IAND(seeds32(3),-16_i_32) == 0) seeds32 (3) = i_seeds (3) - 1023_i_32
IF (IAND(seeds32(4),-32_i_32) == 0) seeds32 (4) = i_seeds (4) - 1023_i_32
END SUBROUTINE InitializeSeeds32
SUBROUTINE InitializeSeeds64 (i_seeds)
IMPLICIT NONE
INTEGER (KIND=i_64), DIMENSION (5), INTENT (IN) :: i_seeds
seeds64 = i_seeds
IF (IAND(seeds64(1),-2_i_64) == 0) seeds64 (1) = i_seeds (1) - 8388607_i_64
IF (IAND(seeds64(2),-512_i_64) == 0) seeds64 (2) = i_seeds (2) - 8388607_i_64
IF (IAND(seeds64(3),-4096_i_64) == 0) seeds64 (3) = i_seeds (3) - 8388607_i_64
IF (IAND(seeds64(4),-131072_i_64) == 0) seeds64 (4) = i_seeds (4) - 8388607_i_64
IF (IAND(seeds64(5),-8388608_i_64) == 0) seeds64 (5) = i_seeds (5) - 8388607_i_64
END SUBROUTINE InitializeSeeds64
SUBROUTINE RandomSeeds (seed)
IMPLICIT NONE
INTEGER, INTENT (IN) :: seed
INTEGER :: i
INTEGER :: n_seeds
INTEGER, DIMENSION (:), ALLOCATABLE :: intrinsic_seeds
INTEGER, DIMENSION (5) :: i_seeds
REAL, DIMENSION (5) :: r_seeds
!IF (RNG_algorithm == 1) THEN
!CALL init_genrand (seed)
!ELSE
CALL RANDOM_SEED (SIZE=n_seeds)
ALLOCATE (intrinsic_seeds(n_seeds))
IF (n_seeds > 0) THEN
intrinsic_seeds (1) = Abs (seed)
DO i = 2, n_seeds
intrinsic_seeds (i) = Mod (8121*intrinsic_seeds(i-1)+28411, 134456)
END DO
!END IF
CALL RANDOM_SEED (PUT=intrinsic_seeds)
DEALLOCATE (intrinsic_seeds)
CALL RANDOM_NUMBER (r_seeds)
i_seeds = Int (HUGE(1_i_32)*r_seeds)
CALL InitializeSeeds32 (Int(i_seeds(1:4), i_32))
CALL InitializeSeeds64 (Int(i_seeds, i_64))
END IF
END SUBROUTINE RandomSeeds
SUBROUTINE UnpredictableSeeds (seed)
IMPLICIT NONE
INTEGER, INTENT (OUT), OPTIONAL :: seed
INTEGER :: times (8)
INTEGER :: unpredictable_seed
CALL DATE_AND_TIME (VALUES=times)
unpredictable_seed = times (6) * (times(2)*times(3)+times(5)) + times (7) * times (5) + times (8) * times (6)
CALL RandomSeeds (unpredictable_seed)
IF (PRESENT(seed)) seed = unpredictable_seed
END SUBROUTINE UnpredictableSeeds
SUBROUTINE UniformScalar_i_sp (number, range)
IMPLICIT NONE
INTEGER (KIND=i_sp), INTENT (OUT) :: number
INTEGER (KIND=i_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
INTEGER (KIND=i_32) :: temp_uniform
!IF (RNG_algorithm == 1) THEN
!CALL next_genrand (number)
!ELSE
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(1), 6), seeds32(1)),-13)
seeds32 (1) = IEOR (ISHFT(IAND(seeds32(1),-2_i_32), 18), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(2), 2), seeds32(2)),-27)
seeds32 (2) = IEOR (ISHFT(IAND(seeds32(2),-8_i_32), 2), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(3), 13), seeds32(3)),-21)
seeds32 (3) = IEOR (ISHFT(IAND(seeds32(3),-16_i_32), 7), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(4), 3), seeds32(4)),-12)
seeds32 (4) = IEOR (ISHFT(IAND(seeds32(4),-128_i_32), 13), temp_uniform)
number = IEOR (IEOR(IEOR(seeds32(1), seeds32(2)), seeds32(3)), seeds32(4))
!END IF
IF (PRESENT(range)) THEN
IF (.FALSE.) THEN
number = Abs (number) / (HUGE(1_i_32)/(range(2)-range(1)+1)) + range (1)
ELSE
number = Int (0.5_r_sp*(1.0_r_sp+REAL(number, r_sp)*normalization_r_sp)*(range(2)-range(1)+1), i_sp) + range (1)
END IF
number = Max (range(1), Min(number, range(2)))
END IF
END SUBROUTINE
SUBROUTINE UniformScalar_i_dp (number, range)
IMPLICIT NONE
INTEGER (KIND=i_dp), INTENT (OUT) :: number
INTEGER (KIND=i_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
INTEGER (KIND=i_64) :: temp_uniform
temp_uniform = ISHFT (IEOR(ISHFT(seeds64(1), 1), seeds64(1)),-53)
seeds64 (1) = IEOR (ISHFT(IAND(seeds64(1),-2_i_64), 10), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds64(2), 24), seeds64(2)),-50)
seeds64 (2) = IEOR (ISHFT(IAND(seeds64(2),-512_i_64), 5), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds64(3), 3), seeds64(3)),-23)
seeds64 (3) = IEOR (ISHFT(IAND(seeds64(3),-4096_i_64), 29), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds64(4), 5), seeds64(4)),-24)
seeds64 (4) = IEOR (ISHFT(IAND(seeds64(4),-131072_i_64), 23), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds64(5), 3), seeds64(4)),-33)
seeds64 (5) = IEOR (ISHFT(IAND(seeds64(5),-8388608_i_64), 8), temp_uniform)
number = IEOR (IEOR(IEOR(IEOR(seeds64(1), seeds64(2)), seeds64(3)), seeds64(4)), seeds64(5))
IF (PRESENT(range)) THEN
IF (.FALSE.) THEN
number = Abs (number) / (HUGE(1_i_64)/(range(2)-range(1)+1)) + range (1)
ELSE
number = Int (0.5_r_sp*(1.0_r_sp+REAL(number, r_sp)*normalization_r_sp)*(range(2)-range(1)+1), i_dp) + range (1)
END IF
number = Max (range(1), Min(number, range(2)))
END IF
END SUBROUTINE
SUBROUTINE RandomBitScalar_i_sp (number)
IMPLICIT NONE
INTEGER (KIND=i_sp), INTENT (OUT) :: number
INTEGER (KIND=i_32) :: temp_uniform
!IF (RNG_algorithm == 1) THEN
!CALL next_genrand (number)
!ELSE
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(1), 6), seeds32(1)),-13)
seeds32 (1) = IEOR (ISHFT(IAND(seeds32(1),-2_i_32), 18), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(2), 2), seeds32(2)),-27)
seeds32 (2) = IEOR (ISHFT(IAND(seeds32(2),-8_i_32), 2), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(3), 13), seeds32(3)),-21)
seeds32 (3) = IEOR (ISHFT(IAND(seeds32(3),-16_i_32), 7), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(4), 3), seeds32(4)),-12)
seeds32 (4) = IEOR (ISHFT(IAND(seeds32(4),-128_i_32), 13), temp_uniform)
number = IEOR (IEOR(IEOR(seeds32(1), seeds32(2)), seeds32(3)), seeds32(4))
!END IF
END SUBROUTINE
SUBROUTINE RandomBitScalar_i_dp (number)
IMPLICIT NONE
INTEGER (KIND=i_dp), INTENT (OUT) :: number
INTEGER (KIND=i_64) :: temp_uniform
temp_uniform = ISHFT (IEOR(ISHFT(seeds64(1), 1), seeds64(1)),-53)
seeds64 (1) = IEOR (ISHFT(IAND(seeds64(1),-2_i_64), 10), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds64(2), 24), seeds64(2)),-50)
seeds64 (2) = IEOR (ISHFT(IAND(seeds64(2),-512_i_64), 5), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds64(3), 3), seeds64(3)),-23)
seeds64 (3) = IEOR (ISHFT(IAND(seeds64(3),-4096_i_64), 29), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds64(4), 5), seeds64(4)),-24)
seeds64 (4) = IEOR (ISHFT(IAND(seeds64(4),-131072_i_64), 23), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds64(5), 3), seeds64(4)),-33)
seeds64 (5) = IEOR (ISHFT(IAND(seeds64(5),-8388608_i_64), 8), temp_uniform)
number = IEOR (IEOR(IEOR(IEOR(seeds64(1), seeds64(2)), seeds64(3)), seeds64(4)), seeds64(5))
END SUBROUTINE
SUBROUTINE RandomBitArray_i_sp (array, n_elements, method)
IMPLICIT NONE
INTEGER (KIND=i_wp), INTENT (IN) :: n_elements
INTEGER (KIND=i_sp), DIMENSION (n_elements), INTENT (OUT) :: array
INTEGER, INTENT (IN), OPTIONAL :: method
INTEGER (KIND=i_wp) :: element
INTEGER (KIND=i_32) :: temp_uniform
DO element = 1, n_elements
!IF (RNG_algorithm == 1) THEN
!CALL next_genrand (array(element))
!ELSE
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(1), 6), seeds32(1)),-13)
seeds32 (1) = IEOR (ISHFT(IAND(seeds32(1),-2_i_32), 18), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(2), 2), seeds32(2)),-27)
seeds32 (2) = IEOR (ISHFT(IAND(seeds32(2),-8_i_32), 2), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(3), 13), seeds32(3)),-21)
seeds32 (3) = IEOR (ISHFT(IAND(seeds32(3),-16_i_32), 7), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(4), 3), seeds32(4)),-12)
seeds32 (4) = IEOR (ISHFT(IAND(seeds32(4),-128_i_32), 13), temp_uniform)
array (element) = IEOR (IEOR(IEOR(seeds32(1), seeds32(2)), seeds32(3)), seeds32(4))
!END IF
END DO
END SUBROUTINE
SUBROUTINE RandomBitArray_i_dp (array, n_elements, method)
IMPLICIT NONE
INTEGER (KIND=i_wp), INTENT (IN) :: n_elements
INTEGER (KIND=i_dp), DIMENSION (n_elements), INTENT (OUT) :: array
INTEGER, INTENT (IN), OPTIONAL :: method
INTEGER (KIND=i_wp) :: element
INTEGER (KIND=i_64) :: temp_uniform
DO element = 1, n_elements
temp_uniform = ISHFT (IEOR(ISHFT(seeds64(1), 1), seeds64(1)),-53)
seeds64 (1) = IEOR (ISHFT(IAND(seeds64(1),-2_i_64), 10), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds64(2), 24), seeds64(2)),-50)
seeds64 (2) = IEOR (ISHFT(IAND(seeds64(2),-512_i_64), 5), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds64(3), 3), seeds64(3)),-23)
seeds64 (3) = IEOR (ISHFT(IAND(seeds64(3),-4096_i_64), 29), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds64(4), 5), seeds64(4)),-24)
seeds64 (4) = IEOR (ISHFT(IAND(seeds64(4),-131072_i_64), 23), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds64(5), 3), seeds64(4)),-33)
seeds64 (5) = IEOR (ISHFT(IAND(seeds64(5),-8388608_i_64), 8), temp_uniform)
array (element) = IEOR (IEOR(IEOR(IEOR(seeds64(1), seeds64(2)), seeds64(3)), seeds64(4)), seeds64(5))
END DO
END SUBROUTINE
SUBROUTINE UniformArrayContiguous_1_i_sp (array, n_extent1, range)
IMPLICIT NONE
INTEGER (KIND=i_wp), INTENT (IN) :: n_extent1
INTEGER (KIND=i_sp), DIMENSION (n_extent1), INTENT (OUT) :: array
INTEGER (KIND=i_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
CALL RandomBitArray_i_sp (array, Int(SIZE(array), KIND=i_wp))
IF (PRESENT(range)) THEN
IF (.FALSE.) THEN
array = Abs (array) / (HUGE(1_i_32)/(range(2)-range(1)+1)) + range (1)
ELSE
array = Int (0.5_r_sp*(1.0_r_sp+REAL(array, r_sp)*normalization_r_sp)*(range(2)-range(1)+1), i_sp) + range (1)
END IF
array = Max (range(1), Min(array, range(2)))
END IF
END SUBROUTINE
SUBROUTINE UniformArrayContiguous_1_i_dp (array, n_extent1, range)
IMPLICIT NONE
INTEGER (KIND=i_wp), INTENT (IN) :: n_extent1
INTEGER (KIND=i_dp), DIMENSION (n_extent1), INTENT (OUT) :: array
INTEGER (KIND=i_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
CALL RandomBitArray_i_dp (array, Int(SIZE(array), KIND=i_wp))
IF (PRESENT(range)) THEN
IF (.FALSE.) THEN
array = Abs (array) / (HUGE(1_i_64)/(range(2)-range(1)+1)) + range (1)
ELSE
array = Int (0.5_r_sp*(1.0_r_sp+REAL(array, r_sp)*normalization_r_sp)*(range(2)-range(1)+1), i_dp) + range (1)
END IF
array = Max (range(1), Min(array, range(2)))
END IF
END SUBROUTINE
SUBROUTINE UniformArrayContiguous_2_i_sp (array, n_extent1, n_extent2, range)
IMPLICIT NONE
INTEGER (KIND=i_wp), INTENT (IN) :: n_extent1, n_extent2
INTEGER (KIND=i_sp), DIMENSION (n_extent1, n_extent2), INTENT (OUT) :: array
INTEGER (KIND=i_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
CALL RandomBitArray_i_sp (array, Int(SIZE(array), KIND=i_wp))
IF (PRESENT(range)) THEN
IF (.FALSE.) THEN
array = Abs (array) / (HUGE(1_i_32)/(range(2)-range(1)+1)) + range (1)
ELSE
array = Int (0.5_r_sp*(1.0_r_sp+REAL(array, r_sp)*normalization_r_sp)*(range(2)-range(1)+1), i_sp) + range (1)
END IF
array = Max (range(1), Min(array, range(2)))
END IF
END SUBROUTINE
SUBROUTINE UniformArrayContiguous_2_i_dp (array, n_extent1, n_extent2, range)
IMPLICIT NONE
INTEGER (KIND=i_wp), INTENT (IN) :: n_extent1, n_extent2
INTEGER (KIND=i_dp), DIMENSION (n_extent1, n_extent2), INTENT (OUT) :: array
INTEGER (KIND=i_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
CALL RandomBitArray_i_dp (array, Int(SIZE(array), KIND=i_wp))
IF (PRESENT(range)) THEN
IF (.FALSE.) THEN
array = Abs (array) / (HUGE(1_i_64)/(range(2)-range(1)+1)) + range (1)
ELSE
array = Int (0.5_r_sp*(1.0_r_sp+REAL(array, r_sp)*normalization_r_sp)*(range(2)-range(1)+1), i_dp) + range (1)
END IF
array = Max (range(1), Min(array, range(2)))
END IF
END SUBROUTINE
SUBROUTINE UniformArrayContiguous_3_i_sp (array, n_extent1, n_extent2, n_extent3, range)
IMPLICIT NONE
INTEGER (KIND=i_wp), INTENT (IN) :: n_extent1, n_extent2, n_extent3
INTEGER (KIND=i_sp), DIMENSION (n_extent1, n_extent2, n_extent3), INTENT (OUT) :: array
INTEGER (KIND=i_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
CALL RandomBitArray_i_sp (array, Int(SIZE(array), KIND=i_wp))
IF (PRESENT(range)) THEN
IF (.FALSE.) THEN
array = Abs (array) / (HUGE(1_i_32)/(range(2)-range(1)+1)) + range (1)
ELSE
array = Int (0.5_r_sp*(1.0_r_sp+REAL(array, r_sp)*normalization_r_sp)*(range(2)-range(1)+1), i_sp) + range (1)
END IF
array = Max (range(1), Min(array, range(2)))
END IF
END SUBROUTINE
SUBROUTINE UniformArrayContiguous_3_i_dp (array, n_extent1, n_extent2, n_extent3, range)
IMPLICIT NONE
INTEGER (KIND=i_wp), INTENT (IN) :: n_extent1, n_extent2, n_extent3
INTEGER (KIND=i_dp), DIMENSION (n_extent1, n_extent2, n_extent3), INTENT (OUT) :: array
INTEGER (KIND=i_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
CALL RandomBitArray_i_dp (array, Int(SIZE(array), KIND=i_wp))
IF (PRESENT(range)) THEN
IF (.FALSE.) THEN
array = Abs (array) / (HUGE(1_i_64)/(range(2)-range(1)+1)) + range (1)
ELSE
array = Int (0.5_r_sp*(1.0_r_sp+REAL(array, r_sp)*normalization_r_sp)*(range(2)-range(1)+1), i_dp) + range (1)
END IF
array = Max (range(1), Min(array, range(2)))
END IF
END SUBROUTINE
SUBROUTINE UniformArrayLoops_1_i_sp (array, range)
IMPLICIT NONE
INTEGER (KIND=i_sp), DIMENSION (:), INTENT (OUT) :: array
INTEGER (KIND=i_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
INTEGER (KIND=i_wp) :: i1
DO i1 = Int (LBOUND(array, DIM=1), KIND=i_wp), Int (UBOUND(array, DIM=1), KIND=i_wp)
CALL RandomBitScalar_i_sp (array(i1))
END DO
IF (PRESENT(range)) THEN
IF (.FALSE.) THEN
array = Abs (array) / (HUGE(1_i_32)/(range(2)-range(1)+1)) + range (1)
ELSE
array = Int (0.5_r_sp*(1.0_r_sp+REAL(array, r_sp)*normalization_r_sp)*(range(2)-range(1)+1), i_sp) + range (1)
END IF
array = Max (range(1), Min(array, range(2)))
END IF
END SUBROUTINE
SUBROUTINE UniformArrayLoops_1_i_dp (array, range)
IMPLICIT NONE
INTEGER (KIND=i_dp), DIMENSION (:), INTENT (OUT) :: array
INTEGER (KIND=i_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
INTEGER (KIND=i_wp) :: i1
DO i1 = Int (LBOUND(array, DIM=1), KIND=i_wp), Int (UBOUND(array, DIM=1), KIND=i_wp)
CALL RandomBitScalar_i_dp (array(i1))
END DO
IF (PRESENT(range)) THEN
IF (.FALSE.) THEN
array = Abs (array) / (HUGE(1_i_64)/(range(2)-range(1)+1)) + range (1)
ELSE
array = Int (0.5_r_sp*(1.0_r_sp+REAL(array, r_sp)*normalization_r_sp)*(range(2)-range(1)+1), i_dp) + range (1)
END IF
array = Max (range(1), Min(array, range(2)))
END IF
END SUBROUTINE
SUBROUTINE UniformArrayLoops_2_i_sp (array, range)
IMPLICIT NONE
INTEGER (KIND=i_sp), DIMENSION (:, :), INTENT (OUT) :: array
INTEGER (KIND=i_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
INTEGER (KIND=i_wp) :: i1, i2
DO i2 = Int (LBOUND(array, DIM=2), KIND=i_wp), Int (UBOUND(array, DIM=2), KIND=i_wp)
DO i1 = Int (LBOUND(array, DIM=1), KIND=i_wp), Int (UBOUND(array, DIM=1), KIND=i_wp)
CALL RandomBitScalar_i_sp (array(i1, i2))
END DO
END DO
IF (PRESENT(range)) THEN
IF (.FALSE.) THEN
array = Abs (array) / (HUGE(1_i_32)/(range(2)-range(1)+1)) + range (1)
ELSE
array = Int (0.5_r_sp*(1.0_r_sp+REAL(array, r_sp)*normalization_r_sp)*(range(2)-range(1)+1), i_sp) + range (1)
END IF
array = Max (range(1), Min(array, range(2)))
END IF
END SUBROUTINE
SUBROUTINE UniformArrayLoops_2_i_dp (array, range)
IMPLICIT NONE
INTEGER (KIND=i_dp), DIMENSION (:, :), INTENT (OUT) :: array
INTEGER (KIND=i_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
INTEGER (KIND=i_wp) :: i1, i2
DO i2 = Int (LBOUND(array, DIM=2), KIND=i_wp), Int (UBOUND(array, DIM=2), KIND=i_wp)
DO i1 = Int (LBOUND(array, DIM=1), KIND=i_wp), Int (UBOUND(array, DIM=1), KIND=i_wp)
CALL RandomBitScalar_i_dp (array(i1, i2))
END DO
END DO
IF (PRESENT(range)) THEN
IF (.FALSE.) THEN
array = Abs (array) / (HUGE(1_i_64)/(range(2)-range(1)+1)) + range (1)
ELSE
array = Int (0.5_r_sp*(1.0_r_sp+REAL(array, r_sp)*normalization_r_sp)*(range(2)-range(1)+1), i_dp) + range (1)
END IF
array = Max (range(1), Min(array, range(2)))
END IF
END SUBROUTINE
SUBROUTINE UniformArrayLoops_3_i_sp (array, range)
IMPLICIT NONE
INTEGER (KIND=i_sp), DIMENSION (:, :, :), INTENT (OUT) :: array
INTEGER (KIND=i_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
INTEGER (KIND=i_wp) :: i1, i2, i3
DO i3 = Int (LBOUND(array, DIM=3), KIND=i_wp), Int (UBOUND(array, DIM=3), KIND=i_wp)
DO i2 = Int (LBOUND(array, DIM=2), KIND=i_wp), Int (UBOUND(array, DIM=2), KIND=i_wp)
DO i1 = Int (LBOUND(array, DIM=1), KIND=i_wp), Int (UBOUND(array, DIM=1), KIND=i_wp)
CALL RandomBitScalar_i_sp (array(i1, i2, i3))
END DO
END DO
END DO
IF (PRESENT(range)) THEN
IF (.FALSE.) THEN
array = Abs (array) / (HUGE(1_i_32)/(range(2)-range(1)+1)) + range (1)
ELSE
array = Int (0.5_r_sp*(1.0_r_sp+REAL(array, r_sp)*normalization_r_sp)*(range(2)-range(1)+1), i_sp) + range (1)
END IF
array = Max (range(1), Min(array, range(2)))
END IF
END SUBROUTINE
SUBROUTINE UniformArrayLoops_3_i_dp (array, range)
IMPLICIT NONE
INTEGER (KIND=i_dp), DIMENSION (:, :, :), INTENT (OUT) :: array
INTEGER (KIND=i_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
INTEGER (KIND=i_wp) :: i1, i2, i3
DO i3 = Int (LBOUND(array, DIM=3), KIND=i_wp), Int (UBOUND(array, DIM=3), KIND=i_wp)
DO i2 = Int (LBOUND(array, DIM=2), KIND=i_wp), Int (UBOUND(array, DIM=2), KIND=i_wp)
DO i1 = Int (LBOUND(array, DIM=1), KIND=i_wp), Int (UBOUND(array, DIM=1), KIND=i_wp)
CALL RandomBitScalar_i_dp (array(i1, i2, i3))
END DO
END DO
END DO
IF (PRESENT(range)) THEN
IF (.FALSE.) THEN
array = Abs (array) / (HUGE(1_i_64)/(range(2)-range(1)+1)) + range (1)
ELSE
array = Int (0.5_r_sp*(1.0_r_sp+REAL(array, r_sp)*normalization_r_sp)*(range(2)-range(1)+1), i_dp) + range (1)
END IF
array = Max (range(1), Min(array, range(2)))
END IF
END SUBROUTINE
SUBROUTINE UniformScalar_r_sp (number, range)
IMPLICIT NONE
REAL (KIND=r_sp), INTENT (OUT) :: number
REAL (KIND=r_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
INTEGER (KIND=i_32) :: temp_uniform
INTEGER (KIND=i_32) :: i_number
!IF (RNG_algorithm == 1) THEN
!CALL next_genrand (i_number)
!ELSE
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(1), 6), seeds32(1)),-13)
seeds32 (1) = IEOR (ISHFT(IAND(seeds32(1),-2_i_32), 18), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(2), 2), seeds32(2)),-27)
seeds32 (2) = IEOR (ISHFT(IAND(seeds32(2),-8_i_32), 2), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(3), 13), seeds32(3)),-21)
seeds32 (3) = IEOR (ISHFT(IAND(seeds32(3),-16_i_32), 7), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(4), 3), seeds32(4)),-12)
seeds32 (4) = IEOR (ISHFT(IAND(seeds32(4),-128_i_32), 13), temp_uniform)
i_number = IEOR (IEOR(IEOR(seeds32(1), seeds32(2)), seeds32(3)), seeds32(4))
!END IF
IF (PRESENT(range)) THEN
number = 0.5_r_sp * ((range(2)-range(1))*REAL(i_number, r_sp)*normalization_r_sp+(range(1)+range(2)))
ELSE
number = 0.5_r_sp * (REAL(i_number, r_sp)*normalization_r_sp+1.0_r_sp)
END IF
IF (.FALSE.) number = Max (number, TINY(number))
END SUBROUTINE
SUBROUTINE UniformScalar_r_dp (number, range)
IMPLICIT NONE
REAL (KIND=r_dp), INTENT (OUT) :: number
REAL (KIND=r_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
INTEGER (KIND=i_32) :: temp_uniform
INTEGER (KIND=i_32) :: i_number
!IF (RNG_algorithm == 1) THEN
!CALL next_genrand (i_number)
!ELSE
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(1), 6), seeds32(1)),-13)
seeds32 (1) = IEOR (ISHFT(IAND(seeds32(1),-2_i_32), 18), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(2), 2), seeds32(2)),-27)
seeds32 (2) = IEOR (ISHFT(IAND(seeds32(2),-8_i_32), 2), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(3), 13), seeds32(3)),-21)
seeds32 (3) = IEOR (ISHFT(IAND(seeds32(3),-16_i_32), 7), temp_uniform)
temp_uniform = ISHFT (IEOR(ISHFT(seeds32(4), 3), seeds32(4)),-12)
seeds32 (4) = IEOR (ISHFT(IAND(seeds32(4),-128_i_32), 13), temp_uniform)
i_number = IEOR (IEOR(IEOR(seeds32(1), seeds32(2)), seeds32(3)), seeds32(4))
!END IF
IF (PRESENT(range)) THEN
number = 0.5_r_dp * ((range(2)-range(1))*REAL(i_number, r_dp)*normalization_r_dp+(range(1)+range(2)))
ELSE
number = 0.5_r_dp * (REAL(i_number, r_dp)*normalization_r_dp+1.0_r_dp)
END IF
IF (.FALSE.) number = Max (number, TINY(number))
END SUBROUTINE
SUBROUTINE RandomBitArray_r_sp (array, n_elements, range)
IMPLICIT NONE
INTEGER (KIND=i_wp), INTENT (IN) :: n_elements
REAL (KIND=r_sp), DIMENSION (n_elements), INTENT (OUT) :: array
REAL (KIND=r_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
INTEGER (KIND=i_wp) :: element
DO element = 1, n_elements
CALL UniformScalar_r_sp (array(element), range)
END DO
END SUBROUTINE
SUBROUTINE RandomBitArray_r_dp (array, n_elements, range)
IMPLICIT NONE
INTEGER (KIND=i_wp), INTENT (IN) :: n_elements
REAL (KIND=r_dp), DIMENSION (n_elements), INTENT (OUT) :: array
REAL (KIND=r_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
INTEGER (KIND=i_wp) :: element
DO element = 1, n_elements
CALL UniformScalar_r_dp (array(element), range)
END DO
END SUBROUTINE
SUBROUTINE UniformArrayContiguous_1_r_sp (array, n_extent1, range)
IMPLICIT NONE
INTEGER (KIND=i_wp), INTENT (IN) :: n_extent1
REAL (KIND=r_sp), DIMENSION (n_extent1), INTENT (OUT) :: array
REAL (KIND=r_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
CALL RandomBitArray_r_sp (array, Int(SIZE(array), KIND=i_wp))
END SUBROUTINE
SUBROUTINE UniformArrayContiguous_1_r_dp (array, n_extent1, range)
IMPLICIT NONE
INTEGER (KIND=i_wp), INTENT (IN) :: n_extent1
REAL (KIND=r_dp), DIMENSION (n_extent1), INTENT (OUT) :: array
REAL (KIND=r_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
CALL RandomBitArray_r_dp (array, Int(SIZE(array), KIND=i_wp))
END SUBROUTINE
SUBROUTINE UniformArrayContiguous_2_r_sp (array, n_extent1, n_extent2, range)
IMPLICIT NONE
INTEGER (KIND=i_wp), INTENT (IN) :: n_extent1, n_extent2
REAL (KIND=r_sp), DIMENSION (n_extent1, n_extent2), INTENT (OUT) :: array
REAL (KIND=r_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
CALL RandomBitArray_r_sp (array, Int(SIZE(array), KIND=i_wp))
END SUBROUTINE
SUBROUTINE UniformArrayContiguous_2_r_dp (array, n_extent1, n_extent2, range)
IMPLICIT NONE
INTEGER (KIND=i_wp), INTENT (IN) :: n_extent1, n_extent2
REAL (KIND=r_dp), DIMENSION (n_extent1, n_extent2), INTENT (OUT) :: array
REAL (KIND=r_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
CALL RandomBitArray_r_dp (array, Int(SIZE(array), KIND=i_wp))
END SUBROUTINE
SUBROUTINE UniformArrayContiguous_3_r_sp (array, n_extent1, n_extent2, n_extent3, range)
IMPLICIT NONE
INTEGER (KIND=i_wp), INTENT (IN) :: n_extent1, n_extent2, n_extent3
REAL (KIND=r_sp), DIMENSION (n_extent1, n_extent2, n_extent3), INTENT (OUT) :: array
REAL (KIND=r_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
CALL RandomBitArray_r_sp (array, Int(SIZE(array), KIND=i_wp))
END SUBROUTINE
SUBROUTINE UniformArrayContiguous_3_r_dp (array, n_extent1, n_extent2, n_extent3, range)
IMPLICIT NONE
INTEGER (KIND=i_wp), INTENT (IN) :: n_extent1, n_extent2, n_extent3
REAL (KIND=r_dp), DIMENSION (n_extent1, n_extent2, n_extent3), INTENT (OUT) :: array
REAL (KIND=r_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
CALL RandomBitArray_r_dp (array, Int(SIZE(array), KIND=i_wp))
END SUBROUTINE
SUBROUTINE UniformArrayLoops_1_r_sp (array, range)
IMPLICIT NONE
REAL (KIND=r_sp), DIMENSION (:), INTENT (OUT) :: array
REAL (KIND=r_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
INTEGER (KIND=i_wp) :: i1
DO i1 = Int (LBOUND(array, DIM=1), KIND=i_wp), Int (UBOUND(array, DIM=1), KIND=i_wp)
CALL UniformScalar_r_sp (array(i1), range)
END DO
END SUBROUTINE
SUBROUTINE UniformArrayLoops_1_r_dp (array, range)
IMPLICIT NONE
REAL (KIND=r_dp), DIMENSION (:), INTENT (OUT) :: array
REAL (KIND=r_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
INTEGER (KIND=i_wp) :: i1
DO i1 = Int (LBOUND(array, DIM=1), KIND=i_wp), Int (UBOUND(array, DIM=1), KIND=i_wp)
CALL UniformScalar_r_dp (array(i1), range)
END DO
END SUBROUTINE
SUBROUTINE UniformArrayLoops_2_r_sp (array, range)
IMPLICIT NONE
REAL (KIND=r_sp), DIMENSION (:, :), INTENT (OUT) :: array
REAL (KIND=r_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
INTEGER (KIND=i_wp) :: i1, i2
DO i2 = Int (LBOUND(array, DIM=2), KIND=i_wp), Int (UBOUND(array, DIM=2), KIND=i_wp)
DO i1 = Int (LBOUND(array, DIM=1), KIND=i_wp), Int (UBOUND(array, DIM=1), KIND=i_wp)
CALL UniformScalar_r_sp (array(i1, i2), range)
END DO
END DO
END SUBROUTINE
SUBROUTINE UniformArrayLoops_2_r_dp (array, range)
IMPLICIT NONE
REAL (KIND=r_dp), DIMENSION (:, :), INTENT (OUT) :: array
REAL (KIND=r_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
INTEGER (KIND=i_wp) :: i1, i2
DO i2 = Int (LBOUND(array, DIM=2), KIND=i_wp), Int (UBOUND(array, DIM=2), KIND=i_wp)
DO i1 = Int (LBOUND(array, DIM=1), KIND=i_wp), Int (UBOUND(array, DIM=1), KIND=i_wp)
CALL UniformScalar_r_dp (array(i1, i2), range)
END DO
END DO
END SUBROUTINE
SUBROUTINE UniformArrayLoops_3_r_sp (array, range)
IMPLICIT NONE
REAL (KIND=r_sp), DIMENSION (:, :, :), INTENT (OUT) :: array
REAL (KIND=r_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
INTEGER (KIND=i_wp) :: i1, i2, i3
DO i3 = Int (LBOUND(array, DIM=3), KIND=i_wp), Int (UBOUND(array, DIM=3), KIND=i_wp)
DO i2 = Int (LBOUND(array, DIM=2), KIND=i_wp), Int (UBOUND(array, DIM=2), KIND=i_wp)
DO i1 = Int (LBOUND(array, DIM=1), KIND=i_wp), Int (UBOUND(array, DIM=1), KIND=i_wp)
CALL UniformScalar_r_sp (array(i1, i2, i3), range)
END DO
END DO
END DO
END SUBROUTINE
SUBROUTINE UniformArrayLoops_3_r_dp (array, range)
IMPLICIT NONE
REAL (KIND=r_dp), DIMENSION (:, :, :), INTENT (OUT) :: array
REAL (KIND=r_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: range
INTEGER (KIND=i_wp) :: i1, i2, i3
DO i3 = Int (LBOUND(array, DIM=3), KIND=i_wp), Int (UBOUND(array, DIM=3), KIND=i_wp)
DO i2 = Int (LBOUND(array, DIM=2), KIND=i_wp), Int (UBOUND(array, DIM=2), KIND=i_wp)
DO i1 = Int (LBOUND(array, DIM=1), KIND=i_wp), Int (UBOUND(array, DIM=1), KIND=i_wp)
CALL UniformScalar_r_dp (array(i1, i2, i3), range)
END DO
END DO
END DO
END SUBROUTINE
SUBROUTINE NormalScalar_r_sp (number, mean_std)
IMPLICIT NONE
REAL (KIND=r_sp), INTENT (OUT) :: number
REAL (KIND=r_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: mean_std
REAL (KIND=r_sp), PARAMETER :: s = 0.449871_r_sp, t = - 0.386595_r_sp, a = 0.19600_r_sp, b = 0.25472_r_sp, r1 = 0.27597_r_sp, r2 &
& = 0.27846_r_sp, r = 1.7156_r_sp
REAL (KIND=r_sp) :: u, v, x, y, q
DO
CALL UniformScalar_r_sp (u)
u = Max (u, TINY(u))
CALL UniformScalar_r_sp (v)
v = r * (v-0.5_r_sp)
x = u - s
y = Abs (v) - t
q = x ** 2 + y * (a*y-b*x)
IF (q < r1) EXIT
IF (q > r2) CYCLE
IF (v**2 <-4.0_r_sp*Log(u)*(u**2)) EXIT
END DO
number = v / u
IF (PRESENT(mean_std)) THEN
number = mean_std (2) * number + mean_std (1)
END IF
END SUBROUTINE
SUBROUTINE NormalScalar_r_dp (number, mean_std)
IMPLICIT NONE
REAL (KIND=r_dp), INTENT (OUT) :: number
REAL (KIND=r_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: mean_std
REAL (KIND=r_dp), PARAMETER :: s = 0.449871_r_dp, t = - 0.386595_r_dp, a = 0.19600_r_dp, b = 0.25472_r_dp, r1 = 0.27597_r_dp, r2 &
& = 0.27846_r_dp, r = 1.7156_r_dp
REAL (KIND=r_dp) :: u, v, x, y, q
DO
CALL UniformScalar_r_dp (u)
u = Max (u, TINY(u))
CALL UniformScalar_r_dp (v)
v = r * (v-0.5_r_dp)
x = u - s
y = Abs (v) - t
q = x ** 2 + y * (a*y-b*x)
IF (q < r1) EXIT
IF (q > r2) CYCLE
IF (v**2 <-4.0_r_dp*Log(u)*(u**2)) EXIT
END DO
number = v / u
IF (PRESENT(mean_std)) THEN
number = mean_std (2) * number + mean_std (1)
END IF
END SUBROUTINE
SUBROUTINE NormalArrayLoops_1_r_sp (array, mean_std)
IMPLICIT NONE
REAL (KIND=r_sp), DIMENSION (:), INTENT (OUT) :: array
REAL (KIND=r_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: mean_std
INTEGER (KIND=i_wp) :: i1
DO i1 = Int (LBOUND(array, DIM=1), KIND=i_wp), Int (UBOUND(array, DIM=1), KIND=i_wp)
CALL NormalScalar_r_sp (array(i1), mean_std)
END DO
END SUBROUTINE
SUBROUTINE NormalArrayLoops_1_r_dp (array, mean_std)
IMPLICIT NONE
REAL (KIND=r_dp), DIMENSION (:), INTENT (OUT) :: array
REAL (KIND=r_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: mean_std
INTEGER (KIND=i_wp) :: i1
DO i1 = Int (LBOUND(array, DIM=1), KIND=i_wp), Int (UBOUND(array, DIM=1), KIND=i_wp)
CALL NormalScalar_r_dp (array(i1), mean_std)
END DO
END SUBROUTINE
SUBROUTINE NormalArrayLoops_2_r_sp (array, mean_std)
IMPLICIT NONE
REAL (KIND=r_sp), DIMENSION (:, :), INTENT (OUT) :: array
REAL (KIND=r_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: mean_std
INTEGER (KIND=i_wp) :: i1, i2
DO i2 = Int (LBOUND(array, DIM=2), KIND=i_wp), Int (UBOUND(array, DIM=2), KIND=i_wp)
DO i1 = Int (LBOUND(array, DIM=1), KIND=i_wp), Int (UBOUND(array, DIM=1), KIND=i_wp)
CALL NormalScalar_r_sp (array(i1, i2), mean_std)
END DO
END DO
END SUBROUTINE
SUBROUTINE NormalArrayLoops_2_r_dp (array, mean_std)
IMPLICIT NONE
REAL (KIND=r_dp), DIMENSION (:, :), INTENT (OUT) :: array
REAL (KIND=r_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: mean_std
INTEGER (KIND=i_wp) :: i1, i2
DO i2 = Int (LBOUND(array, DIM=2), KIND=i_wp), Int (UBOUND(array, DIM=2), KIND=i_wp)
DO i1 = Int (LBOUND(array, DIM=1), KIND=i_wp), Int (UBOUND(array, DIM=1), KIND=i_wp)
CALL NormalScalar_r_dp (array(i1, i2), mean_std)
END DO
END DO
END SUBROUTINE
SUBROUTINE NormalArrayLoops_3_r_sp (array, mean_std)
IMPLICIT NONE
REAL (KIND=r_sp), DIMENSION (:, :, :), INTENT (OUT) :: array
REAL (KIND=r_sp), DIMENSION (2), INTENT (IN), OPTIONAL :: mean_std
INTEGER (KIND=i_wp) :: i1, i2, i3
DO i3 = Int (LBOUND(array, DIM=3), KIND=i_wp), Int (UBOUND(array, DIM=3), KIND=i_wp)
DO i2 = Int (LBOUND(array, DIM=2), KIND=i_wp), Int (UBOUND(array, DIM=2), KIND=i_wp)
DO i1 = Int (LBOUND(array, DIM=1), KIND=i_wp), Int (UBOUND(array, DIM=1), KIND=i_wp)
CALL NormalScalar_r_sp (array(i1, i2, i3), mean_std)
END DO
END DO
END DO
END SUBROUTINE
SUBROUTINE NormalArrayLoops_3_r_dp (array, mean_std)
IMPLICIT NONE
REAL (KIND=r_dp), DIMENSION (:, :, :), INTENT (OUT) :: array
REAL (KIND=r_dp), DIMENSION (2), INTENT (IN), OPTIONAL :: mean_std
INTEGER (KIND=i_wp) :: i1, i2, i3
DO i3 = Int (LBOUND(array, DIM=3), KIND=i_wp), Int (UBOUND(array, DIM=3), KIND=i_wp)
DO i2 = Int (LBOUND(array, DIM=2), KIND=i_wp), Int (UBOUND(array, DIM=2), KIND=i_wp)
DO i1 = Int (LBOUND(array, DIM=1), KIND=i_wp), Int (UBOUND(array, DIM=1), KIND=i_wp)
CALL NormalScalar_r_dp (array(i1, i2, i3), mean_std)
END DO
END DO
END DO
END SUBROUTINE
SUBROUTINE SaveRandomSeeds (rng, state)
INTEGER (KIND=i_wp), INTENT (IN), OPTIONAL :: rng
TYPE (Random_Real_Seeds), INTENT (OUT), OPTIONAL :: state
IF (PRESENT(rng)) rng_seeds(rng)% seeds = seeds32
IF (PRESENT(state)) state% seeds = seeds32
END SUBROUTINE
SUBROUTINE ResetRandomSeeds (rng, state)
INTEGER (KIND=i_wp), INTENT (IN), OPTIONAL :: rng
TYPE (Random_Real_Seeds), INTENT (OUT), OPTIONAL :: state
IF (PRESENT(rng)) rng_seeds(rng)% seeds = 0
IF (PRESENT(state)) state% seeds = 0
END SUBROUTINE
SUBROUTINE RestoreRandomSeeds (rng, state)
INTEGER (KIND=i_wp), INTENT (IN), OPTIONAL :: rng
TYPE (Random_Real_Seeds), INTENT (IN), OPTIONAL :: state
IF (PRESENT(rng)) THEN
IF (ALL(rng_seeds(rng)% seeds == 0)) THEN
WRITE (0,*) "Trying to reset seeds to invalid state rng=", rng
STOP
END IF
seeds32 = rng_seeds(rng)% seeds
END IF
IF (PRESENT(state)) seeds32 = state% seeds
END SUBROUTINE
END MODULE Random_Numbers
