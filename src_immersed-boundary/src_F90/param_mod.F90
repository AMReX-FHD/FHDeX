module param

  use amrex_fort_module, only : rt => amrex_real

! Parameters limiting user-specifed input.
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  ! Maximum number of items for specifying initial conditions
  integer, parameter :: dim_ic = 500
  ! Maximum number of items for specifying boundary conditions
  integer, parameter :: dim_bc = 500
  ! Maximum number of items for specifying point sources
  integer, parameter :: dim_ps = 5000
  ! Maximum number of solids phases
  integer, parameter :: dim_m = 10
  ! Maximum number of gas species
  integer, parameter :: dim_n_g = 100
  ! Maximum number of solids species per phase.
  integer, parameter :: dim_n_s = 100
  ! Maximum number of user-defined output files
  integer, parameter :: dim_usr = 128

  ! Number of Equation types:
  !  1) Gas pressure
  !  2) Gas and solids U-Momentum equation
  !  3) Gas and solids V-Momentum equation
  !  4) Gas and solids W-Momentum equation
  integer, parameter :: dim_eqs = 4

! Parameters describing problem size: (set from user input)
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

! Parameters for testing if user input was specifed.
      real(rt), parameter :: undefined = 9.87654321d31
      integer,      parameter :: undefined_i = 987654321
      character,    parameter :: undefined_c = ' '

! Cutoffs for large and small numbers
      real(rt), parameter :: large_number = 1.0d32
      real(rt), parameter :: small_number = 1.0d-15

! Common parameter constants
      real(rt), parameter :: zero = 0.0d0
      real(rt), parameter :: half = 0.5d0
      real(rt), parameter :: one  = 1.0d0
      real(rt), parameter :: two  = 2.0d0

      real(rt), parameter :: my_huge  = 1.d200

      interface is_defined
         module procedure is_defined_db
         module procedure is_defined_i
      end interface is_defined

      interface is_undefined
         module procedure is_undefined_db
         module procedure is_undefined_i
      end interface is_undefined

   contains

      pure logical function is_defined_db(x)
         real(rt), intent(in) :: x
         is_defined_db = .not.equal(x, undefined)
      end function is_defined_db

      pure logical function is_defined_i(x)
         integer, intent(in) :: x
         is_defined_i = (x /= undefined_i)
      end function is_defined_i

      pure logical function is_undefined_db(x)
         real(rt), intent(in) :: x
         is_undefined_db = equal(x, undefined)
      end function is_undefined_db

      pure logical function is_undefined_i(x)
         integer, intent(in) :: x
         is_undefined_i = (x == undefined_i)
      end function is_undefined_i

      pure logical function equal(x, y)
         real(rt), intent(in) :: x, y
         equal = (abs(x-y) < epsilon(x))
      end function equal

      END MODULE param
