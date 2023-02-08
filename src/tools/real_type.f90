!> 
!! <B> REAL NUMBERS PRECISION IN CHORAL: </B> selects simple/double/quad 
!!
!! Real kind definition:
!! - \ref real_type::sp "SP" = simple precision
!! - \ref real_type::dp "DP" = double precision
!! - \ref real_type::qp "QP" = quad   precision
!!
!! Real numbers in CHORAL have the type 
!! \code{f90}  real(kind=RP) \endcode
!!\li \ref real_type::rp "RP" = real kind used in all the library
!!\li RP is determined relativelly to the preprocessor variable 'RPC'
!!\li That variable 'RPC' is set in 'choral/CONFIG.in'
!!\li Depending on the preprocessor variable 'RPC', 
!!    - 'RP' is set to 'SP', 'DP' or 'QP'
!!    - \ref real_type::real_tol "REAL_TOL" 
!!      is set to a given small value
!!..
!!
!! <b> Tools on real numbers </b>
!!\li  \ref real_type::real_equal "equal"= 
!!     check equality of real numbers up to 
!!     a tolerance of \ref real_type::real_tol "REAL_TOL",
!!\li  \ref real_type::re "re"      = 
!!     conversion of integers and rationals to real,
!!\li  r_1d  = derived type = one-dimensional real array:
!!             this type allows to manage 'arraus of arrays'!
!!\li  r_2d  = same as r_1d in dimension 2
!! 
!! @author Charles Pierre, april 2019.
!>

module real_type
  
  implicit none

  private

  !! CONSTANTS
  public :: REAL_TOL,cal_TOL,RP                !! tested

  !! DERIVED TYPES
  !public :: r_1d, r_2d

  !! ROUTINES
  public :: equal                   !! tested
  public :: re                      !! tested
  !public :: allocMem, freeMem, print


  !       %----------------------------------------%
  !       |                                        |
  !       |          CONSTANTS                     |
  !       |                                        |
  !       %----------------------------------------%

  !> simple precision for real numbers
  integer, parameter :: SP = selected_real_kind(6, 37)

  !> double precision for real numbers
  integer, parameter :: DP = selected_real_kind(12) 

  !> triple precision for real numbers
  integer, parameter :: Tp = selected_real_kind(17) 

  !> quadruple precision for real numbers
  integer, parameter :: QP = selected_real_kind(32) 

  !>  real(kind=RP) = real precision in the code
  !>  REAL_TOL = epsilon to test real equality

  integer , parameter :: RP = selected_real_kind(12) 
  real(RP), parameter :: REAL_TOL = 1E-14_RP
  real(RP), parameter :: cal_TOL = 1E-8_RP


  !       %----------------------------------------%
  !       |                                        |
  !       |          DERIVED TYPE                  |
  !       |                                        |
  !       %----------------------------------------%

 

  !       %----------------------------------------%
  !       |                                        |
  !       |          INTERFACES                    |
  !       |                                        |
  !       %----------------------------------------%



  !> Test real equality
  interface equal 
     module procedure  real_equal
  end interface equal

  !> conversion integers or rational to real
  !> - re(1)   = 1.0_RP           
  !> - re(1,3) = 1.0_RP / 3.0_RP
  interface re
     module procedure  rf, ri
  end interface re

 

contains

  !> Destructor
  !!

  !> <b> Check real equality </b>
  !>
  !> Test if |x - y| < \ref real_type::real_tol "REAL_TOL"
  !>
  function real_equal(u, v) result(res)
    logical :: res
    real(RP), intent(in) :: u, v

    res = abs(u-v) < REAL_TOL

  end function real_equal


  !> Conversion of rational fraction n/d --> real(RP)
  function rf(n, d) 
    integer, intent(in) :: n, d
    real(RP)          :: rf

    rf = real(n, RP)/real(d, RP)

  end function rf

  !> Conversion integer --> real(RP)
  function ri(n) 
    integer, intent(in) :: n
    real(RP)          :: ri

    ri = real(n, RP)

  end function ri


end module real_type
