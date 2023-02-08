!>
!!
!! <B> OPENMP OPERATIONS ON 1-DIMENSIONAL REAL ARRAYS </B>
!!
!!    
!!
!! @author Charles Pierre
!>

module R1d_mod

  use real_type, only: RP, equal
  use basic_tools
  !$ use OMP_LIB

  implicit none
  private

  !       %----------------------------------------%
  !       |                                        |
  !       |          PUBLIC DATA                   |
  !       |                                        |
  !       %----------------------------------------%
  ! GENERIC

  public :: axpy         !! TESTED    
  public :: xpay         !! TESTED         
  public :: scalProd     !! TESTED
  public :: norm_2       !! TESTED
     

  !       %----------------------------------------%
  !       |                                        |
  !       |       GENERIc SUBROUTINES              |
  !       |                                        |
  !       %----------------------------------------%
  !> x = a*x + y 
  interface axpy
     module procedure R1d_axpy      !! TESTED
  end interface axpy

  !> x = x + b*y  // OR // z = x + b*y 
  interface xpay
     module procedure xpay_1      !! TESTED
     module procedure xpay_2      !! TESTED
  end interface xpay

  !> scalar product
  interface scalProd
     module procedure R1d_scalProd   !! TESTED
  end interface scalProd

  !> L2 vector norm
  interface norm_2
     module procedure R1d_norm_2     !! TESTED
  end interface norm_2


contains

  !> x(:) := a*x(:) + y(:)
  !>
  subroutine R1d_axpy(a, x, y)
    real(RP), dimension(:), intent(inout) :: x
    real(RP)              , intent(in)    :: a
    real(RP), dimension(:), intent(in)    :: y

    integer :: ii, n2


    if ( size(y,1) < size(x,1)) call quit(&
         &"R1d_mod: R1d_axpy: 1" )


    n2 = size(x,1)
    !$OMP PARALLEL 
    !$OMP DO
    do ii=1, n2
       x(ii) = a*x(ii) + y(ii) 
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine R1d_axpy

  !> x(:) := x(:) + a*y(:)
  !>
  subroutine xpay_1(x, a, y)
    real(RP), dimension(:), intent(inout) :: x
    real(RP)              , intent(in)    :: a
    real(RP), dimension(:), intent(in)    :: y

    integer :: ii, n2


    if ( size(y,1) < size(x,1)) call quit(&
         &"R1d_mod: xpay_1: 1" )


    n2 = size(x,1)
    if (equal(a,1._RP)) then
       !$OMP PARALLEL 
       !$OMP DO
       do ii=1, n2
          x(ii) = x(ii) + y(ii) 
       end do
       !$OMP END DO
       !$OMP END PARALLEL
    else
       !$OMP PARALLEL 
       !$OMP DO
       do ii=1, n2
          x(ii) = x(ii) + a*y(ii) 
       end do
       !$OMP END DO
       !$OMP END PARALLEL
    end if

  end subroutine xpay_1


  !> z(:) := x(:) + a*y(:)
  !>
  subroutine xpay_2(z, x, a, y)
    real(RP), dimension(:), intent(out) :: z
    real(RP)              , intent(in)  :: a
    real(RP), dimension(:), intent(in)  :: x, y

    integer :: ii, n2


    if ( size(x,1) < size(z,1)) call quit(&
         &"R1d_mod: xpay_2: 1" )
    if ( size(y,1) < size(z,1)) call quit(&
         &"R1d_mod: xpay_2: 2" )


    n2 = size(z,1)
    if (equal(a,1._RP)) then
       !$OMP PARALLEL 
       !$OMP DO
       do ii=1, n2
          z(ii) = x(ii) + y(ii) 
       end do
       !$OMP END DO
       !$OMP END PARALLEL

    else
       !$OMP PARALLEL 
       !$OMP DO
       do ii=1, n2
          z(ii) = x(ii) + a*y(ii) 
       end do
       !$OMP END DO
       !$OMP END PARALLEL

    end if

  end subroutine xpay_2


  !> scalar product
  !>
  function R1d_scalProd(x, y)  result(res)
    real(RP), dimension(:), intent(in) :: x, y
    real(RP)                           :: res
    integer :: ii, n2
    

    if ( size(x,1) /= size(y,1)) call quit(&
         &"R1d_mod: R1d_scalProd: 1" )


    res = 0._RP
    n2 = size(x,1)
    !$OMP PARALLEL 
    !$OMP DO reduction(+:res) 
    do ii=1, n2
       res = res + x(ii) * y(ii)
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end function R1d_scalProd

  !> Euclidian norm
  !>
  function R1d_norm_2(x)  result(res)
    real(RP), dimension(:), intent(in) :: x
    real(RP)                           :: res

    res = scalProd(x,x)
    res = sqrt(res)

  end function R1d_norm_2

end module R1d_mod
