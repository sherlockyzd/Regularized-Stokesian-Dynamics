!>
!! 
!!   <B>  LINEAR ALGEBRA TOOLS  </B> 
!!
!! @author Charles Pierre
!!
!>

module algebra_lin

  use real_type, only: RP
  use basic_tools

  implicit none  
  private

  !       %----------------------------------------%
  !       |                                        |
  !       |          PUBLIC DATA                   |
  !       |                                        |
  !       %----------------------------------------%
  !!
  !! LINEAR ALGEBRA
  !!
  !! GENERIC
  public :: matVecProd         ! tested 
  !! NON GENERIC

  public :: invtrisup          ! tested  



  !       %----------------------------------------%
  !       |                                        |
  !       |       GENERIc SUBROUTINES              |
  !       |                                        |
  !       %----------------------------------------%


  interface matVecProd
     module procedure algebra_lin_matVecProd   !! TESTED
  end interface matVecProd

contains

  subroutine algebra_lin_matVecProd(y,A,X)
    real(RP), dimension(:)  , intent(in) :: X
    real(RP), dimension(:)  , intent(out):: y
    real(RP), dimension(size(y,1), size(x,1)), intent(in) :: A

    integer :: ii


    if ( size(A,2)/= size(x,1) )       &
         &  call quit( &
            &  'algebra_lin: algebra_lin_matVecProd: 1' )
    if ( size(A,1)/= size(y,1) )       &
         &  call quit( &
            &  'algebra_lin: algebra_lin_matVecProd: 2' )

    
    y = x(1)*A(:,1)
    do ii=2, size(x,1)
       y = y + x(ii)*A(:,ii)
    end do

  end subroutine algebra_lin_matVecProd  


  !> Calcul de x=M^{-1}vec avec M triangulaire superieure
  !>
  subroutine invtrisup(res,M,vec)

    real(RP), dimension(:,:), intent(in)  :: M
    real(RP), dimension(:)  , intent(in)  :: vec
    real(RP), dimension(:)  , intent(out) :: res

    integer :: ii,ji,ni


    if ( size(M,1)/= size(M,2) )   call quit( &
            &  'algebra_lin: invtrisup: 1' )
    if ( size(M,1)/= size(vec,1) ) call quit( &
            &  'algebra_lin: invtrisup: 2' )
    if ( size(M,1)/= size(res,1) ) call quit( &
            &  'algebra_lin: invtrisup: 3' )


    ni=size(M,1)
    res(ni)=vec(ni)/M(ni,ni)

    do ii=ni-1,1,-1
       res(ii)=vec(ii)
       do ji=ii+1,ni
          res(ii)=res(ii)-M(ii,ji)*res(ji)
       end do
       res(ii)=res(ii)/M(ii,ii)
    end do

  end subroutine invtrisup



end module algebra_lin


