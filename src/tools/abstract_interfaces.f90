!>
!!
!!<B>   DEFINITION OF ABSTRACT_INTERFACES FOR THE LIBRARY CHORAL </B>
!!
!! When a function or a subroutine is passed as an argument
!! or used inside a derived type: 
!!<br> it requires the definition of an abstract interface.
!!
!! Abstract interfaces used in the library choral 
!! are defined in this module.
!! 
!! @author Charles Pierre, april 2019
!!
!> 

module abstract_interfaces

  use real_type, only: RP
  
  implicit none
  private

  public :: RnToRn


  abstract interface


     !> Abstract interface:  \f$ f:~ R^n \mapsto R^n \f$  
     subroutine RnToRn(Fy, y)
       import :: RP
       real(RP), dimension(:), intent(out) :: Fy
       real(RP), dimension(:), intent(in)  :: y
     end subroutine RnToRn
     
  end interface

end module abstract_interfaces
