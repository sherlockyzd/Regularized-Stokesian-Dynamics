!>
!!
!!<B> BASIC TOOLS </B> 
!!
!! @author Charles Pierre
!> 

module basic_tools 

  implicit none
  private

  public :: quit, warning


contains

  !> Stop program execution, display an error messahe
  !> 
  !>@param[in] message  error message
  !> 
  subroutine quit(message)
    character(len=*), intent(in) :: message
    
    write(*,*)
    write(*,*) "basic_tools     : quit"
    write(*,*) "ERROR = ", trim(message)
    stop -1

  end subroutine quit
  
  !>Warning message
  !>
  !>  - display a warning message if verb>0
  !>  - stop program execution if the variable NOWARN 
  !>    has been set to 1 in choral/CONFIG.in
  !>
  !>@param[in] message  warning message
  !>@param[in] verb     verbosity
  !> 
  subroutine warning(message, verb)
    character(len=*), intent(in) :: message
    integer         , intent(in) :: verb
    
#ifdef NOWARN
    write(*,*) "basic_tools     : "
    write(*,*) "WARNING: ", trim(message)
    stop -1

#else
    if (verb>0) then
       write(*,*) "basic_tools     : "
       write(*,*) "WARNING: ", trim(message)       
    end if
#endif

  end subroutine warning

end module basic_tools
