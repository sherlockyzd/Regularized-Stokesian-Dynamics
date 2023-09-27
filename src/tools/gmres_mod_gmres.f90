!!
!!  TEST ON gmres_mod_gmres
!!
!!                


module gmres_mod_gmres

  use basic_tools 
  use real_type
  use algebra_lin
  use gmres_mod

  implicit none
  private


  public::simplegmres

contains


  subroutine simplegmres(xout,n,rhs,A,iter,nb,res)
    integer , intent(in)   :: n
    real(RP), dimension(n,n),intent(in) :: A
    real(RP), dimension(n),intent(in) :: rhs
    real(RP), dimension(n),intent(inout) :: xout
    integer , intent(out)   :: iter, nb
    real(RP), intent(out)   :: res


    real(RP), dimension(n)   ::  sol2
    real(RP) ::  tol
    integer:: itermax


    !xout = 0.0_RP
    itermax=200
    !! searching sol2 solution of A*x = rhs
    tol = REAL_TOL * 100.0_RP
    call gmres(xout, iter, nb, res, rhs, A_prod, tol, itermax, 2, 0)
       
    !! discrepency between the exact solution 
    !! and the computed solution


    !print*, "gmres_mod_gmres_test = ok"
    print*, "iters=============", iter
    print*, "final residual====", res

    contains

      subroutine A_prod(y,x)
        real(RP), dimension(:), intent(out) :: y
        real(RP), dimension(:), intent(in)  :: x

        call matVecProd(y,A,X)

      end subroutine A_prod

    end subroutine simplegmres
  
end module gmres_mod_gmres
