!>
!!
!! <B> GMRES LINEAR SOLVER  </B>
!!
!! Source = Youssef SAAD, 'Iterative methods for sparse linear system'
!! https://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf
!>

module gmres_mod

  use real_type
  use basic_tools
  use abstract_interfaces, only: RnToRn
  use algebra_lin, only: invtrisup
  use R1d_mod, only: norm_2, axpy, scalProd, xpay

  implicit none  
  private 
  
  public :: gmres   !! TESTED

contains


  !> GMRES (no preconditioning)
  !> 
  !>
  !> for the linear system \f$ Ax = b\f$
  !>
  !> INPUT/OUTPUT :
  !> \li            x = initial guess / solution
  !>
  !> OUTPUT :
  !> \li            res  = final residual
  !> \li            iter = number of performed iterations
  !> \li            iter = -1 = resolution failure     
  !>
  !> INPUT :
  !> \li            b     = RHS
  !> \li            A     = \f$ A~:~~~x \mapsto A x \f$
  !>                        matrix/vector product (procedural)
  !> \li            tol   = tolerance
  !> \li            itMax = maximal iteration-number
  !> \li            rst   = restart number
  !> \li            verb  = verbosity
  !>
  
  subroutine gmres(x, iter, nbPrd, res, &
       & b, A, tol, itmax, rst, verb)

    real(RP), dimension(:), intent(inout) :: x
    integer               , intent(out)   :: iter, nbPrd
    real(RP)              , intent(out)   :: res
    procedure(RnToRn)                     :: A
    real(RP), dimension(:), intent(in)    :: b
    real(RP)              , intent(in)    :: tol
    integer               , intent(in)    :: itmax, rst, verb

    real(RP), dimension(rst+1, size(x,1) ) :: V
    real(RP), dimension(rst+1, rst       ) :: H

    real(RP), dimension(size(x,1)) :: r, w
    real(RP), dimension(rst      ) :: sn, cs, y
    real(RP), dimension(rst +1   ) :: s

    real(RP) :: nb2, nr2, temp
    integer  :: nn, ii, kk

    if (verb>1) write(*,*) 'gmres_mod       : gmres'

    iter  = 0
    nbPrd = 0


    nb2 = norm_2(b)
    nn = size(x,1)
    if (nb2/re(nn)<1E-8_RP) nb2=sqrt(re(nn))
    nb2 = 1._RP/nb2

    call A(r, x)
    nbPrd = 1
    r   = b-r
    nr2 = norm_2(r)
    res = nr2 * nb2

    if (verb>2) write(*,*)'  iter', iter,' residual', res
    if (res<tol) return

    V  = 0._RP
    H  = 0._RP
    cs = 0._RP
    sn = 0._RP

    do iter=1, itmax

       if (verb>2) write(*,*)'  iter', iter,' residual', res

       s     = 0._RP
       s(1)  = nr2

       V(1,:)  = r / nr2
      
       do ii = 1, rst  ! orthonormal basis using Gram-Schmidt

          call A(w, V(ii,:))
          nbPrd = nbPrd + 1

          do kk = 1, ii
             H(kk,ii)= ScalProd( w, V(kk,:))
             w = w - H(kk,ii)*V(kk,:)
          end do

          H(ii+1,ii) = norm_2(w)

          V(ii+1, : ) = w / H(ii+1,ii)

          ! apply Givens rotation
          !
          do kk = 1, ii-1 
             temp       =  cs(kk)*H(kk,ii) + sn(kk)*H(kk+1,ii)
             H(kk+1,ii) = -sn(kk)*H(kk,ii) + cs(kk)*H(kk+1,ii)
             H(kk,ii)   = temp
          end do

          call grotmat(cs(ii), sn(ii), H(ii,ii), H(ii+1,ii) )

          temp    =  cs(ii)*s(ii)             ! approximate residual norm
          s(ii+1) = -sn(ii)*s(ii)
          s(ii)   = temp

          H(ii,ii)   = cs(ii)*H(ii,ii) + sn(ii)*H(ii+1,ii)
          H(ii+1,ii) = 0._RP

          res  = abs(s(ii+1)) * nb2

          if (verb>2) write(*,*)'  iter', iter,' residual', res

          if ( res < tol ) then

             call invtrisup(y(1:ii), H(1:ii,1:ii), s(1:ii))
             do kk=1, ii
                x    = x + V(kk, : )*y(kk)
             end do

             return

          end if

       end do

       ! update approximation
       !
       call invtrisup(y(1:rst), H(1:rst,1:rst), s(1:rst))
       do kk=1, rst
          x    = x + V(kk, : )*y(kk)
       end do

       ! compute residual
       !
       call  A(r, x)                  
       nbPrd = nbPrd + 1
       r     = b-r
       nr2   = norm_2(r)
       res   = nr2 * nb2
       if (verb>2) write(*,*)'  iter', iter,' residual', res

       if ( res < tol ) return

    end do

    ! Convergence failure flag
    call warning("gmres_mod: gmres: not converged", verb-1)
    iter = -1

  end subroutine gmres

  !> Matrice de rotation de Givens
  subroutine grotmat(cs, sn, a, b)

    real(RP), intent(in)  :: a,b
    real(RP), intent(out) :: cs,sn

    real(RP)              :: tmp

    if (abs(b)<REAL_TOL) then
       cs = 1._RP
       sn = 0._RP

    else if (abs(b)>abs(a)) then
       tmp = a/b
       sn  = 1._RP / sqrt( 1._RP + tmp**2)
       cs  = tmp*sn

    else
       tmp = b/a
       cs  = 1._RP / sqrt( 1._RP + tmp**2)
       sn  = tmp*cs

    end if

  end subroutine grotmat

end module gmres_mod

