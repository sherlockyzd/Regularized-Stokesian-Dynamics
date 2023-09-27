module conjugate_gradient_method
use omp_lib
implicit none

PRIVATE
public:: conjugate_gradient

contains


SUBROUTINE conjugate_gradient(N,A,b,x)

!character(len=32) :: arg, mode
integer,intent(in)::  N
real*8, intent(in)::  A(N,N),b(N)
real*8, intent(inout)::  x(N)

integer :: i,j,k,it,NUMPROC 
real*8 :: Ap(N),r(N,2), p(N,2)
real*8 :: alpha, beta, rr, rr2,  bnorm, pAp
real*8, parameter :: eps = 1.0E-9
!real*8 :: t1, t2, t3

r = 0.0
p = 0.0
r(:,1) = b
p(:,1) = r(:,1)
rr = dotprod(b,b)
bnorm = sqrt(rr)

! method
it = 0
do while (sqrt(rr) / bnorm > eps)
    Ap = multmv(A,p(:,1))
    pAp =  dotprod(Ap,p(:,1))
    alpha = rr / pAp   
    r(:,2) = r(:,1) - alpha * Ap 
    rr2 = dotprod(r(:,2),r(:,2)) 
    beta = rr2 / rr
    p(:,2) = r(:,2) + beta * p(:,1) 
    x = x + alpha * p(:,1)
    rr = rr2
    p(:,1) = p(:,2)
    r(:,1) = r(:,2)
    it = it + 1
enddo

print "('||r||/||b|| = ',ES10.2, ' on ',i5 ' iters.')", sqrt(rr) / bnorm, it
print '("CHECK of solve (method): ",L)',check_solve(A,x,b,eps)
end SUBROUTINE conjugate_gradient

! Dot product
real*8 function dotprod(u, v)
real*8, dimension(:), intent(in) :: u, v
integer :: i
dotprod = 0.0
if(size(u) /= size(v))return
!$OMP PARALLEL DO REDUCTION(+:DOTPROD)
do i=1,size(u)
    dotprod = dotprod + u(i) * v(i)
enddo
end function dotprod

! Matrix multiplication
function multmat(A,B)
real*8, dimension(:,:), intent(in) :: A, B
real*8, dimension(size(A, dim = 1),size(B, dim = 2))  :: multmat
real*8, dimension(size(A, dim = 2)) :: row
integer :: i,j,k,m,n,p

multmat = 0.
n = size(A, dim = 1)
m = size(B, dim = 2)
p = size(A, dim = 2)
if(size(A, dim = 2) /= size(B, dim = 1))return

!$OMP PARALLEL DO PRIVATE(i,j,k,row) SHARED(A,B,multmat,m,n,p)
do i = 1, n
    row(:) = A(i,:)
    do j = 1, m
    do k = 1, p
        multmat(i,j) = multmat(i,j) + row(k) * B(k,j)
    enddo
    enddo
enddo
end function multmat

! Matrix on vector multiplication
function multmv(A,v)
real*8, dimension(:,:), intent(in) :: A
real*8, dimension(:), intent(in) :: v
real*8, dimension(size(A, dim = 1))  :: multmv
integer :: i,j,m,n

multmv = 0.
n = size(A, dim = 1)
m = size(v)
if(size(A, dim = 2) /= m)return

!$OMP PARALLEL DO PRIVATE(i,j) SHARED(A,v,multmv,m,n)
do i = 1, n
    do j = 1, m
        multmv(i) = multmv(i) + A(i,j) * v(j)
    enddo
enddo
end function multmv

! Check the solution of the linear system
logical function check_solve(A,x,b,eps)
real*8, dimension(:,:), intent(in) :: A
real*8, dimension(:), intent(in) :: x, b
real*8, intent(in), optional :: eps
real*8 :: e = 1.0E-5
real*8 :: summ
integer :: i,j,k,n
logical :: flag
 check_solve = .true.
n = size(A, dim = 1)
if (check_solve) then
    !$OMP PARALLEL DO PRIVATE(i,j,summ) SHARED(a,x,b,check_solve,e,n)
    do i = 1, n
        summ = b(i)
        do j = 1, n
            summ = summ - A(i,j) * x(j)
        enddo
        if (abs(summ) > e) then
            !$OMP CRITICAL
             check_solve = .false.
            !$OMP END CRITICAL
        endif
    enddo
    !$OMP END PARALLEL DO
endif
end function check_solve

end module conjugate_gradient_method
