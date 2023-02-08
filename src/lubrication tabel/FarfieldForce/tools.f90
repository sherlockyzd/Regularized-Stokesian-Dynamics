
!****************************************

      MODULE TENSORS  
      implicit none  
      real*8,parameter :: PAI=3.141592653589793238462643383279502884_8
      real*8,parameter :: sqrtPI=1.7724538509055160272981674833411_8
      real*8,parameter :: pai2=9.8696044010893586188344909998762_8
      real*8,parameter :: PIP5=0.564189583547756286948079451560772585844_8
      REAL*8 :: KDirac(3,3),kd(3,3),EPS(3,3,3),Y2(5,3,3),Y21(5,3,3),Y22(5,3,3)
      REAL*8 :: II(3,3,3,3),L1(3,3),L2(3,3),L3(3,3),L4(3,3,3),L5(3,3,3)
      real*8 :: L6(3,3,3),L7(3,3,3,3),L8(3,3,3,3),L9(3,3,3,3)
      REAL*8 :: L4_tr(5,3),L5_tr(5,3),L6_tr(5,3),L7_tr(5,5),L8_tr(5,5),L9_tr(5,5)
      REAL*8 :: X11A,X12A,Y11A,Y12A,Y11B,Y12B,X11C,X12C,Y11C,Y12C,X11G,X12G,Y11G,Y12G
      real*8 :: Y11H,Y12H,X11M,X12M,Y11M,Y12M,Z11M,Z12M
      END MODULE TENSORS

      SUBROUTINE INIT_KDirac()
      USE TENSORS
      IMPLICIT NONE

      KDirac=0.D0

      KDirac(1,1)=1.D0
      KDirac(2,2)=1.D0
      KDirac(3,3)=1.D0
      KD=KDirac

      RETURN
      END

!***********************************************************
!***********************************************************
!***********************************************************

!  lapack based matrix inversion

      SUBROUTINE MATREV(A,NN)
      IMPLICIT NONE
      INTEGER,intent(in):: NN
      REAL*8,intent(inout):: A(NN,NN)
      INTEGER INFO,I,J
      CHARACTER*1 UPLO
      PARAMETER (UPLO = 'U')

      CALL DPOTRF(UPLO,NN,A,NN,INFO)
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, the leading minor of order k is not
!               positive definite, and the factorization could not be
!               completed.
      IF(INFO.GT.0) THEN
      WRITE(*,*) 'k=',INFO
      WRITE(*,*) 'The leading minor of order k is not'
      WRITE(*,*) 'positive definite, and the factorization could not be'
      WRITE(*,*) 'completed.'
      ENDIF
      CALL DPOTRI(UPLO,NN,A,NN,INFO)

      DO I=2,NN
       DO J=1,I-1
        A(I,J)=A(J,I)
       ENDDO
      ENDDO

      RETURN
      END


    SUBROUTINE CALC_EPSR(EPSR,RW)
    IMPLICIT NONE
    REAL*8,intent(out):: EPSR(3,3)
    REAL*8,intent(in):: RW(3)

    EPSR=0.D0

    EPSR(1,2)= RW(3)
    EPSR(2,3)= RW(1)
    EPSR(3,1)= RW(2)
    
    EPSR(2,1)=-RW(3)
    EPSR(3,2)=-RW(1)
    EPSR(1,3)=-RW(2)

    RETURN
    END

      SUBROUTINE INIT_EPS()
      USE TENSORS
      IMPLICIT NONE

      EPS=0.D0

      EPS(1,2,3) = 1.D0
      EPS(1,3,2) = -1.D0
      EPS(2,1,3) = -1.D0
      EPS(2,3,1) = 1.D0
      EPS(3,1,2) = 1.D0
      EPS(3,2,1) = -1.D0

      RETURN
      END

      SUBROUTINE INIT_Y2()
      USE TENSORS
      IMPLICIT NONE
      Y2=0.D0

      Y2(1,1,2)=1.D0/SQRT(2.D0)
      Y2(1,2,1)=1.D0/SQRT(2.D0)

      Y2(2,2,3)=-1.D0/SQRT(2.D0)
      Y2(2,3,2)=-1.D0/SQRT(2.D0)

      Y2(3,1,1)=-1.D0/SQRT(6.D0)
      Y2(3,2,2)=-1.D0/SQRT(6.D0)
      Y2(3,3,3)=SQRT(2.D0/3.D0)

      Y2(4,1,3)=-1.D0/SQRT(2.D0)
      Y2(4,3,1)=-1.D0/SQRT(2.D0)

      Y2(5,1,1)=1.D0/SQRT(2.D0)
      Y2(5,2,2)=-1.D0/SQRT(2.D0)

      Y21=Y2
      Y22=Y2
      RETURN
      END

      SUBROUTINE INIT_II()
      USE TENSORS
      IMPLICIT NONE
      integer i,j,l,m

      II=0.D0
      do i=1,3
       do j=1,3
        do l=1,3
         do m=1,3

         II(i,j,l,m)=(3.0_8*Kdirac(i,m)*Kdirac(j,l) &
           +3.0_8*Kdirac(i,l)*Kdirac(j,m)-2.0_8*Kdirac(i,j)*Kdirac(l,m))/6.0_8
         enddo
        enddo
       enddo
      enddo

      RETURN
      END

      SUBROUTINE mulT3T2_front(t3,t2,v)
      IMPLICIT NONE
      REAL*8,intent(in):: t3(3,3,3),t2(3,3)
      real(8),intent(out)::v(3)
      INTEGER I,J

      v=0.D0

      DO I=1,3
       DO J=1,3
          v = v + t3(I,J,1:3)*t2(I,J)
       ENDDO
      ENDDO

      RETURN
      END


      subroutine RTD_front_transY2(C1,C1_CART)
      use tensors,only:Y21
      IMPLICIT NONE
      REAL*8,intent(out):: C1(5,3)
      REAL*8,intent(in):: C1_CART(3,3,3)          ! radiu      
      integer I,j,k,l   
      !real*8 V(3)
      C1=0.0_8
      do k=1,5
          do l=1,3
          do I=1,3
              do J=1,3
                C1(k,l) = C1(k,l) + C1_cart(I,J,l)*Y21(k,I,J)
              enddo
          enddo
          enddo
      enddo
      end
   



      SUBROUTINE mulT3T2_back(t3,t2,v)
      IMPLICIT NONE
      REAL*8,intent(in):: t3(3,3,3),t2(3,3)
      real(8),intent(out)::v(3)
      INTEGER I,J

      v=0.D0

      DO I=1,3
       DO J=1,3
          v = v + t3(1:3,I,J)*t2(I,J)
       ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE mulT2aT4T2b(t2a,t4,t2b,a)
      IMPLICIT NONE
      REAL*8,intent(out):: a
      REAL*8,intent(in):: t4(3,3,3,3),t2a(3,3),t2b(3,3)
      INTEGER I,J,K,L

      a=0.D0

      DO I=1,3
       DO J=1,3
        DO K=1,3
         DO L=1,3
          a = a + t4(I,J,K,L)*t2a(I,J)*t2b(K,L)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END

      subroutine DD_transY2(B1,B1_CART)
      use tensors
      IMPLICIT NONE
      REAL*8,intent(in):: B1_CART(3,3,3,3)
      REAL*8,intent(out):: B1(5,5)
      real*8 a,t2a(3,3),t2b(3,3)
      integer I,j

      DO I=1,5
       DO J=1,5
        t2a=Y21(I,1:3,1:3)
        t2b=Y22(J,1:3,1:3)
        CALL mulT2aT4T2b(t2a,B1_CART,t2b,a) 
        B1(I,J)=a
       ENDDO
      ENDDO
     end subroutine DD_transY2

      SUBROUTINE QQ_TR(B1,B1_CART)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: B1(5,5)
      REAL*8,intent(in):: B1_CART(3,3,3,3)       

      INTEGER I,J
      real*8 a,t2a(3,3),t2b(3,3)

      DO I=1,5
       DO J=1,5
        t2a=Y21(I,1:3,1:3)
        t2b=Y22(J,1:3,1:3)
        CALL mulT2aT4T2b(t2a,B1_CART,t2b,a) 
        B1(I,J)=a
       ENDDO
      ENDDO
      end SUBROUTINE QQ_TR

      SUBROUTINE GH_TR(C1,GC1_CART)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: C1(3,5)
      REAL*8,intent(in):: GC1_CART(3,3,3) 

      !if(front)then
       !call RTD_front_transY2(C1,GC1_CART)
       !call RTD_front_transY2(C2,HC1_CART)
      !else
       call RTD_back_transY2(C1,GC1_CART)
      !endif

      END SUBROUTINE GH_TR

      subroutine RTD_back_transY2(C1,C1_CART)
      use tensors,only:Y22
      IMPLICIT NONE
      REAL*8,intent(out):: C1(3,5)
      REAL*8,intent(in):: C1_CART(3,3,3)          ! radiu 
      integer I  
      real*8 V(3),t2b(3,3)

      DO I=1,5
       t2b=Y22(I,1:3,1:3)
       CALL mulT3T2_back(C1_CART,t2b,V)
       C1(1:3,I)=V
      ENDDO
     end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*********************************

      SUBROUTINE Init_Jg(rn,dist,Jg)
      !use LATTICE_BASE,only:lammda
      use tensors, only:kdirac
      IMPLICIT NONE
      real*8,intent(in):: rn(3),dist
      REAL*8,intent(out):: Jg(3,3)

      integer i,j

      forall(i=1:3,j=1:3)
         Jg(i,j)=kdirac(i,j)/dist+rn(i)*rn(j)/(dist**3.0_8)
      end forall
      end SUBROUTINE Init_Jg


      SUBROUTINE Init_D_Jg(rn,dist,D_Jg)
      !use LATTICE_BASE,only:lammda
      use tensors, only:kdirac
      IMPLICIT NONE
      real*8,intent(in):: rn(3),dist
      REAL*8,intent(out):: D_Jg(3,3,3)

      integer i,j,l

      forall(i=1:3,j=1:3,l=1:3)
         D_Jg(i,j,l)=(-kdirac(i,j)*rn(l)+kdirac(i,l)*rn(j)+kdirac(j,l)*rn(i))/(dist**3.0_8) &
                   -3.0_8*rn(i)*rn(j)*rn(l)/(dist**5.0_8)
      end forall
      end SUBROUTINE Init_D_Jg


      SUBROUTINE Init_DD_Jg(rn,dist,DD_Jg)
      !use LATTICE_BASE,only:lammda
      use tensors, only:kdirac
      IMPLICIT NONE
      real*8,intent(in):: rn(3),dist
      REAL*8,intent(out):: DD_Jg(3,3,3,3)

      integer i,j,m,l

      forall(i=1:3,j=1:3,l=1:3,m=1:3)
         DD_Jg(i,j,m,l)=(-kdirac(i,j)*kdirac(l,m)+kdirac(i,l)*kdirac(j,m) &
                  +kdirac(j,l)*kdirac(i,m))/(dist**3.0_8) &
                  +15.0_8*rn(i)*rn(j)*rn(l)*rn(m)/(dist**7.0_8) &
                  -3.0_8*(-kdirac(i,j)*rn(l)*rn(m)+kdirac(i,l)*rn(j)*rn(m) &
                  +kdirac(j,l)*rn(i)*rn(m)+kdirac(i,m)*rn(j)*rn(l) &
                  +kdirac(j,m)*rn(i)*rn(l)+kdirac(l,m)*rn(i)*rn(j))/(dist**5.0_8)
      end forall
      end SUBROUTINE Init_DD_Jg

      SUBROUTINE Init_D2_Jg(rn,dist,D2_Jg)
      !use LATTICE_BASE,only:lammda
      use tensors, only:kdirac
      IMPLICIT NONE
      real*8,intent(in):: rn(3),dist
      REAL*8,intent(out):: D2_Jg(3,3)

      integer i,j

      forall(i=1:3,j=1:3)
         D2_Jg(i,j)=2.0_8*kdirac(i,j)/(dist**3.0_8) &
                   -6.0_8*rn(i)*rn(j)/(dist**5.0_8)
      end forall
      end SUBROUTINE Init_D2_Jg

      SUBROUTINE Init_DD2_Jg(rn,dist,DD2_Jg)
      !use LATTICE_BASE,only:lammda
      use tensors, only:kdirac
      IMPLICIT NONE
      real*8,intent(in):: rn(3),dist
      REAL*8,intent(out):: DD2_Jg(3,3,3)

      integer i,j,k

      forall(i=1:3,j=1:3,k=1:3)
         DD2_Jg(i,j,k)=-6.0_8*(kdirac(i,j)*rn(k)+kdirac(i,k)*rn(j)+kdirac(j,k)*rn(i))/(dist**5.0_8) &
                   +30.0_8*rn(i)*rn(j)*rn(k)/(dist**7.0_8)
      end forall
      end SUBROUTINE Init_DD2_Jg

      SUBROUTINE Init_DDD2_Jg(rn,dist,DDD2_Jg)
      !use LATTICE_BASE,only:lammda
      use tensors, only:kdirac
      IMPLICIT NONE
      real*8,intent(in):: rn(3),dist
      REAL*8,intent(out):: DDD2_Jg(3,3,3,3)

      integer i,j,k,l

      forall(i=1:3,j=1:3,l=1:3,k=1:3)
         DDD2_Jg(i,j,k,l)=-6.0_8*(kdirac(i,j)*kdirac(k,l)+kdirac(i,l)*kdirac(j,k) &
                  +kdirac(j,l)*kdirac(i,k))/(dist**5.0_8) &
                  -210.0_8*rn(i)*rn(j)*rn(l)*rn(k)/(dist**9.0_8) &
                  +30.0_8*(-kdirac(i,j)*rn(l)*rn(k)+kdirac(i,l)*rn(j)*rn(k) &
                  +kdirac(j,l)*rn(i)*rn(k)+kdirac(i,k)*rn(j)*rn(l) &
                  +kdirac(j,k)*rn(i)*rn(l)+kdirac(l,k)*rn(i)*rn(j))/(dist**7.0_8)
      end forall
      end SUBROUTINE Init_DDD2_Jg   


      SUBROUTINE Init_Rg(D_Jg,Rg)
      !use LATTICE_BASE,only:lammda
      use tensors, only:EPS
      IMPLICIT NONE
      real*8,intent(in):: D_Jg(3,3,3)
      REAL*8,intent(out):: Rg(3,3)

      integer i,j,k,l

      Rg=0.0_8
      do k=1,3
       do l=1,3
        forall(i=1:3,j=1:3)
         Rg(i,j)=Rg(i,j)-0.5_8*EPS(j,k,l)*D_Jg(i,l,k)
        end forall
       enddo
      enddo

      end SUBROUTINE Init_Rg

      SUBROUTINE Init_D_Rg(DD_Jg,D_Rg)
      !use LATTICE_BASE,only:lammda
      use tensors, only:EPS
      IMPLICIT NONE
      real*8,intent(in):: DD_Jg(3,3,3,3)
      REAL*8,intent(out):: D_Rg(3,3,3)

      integer i,j,l,m,n

      D_Rg=0.0_8
      do m=1,3
       do n=1,3
        forall(i=1:3,j=1:3,l=1:3)
         D_Rg(i,j,l)=D_Rg(i,j,l)-0.5_8*EPS(j,m,n)*DD_Jg(i,n,l,m)
        end forall
       enddo
      enddo

      end SUBROUTINE Init_D_Rg

      SUBROUTINE Init_D2_Rg(DD2_Jg,D2_Rg)
      !use LATTICE_BASE,only:lammda
      use tensors, only:EPS
      IMPLICIT NONE
      real*8,intent(in):: DD2_Jg(3,3,3)
      REAL*8,intent(out):: D2_Rg(3,3)

      integer i,j,k,l

      D2_Rg=0.0_8
      do k=1,3
       do l=1,3
        forall(i=1:3,j=1:3)
         D2_Rg(i,j)=D2_Rg(i,j)-0.5_8*EPS(j,k,l)*DD2_Jg(i,l,k)
        end forall
       enddo
      enddo

      end SUBROUTINE Init_D2_Rg


      SUBROUTINE Init_Kg(D_Jg,Kg)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:EPS
      IMPLICIT NONE
      real*8,intent(in):: D_Jg(3,3,3)
      REAL*8,intent(out):: Kg(3,3,3)

      integer i,j,k

      forall(i=1:3,j=1:3,k=1:3)

         Kg(i,j,k)=0.5_8*(D_Jg(i,j,k)+D_Jg(i,k,j))

      end forall
      end SUBROUTINE Init_Kg

      SUBROUTINE Init_D_Kg(DD_Jg,D_Kg)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:EPS
      IMPLICIT NONE
      real*8,intent(in):: DD_Jg(3,3,3,3)
      REAL*8,intent(out):: D_Kg(3,3,3,3)

      integer i,j,k,l

      forall(i=1:3,j=1:3,k=1:3,l=1:3)

         D_Kg(i,j,k,l)=0.5_8*(DD_Jg(i,j,l,k)+DD_Jg(i,k,l,j))

      end forall

      end SUBROUTINE Init_D_Kg

      SUBROUTINE Init_D2_Kg(DD2_Jg,D2_Kg)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:EPS
      IMPLICIT NONE
      real*8,intent(in):: DD2_Jg(3,3,3)
      REAL*8,intent(out):: D2_Kg(3,3,3)

      !integer i,j,k

         D2_Kg=DD2_Jg

      end SUBROUTINE Init_D2_Kg

      SUBROUTINE Init_DD2_Kg(DDD2_Jg,DD2_Kg)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:EPS
      IMPLICIT NONE
      real*8,intent(in):: DDD2_Jg(3,3,3,3)
      REAL*8,intent(out):: DD2_Kg(3,3,3,3)

      !integer i,j,k

         DD2_Kg=DDD2_Jg

      end SUBROUTINE Init_DD2_Kg
!******************************************************

      SUBROUTINE Init_L(ei)
      use tensors
      IMPLICIT NONE
      real*8,intent(in):: ei(3)
      !REAL*8,intent(out):: DD2_Kg(3,3,3,3)
      integer:: i,j,k,l

      forall(i=1:3,j=1:3)
        L1(i,j)=ei(i)*ei(j)
      end forall

      L2=kD-L1

      call CALC_EPSR(L3,ei)

      L4=0.0_8

      forall(i=1:3,j=1:3,k=1:3)
         L4(i,j,k)= (L1(i,j)-1/3.0_8*kD(i,j))*ei(k)
      end forall

      L5=0.0_8
      forall(i=1:3,j=1:3,k=1:3)
          L5(i,j,k)=ei(i)*kD(j,k)+ei(j)*kD(i,k)-2.0_8*ei(i)*ei(j)*ei(k)
      end forall

      write(*,*) 'L1'
      do i=1,3
        write(*,*)L1(i,:)
      enddo
      write(*,*) 'L2'
      do i=1,3
        write(*,*)L2(i,:)
      enddo
      write(*,*) 'L3'
      do i=1,3
        write(*,*)L3(i,:)
      enddo
      call RTD_front_transY2(L4_tr,L4)
      call RTD_front_transY2(L5_tr,L5)

      L6=0.0_8
      do i=1,3
          do j=1,3
            do k=1,3
                do l=1,3
                L6(i,j,k)=L6(i,j,k)+EPS(i,k,l)*ei(l)*ei(j)+EPS(j,k,l)*ei(l)*ei(i)
                enddo
            enddo
          enddo
      enddo
      call RTD_front_transY2(L6_tr,L6)

      L7=0.0_8
      forall(i=1:3,j=1:3,k=1:3,l=1:3)
        L7(i,j,k,l)=1.5_8*(L1(i,j)-1/3.0_8*kD(i,j))*(L1(k,l)-1/3.0_8*kD(k,l))
      end forall
      call QQ_TR(L7_tr,L7)



      L8=0.0_8
      forall(i=1:3,j=1:3,k=1:3,l=1:3)
        L8(i,j,k,l)=0.5_8*(L1(i,k)*kD(j,l)+L1(j,k)*kD(i,l)+L1(i,l)*kD(j,k)+ &
            & L1(j,l)*kD(i,k)-4.0_8*L1(i,j)*L1(k,l))
      end forall
      call QQ_TR(L8_tr,L8)


      L9=0.0_8
      forall(i=1:3,j=1:3,k=1:3,l=1:3)
          L9(i,j,k,l)=0.5_8*(kD(i,k)*kD(j,l)+kD(j,k)*kD(i,l)-kD(i,j)*kD(k,l) &
              & +L1(i,j)*kD(k,l)+L1(k,l)*kD(i,j)-L1(i,k)*kD(j,l)-L1(j,k)*kD(i,l) &
              & -L1(i,l)*kD(j,k) -L1(j,l)*kD(i,k)+L1(i,j)*L1(k,l))
      end forall
      call QQ_TR(L9_tr,L9)

      write(*,*) 'L4_tr'
      do i=1,5
        write(*,*)L4_tr(i,:)
      enddo

      write(*,*) 'L5_tr'
      do i=1,5
        write(*,*)L5_tr(i,:)
      enddo

      write(*,*) 'L6_tr'
      do i=1,5
        write(*,*)L6_tr(i,:)
      enddo

      write(*,*) 'L7_tr'
      do i=1,5
        write(*,*)L7_tr(i,:)
      enddo

      write(*,*) 'L8_tr'
      do i=1,5
        write(*,*)L8_tr(i,:)
      enddo

      write(*,*) 'L9_tr'
      do i=1,5
        write(*,*)L9_tr(i,:)
      enddo
      end SUBROUTINE Init_L