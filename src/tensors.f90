
      SUBROUTINE INIT_KDirac()
      USE TENSORS
      use prutil,only:kd
      IMPLICIT NONE

      KDirac=0.D0

      KDirac(1,1)=1.D0
      KDirac(2,2)=1.D0
      KDirac(3,3)=1.D0
      kd=kdirac
      RETURN
      END

!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE INIT_EPS()
      USE TENSORS
      use prutil,only:per
      IMPLICIT NONE

      EPS=0.D0

      EPS(1,2,3) = 1.D0
      EPS(1,3,2) = -1.D0
      EPS(2,1,3) = -1.D0
      EPS(2,3,1) = 1.D0
      EPS(3,1,2) = 1.D0
      EPS(3,2,1) = -1.D0

      per=eps

      RETURN
      END

!***********************************************************
!***********************************************************
!***********************************************************

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
     
      SUBROUTINE INIT_Y20()
      USE TENSORS
      IMPLICIT NONE
      Y2=0.D0

      Y2(1,1,1)=(SQRT(3.D0)+1)/2.D0
      Y2(1,2,2)=(SQRT(3.D0)-1)/2.D0

      Y2(2,1,2)=SQRT(2.D0)
      !Y2(2,3,2)=-1.D0/SQRT(2.D0)

      Y2(3,1,1)=(SQRT(3.D0)-1)/2.D0
      Y2(3,2,2)=(SQRT(3.D0)+1)/2.D0
      !Y2(3,3,3)=SQRT(2.D0/3.D0)

      Y2(4,1,3)=SQRT(2.D0)
      !Y2(4,3,1)=-1.D0/SQRT(2.D0)

      Y2(5,2,3)=SQRT(2.D0)
      !Y2(5,2,2)=-1.D0/SQRT(2.D0)

      Y21=Y2
      Y22=Y2
      RETURN
      END

      SUBROUTINE INIT_Y21()
      USE TENSORS
      IMPLICIT NONE
      Y21=0.D0

      Y21(1,1,1)=1.D0
      Y21(1,3,3)=-1.D0

      Y21(2,1,2)=2.D0

      Y21(5,2,2)=1.D0
      Y21(5,3,3)=-1.D0

      Y21(3,1,3)=2.D0
      Y21(4,2,3)=2.D0
!******************
      Y2=Y21
!******************

      Y22=0.D0

      Y22(1,1,1)=1.D0
      Y22(1,3,3)=-1.D0

      Y22(2,2,1)=2.D0

      Y22(5,2,2)=1.D0
      Y22(5,3,3)=-1.D0

      Y22(3,3,1)=2.D0
      Y22(4,3,2)=2.D0
      RETURN
      END
!***********************************************************
!***********************************************************
!***********************************************************

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


      Iunit(3,3,3,3)=0.D0

      Iunit(1,1,1,1)=1.d0
      Iunit(2,2,2,2)=1.d0
      Iunit(3,3,3,3)=1.d0


      IIiso=0.D0
      do i=1,3
       do j=1,3
        do l=1,3
         do m=1,3

         IIiso(i,j,l,m)=(Kdirac(i,m)*Kdirac(j,l) &
           +Kdirac(i,l)*Kdirac(j,m)+Kdirac(i,j)*Kdirac(l,m))
         enddo
        enddo
       enddo
      enddo

      RETURN
      END



     
!***********************************************************
!***********************************************************
!***********************************************************


      SUBROUTINE CALC_RR(RR,RW)
      IMPLICIT NONE
      REAL*8,intent(out):: RR(3,3)
      REAL*8,intent(in)::RW(3)
      INTEGER I,J

      DO I=1,3
       DO J=1,3
        RR(I,J)=RW(I)*RW(J)
       ENDDO
      ENDDO

      RETURN
      END

!***********************************************************
!***********************************************************
!***********************************************************

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

!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE CALC_UR(t3UR,RW)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: t3UR(3,3,3)
      REAL*8,intent(in):: RW(3)
      INTEGER I

      t3UR=0.D0

      DO I=1,3
       t3UR(I,I,1:3)= RW
      ENDDO

      RETURN
      END


!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE CALC_RRR(t3RRR,RW)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: t3RRR(3,3,3)
      REAL*8,intent(in):: RW(3)
      INTEGER I,J,K

      t3RRR=0.D0

      DO I=1,3
       DO J=1,3
        DO K=1,3
         t3RRR(I,J,K)= RW(I)*RW(J)*RW(K)
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END


!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE CALC_EPSRR(t3EPSRR,RW)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: t3EPSRR(3,3,3)
      REAL*8,intent(in)::RW(3)
      
      REAL*8::EPSR(3,3)
      INTEGER I,J,K

      t3EPSRR=0.D0

      CALL CALC_EPSR(EPSR,RW)

      DO I=1,3
       DO J=1,3
        DO K=1,3
         t3EPSRR(I,J,K) = EPSR(I,J)*RW(K)
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE CALC_EPSRRlim(t3EPSRR,RW)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: t3EPSRR(3,3,3)
      REAL*8,intent(in)::RW(3)
      
      REAL*8::EPSR(3,3)
      INTEGER I,J,K

      t3EPSRR=0.D0

      CALL CALC_EPSR(EPSR,RW)

      DO I=1,3
       DO J=1,3
        DO K=1,3
         t3EPSRR(I,J,K) = EPSR(J,I)*RW(K)
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE CALC_EPSRRmil(t3EPSRR,RW)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: t3EPSRR(3,3,3)
      REAL*8,intent(in)::RW(3)
      
      REAL*8::EPSR(3,3)
      INTEGER I,J,K

      t3EPSRR=0.D0

      CALL CALC_EPSR(EPSR,RW)

      DO I=1,3
       DO J=1,3
        DO K=1,3
         t3EPSRR(I,J,K) = EPSR(K,I)*RW(J)
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END



!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE CALC_t4D0(t4D0,RR)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(inout):: t4D0(3,3,3,3)
      REAL*8,intent(in)::RR(3,3)
      REAL*8 TMP(3,3)

      TMP = RR-KDirac/3.D0

      CALL t2at2b_t4(t4D0,TMP,TMP)
      t4D0 = 3.D0/2.D0*t4D0

      RETURN
      END

!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE CALC_t4D1(t4D1,RR)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(inout):: t4D1(3,3,3,3)
      REAL*8,intent(in):: RR(3,3)
      real(8)::TMP(3,3,3,3)

      t4D1=0.D0

!     DELTAjlRiRk
      CALL t2at2b_t4(TMP,KDirac,RR)
      CALL t4trans(TMP,1,4)
      t4D1=t4D1+TMP
!     DELTAikRjRl
      CALL t2at2b_t4(TMP,KDirac,RR)
      CALL t4trans(TMP,2,3)
      t4D1=t4D1+TMP
!     DELTAjkRiRl
      CALL t2at2b_t4(TMP,KDirac,RR)
      CALL t4trans(TMP,1,3)
      t4D1=t4D1+TMP
!     DELTAilRjRk
      CALL t2at2b_t4(TMP,KDirac,RR)
      CALL t4trans(TMP,2,4)
      t4D1=t4D1+TMP
!     RiRjRkRl
      CALL t2at2b_t4(TMP,RR,RR)
      t4D1=t4D1 -4.D0*TMP

      t4D1=0.5D0*t4D1

      RETURN
      END

!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE CALC_t4D2(t4D2,RR)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(inout):: t4D2(3,3,3,3)
      REAL*8,intent(in):: RR(3,3)
      REAL*8 TMP(3,3,3,3)

      t4D2=0.D0

!     DELTAikDELTAjl
      CALL t2at2b_t4(TMP,KDirac,KDirac)
      CALL t4trans(TMP,2,3)
      t4D2=t4D2+TMP
!     DELTAjkDELTAil
      CALL t2at2b_t4(TMP,KDirac,KDirac)
      CALL t4trans(TMP,1,3)
      t4D2=t4D2+TMP
!     DELTAijDELTAkl
      CALL t2at2b_t4(TMP,KDirac,KDirac)
      t4D2=t4D2-TMP
!     DELTAklRiRj
      CALL t2at2b_t4(TMP,RR,KDirac)
      t4D2=t4D2+TMP
!     DELTAijRkRl
      CALL t2at2b_t4(TMP,KDirac,RR)
      t4D2=t4D2+TMP
!     RiRjRkRl
      CALL t2at2b_t4(TMP,RR,RR)
      t4D2=t4D2-3.D0*TMP

      t4D2 = 0.5D0*t4D2

      RETURN
      END

!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE t2at2b_t4(t4,t2a,t2b)
      IMPLICIT NONE
      REAL*8,intent(out):: t4(3,3,3,3)
      real(8),intent(in)::t2a(3,3),t2b(3,3)
      INTEGER I,J,K,L

      t4=0.D0

      DO I=1,3
       DO J=1,3
        DO K=1,3
         DO L=1,3
          t4(I,J,K,L)=t2a(I,J)*t2b(K,L)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END

!***********************************************************
!***********************************************************
!***********************************************************

!     transposes two indices of 3x3x3x3 matrix i1,i2, rule i1<i2

      SUBROUTINE t4trans(t4,i1,i2)
      IMPLICIT NONE
      REAL*8,intent(inout):: t4(3,3,3,3)
      INTEGER,intent(in)::i1,i2

      INTEGER I,J,K,L
      REAL*8 t4old(3,3,3,3)
      t4old=t4 
     
      IF ((i1.EQ.1).AND.(i2.EQ.2)) THEN
       DO I=1,3
        DO J=1,3
         DO K=1,3
          DO L=1,3
           t4(I,J,K,L)=t4old(J,I,K,L)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ELSE IF ((i1.EQ.1).AND.(i2.EQ.3)) THEN
       DO I=1,3
        DO J=1,3
         DO K=1,3
          DO L=1,3
           t4(I,J,K,L)=t4old(K,J,I,L)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ELSE IF ((i1.EQ.1).AND.(i2.EQ.4)) THEN
       DO I=1,3
        DO J=1,3
         DO K=1,3
          DO L=1,3
           t4(I,J,K,L)=t4old(L,J,K,I)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ELSE IF ((i1.EQ.2).AND.(i2.EQ.3)) THEN
       DO I=1,3
        DO J=1,3
         DO K=1,3
          DO L=1,3
           t4(I,J,K,L)=t4old(I,K,J,L)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ELSE IF ((i1.EQ.2).AND.(i2.EQ.4)) THEN
       DO I=1,3
        DO J=1,3
         DO K=1,3
          DO L=1,3
           t4(I,J,K,L)=t4old(I,L,K,J)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ELSE IF ((i1.EQ.3).AND.(i2.EQ.4)) THEN
       DO I=1,3
        DO J=1,3
         DO K=1,3
          DO L=1,3
           t4(I,J,K,L)=t4old(I,J,L,K)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDIF
 
      RETURN
      END


!     transposes two indices of 3x3x3 matrix i1,i2, rule i1<i2

      SUBROUTINE t3trans(t3,i1,i2)
      IMPLICIT NONE
      REAL*8,intent(inout):: t3(3,3,3)
      INTEGER,intent(in)::i1,i2

      INTEGER I,J,K
      REAL*8 t3old(3,3,3)
      t3old=t3
     
      IF ((i1.EQ.1).AND.(i2.EQ.2)) THEN
       DO I=1,3
        DO J=1,3
         DO K=1,3
    
           t3(I,J,K)=t3old(J,I,K)
   
         ENDDO
        ENDDO
       ENDDO
      ELSE IF ((i1.EQ.1).AND.(i2.EQ.3)) THEN
       DO I=1,3
        DO J=1,3
         DO K=1,3
      
           t3(I,J,K)=t3old(K,J,I)
   
         ENDDO
        ENDDO
       ENDDO

      ELSE IF ((i1.EQ.2).AND.(i2.EQ.3)) THEN
       DO I=1,3
        DO J=1,3
         DO K=1,3
  
           t3(I,J,K)=t3old(I,K,J)
   
         ENDDO
        ENDDO
       ENDDO
      ENDIF
 
      RETURN
      END

!***********************************************************
!***********************************************************
!***********************************************************

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

!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE mulT2aT2T2b(t2a,t2,a)
      IMPLICIT NONE
      REAL*8,intent(in):: t2(3,3),t2a(3,3)
      REAL*8,intent(out)::a
      INTEGER I,J

      a=0.D0

      DO I=1,3
       DO J=1,3
          a = a + t2(I,J)*t2a(I,J)
       ENDDO
      ENDDO
      RETURN
      END
      
!***********************************************************
!***********************************************************
!***********************************************************

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

!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE EIGEN(A,D,NN)
      IMPLICIT NONE
      INTEGER,intent(in):: NN
      REAL*8,intent(in):: A(NN,NN)
      REAL*8,intent(out):: D(NN)
      
      REAL*8 E(NN-1),TAU(NN-1),WORK(2*NN)
      INTEGER INFO

      CALL DSYTRD('U', NN, A, NN, D, E, TAU, WORK, 2*NN, INFO )
      IF(INFO.NE.0) WRITE(*,*) 'EIGEN ERROR'

      CALL DORGTR('U', NN, A, NN, TAU, WORK, 2*NN, INFO )
      IF(INFO.NE.0) WRITE(*,*) 'EIGEN ERROR'

      CALL DSTEQR('V', NN, D, E, A, NN, WORK, INFO )
      IF(INFO.NE.0) WRITE(*,*) 'EIGEN ERROR'

      RETURN
      END
      
!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE MAT_INV3(LI,LR,DET_LR)
      IMPLICIT NONE
      REAL*8,intent(in)::LR(3,3)
      REAL*8,intent(out)::DET_LR,LI(3,3)
      
      det_LR=LR(1,1)*LR(2,2)*LR(3,3)+LR(1,2)*LR(2,3)*LR(3,1)+ &
             LR(1,3)*LR(2,1)*LR(3,2)-LR(1,3)*LR(2,2)*LR(3,1)-   &    
             LR(2,3)*LR(3,2)*LR(1,1)-LR(3,3)*LR(1,2)*LR(2,1)

      LI(1,1)=(1.D0/det_LR)*(LR(2,2)*LR(3,3)-LR(2,3)*LR(3,2))
      LI(1,2)=(1.D0/det_LR)*(LR(1,3)*LR(3,2)-LR(1,2)*LR(3,3))
      LI(1,3)=(1.D0/det_LR)*(LR(1,2)*LR(2,3)-LR(1,3)*LR(2,2))

      LI(2,1)=(1.D0/det_LR)*(LR(2,3)*LR(3,1)-LR(2,1)*LR(3,3))
      LI(2,2)=(1.D0/det_LR)*(LR(1,1)*LR(3,3)-LR(1,3)*LR(3,1))
      LI(2,3)=(1.D0/det_LR)*(LR(1,3)*LR(2,1)-LR(1,1)*LR(2,3))

      LI(3,1)=(1.D0/det_LR)*(LR(2,1)*LR(3,2)-LR(2,2)*LR(3,1))
      LI(3,2)=(1.D0/det_LR)*(LR(1,2)*LR(3,1)-LR(1,1)*LR(3,2))
      LI(3,3)=(1.D0/det_LR)*(LR(1,1)*LR(2,2)-LR(1,2)*LR(2,1))
      
      RETURN
      END


      subroutine DD_transY2_philip(XM12,XM12T)
      IMPLICIT NONE
      REAL*8,intent(out):: XM12(5,5)
      REAL*8,intent(in):: XM12T(3,3,3,3)          ! radiu      
 
      INTEGER I,J
      XM12(1,1) = XM12T(1,1,1,1)-2.00*XM12T(1,1,3,3)+XM12T(3,3,3,3)
      XM12(1,2) = 2.00*(XM12T(1,1,2,1)-XM12T(1,2,3,3))
      XM12(1,3) = 2.00*(XM12T(1,1,3,1)-XM12T(1,3,3,3))
      XM12(1,4) = 2.00*(XM12T(1,1,3,2)-XM12T(2,3,3,3))
      XM12(1,5) = XM12T(1,1,2,2)-XM12T(1,1,3,3)-XM12T(2,2,3,3) &
       +XM12T(3,3,3,3)
      XM12(2,2) = 4.00*XM12T(1,2,2,1)
      XM12(2,3) = 4.00*XM12T(1,2,3,1)
      XM12(2,4) = 4.00*XM12T(1,2,3,2)
      XM12(2,5) = 2.00*(XM12T(1,2,2,2)-XM12T(1,2,3,3))
      XM12(3,3) = 4.00*(XM12T(1,3,3,1))
      XM12(3,4) = 4.00*XM12T(1,3,3,2)
      XM12(3,5) = 2.00*(XM12T(2,2,3,1)-XM12T(1,3,3,3))
      XM12(4,4) = 4.00*XM12T(2,3,3,2)
      XM12(4,5) = 2.00*(XM12T(2,2,3,2)-XM12T(2,3,3,3))
      XM12(5,5) = XM12T(2,2,2,2) - 2.00*XM12T(2,2,3,3)+XM12T(3,3,3,3)

      DO I = 1,5
       DO J = I+1,5
       XM12(J,I) = XM12(I,J)
       enddo
      enddo    

     end



      subroutine RTD_transY2_philip(H12,H12T)
      IMPLICIT NONE
      REAL*8,intent(out):: H12(3,5)
      REAL*8,intent(in):: H12T(3,3,3)          ! radiu      

      H12(1,1) = H12T(1,1,1) - H12T(3,3,1)
      H12(1,2) = 2.00*H12T(1,2,1)
      H12(1,3) = 2.00*H12T(1,3,1)
      H12(1,4) = 2.00*H12T(2,3,1)
      H12(1,5) = H12T(2,2,1)-H12T(3,3,1)
      H12(2,1) = H12T(1,1,2)-H12T(3,3,2)
      H12(2,2) = 2.00*H12T(1,2,2)
      H12(2,3) = 2.00*H12T(1,3,2)
      H12(2,4) = 2.00*H12T(2,3,2)
      H12(2,5) = H12T(2,2,2)-H12T(3,3,2)
      H12(3,1) = H12T(1,1,3) - H12T(3,3,3)
      H12(3,2) = 2.00*H12T(1,2,3)
      H12(3,3) = 2.00*H12T(1,3,3)
      H12(3,4) = 2.00*H12T(2,3,3)
      H12(3,5) = H12T(2,2,3) - H12T(3,3,3)
     end


      subroutine Logtoreal(Ireal,L)
      IMPLICIT NONE
      REAL*8,intent(out):: Ireal
      logical,intent(in):: L          ! radiu      
      if(L) then
      Ireal=1
      else 
      Ireal=0
      endif
      end

      subroutine RTD_CSDtransY2(g,gtensor)
      use prutil,only:cp
      IMPLICIT NONE
      REAL*8,intent(out):: g(5,3)
      REAL*8,intent(in):: gtensor(3,3,3)          ! radiu  
        g(1,1)=gtensor(1,1,1)-gtensor(3,3,1)
        g(1,2)=gtensor(1,1,2)-gtensor(3,3,2)
        g(1,3)=gtensor(1,1,3)-gtensor(3,3,3)
        g(2,1)=2.0_cp*gtensor(1,2,1)
        g(2,2)=2.0_cp*gtensor(1,2,2)
        g(2,3)=2.0_cp*gtensor(1,2,3)
        g(3,1)=2.0_cp*gtensor(1,3,1)
        g(3,2)=2.0_cp*gtensor(1,3,2)
        g(3,3)=2.0_cp*gtensor(1,3,3)
        g(4,1)=2.0_cp*gtensor(2,3,1)
        g(4,2)=2.0_cp*gtensor(2,3,2)
        g(4,3)=2.0_cp*gtensor(2,3,3)
        g(5,1)=gtensor(2,2,1)-gtensor(3,3,1)
        g(5,2)=gtensor(2,2,2)-gtensor(3,3,2)
        g(5,3)=gtensor(2,2,3)-gtensor(3,3,3)
      end
        

      subroutine DD_CSDtransY2(m,mtensor)
      use prutil,only:cp
      IMPLICIT NONE
      REAL*8,intent(out):: m(5,5)
      REAL*8,intent(in):: mtensor(3,3,3,3)          ! radiu            
        ! line 1 of m 5*5
        m(1,1) = &
        & mtensor(1,1,1,1)-mtensor(1,1,3,3)-mtensor(3,3,1,1)+mtensor(3,3,3,3) 
        m(1,2) = &
        & mtensor(1,1,1,2)+mtensor(1,1,2,1)-mtensor(3,3,1,2)-mtensor(3,3,2,1) 
        m(1,3) = &
        & mtensor(1,1,1,3)+mtensor(1,1,3,1)-mtensor(3,3,1,3)-mtensor(3,3,3,1) 
        m(1,4) = &
        & mtensor(1,1,2,3)+mtensor(1,1,3,2)-mtensor(3,3,2,3)-mtensor(3,3,3,2) 
        m(1,5) = &
        & mtensor(1,1,2,2)-mtensor(1,1,3,3)-mtensor(3,3,2,2)+mtensor(3,3,3,3) 
        
        ! line 2 of m 5*5
        m(2,1) = 2.0_cp*(mtensor(1,2,1,1)-mtensor(1,2,3,3))
        m(2,2) = 2.0_cp*(mtensor(1,2,1,2)+mtensor(1,2,2,1))
        m(2,3) = 2.0_cp*(mtensor(1,2,1,3)+mtensor(1,2,3,1))
        m(2,4) = 2.0_cp*(mtensor(1,2,2,3)+mtensor(1,2,3,2))
        m(2,5) = 2.0_cp*(mtensor(1,2,2,2)-mtensor(1,2,3,3))
        
        ! line 3 of m 5*5
        m(3,1) = 2.0_cp*(mtensor(1,3,1,1)-mtensor(1,3,3,3))
        m(3,2) = 2.0_cp*(mtensor(1,3,1,2)+mtensor(1,3,2,1))
        m(3,3) = 2.0_cp*(mtensor(1,3,1,3)+mtensor(1,3,3,1))
        m(3,4) = 2.0_cp*(mtensor(1,3,2,3)+mtensor(1,3,3,2))
        m(3,5) = 2.0_cp*(mtensor(1,3,2,2)-mtensor(1,3,3,3))
        
        ! line 4 of m 5*5
        m(4,1) = 2.0_cp*(mtensor(2,3,1,1)-mtensor(2,3,3,3))
        m(4,2) = 2.0_cp*(mtensor(2,3,1,2)+mtensor(2,3,2,1))
        m(4,3) = 2.0_cp*(mtensor(2,3,1,3)+mtensor(2,3,3,1))
        m(4,4) = 2.0_cp*(mtensor(2,3,2,3)+mtensor(2,3,3,2))
        m(4,5) = 2.0_cp*(mtensor(2,3,2,2)-mtensor(2,3,3,3))
        
        ! line 5 of m 5*5
        m(5,1) = &
        & mtensor(2,2,1,1)-mtensor(2,2,3,3)-mtensor(3,3,1,1)+mtensor(3,3,3,3) 
        m(5,2) = &
        & mtensor(2,2,1,2)+mtensor(2,2,2,1)-mtensor(3,3,1,2)-mtensor(3,3,2,1) 
        m(5,3) = &
        & mtensor(2,2,1,3)+mtensor(2,2,3,1)-mtensor(3,3,1,3)-mtensor(3,3,3,1) 
        m(5,4) = &
        & mtensor(2,2,2,3)+mtensor(2,2,3,2)-mtensor(3,3,2,3)-mtensor(3,3,3,2) 
        m(5,5) = &
        & mtensor(2,2,2,2)-mtensor(2,2,3,3)-mtensor(3,3,2,2)+mtensor(3,3,3,3) 
     end


      subroutine EIN_CSDtransY2(einf_sp,einf_bg,NN)
      use prutil,only:cp
      IMPLICIT NONE
      REAL*8,intent(in):: einf_bg(3,3)
      integer,intent(in)::NN          ! radiu 
      REAL*8,intent(out):: einf_sp(5*NN)
      integer I
      Do i=1,NN 
        einf_sp(5*(i-1)+1) = -einf_bg(3,3)+einf_bg(1,1)
        einf_sp(5*(i-1)+2) = 2.0_cp*einf_bg(1,2)
        einf_sp(5*(i-1)+3) = 2.0_cp*einf_bg(1,3)
        einf_sp(5*(i-1)+4) = 2.0_cp*einf_bg(2,3)
        einf_sp(5*(i-1)+5) = -einf_bg(3,3)+einf_bg(2,2)
     enddo
     end




     subroutine EI_transY2(einf_sp,einf_bg)
      use tensors,only:Y21
      IMPLICIT NONE
      REAL*8,intent(in):: einf_bg(3,3)
      REAL*8,intent(out):: einf_sp(5)
      integer I
      real*8 atemp

      DO I=1,5
        call mulT2aT2T2b(Y21(I,1:3,1:3),einf_bg,aTemp)
        einf_sp(I)=aTemp
      ENDDO

      end

      subroutine EIN_transY2(EIN,EIJ,NN)
      use prutil,only:cp
      IMPLICIT NONE
      REAL*8,intent(in):: EIJ(3,3)
      integer,intent(in)::NN          ! radiu 
      REAL*8,intent(out):: EIN(5*NN)
      integer I
      real*8 EI(5)

      call EI_transY2(EI,EIJ)

      Do i=1,NN 
        EIN(5*(i-1)+1:5*i) = EI(1:5)
      enddo
     end

        !CALL mulT2aT4T2b(Y2(I,1:3,1:3),B1_CART,Y2(J,1:3,1:3),a) 
        !B1(I,J)=a 
       ! call mulT2aT2T2b(Y21(I,1:3,1:3),Eij,aTemp)
       ! EI(I)=aTemp

      !WRITE(*,*) "EI=", EI

     ! do J=1,5*NN
     !   i=mod(J,5)
     !   if(i.eq.0) i=5
     !   !WRITE(*,*) "i=", i
     !   EIN(J)=EI(i)
     ! enddo


      subroutine DD_transY2(B1,B1_CART)
      use tensors
      IMPLICIT NONE
      REAL*8,intent(in):: B1_CART(3,3,3,3)
      REAL*8,intent(out):: B1(5,5)
      real*8 a
      integer I,j
      DO I=1,5
       DO J=1,5
        CALL mulT2aT4T2b(Y21(I,1:3,1:3),B1_CART,Y22(J,1:3,1:3),a) 
        B1(I,J)=a
       ENDDO
      ENDDO

     end

      subroutine RTD_back_transY2(C1,C1_CART)
      use tensors,only:Y22
      IMPLICIT NONE
      REAL*8,intent(out):: C1(3,5)
      REAL*8,intent(in):: C1_CART(3,3,3)          ! radiu 
      integer I  
      real*8 V(3)   

      DO I=1,5
       CALL mulT3T2_back(C1_CART,Y22(I,1:3,1:3),V)
       C1(1:3,I)=V
      ENDDO
     end

      subroutine RTD_front_transY2(C1,C1_CART)
      use tensors,only:Y21
      IMPLICIT NONE
      REAL*8,intent(out):: C1(5,3)
      REAL*8,intent(in):: C1_CART(3,3,3)          ! radiu      
      integer I    
      real*8 V(3)

      DO I=1,5
       CALL mulT3T2_front(C1_CART,Y21(I,1:3,1:3),V)
       C1(I,1:3)=V
      ENDDO
     end

      SUBROUTINE init_linspace(a,b,n,x)
      real*8,intent(in)::a,b
      integer,intent(in)::n
      real*8,intent(out)::x(n+1)

      integer::i

      x=(/(i,i=0,n)/)*1.0d0/n*(b-a)+a
      
      end SUBROUTINE init_linspace
     !****************************************************
!*******************************************************

      SUBROUTINE QQ_TR(B1,B1_CART)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: B1(5,5)
      REAL*8,intent(in):: B1_CART(3,3,3,3)       

      INTEGER I,J
      real*8 a

      DO I=1,5
       DO J=1,5
        CALL mulT2aT4T2b(Y21(I,1:3,1:3),B1_CART,Y22(J,1:3,1:3),a) 
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


      subroutine subToAll_particle(NN,Nswimer,Nfilament,U_pos_swimer,U_pos_filament,U_pos)
      IMPLICIT NONE
      integer,intent(in)::NN,Nswimer,Nfilament
      real*8,intent(in)::U_pos_swimer(6*Nswimer),U_pos_filament(6*Nfilament)
      real*8,intent(inout)::U_pos(6*NN)

      U_pos(1:3*Nswimer)=U_pos_swimer(1:3*Nswimer)
      U_pos(3*Nswimer+1:3*Nswimer+3*Nfilament)=U_pos_filament(1:3*Nfilament)
      U_pos(3*NN+1:3*NN+3*Nswimer)=U_pos_swimer(3*Nswimer+1:6*Nswimer)
      U_pos(3*NN+3*Nswimer+1:3*NN+3*Nswimer+3*Nfilament)=U_pos_filament(3*Nfilament+1:6*Nfilament)
      end subroutine subToAll_particle

      subroutine AllTosub_particle(NN,Nswimer,Nfilament,U_pos_swimer,U_pos_filament,U_pos)
      IMPLICIT NONE
      integer,intent(in)::NN,Nswimer,Nfilament
      real*8,intent(out)::U_pos_swimer(6*Nswimer),U_pos_filament(6*Nfilament)
      real*8,intent(in)::U_pos(6*NN)

      U_pos_swimer(1:3*Nswimer)=U_pos(1:3*Nswimer)
      U_pos_filament(1:3*Nfilament)=U_pos(3*Nswimer+1:3*Nswimer+3*Nfilament)
      U_pos_swimer(3*Nswimer+1:6*Nswimer)=U_pos(3*NN+1:3*NN+3*Nswimer)
      U_pos_filament(3*Nfilament+1:6*Nfilament)=U_pos(3*NN+3*Nswimer+1:3*NN+3*Nswimer+3*Nfilament)
      end subroutine AllTosub_particle
