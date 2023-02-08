! Procedures here are most computationally costly.
! For better performance it is advised to replace them with parallel
! implementations.


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


!***********************************************************
!***********************************************************
!***********************************************************

! lapack based Cholesky decomposition

      SUBROUTINE CHOLESKY(B,NN)
      IMPLICIT NONE
      INTEGER,intent(in):: NN
      REAL*8,intent(inout):: B(NN,NN)
      CHARACTER*1 UPLO
      PARAMETER (UPLO = 'L')
      INTEGER I,J,INFO

      CALL DPOTRF(UPLO,NN,B,NN,INFO)
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
      STOP
      ENDIF

!     RETURN     !!!  The strictly upper triangular part of A is not referenced.
      DO I=1,NN-1
      DO J=I+1,NN
       B(I,J)=0.D0
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CHOLESKY

!**************************************************
!******************************************************************





      SUBROUTINE M_Tr_front(APP,APQ,AQQ,A_Tr_front,A_Tr_back)
      USE SIZE,only:NN           ! NN 
      !use LATTICE_BASE,only:lammda
      !use tensors, only:EPS
      IMPLICIT NONE
      real*8,intent(in):: APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8,intent(out):: A_Tr_front(11*NN,11*NN),A_Tr_back(11*NN,11*NN)

      integer i

         A_Tr_front=0.0_8
         A_Tr_back=0.0_8

           forall(i=1:6*NN)
            A_Tr_front(i,i)=1.0_8
           end forall
           A_Tr_front(1:6*NN,6*NN+1:11*NN)=-APQ
           A_Tr_front(6*NN+1:11*NN,6*NN+1:11*NN)=-AQQ

           forall(i=1:5*NN)
            A_Tr_back(6*NN+i,6*NN+i)=-1.0_8
           end forall
           A_Tr_back(1:6*NN,1:6*NN)=APP
           A_Tr_back(6*NN+1:11*NN,1:6*NN)=transpose(APQ)
      end SUBROUTINE M_Tr_front


      SUBROUTINE Part_to_Whole(APP,APQ,AQQ,A_grmb)
      USE SIZE,only:NN           ! NN 
      IMPLICIT NONE
      real*8,intent(in):: APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8,intent(out):: A_grmb(11*NN,11*NN)

      A_grmb(1:6*NN,1:6*NN)=APP
      A_grmb(6*NN+1:11*NN,6*NN+1:11*NN)=AQQ
      A_grmb(1:6*NN,6*NN+1:11*NN)=APQ
      A_grmb(6*NN+1:11*NN,1:6*NN)=transpose(APQ)
      end SUBROUTINE Part_to_Whole

      SUBROUTINE Whole_to_Part(APP,APQ,AQQ,A_grmb)
      USE SIZE,only:NN           ! NN 
      IMPLICIT NONE
      real*8,intent(out):: APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8,intent(in):: A_grmb(11*NN,11*NN)

      APP=A_grmb(1:6*NN,1:6*NN)
      AQQ=A_grmb(6*NN+1:11*NN,6*NN+1:11*NN)
      APQ=A_grmb(1:6*NN,6*NN+1:11*NN)

      end SUBROUTINE Whole_to_Part
!*********************************************

      subroutine arrangResist(NN,rfu,rfe,rse,RPP0,RPQ0,RQQ0,RPP,RPQ,RQQ)
      implicit none
      INTEGER,intent(in):: NN
      REAL*8,intent(in):: rfu(6*NN,6*NN),rfe(6*NN,5*NN),rse(5*NN,5*NN)
      real(8),intent(in)::RPP0(6*NN,6*NN),RPQ0(6*NN,5*NN),RQQ0(5*NN,5*NN)
      real(8),intent(out)::RPP(6*NN,6*NN),RPQ(6*NN,5*NN),RQQ(5*NN,5*NN)

      RPP=RPP0+rfu
      RPQ=RPQ0+rfe
      RQQ=RQQ0+rse
      end

      subroutine averageU(RPP,U_ave,NN)
      use tensors
      use SYS_property,only:mu_f
      use CONFIG,only:radii_test
      INTEGER,intent(in):: NN
      REAL*8,intent(in):: RPP(6*NN,6*NN)
      real(8),intent(out)::U_ave
      
      REAL*8 UU(6*NN),F(6*NN)
      integer I

      U_ave=0.0_8
      F=0.0_8
      do i=1,NN
         F(3*(i-1)+1)=1.0_8*(6.0_8*pai*mu_f*radii_test)
      enddo

      UU=matmul(RPP,F)
      do i=1,NN
        U_ave=U_ave+UU(3*(i-1)+1)
      enddo

       U_ave=U_ave/NN
      end
!************************************

      subroutine Mualphabeta(mu,alpha,beta,NN)
      use tensors
      use SYS_property,only:mu_f
      use CONFIG,only:radii_test
      implicit none
      INTEGER,intent(in):: NN
      REAL*8,intent(in):: mu(5*NN,5*NN)
      real(8),intent(out)::alpha,beta
      
      REAL*8 mu_ave(5,5),muTemp(5,5)
      integer I,J
  
       mu_ave=0.0_8
        do I=1,NN
         do J=1,NN
         muTemp=mu(5*(I-1)+1:5*I,5*(J-1)+1:5*J)
         mu_ave=mu_ave+muTemp
         enddo
        enddo

      mu_ave=-mu_ave/NN/((20.0_8/3.0_8)*mu_f*pai*radii_test**3)
      alpha=mu_ave(3,3)-1
      beta=mu_ave(1,1)-1

#ifdef mucal
      real*8 DijDkl(3,3,3,3),A(3,3,3,3),B(3,3,3,3),aTemp,A_Y2(5,5),B_Y2(5,5)

      call t2at2b_t4(DijDkl,Kdirac,Kdirac)
      A=Iunit-1.0_8/3.0_8*DijDkl
      B=II


     Call DD_transY2(A_Y2,A)
     Call DD_transY2(B_Y2,B)

     write(*,*) "ok1--------------------------"
     CALL MATREV(A_Y2,5)
     CALL MATREV(B_Y2,5)
     write(*,*) "ok2--------------------------"
     !A_Y2=B_Y2+A_Y2
     !CALL MATREV(A_Y2,5*NN)

      write(*,*) 'A_Y2=',A_Y2(1,:)
      write(*,*) 'A_Y2=',A_Y2(2,:)
      write(*,*) 'A_Y2=',A_Y2(3,:)
      write(*,*) 'A_Y2=',A_Y2(4,:)
      write(*,*) 'A_Y2=',A_Y2(5,:)
      write(*,*) 'B_Y2=',b_Y2(1,:)
      write(*,*) 'B_Y2=',b_Y2(2,:)
      write(*,*) 'B_Y2=',b_Y2(3,:)
      write(*,*) 'B_Y2=',b_Y2(4,:)
      write(*,*) 'B_Y2=',b_Y2(5,:)
     !write(*,*) 'mulAA,mulAB,mulBB',mulAA,mulAB,mulBB
      write(*,*) 'mu_ave=',mu_ave(1,:)
      write(*,*) 'mu_ave=',mu_ave(2,:)
      write(*,*) 'mu_ave=',mu_ave(3,:)
      write(*,*) 'mu_ave=',mu_ave(4,:)
      write(*,*) 'mu_ave=',mu_ave(5,:)
#endif

      end

!***********************************************

      subroutine Kpermeability(KperTT,KperRR,RPP,NN)
      use tensors
      use SYS_property,only:mu_f
      use CONFIG,only:radii_test
      implicit none
      INTEGER,intent(in):: NN
      REAL*8,intent(in):: RPP(6*NN,6*NN)
      real(8),intent(out)::KperTT,KperRR
      
      REAL*8 RTT(3*NN,3*NN),RRR(3*NN,3*NN),T(3,3),R(3,3)
      integer I,J,m,n
  
      
      RTT=RPP(1:3*NN,1:3*NN)
      RRR=RPP(3*NN+1:6*NN,3*NN+1:6*NN)
      T=0.0_8
      R=0.0_8


      do I=1,NN
       do J=1,NN
         forall(m=1:3,n=1:3)
           T(m,n)=T(m,n)+RTT(m+3*(I-1),n+3*(J-1))
           R(m,n)=R(m,n)+RRR(m+3*(I-1),n+3*(J-1))
         end forall
       enddo
      enddo

      T=T/NN
      R=R/NN

   ! write(*,*) 'Rot=',R(1,:)
   ! write(*,*) 'Rot=',R(2,:)
   ! write(*,*) 'Rot=',R(3,:)

    KperTT=0.d0
    KperRR=0.d0
      DO I=1,3
         KperTT=KperTT+T(i,i)
         KperRR=KperRR+R(i,i)
      ENDDO

      KperTT=KperTT/3.d0/(6.0_8*pai*mu_f*radii_test)
      KperRR=KperRR/3.d0/(8.0_8*pai*mu_f*radii_test**3)

      end

!*****************************************************************
!**********************************************************
!**********************************************************
!**********************************************************


      SUBROUTINE MOB_TO_FRI(APP,APQ,AQQ,NN)
      IMPLICIT NONE
      INTEGER,intent(in):: NN
      REAL*8,intent(in):: APP(6*NN,6*NN)
      REAL*8,intent(inout):: APQ(6*NN,5*NN),AQQ(5*NN,5*NN)

      REAL*8 AQP(5*NN,6*NN)

      CALL MATREV(APP,6*NN)

      AQP = TRANSPOSE(APQ)
      APQ = MATMUL(APP,APQ)
      AQQ = AQQ + MATMUL(AQP,APQ)

      RETURN
      END


!***********************************************************
!***********************************************************
!***********************************************************


      SUBROUTINE FRI_TO_MOB(APP,APQ,AQQ,NN)
      IMPLICIT NONE
      INTEGER,intent(in):: NN
      REAL*8,intent(in):: APP(6*NN,6*NN)
      REAL*8,intent(inout):: APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8 AQP(5*NN,6*NN)

      CALL MATREV(APP,6*NN)
      AQP = TRANSPOSE(APQ)
      APQ = MATMUL(APP,APQ)
      AQQ = AQQ - MATMUL(AQP,APQ)

      RETURN
      END


!***********************************************************
!***********************************************************
!***********************************************************


      SUBROUTINE FRI_TO_MOB_RED(APP,APQ,AQQ,NN)
      IMPLICIT NONE
      INTEGER,intent(in):: NN
      REAL*8,intent(in):: APP(6*NN,6*NN)
      REAL*8,intent(inout):: AQQ(5*NN,5*NN),APQ(6*NN,5*NN)
      REAL*8 AQP(5*NN,6*NN)

      AQP = TRANSPOSE(APQ)
      APQ = MATMUL(APP,APQ)
      AQQ = AQQ - MATMUL(AQP,APQ)

      RETURN
      END


!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE INVFRI_TO_FRI(APP,APQ,AQQ,NN)
      IMPLICIT NONE
      INTEGER ,intent(in)::NN
      REAL*8 ,intent(inout)::APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8 AP(11*NN,11*NN)

      AP(11*NN,11*NN) = 0.D0

      AP(1:6*NN,1:6*NN) = APP
      AP(1:6*NN,6*NN+1:11*NN) = APQ 
      AP(6*NN+1:11*NN,1:6*NN) = TRANSPOSE(APQ)
      AP(6*NN+1:11*NN,6*NN+1:11*NN) = AQQ 

      CALL MATREV(APP,6*NN)
      CALL MATREV(AP,11*NN)

      APP = AP(1:6*NN,1:6*NN)
      APQ = AP(1:6*NN,6*NN+1:11*NN)
      AQQ = AP(6*NN+1:11*NN,6*NN+1:11*NN)

      RETURN
      END


!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE APQtoAQP(TPQ,TQP,NN)
      IMPLICIT NONE
      integer,intent(in):: NN
      REAL*8,intent(out):: TQP(5*NN,6*NN)
      REAL*8,intent(in):: TPQ(6*NN,5*NN)
      
      integer::i,j

      TQP = 0.D0

      DO I=1,6*NN
       DO J=1,5*NN
        TQP(J,I) = -TPQ(I,J)
       ENDDO
      ENDDO

      RETURN
      END

!******************************************

     subroutine MobToResist_FTS(APP,APQ,AQQ,NN,RPP,RPQ,RQQ)
      implicit none
      INTEGER,intent(in):: NN
      REAL*8,intent(in):: APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      real(8),intent(out)::RPP(6*NN,6*NN),RPQ(6*NN,5*NN),RQQ(5*NN,5*NN)
      
      REAL*8 Mtemp(11*NN,11*NN)

      Mtemp(1:6*NN,1:6*NN) =APP
      Mtemp(1:6*NN,6*NN+1:11*NN)=APQ
      Mtemp(6*NN+1:11*NN,1:6*NN)=transpose(APQ)
      Mtemp(6*NN+1:11*NN,6*NN+1:11*NN)=AQQ
      write(*,*) 'ok0-------------------'
      CALL MATREV(Mtemp,11*NN)
      write(*,*) 'ok2-------------------'
      RPP=Mtemp(1:6*NN,1:6*NN)
      RPQ=Mtemp(1:6*NN,6*NN+1:11*NN)
      RQQ=Mtemp(6*NN+1:11*NN,6*NN+1:11*NN)
      !RQQ=AQQ
      !CALL MATREV(RQQ,5*NN)

      end

      subroutine MobToResist_FT(APP,NN,RPP)
      implicit none
      INTEGER,intent(in):: NN
      REAL*8,intent(in):: APP(6*NN,6*NN)
      real(8),intent(out)::RPP(6*NN,6*NN)
      

      RPP=APP
      CALL MATREV(RPP,6*NN)

      end



      subroutine arrangement(RPP,RPQ,RQQ,mu,NN)
      implicit none
      INTEGER,intent(in):: NN
      REAL*8,intent(in):: RPP(6*NN,6*NN),RPQ(6*NN,5*NN),RQQ(5*NN,5*NN)
      real(8),intent(out)::mu(5*NN,5*NN)

      real*8 RPPINV(6*NN,6*NN)!,RQQTemp(5*NN,5*NN),muTemp(5,5)

      RPPINV=RPP
      CALL MATREV(RPPINV,6*NN)

      mu=matmul(matmul(transpose(RPQ),RPPINV),RPQ)-RQQ
      !mu=-RQQ
      end
!**********************************************