
! EVALUATES THE 11NNx11NN ROTNE-PRAGER MOBILITY MATRIX.
! DISPLAYS BLOKS:
! PP  6NNx6NN
! PQ  6NNx5NN
! QQ  5NNx5NN

!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE GRPERY_MOB(APP,APQ,AQQ,CONF,RADII,NN)
      use LATTICE_BASE,only:muCalculate
      IMPLICIT NONE
      INTEGER,intent(in):: NN
      real(8),intent(out)::APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8,intent(in):: CONF(3,NN),RADII(NN)

      real(8)::ARR(6*NN,6*NN)
      
      CALL GRPERY_INV_FRI(APP,APQ,AQQ,CONF,RADII,NN)

       if( muCalculate) then
        ARR = APP
        CALL INVFRI_TO_FRI(ARR,APQ,AQQ,NN)
        CALL FRI_TO_MOB_RED(APP,APQ,AQQ,NN)
       endif
      END

!***********************************************************
!***********************************************************
!***********************************************************


!******************************************************

      SUBROUTINE GRPERY_INV_FRI(APP,APQ,AQQ,p_pos,RADII,NN)
      USE TENSORS
      use stokesian_Green
      use SYS_property,only:mu_f
      use tensors,only:pai
      IMPLICIT NONE
      INTEGER,intent(in):: NN
      REAL*8,intent(inout):: APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8,intent(in)::p_pos(3,NN),RADII(NN)

      REAL*8 aaI,aaJ,grmobmx(11*NN,11*NN)
      REAL*8 C1PP(3,3),C1QQ(5,5)
      REAL*8 R(3),DMIN
      REAL*8 TRP(3,3)
      REAL*8 TTP(3,3),RRP(3,3),QQP_TR(5,5),GPQ_TR(3,5),HPQ_TR(3,5)
      INTEGER alpha,beta,I,J
      real*8 a(NN,NN,3,3),b(NN,NN,3,3),c(NN,NN,3,3)
      real*8 g(NN,NN,3,5),h(NN,NN,3,5),m(NN,NN,5,5)

      !CALL CALC_LATTICE_INVERSE(LATTICE,EWS_ERR)
      APP=0.D0
      APQ=0.D0
      AQQ=0.D0
      a=0.0_8
      b=0.0_8
      g=0.0_8
      c=0.0_8
      h=0.0_8
      m=0.0_8
!**************************************************************
!**************************************************************
!****************self interaction******************************
  self: DO alpha=1,NN
        beta=alpha
       CALL ROTNE_PRAGER_TT_SELF(C1PP,RADII(alpha))
        a(alpha,beta,1:3,1:3)=C1PP
        !b(alpha,beta,i,j)=0.0_cp
       CALL ROTNE_PRAGER_RR_SELF(C1PP,RADII(alpha))
        c(alpha,beta,1:3,1:3)=C1PP

       CALL ROTNE_PRAGER_DD_SELF_Y2(C1QQ,RADII(alpha))         
        m(alpha,beta,1:5,1:5) = C1QQ

      ENDDO self
!****************end self interaction**************************
!**************************************************************
!**************************************************************

!**************************************************************
!**************************************************************
!****************pair interaction******************************



! (U,Omega,-E) = - GrMob * (F,T,S) fluid force on particle
! (F,T,S) = -GrRe * (U,Omega,-E) invert
! E vector conversion
! EV1=E11-E33, EV2=2E12, EV3=2E13, EV4=2E23, EV5=E22-E33
! in the leftside, it should be -E, so the elements are -E11+E33, etc

    pairwisealpha: do alpha=1,NN
         pairwisebeta: do beta=alpha+1,NN
            aaI=RADII(alpha)
            aaJ=RADII(beta)
            R=p_pos(1:3,alpha)-p_pos(1:3,beta)
            
            
            DMIN=SQRT( SUM(R**2))
            !CALL PER_SKEW_CORRECTION(R,DMIN)
            IF(DMIN<=(aaI+aaJ)) THEN
              call YAMAKAWA_IJ(TTP,RRP,TRP,QQP_TR,GPQ_TR,HPQ_TR,R,aaI,aaJ)
            else
              call ROTNE_PRAGER_IJ(TTP,RRP,TRP,QQP_TR,GPQ_TR,HPQ_TR,R,aaI,aaJ)
            ENDIF



    !*************************************************************
           a(alpha,beta,1:3,1:3)=TTP
           b(alpha,beta,1:3,1:3)=TRP
           c(alpha,beta,1:3,1:3)=RRP
           
           forall(i=1:3,j=1:5)
           g(alpha,beta,i,j)=GPQ_TR(i,j)
           h(alpha,beta,i,j)=HPQ_TR(i,j)
           end forall

           m(alpha,beta,1:5,1:5)=QQP_TR

    !*************************************************************
            forall(i=1:3,j=1:3)
            a(beta,alpha,i,j)=a(alpha,beta,i,j)
            b(beta,alpha,i,j)=-b(alpha,beta,i,j)
            c(beta,alpha,i,j)=c(alpha,beta,i,j)
            end forall

            forall(i=1:3,j=1:5)
            g(beta,alpha,i,j)=-g(alpha,beta,i,j)
            h(beta,alpha,i,j)=h(alpha,beta,i,j)
            end forall

            forall(i=1:5,j=1:5)
            m(beta,alpha,i,j)=m(alpha,beta,i,j)
            end forall

            end do pairwisebeta
    end do pairwisealpha


    do alpha=1,NN
      do beta=1,NN
            do i=1,3
              do j=1,3
              !fill a,b,c
              grmobmx(3*(alpha-1)+i,3*(beta-1)+j)=a(alpha,beta,i,j)
              grmobmx(3*(alpha-1)+i,3*NN+3*(beta-1)+j)=b(alpha,beta,i,j)
              grmobmx(3*NN+3*(alpha-1)+i,3*NN+3*(beta-1)+j)=c(alpha,beta,i,j)
              !fill btilda by symmetry
              grmobmx(3*NN+3*(beta-1)+j,3*(alpha-1)+i)=b(alpha,beta,i,j)
              end do
            end do
            do i=1,3
              do j=1,5
              !fill g,h
              grmobmx(3*(alpha-1)+i,6*NN+5*(beta-1)+j)=g(alpha,beta,i,j)
              grmobmx(3*NN+3*(alpha-1)+i,6*NN+5*(beta-1)+j)=h(alpha,beta,i,j)
              !fill gtilda,htilda by symmetry
              grmobmx(6*NN+5*(beta-1)+j,3*(alpha-1)+i)=g(alpha,beta,i,j)
              grmobmx(6*NN+5*(beta-1)+j,3*NN+3*(alpha-1)+i)=h(alpha,beta,i,j)

              end do
            end do
            do i=1,5
              do j=1,5
              !fill in m
              grmobmx(6*NN+5*(alpha-1)+i,6*NN+5*(beta-1)+j)=m(alpha,beta,i,j)
              end do
            end do
      end do
    end do

    !open(unit=101,file='grmobmx',status='replace')
    !write(101,*) grmobmx
    !close(unit=101)
          grmobmx=grmobmx/(pai*mu_f)
          APP=grmobmx(1:6*NN,1:6*NN)
          APQ=grmobmx(1:6*NN,6*NN+1:11*NN)
          AQQ=grmobmx(6*NN+1:11*NN,6*NN+1:11*NN)

  END




!***********************************************************

      SUBROUTINE ROTNE_PRAGER_TT_SELF(A1,aaI)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: A1(3,3)
      real(8),intent(in)::aaI        ! radius

      A1=1.D0/(6.0_8*aai)*KDirac

      END


      SUBROUTINE ROTNE_PRAGER_RR_SELF(A1,aai)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: A1(3,3)
      REAL*8,intent(in)::aai
      A1=1.0_8/(8.0_8*aai**3)*KDirac

      END

      SUBROUTINE ROTNE_PRAGER_DD_SELF_Y2(B1,aaI)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: B1(5,5)
      REAL*8,intent(in):: aaI       ! radiu      
     
      REAL*8 B1_CART(3,3,3,3)


      CALL ROTNE_PRAGER_DD_SELF(B1_CART,aaI)

      Call DD_transY2(B1,B1_CART)

      END

      SUBROUTINE ROTNE_PRAGER_DD_SELF(B1,aaI)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: B1(3,3,3,3)
      REAL*8,intent(in):: aaI          ! radiu  


      B1=3.0_8/(20.0_8*aai**3)*II

      END






!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE YAMAKAWA_TT_IJ(A1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: A1(3,3)
      REAL*8,intent(in):: aaI,aaJ,R(3)           ! radius
      
      REAL*8 DIST,DIST2,DIST3,RW(3),RR(3,3)
      REAL*8 M0TT          ! tt single sphere mobility 
      PARAMETER(M0TT=1.D0/6.D0)
      REAL*8 PRE,PU,PRR,aaIaaJ

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST

      CALL CALC_RR(RR,RW)

      IF (aaI.LE.aaJ) THEN
       IF (DIST.LE.(aaJ-aaI)) THEN
        A1=M0TT/aaJ*KDirac
       ELSE
        DIST2=DIST**2
        DIST3=DIST**3
        aaIaaJ=aaI-aaJ
        PRE=1.D0/(6*32*aaI*aaJ)
        PU=16*(aaI+aaJ)-(aaIaaJ**2+3*DIST2)**2/DIST3
        PRR=3*(aaIaaJ**2 - DIST2)**2/DIST3
        A1=PRE*(PU*KDirac+PRR*RR)
       ENDIF
      ELSE
       IF (DIST.LE.(aaI-aaJ)) THEN
        A1=M0TT/aaI*KDirac
       ELSE
        DIST2=DIST**2
        DIST3=DIST**3
        aaIaaJ=aaI-aaJ
        PRE=1.D0/(6*32*aaI*aaJ)
        PU=16*(aaI+aaJ)-(aaIaaJ**2+3*DIST2)**2/DIST3
        PRR=3*(aaIaaJ**2 - DIST2)**2/DIST3
        A1=PRE*(PU*KDirac+PRR*RR)
       ENDIF       
      ENDIF

      END
!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE YAMAKAWA_RR_IJ(A1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: A1(3,3)
      REAL*8,intent(in):: aaI,aaJ ,R(3)          ! radiu
      
      REAL*8 DIST,DIST2,DIST3,RW(3),RR(3,3)
      REAL*8 M0RR          ! rr single sphere mobility 
      PARAMETER(M0RR=1.D0/8.D0)
      REAL*8 PRE,PU,PRR,aaIaaJ,aaI2,aaJ2,aaI3,aaJ3

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST

      CALL CALC_RR(RR,RW)

      IF (aaI.LE.aaJ) THEN
       IF (DIST.LE.(aaJ-aaI)) THEN
        A1=M0RR/aaJ**3*KDirac
       ELSE
        DIST2=DIST**2
        DIST3=DIST**3
        aaIaaJ=aaI-aaJ
        aaI2=aaI**2
        aaJ2=aaJ**2
        aaI3=aaI**3
        aaJ3=aaJ**3
        PRE=1.D0/(8*64*aaI3*aaJ3)

        PU=  5*DIST3 - 27*DIST*(aaI2+aaJ2)  &
          + 32*(aaI3+aaJ3) &
          - 9*(aaI2-aaJ2)**2/DIST  &
          - aaIaaJ**4*(aaI2 + 4*aaI*aaJ + aaJ2)/DIST3 

        PRR=3*(aaIaaJ**2 - DIST2)**2 &
         *((aaI2 + 4*aaI*aaJ + aaJ2)-DIST2)/DIST3

        A1=PRE*(PU*KDirac+PRR*RR)
       ENDIF
      ELSE
       IF (DIST.LE.(aaI-aaJ)) THEN
        A1=M0RR/aaI**3*KDirac
       ELSE
        DIST2=DIST**2
        DIST3=DIST**3
        aaIaaJ=aaI-aaJ
        aaI2=aaI**2
        aaJ2=aaJ**2
!-----------------------------------------------------------
        aaI3=aaI**3
        aaJ3=aaJ**3
!-----------------------------------------------------------
        PRE=1.D0/(8*64*(aaI**3)*(aaJ**3))

        PU=  5*DIST3 - 27*DIST*(aaI2+aaJ2)  &
          + 32*(aaI3+aaJ3) &
          - 9*(aaI2-aaJ2)**2/DIST  &
          - aaIaaJ**4*(aaI2 + 4*aaI*aaJ + aaJ2)/DIST3

        PRR=3*(aaIaaJ**2 - DIST2)**2 &
         *((aaI2 + 4*aaI*aaJ + aaJ2)-DIST2)/DIST3

        A1=PRE*(PU*KDirac+PRR*RR)
       ENDIF       
      ENDIF

      END
!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE YAMAKAWA_TD_IJ_Y2(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: C1(3,5)
      REAL*8,intent(in):: aaI,aaJ ,R(3)          ! radiu 

      REAL*8 C1_CART(3,3,3)


      CALL YAMAKAWA_TD_IJ(C1_CART,R,aaI,aaJ)


      !if(front)then
      !  call RTD_front_transY2(C1,C1_CART)
      !else
        call RTD_back_transY2(C1,C1_CART)
      !endif

      END
!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE YAMAKAWA_RT_IJ(A1,R,aaI,aaJ)
      IMPLICIT NONE
      REAL*8,intent(out):: A1(3,3)
      REAL*8,intent(in):: aaI,aaJ ,R(3)          ! radiu
      
      REAL*8 DIST,RW(3),EPSR(3,3)
      REAL*8 PRE,PEPSR

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST

      CALL CALC_EPSR(EPSR,RW)

      IF (aaI.LE.aaJ) THEN
       IF (DIST.LE.(aaJ-aaI)) THEN
        A1=0.D0
       ELSE
        PRE=1.D0/(16*8*aaI**3*aaJ)
        PEPSR=((aaI-aaJ)+DIST)**2 &
        *(aaJ**2+2*aaJ*(aaI+DIST)-3*(aaI-DIST)**2) &
        /DIST**2
        A1=PRE*PEPSR*EPSR
       ENDIF
      ELSE
       IF (DIST.LE.(aaI-aaJ)) THEN
        A1=DIST/(8*(aaI**3))*EPSR
       ELSE
        PRE=1.D0/(16*8*aaI**3*aaJ)
        PEPSR=((aaI-aaJ)+DIST)**2 &
        *(aaJ**2+2*aaJ*(aaI+DIST)-3*(aaI-DIST)**2) &
        /DIST**2
        A1=PRE*PEPSR*EPSR
       ENDIF       
      ENDIF
      
      RETURN
      END

!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE YAMAKAWA_RD_IJ(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: C1(3,3,3)
      REAL*8,intent(in):: aaI,aaJ ,R(3)          ! radiu 

      REAL*8 t3EPSRR(3,3,3)
      REAL*8 DIST,DIST2,RW(3)

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST2=DIST**2

      CALL CALC_EPSRR(t3EPSRR,RW)

      C1=0.D0

      IF (aaI.LE.aaJ) THEN
       IF (DIST.LE.(aaJ-aaI)) THEN
        C1=0.D0
       ELSE

        C1 = -3.D0/(256.D0*(aaI**3)*(aaJ**3)*DIST**3) &
             *((aaI-aaJ)**2 - DIST2)**2 &
             *((aaI**2 + 4*aaI*aaJ + aaJ**2)-DIST2) &
             *t3EPSRR

       ENDIF
      ELSE
       IF (DIST.LE.(aaI-aaJ)) THEN
        C1=0.D0
       ELSE

        C1 = -3.D0/(256.D0*(aaI**3)*(aaJ**3)*DIST**3) &
             *((aaI-aaJ)**2 - DIST2)**2 &
             *((aaI**2 + 4*aaI*aaJ + aaJ**2)-DIST2) &
             *t3EPSRR

       ENDIF       
      ENDIF

      END

!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE YAMAKAWA_RD_IJ_Y2(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: C1(3,5)
      REAL*8,intent(in):: aaI,aaJ ,R(3)          ! radiu 

      REAL*8 C1_CART(3,3,3)


      CALL YAMAKAWA_RD_IJ(C1_CART,R,aaI,aaJ)

      !if(front)then
       ! call RTD_front_transY2(C1,C1_CART)
      !else
        call RTD_back_transY2(C1,C1_CART)
      !endif

      END
!***********************************************************
!***********************************************************
!***********************************************************


      SUBROUTINE YAMAKAWA_TD_IJ(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: C1(3,3,3)
      REAL*8,intent(in):: aaI,aaJ ,R(3)          ! radiu 

      REAL*8 t3UR(3,3,3),t3RRR(3,3,3)
      REAL*8 DIST,DIST2,DIST4,RW(3),PRE,PUR,PRRR

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST2=DIST**2
      DIST4=DIST**4

      CALL CALC_UR(t3UR,RW)
      CALL CALC_RRR(t3RRR,RW)

      C1=0.D0

      IF (aaI.LE.aaJ) THEN
       IF (DIST.LE.(aaJ-aaI)) THEN
        C1=-3.D0*DIST/(20.D0*aaJ**3)*t3UR
       ELSE

        PRE = 1.D0/(aaI*aaJ**3)

        PUR = (10.D0*DIST2 - 24.D0*aaI*DIST &
             + 15.D0*(aaI**2-aaJ**2)  &
             - (aaI-aaJ)**5*(aaI+5.D0*aaJ)/DIST4 )/320.D0

        PRRR = ((aaI-aaJ)**2-DIST2)**2 &
              *((aaI-aaJ)*(aaI+5*aaJ)-DIST2) &
              /(128.D0*DIST4)

        C1 = PRE*(PUR*t3UR+PRRR*t3RRR)

       ENDIF
      ELSE
       IF (DIST.LE.(aaI-aaJ)) THEN
        C1=0.D0
       ELSE

        PRE = 1.D0/(aaI*aaJ**3)

        PUR = (10.D0*DIST2 - 24.D0*aaI*DIST &
             + 15.D0*(aaI**2-aaJ**2)  &
             - (aaI-aaJ)**5*(aaI+5.D0*aaJ)/DIST4 )/320.D0

        PRRR = ((aaI-aaJ)**2-DIST2)**2 &
              *((aaI-aaJ)*(aaI+5*aaJ)-DIST2) &
              /(128.D0*DIST4)

        C1 = PRE*(PUR*t3UR+PRRR*t3RRR)

       ENDIF       
      ENDIF

      END

!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE YAMAKAWA_DD_IJ_Y2(B1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: B1(5,5)
      REAL*8,intent(in):: aaI,aaJ ,R(3)          ! radiu      
      
      REAL*8 B1_CART(3,3,3,3)

      CALL YAMAKAWA_DD_IJ(B1_CART,R,aaI,aaJ)

      Call DD_transY2(B1,B1_CART)


      RETURN
      END

!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE YAMAKAWA_DD_IJ(B1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: B1(3,3,3,3)
      REAL*8,intent(in):: aaI,aaJ ,R(3)          ! radiu 

      REAL*8 t4D1(3,3,3,3),t4D2(3,3,3,3),t4D0(3,3,3,3)
      REAL*8 DIST,DIST3,DIST5,RW(3),RR(3,3)
      REAL*8 aaI2,aaJ2,aaI3,aaJ3

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST5=DIST**5
      DIST3=DIST**3
      aaI2=aaI**2
      aaJ2=aaJ**2
      aaI3=aaI**3
      aaJ3=aaJ**3

      CALL CALC_RR(RR,RW)

      CALL CALC_t4D0(t4D0,RR)
      CALL CALC_t4D1(t4D1,RR)
      CALL CALC_t4D2(t4D2,RR)
      t4D2 = t4D2 - t4D1

      B1=0.D0

      IF (aaI.LE.aaJ) THEN
       IF (DIST.LE.(aaJ-aaI)) THEN
        B1= 3.D0/(aaJ**3*20.D0)*(t4D0+t4D1+t4D2)
       ELSE

      B1=B1+3.D0/(1280.D0*aaI3*aaJ3) &
          *( &
          +3.D0*(aaI-aaJ)**6*(aaI2+6.D0*aaI*aaJ+aaJ2)/DIST5 &
          -10.D0*(aaI-aaJ)**4*(aaI2+4.D0*aaI*aaJ+aaJ2)/DIST3 &
          +32.D0*(aaI3+aaJ3) &
          -30.D0*(aaI2+aaJ2)*DIST &
          +5.D0*DIST3 &
           )*t4D0

      B1=B1+3.D0/(1280.D0*aaI3*aaJ3) &
           *( &
           -2.D0*(aaI-aaJ)**6*(aaI2+6.D0*aaI*aaJ+aaJ2)/DIST5 &
           +5.D0*(aaI-aaJ)**4*(aaI2+4.D0*aaI*aaJ+aaJ2)/DIST3 &
           -15.D0*(aaI2-aaJ2)**2/DIST &
           +32.D0*(aaI3+aaJ3) &
           -25.D0*(aaI2+aaJ2)*DIST &
           +5.D0*DIST3 &
           )*t4D1

      B1=B1+3.D0/(2560.D0*aaI3*aaJ3) &
           *( &
           +(aaI-aaJ)**6*(aaI2+6.D0*aaI*aaJ+aaJ2)/DIST5 &
           -30.D0*(aaI2-aaJ2)**2/DIST &
           +64.D0*(aaI3+aaJ3) &
           -40.D0*(aaI2+aaJ2)*DIST &
           +5.D0*DIST3 &
           )*t4D2

       ENDIF
      ELSE
       IF (DIST.LE.(aaI-aaJ)) THEN
        B1= 3.D0/(aaI**3*20.D0)*(t4D0+t4D1+t4D2)
       ELSE

      B1=B1+3.D0/(1280.D0*aaI3*aaJ3) &
           *( &
           +3.D0*(aaI-aaJ)**6*(aaI2+6.D0*aaI*aaJ+aaJ2)/DIST5 &
           -10.D0*(aaI-aaJ)**4*(aaI2+4.D0*aaI*aaJ+aaJ2)/DIST3 &
           +32.D0*(aaI3+aaJ3) &
           -30.D0*(aaI2+aaJ2)*DIST &
           +5.D0*DIST3 &
            )*t4D0

      B1=B1+3.D0/(1280.D0*aaI3*aaJ3) &
           *( &
           -2.D0*(aaI-aaJ)**6*(aaI2+6.D0*aaI*aaJ+aaJ2)/DIST5 &
           +5.D0*(aaI-aaJ)**4*(aaI2+4.D0*aaI*aaJ+aaJ2)/DIST3 &
           -15.D0*(aaI2-aaJ2)**2/DIST &
           +32.D0*(aaI3+aaJ3) &
           -25.D0*(aaI2+aaJ2)*DIST &
           +5.D0*DIST3 &
           )*t4D1

      B1=B1+3.D0/(2560.D0*aaI3*aaJ3) &
           *( &
           +(aaI-aaJ)**6*(aaI2+6.D0*aaI*aaJ+aaJ2)/DIST5 &
           -30.D0*(aaI2-aaJ2)**2/DIST &
           +64.D0*(aaI3+aaJ3) &
           -40.D0*(aaI2+aaJ2)*DIST &
           +5.D0*DIST3 &
           )*t4D2

       ENDIF       
      ENDIF


      END
!**********************************************************


     SUBROUTINE YAMAKAWA_IJ(TTPY,RRPY,TRPY,QQPY_TR,GPQY_TR,HPQY_TR,R,aaI,aaJ)
      IMPLICIT NONE
      REAL*8,intent(out):: TTPY(3,3),RRPY(3,3),TRPY(3,3),QQPY_TR(5,5)
      real*8,intent(out):: GPQY_TR(3,5),HPQY_TR(3,5)
      REAL*8,intent(in):: aaI,aaJ,R(3)      ! radius

      real*8 QQPY(3,3,3,3),GPQY(3,3,3),HPQY(3,3,3),RTPY0(3,3)
      CALL YAMAKAWA_TT_IJ(TTPY,R,aaI,aaJ)
      CALL YAMAKAWA_RR_IJ(RRPY,R,aaI,aaJ)
      CALL YAMAKAWA_RT_IJ(RTPY0,-R,aaJ,aaI)
      TRPY=transpose(RTPY0)
      CALL YAMAKAWA_TD_IJ(GPQY,R,aaI,aaJ)
      CALL YAMAKAWA_RD_IJ(HPQY,R,aaI,aaJ)
      CALL YAMAKAWA_DD_IJ(QQPY,R,aaI,aaJ)
      call GH_TR(GPQY_TR,GPQY)
      call GH_TR(HPQY_TR,HPQY)
      call QQ_TR(QQPY_TR,QQPY)

      end SUBROUTINE YAMAKAWA_IJ


     SUBROUTINE overlap_Yamakawa_CORRECTION_IJ(TTPC,RRPC,TRPC,QQPC_TR,GPQC_TR,HPQC_TR,R,aaI,aaJ)
     use stokesian_Green
     IMPLICIT NONE
     REAL*8,intent(in):: aaI,aaJ,R(3)      ! radius
     REAL*8,intent(out):: TTPC(3,3),RRPC(3,3),TRPC(3,3),QQPC_TR(5,5)
     real*8,intent(out):: GPQC_TR(3,5),HPQC_TR(3,5)
     
     REAL*8:: TTPY(3,3),RRPY(3,3),TRPY(3,3),QQPY_TR(5,5)
     real*8:: GPQY_TR(3,5),HPQY_TR(3,5)
     REAL*8:: TTP(3,3),RRP(3,3),TRP(3,3),QQP_TR(5,5)
     real*8:: GPQ_TR(3,5),HPQ_TR(3,5)

     call YAMAKAWA_IJ(TTPY,RRPY,TRPY,QQPY_TR,GPQY_TR,HPQY_TR,R,aaI,aaJ)
     call ROTNE_PRAGER_IJ(TTP,RRP,TRP,QQP_TR,GPQ_TR,HPQ_TR,R,aaI,aaJ)

     TTPC=TTPY-TTP
     RRPC=RRPY-RRP
     TRPC=TRPY-TRP
     QQPC_TR=QQPY_TR-QQP_TR
     GPQC_TR=GPQY_TR-GPQ_TR
     HPQC_TR=HPQY_TR-HPQ_TR

     end SUBROUTINE overlap_Yamakawa_CORRECTION_IJ

     SUBROUTINE overlap_Regular_CORRECTION_IJ(TTPC,RRPC,TRPC,QQPC_TR,GPQC_TR,HPQC_TR,R,aaI,aaJ)
     use stokesian_Green
     use Regularization
     IMPLICIT NONE
     REAL*8,intent(in):: aaI,aaJ,R(3)      ! radius
     REAL*8,intent(out):: TTPC(3,3),RRPC(3,3),TRPC(3,3),QQPC_TR(5,5)
     real*8,intent(out):: GPQC_TR(3,5),HPQC_TR(3,5)
     
     REAL*8:: TTPY(3,3),RRPY(3,3),TRPY(3,3),QQPY_TR(5,5)
     real*8:: GPQY_TR(3,5),HPQY_TR(3,5)
     REAL*8:: TTP(3,3),RRP(3,3),TRP(3,3),QQP_TR(5,5)
     real*8:: GPQ_TR(3,5),HPQ_TR(3,5)

     call Regular_ROTNE_PRAGER_IJ(TTPY,RRPY,TRPY,QQPY_TR,GPQY_TR,HPQY_TR,R,aaI,aaJ)
     call ROTNE_PRAGER_IJ(TTP,RRP,TRP,QQP_TR,GPQ_TR,HPQ_TR,R,aaI,aaJ)

     TTPC=TTPY-TTP
     RRPC=RRPY-RRP
     TRPC=TRPY-TRP
     QQPC_TR=QQPY_TR-QQP_TR
     GPQC_TR=GPQY_TR-GPQ_TR
     HPQC_TR=HPQY_TR-HPQ_TR

     end SUBROUTINE overlap_Regular_CORRECTION_IJ