
! EVALUATES THE 11NNx11NN ROTNE-PRAGER MOBILITY MATRIX.
! DISPLAYS BLOKS:
! PP  6NNx6NN
! PQ  6NNx5NN
! QQ  5NNx5NN

!***********************************************************
!***********************************************************
!**********************************************************
!******************************************************


  Module Wall_interaction
  USE TENSORS
  use stokesian_Green
  use SYS_property,only:mu_f
  use tensors,only:pai
  !use tensors,only:kdirac,EPS,pip5,IIiso,pai,II
  !use LATTICE_BASE,only:lammda,phi,LI,LR,LV,NR,NI
      
   IMPLICIT NONE
   private

   public:: Wall_GRPERY_INV_FRI

   contains

      SUBROUTINE Wall_GRPERY_INV_FRI(APP,p_pos,RADII,NN)
      USE TENSORS
      use stokesian_Green
      use SYS_property,only:mu_f
      use tensors,only:pai
      IMPLICIT NONE
      INTEGER,intent(in):: NN
      REAL*8,intent(inout):: APP(3*NN,3*NN)
      REAL*8,intent(in)::p_pos(3,NN),RADII(NN)

      REAL*8 aaI,aaJ

      REAL*8 R(3),DMIN,RL(3)

      REAL*8 TTP(3,3)
      INTEGER alpha,beta,I,J
      real*8 a(NN,NN,3,3)

      !CALL CALC_LATTICE_INVERSE(LATTICE,EWS_ERR)
      APP=0.D0

      a=0.0_8

!**************************************************************
!**************************************************************
!****************self interaction******************************
  self: DO alpha=1,NN
        if(p_pos(3,alpha)<=10.0) then
          beta=alpha
          aaI=RADII(alpha)
          aaJ=RADII(beta)
          R=0.0_8
          R(3)=2.0_8*p_pos(3,alpha)
          
          
          DMIN=SQRT( SUM(R**2))
          !CALL PER_SKEW_CORRECTION(R,DMIN)
          !IF(DMIN<=(aaI+aaJ)) THEN
          !  call YAMAKAWA_IJ(TTP,RRP,TRP,QQP_TR,GPQ_TR,HPQ_TR,R,aaI,aaJ)
          !else
          call ROTNE_PRAGER_TT_IJ(TTP,R,aaI,aaJ)
            !call ROTNE_PRAGER_IJ(TTP,RRP,TRP,QQP_TR,GPQ_TR,HPQ_TR,R,aaI,aaJ)
              !ENDIF
          call Wall_P_Right_tensor(TTP)
      !*************************************************************
             a(alpha,beta,1:3,1:3)=TTP

         !CALL ROTNE_PRAGER_TT_SELF(C1PP,RADII(alpha))
         ! a(alpha,beta,1:3,1:3)=C1PP
          !b(alpha,beta,i,j)=0.0_cp
         !CALL ROTNE_PRAGER_RR_SELF(C1PP,RADII(alpha))
          !c(alpha,beta,1:3,1:3)=C1PP

         !CALL ROTNE_PRAGER_DD_SELF_Y2(C1QQ,RADII(alpha))         
         ! m(alpha,beta,1:5,1:5) = C1QQ
      ENDIF
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
          if(p_pos(3,alpha)<=10.0.or.p_pos(3,beta)<=10.0)then
            call Wall_P_Left_vector(p_pos(1:3,beta),RL)
            R=p_pos(1:3,alpha)-RL
            
            
            DMIN=SQRT( SUM(R**2))
            !CALL PER_SKEW_CORRECTION(R,DMIN)
            !IF(DMIN<=(aaI+aaJ)) THEN
            !  call YAMAKAWA_IJ(TTP,RRP,TRP,QQP_TR,GPQ_TR,HPQ_TR,R,aaI,aaJ)
            !else
            call ROTNE_PRAGER_TT_IJ(TTP,R,aaI,aaJ)
            !ENDIF
            call Wall_P_Right_tensor(TTP)



    !*************************************************************
           a(alpha,beta,1:3,1:3)=TTP

            forall(i=1:3,j=1:3)
            a(beta,alpha,i,j)=a(alpha,beta,i,j)
            end forall
          endif
          end do pairwisebeta
    end do pairwisealpha


    do alpha=1,NN
      do beta=1,NN
            do i=1,3
              do j=1,3
              !fill a,b,c
              APP(3*(alpha-1)+i,3*(beta-1)+j)=a(alpha,beta,i,j)
              end do
            end do
      end do
    end do


          APP=APP/(pai*mu_f)


  END SUBROUTINE Wall_GRPERY_INV_FRI

    SUBROUTINE Wall_P_Left_vector(R,RL)
      IMPLICIT NONE
      REAL*8,intent(in):: R(3)
      REAL*8,intent(out)::RL(3)
      
      !REAL*8 FLJ(3*NN)    ! STRECHING and BENDING FORCES
      !INTEGER I,J
      RL=R
      RL(3)=-R(3)

    END

    SUBROUTINE Wall_P_Left_tensor(R)
      IMPLICIT NONE
      REAL*8,intent(inout):: R(3,3)
      !REAL*8,intent(out)::RL(3,3)
      
      !REAL*8 FLJ(3*NN)    ! STRECHING and BENDING FORCES
      !INTEGER I,J
      !RL=R
      R(3,:)=-R(3,:)

    END

    SUBROUTINE Wall_P_Right_tensor(R)
      IMPLICIT NONE
      REAL*8,intent(inout):: R(3,3)
      !REAL*8,intent(out)::RR(3,3)
      
      !REAL*8 FLJ(3*NN)    ! STRECHING and BENDING FORCES
      !INTEGER I,J
      !RR=R
      R(:,3)=-R(:,3)

    END

  End Module Wall_interaction