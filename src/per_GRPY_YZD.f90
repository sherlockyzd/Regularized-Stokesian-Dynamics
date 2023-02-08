
! EVALUATES THE 11NNx11NN ****PERIOD**** ROTNE-PRAGER MOBILITY MATRIX.
! DISPLAYS BLOKS:
! PP  6NNx6NN
! PQ  6NNx5NN
! QQ  5NNx5NN

!***********************************************************
!***********************************************************
!***********************************************************
!***************************************************************************

    SUBROUTINE PER_GRPERY_INV_FRI(APP,APQ,AQQ,p_pos,RADII,NN,LATTICE,EWS_ERR)
      USE TENSORS
      use Ewald_summation
      use prutil,only:cp
      use SYS_property,only:mu_f
      use method,only:correction_method
      USE OMP_LIB
      use period_bdy_tools,only:CALC_LATTICE_INVERSE,PER_SKEW_CORRECTION

      IMPLICIT NONE
      INTEGER,intent(in):: NN
      REAL*8,intent(out):: APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8,intent(in)::  p_pos(3,NN),RADII(NN),LATTICE(3,3),EWS_ERR

      REAL*8 aaI,aaJ,grmobmx(11*NN,11*NN)
      REAL*8 TTS(3,3),RRS(3,3),QQS_TR(5,5)
      REAL*8 TTP(3,3),RRP(3,3),QQP_TR(5,5),GPQ_TR(3,5),HPQ_TR(3,5),TRP(3,3)
      real*8 TTPC(3,3),RRPC(3,3),TRPC(3,3),QQPC_TR(5,5)
      real*8 GPQC_TR(3,5),HPQC_TR(3,5)
      REAL*8 R(3),DMIN
      INTEGER alpha,beta,I,J,thread_id
      real*8 a(NN,NN,3,3),b(NN,NN,3,3),c(NN,NN,3,3)
      real*8 g(NN,NN,3,5),h(NN,NN,3,5),m(NN,NN,5,5)


      CALL CALC_LATTICE_INVERSE(LATTICE,EWS_ERR)
      APP=0.D0
      APQ=0.D0
      AQQ=0.D0
      a=0.0_cp
      b=0.0_cp
      g=0.0_cp
      c=0.0_cp
      h=0.0_cp
      m=0.0_cp
!**************************************************************
!**************************************************************
!****************self interaction******************************

!$omp parallel default(NONE)shared(a,c,m,NN,Radii,EWS_ERR)&
!$omp private(alpha,beta,TTS,RRS,QQS_TR,thread_id)
!$omp do schedule(static)
  self: DO alpha=1,NN
    !thread_id=OMP_GET_THREAD_NUM()
    !write(*,*) 'thread_id=',thread_id
        beta=alpha
       CALL PER_ROTNE_PRAGER_SELF(TTS,RRS,QQS_TR,RADII(alpha),EWS_ERR)
!!$omp critical
!!$omp flush(a,c,m)
        a(alpha,beta,:,:)=TTS
        c(alpha,beta,:,:)=RRS
        m(alpha,beta,:,:)=QQS_TR
!!$omp end critical
      ENDDO self
!$omp end do
!$omp end parallel 

!****************end self interaction**************************

!****************pair interaction******************************



! (U,Omega,-E) = - GrMob * (F,T,S) fluid force on particle
! (F,T,S) = -GrRe * (U,Omega,-E) invert
! E vector conversion
! EV1=E11-E33, EV2=2E12, EV3=2E13, EV4=2E23, EV5=E22-E33
! in the leftside, it should be -E, so the elements are -E11+E33, etc


!!!!!!$omp reduction( +:)
!$omp parallel default(NONE)shared(a,b,c,g,h,m,NN,RADII,p_pos,correction_method,EWS_ERR)&
!$omp private(alpha,beta,i,j,TTPC,RRPC,TRPC,QQPC_TR,GPQC_TR,HPQC_TR,TTP,RRP,TRP,QQP_TR,&
!$omp GPQ_TR,HPQ_TR,R,DMIN,aaI,aaJ)
!$omp do schedule(guided)
pairwisealpha: do alpha=1,NN
     pairwisebeta: do beta=alpha+1,NN
        aaI=RADII(alpha)
        aaJ=RADII(beta)
        R=p_pos(1:3,alpha)-p_pos(1:3,beta)
        DMIN=SQRT( SUM(R**2))
        CALL PER_SKEW_CORRECTION(R,DMIN)
        
        CALL PER_ROTNE_PRAGER_IJ(TTP,RRP,TRP,QQP_TR,GPQ_TR,HPQ_TR,R,aaI,aaJ,EWS_ERR)
       
        if(correction_method.ne.0) then
          IF(DMIN<=(aaI+aaJ)) THEN
            if(correction_method.eq.1) then
              call overlap_Yamakawa_CORRECTION_IJ(TTPC,RRPC,TRPC,QQPC_TR,GPQC_TR,HPQC_TR,R,aaI,aaJ) 
              write(*,*) 'make overlap_Yamakawa_CORRECTION'
            elseif(correction_method.eq.2)then
              call overlap_Regular_CORRECTION_IJ(TTPC,RRPC,TRPC,QQPC_TR,GPQC_TR,HPQC_TR,R,aaI,aaJ) 
              write(*,*) 'make overlap_Regularization_CORRECTION'
            endif

           TTP=TTP+TTPC
           RRP=RRP+RRPC
           TRP=TRP+TRPC
           QQP_TR=QQP_TR+QQPC_TR
           GPQ_TR=GPQ_TR+GPQC_TR
           HPQ_TR=HPQ_TR+HPQC_TR

          ENDIF
        endif



!*************************************************************
!!$omp critical
!!$omp flush(a,b,c,g,h,m)
      forall(i=1:3,j=1:3)
       a(alpha,beta,i,j)=TTP(i,j)
       b(alpha,beta,i,j)=TRP(i,j)
       c(alpha,beta,i,j)=RRP(i,j)
      end forall

       forall(i=1:3,j=1:5)
       g(alpha,beta,i,j)=GPQ_TR(i,j)
       h(alpha,beta,i,j)=HPQ_TR(i,j)
       end forall
       forall(i=1:5,j=1:5)
       m(alpha,beta,i,j)=QQP_TR(i,j)
       end forall

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
!!$omp end critical
        end do pairwisebeta
end do pairwisealpha
!$omp end do
!$omp end parallel 



!$omp parallel default(NONE)shared(grmobmx,a,b,c,g,h,m,NN)&
!$omp private(alpha,beta,i,j)
!$OMP DO collapse(2)
do alpha=1,NN
  do beta=1,NN
      forall(i=1:3,j=1:3)
          !fill a,b,c
          grmobmx(3*(alpha-1)+i,3*(beta-1)+j)=a(alpha,beta,i,j)
          grmobmx(3*(alpha-1)+i,3*NN+3*(beta-1)+j)=b(alpha,beta,i,j)
          grmobmx(3*NN+3*(alpha-1)+i,3*NN+3*(beta-1)+j)=c(alpha,beta,i,j)
          !fill btilda by symmetry
          grmobmx(3*NN+3*(beta-1)+j,3*(alpha-1)+i)=b(alpha,beta,i,j)
      end forall
      forall(i=1:3,j=1:5)
          !fill g,h
          grmobmx(3*(alpha-1)+i,6*NN+5*(beta-1)+j)=g(alpha,beta,i,j)
          grmobmx(3*NN+3*(alpha-1)+i,6*NN+5*(beta-1)+j)=h(alpha,beta,i,j)
          !fill gtilda,htilda by symmetry
          grmobmx(6*NN+5*(beta-1)+j,3*(alpha-1)+i)=g(alpha,beta,i,j)
          grmobmx(6*NN+5*(beta-1)+j,3*NN+3*(alpha-1)+i)=h(alpha,beta,i,j)
      end forall
      forall(i=1:5,j=1:5)
          !fill in m
          grmobmx(6*NN+5*(alpha-1)+i,6*NN+5*(beta-1)+j)=m(alpha,beta,i,j)
      end forall
  end do
end do
!$omp end do
!$omp end parallel


!open(unit=101,file='grmobmx',status='replace')
!write(101,*) grmobmx
!close(unit=101)
      
      grmobmx=grmobmx/(pai*mu_f)
      APP=grmobmx(1:6*NN,1:6*NN)
      APQ=grmobmx(1:6*NN,6*NN+1:11*NN)
      AQQ=grmobmx(6*NN+1:11*NN,6*NN+1:11*NN)

  END


!****************end pair interaction**************************
!**************************************************************
!**************************************************************



