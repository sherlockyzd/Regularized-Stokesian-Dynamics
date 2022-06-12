
! EVALUATES THE 11NNx11NN ROTNE-PRAGER MOBILITY MATRIX.
! DISPLAYS BLOKS:
! PP  6NNx6NN
! PQ  6NNx5NN
! QQ  5NNx5NN

!******************************************************

      SUBROUTINE GRPERY_MOB(grmobmx,p_pos,RADII)
      IMPLICIT NONE
      REAL*8,intent(in)::p_pos(3,2),RADII(2)
      REAL*8,intent(out):: grmobmx(22,22)

      REAL*8 aaI,aaJ
      REAL*8 C1PP(3,3),C1QQ(5,5)
      REAL*8 R(3),DMIN
      REAL*8 GPQ(3,3,3),HPQ(3,3,3),TRP(3,3)
      REAL*8 TTP(3,3),RRP(3,3),QQP(3,3,3,3),QQP_TR(5,5),GPQ_TR(3,5),HPQ_TR(3,5)
      INTEGER alpha,beta,I,J,k,l,NN
      real*8 a(2,2,3,3),b(2,2,3,3),c(2,2,3,3)
      real*8 g(2,2,3,5),h(2,2,3,5),m(2,2,5,5)

      !CALL CALC_LATTICE_INVERSE(LATTICE,EWS_ERR)
      a=0.0_8
      b=0.0_8
      g=0.0_8
      c=0.0_8
      h=0.0_8
      m=0.0_8
      NN=2
      
      !grmobmx=0.0_8
!**************************************************************
!**************************************************************
!****************self interaction******************************
  self: DO alpha=1,2
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

    pairwisealpha: do alpha=1,2
         pairwisebeta: do beta=alpha+1,2
            aaI=RADII(alpha)
            aaJ=RADII(beta)
            R=p_pos(1:3,alpha)-p_pos(1:3,beta)
            
            
            DMIN=SQRT( SUM(R**2))

          call ROTNE_PRAGER_IJ(TTP,RRP,TRP,QQP_TR,GPQ_TR,HPQ_TR,R,aaI,aaJ)



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


    do alpha=1,2
      do beta=1,2
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
     
  END SUBROUTINE GRPERY_MOB






!************************************************************
!***********************************************************

      SUBROUTINE ROTNE_PRAGER_IJ(TTP,RRP,TRP,QQP_TR,GPQ_TR,HPQ_TR,RN,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: TTP(3,3),RRP(3,3),TRP(3,3),QQP_TR(5,5)
      real*8,intent(out):: GPQ_TR(3,5),HPQ_TR(3,5)
      REAL*8,intent(in):: aaI,aaJ,RN(3)      ! radius
      
      REAL*8 TTPAA(3,3),TTPAB(3,3)
      REAL*8 RRPAA(3,3),TRPAA(3,3)
      REAL*8 QQPAA(3,3,3,3),QQPAB(3,3,3,3)
      REAL*8 GPQAA(3,3,3),HPQAA(3,3,3)
      real*8 QQP(3,3,3,3),GPQ(3,3,3),HPQ(3,3,3)

      INTEGER i,j,k,l,m
      REAL*8 DIST,nondim_a3,nondim_a
      real*8 Jg(3,3),D_Jg(3,3,3),DD_Jg(3,3,3,3),D2_Jg(3,3),DD2_Jg(3,3,3),DDD2_Jg(3,3,3,3)
      real*8 Rg(3,3),D_Rg(3,3,3),D2_Rg(3,3)
      real*8 Kg(3,3,3),D_Kg(3,3,3,3),D2_Kg(3,3,3),DD2_Kg(3,3,3,3)
      
      nondim_a=1.0_8/(6.0_8*aaI)
      nondim_a3=1.0_8/(6.0_8*aaI*aaI*aaI)
      TTP=0.0_8
      RRP=0.D0
      QQP=0.d0
      GPQ=0.0_8
      HPQ=0.0_8
      TRP =0.D0

     DIST=SQRT( SUM(RN*RN))

     call Init_Jg(rn,dist,Jg)
     call Init_D_Jg(rn,dist,D_Jg)
     call Init_DD_Jg(rn,dist,DD_Jg)
     call Init_D2_Jg(rn,dist,D2_Jg)
     call Init_DD2_Jg(rn,dist,DD2_Jg)
     call Init_DDD2_Jg(rn,dist,DDD2_Jg)

     call Init_Rg(D_Jg,Rg)
     call Init_D_Rg(DD_Jg,D_Rg)
     call Init_D2_Rg(DD2_Jg,D2_Rg)

     call Init_Kg(D_Jg,Kg)
     call Init_D_Kg(DD_Jg,D_Kg)
     call Init_D2_Kg(DD2_Jg,D2_Kg)
     call Init_DD2_Kg(DDD2_Jg,DD2_Kg)

     TTPAA= 0.75_8*aaI*Jg
     TTPAB= 0.75_8*aaI*(aaI*aaI+aaJ*aaJ)/6.0_8*D2_Jg
     TTP = (TTPAA + TTPAB)

     RRPAA=0.0_8
     do k=1,3
      do l=1,3
       forall(i=1:3,j=1:3)
         RRPAA(i,j)=RRPAA(i,j)+EPS(i,k,l)*D_Rg(l,j,k)
       end forall
      enddo
     enddo
     RRP = 3.0_8*aai*aai*aai/8.0_8*RRPAA

     QQPAA=0.0_8
     QQPAB=0.0_8
     forall(i=1:3,j=1:3,k=1:3,l=1:3)
     QQPAA(i,j,k,l)=0.5_8*(D_Kg(i,k,l,j)+D_Kg(j,k,l,i))
     QQPAB(i,j,k,l)=0.5_8*(DD2_Kg(i,k,l,j)+DD2_Kg(j,k,l,i))
     end forall

     QQP=-3.0_8/4.0_8*aai*aai*aai*QQPAA &
        -3.0_8/40.0_8*aai*aai*aai*(aai*aai+aaj*aaj)*QQPAB

     GPQAA=Kg+(aai*aai/6.0_8+0.1_8*aaj*aaj)*D2_Kg

     GPQ=-0.75*aai*GPQAA

     HPQAA=0.0_8
     do l=1,3
      do m=1,3
       forall(i=1:3,j=1:3,k=1:3)
        HPQAA(i,j,k)=HPQAA(i,j,k)+EPS(i,l,m)*(D_Kg(m,j,k,l)+0.1_8*aaj*aaj*DD2_Kg(m,j,k,l))
       end forall
      enddo
     enddo
     HPQ=-3.0_8/8.0_8*aai*aai*aai*HPQAA

     TRP = aai*(3.0_8/4.0_8)*Rg+aai*aai*aai/8.0_8*D2_Rg

      TTP = (TTP)*nondim_a
      RRP = (RRP)*nondim_a3
      QQP = (QQP)*nondim_a3
      GPQ = (GPQ)*nondim_a
      HPQ = (HPQ)*nondim_a3
      TRP = (TRP)*nondim_a
      call GH_TR(GPQ_TR,GPQ)
      call GH_TR(HPQ_TR,HPQ)
      call QQ_TR(QQP_TR,QQP)
      END SUBROUTINE ROTNE_PRAGER_IJ

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
      REAL*8 B1_CART(3,3,3,3),a
      INTEGER I,J
      CALL ROTNE_PRAGER_DD_SELF(B1_CART,aaI)
      Call DD_transY2(B1,B1_CART)
      END

      SUBROUTINE ROTNE_PRAGER_DD_SELF(B1,aaI)
      USE TENSORS
      IMPLICIT NONE
      REAL*8,intent(out):: B1(3,3,3,3)
      REAL*8,intent(in):: aaI          ! radiu  

      integer:: i,j,k,l
      B1=3.0_8/(20.0_8*aai**3)*II

      END

!***********************************************************
!***********************************************************

      SUBROUTINE FarfieldScalar(Resist)
      use TENSORS
      IMPLICIT NONE
      REAL*8,intent(in):: Resist(22,22)

      integer::i,j
      real*8 ::temp,abc(3,3),gh(5,3),mm(5,5)

      abc=Resist(1:3,1:3)
      X11A=abc(1,1)
      Y11A=abc(2,2)
      abc=Resist(1:3,4:6)
      X12A=abc(1,1)
      Y12A=abc(2,2)
      abc=Resist(7:9,1:3)
      Y11B=abc(2,3)
      abc=Resist(7:9,4:6)
      Y12B=abc(2,3)
      abc=Resist(7:9,7:9)
      X11C=abc(1,1)
      Y11C=abc(2,2)
      abc=Resist(7:9,10:12)
      X12C=abc(1,1)
      Y12C=abc(2,2)
      gh=Resist(13:17,1:3)
      X11G=gh(3,1)/L4_tr(3,1)
      write(*,*) 'X11G',X11G,gh(5,1)/L4_tr(5,1)
      Y11G=gh(1,2)/L5_tr(1,2)
      write(*,*) 'Y11G',Y11G,gh(4,3)/L5_tr(4,3)
      gh=Resist(13:17,4:6)
      X12G=gh(3,1)/L4_tr(3,1)
      Y12G=gh(1,2)/L5_tr(1,2)

      gh=Resist(13:17,7:9)
      Y11H=gh(1,3)/L6_tr(1,3)
      !write(*,*) 'Y11G',Y11G,gh(4,3)/L4_tr(4,3)
      gh=Resist(13:17,10:12)
      Y12H=gh(1,3)/L6_tr(1,3)

      mm=Resist(13:17,13:17)
      Y11M=mm(1,1)
      Z11M=mm(2,2)
      X11M=(mm(3,3)-L9_tr(3,3)*Z11M)/L7_tr(3,3)
      write(*,*) 'X11M',X11M,(mm(5,5)-L9_tr(5,5)*Z11M)/L7_tr(5,5)

      mm=Resist(13:17,18:22)
      Y12M=mm(1,1)
      Z12M=mm(2,2)
      X12M=(mm(3,3)-L9_tr(3,3)*Z12M)/L7_tr(3,3)

      end

