   Module stokesian_Green
   use tensors,only:kdirac,EPS
      
   IMPLICIT NONE
   private

   public:: ROTNE_PRAGER_IJ_FTS,ROTNE_PRAGER_IJ_FT,ROTNE_PRAGER_TT_IJ
   contains


      SUBROUTINE ROTNE_PRAGER_IJ_FTS(TTP,RRP,TRP,QQP_TR,GPQ_TR,HPQ_TR,RN,aaI,aaJ)
      !USE TENSORS
      !use stokesian_Green
      IMPLICIT NONE
      REAL*8,intent(out):: TTP(3,3),RRP(3,3),TRP(3,3),QQP_TR(5,5)
      real*8,intent(out):: GPQ_TR(3,5),HPQ_TR(3,5)
      REAL*8,intent(in):: aaI,aaJ,RN(3)      ! radius
      
      REAL*8 TTPAA(3,3),TTPAB(3,3)
      REAL*8 RRPAA(3,3)
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

     DIST=SQRT( SUM(RN**2))

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
      END SUBROUTINE ROTNE_PRAGER_IJ_FTS

      SUBROUTINE ROTNE_PRAGER_IJ_FT(TTP,RRP,TRP,RN,aaI,aaJ)
      !USE TENSORS
      !use stokesian_Green
      IMPLICIT NONE
      REAL*8,intent(out):: TTP(3,3),RRP(3,3),TRP(3,3)
      REAL*8,intent(in):: aaI,aaJ,RN(3)      ! radius
      
      REAL*8 TTPAA(3,3),TTPAB(3,3)
      REAL*8 RRPAA(3,3)


      INTEGER i,j,k,l,m
      REAL*8 DIST,nondim_a3,nondim_a
      real*8 Jg(3,3),D_Jg(3,3,3),DD_Jg(3,3,3,3),D2_Jg(3,3),DD2_Jg(3,3,3)
      real*8 Rg(3,3),D_Rg(3,3,3),D2_Rg(3,3)
      
      nondim_a=1.0_8/(6.0_8*aaI)
      nondim_a3=1.0_8/(6.0_8*aaI*aaI*aaI)
      TTP=0.0_8
      RRP=0.0_8


     DIST=SQRT( SUM(RN**2))

     call Init_Jg(rn,dist,Jg)
     call Init_D_Jg(rn,dist,D_Jg)
     call Init_DD_Jg(rn,dist,DD_Jg)
     call Init_D2_Jg(rn,dist,D2_Jg)
     call Init_DD2_Jg(rn,dist,DD2_Jg)


     call Init_Rg(D_Jg,Rg)
     call Init_D_Rg(DD_Jg,D_Rg)
     call Init_D2_Rg(DD2_Jg,D2_Rg)



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


     
     TRP = aai*(3.0_8/4.0_8)*Rg+aai*aai*aai/8.0_8*D2_Rg

      TTP = (TTP)*nondim_a
      RRP = (RRP)*nondim_a3
      TRP = (TRP)*nondim_a


      END SUBROUTINE ROTNE_PRAGER_IJ_FT








!****************************************************************************************************


      SUBROUTINE ROTNE_PRAGER_TT_IJ(TTP,RN,aaI,aaJ)
      !USE TENSORS
      !use stokesian_Green
      IMPLICIT NONE
      REAL*8,intent(out):: TTP(3,3)
      REAL*8,intent(in):: aaI,aaJ,RN(3)      ! radius
      
      REAL*8 TTPAA(3,3),TTPAB(3,3)

      INTEGER i,j,k,l,m
      REAL*8 DIST,nondim_a
      real*8 Jg(3,3),D2_Jg(3,3)
      
      nondim_a=1.0_8/(6.0_8*aaI)

      TTP=0.0_8


     DIST=SQRT( SUM(RN**2))

     call Init_Jg(rn,dist,Jg)

     call Init_D2_Jg(rn,dist,D2_Jg)


     TTPAA= 0.75_8*aaI*Jg
     TTPAB= 0.75_8*aaI*(aaI*aaI+aaJ*aaJ)/6.0_8*D2_Jg
     TTP = (TTPAA + TTPAB)
     TTP = (TTP)*nondim_a
    END SUBROUTINE ROTNE_PRAGER_TT_IJ
!******************************************************************************************

      SUBROUTINE Init_Jg(rn,dist,Jg)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:kdirac
      !IMPLICIT NONE
      real*8,intent(in):: rn(3),dist
      REAL*8,intent(out):: Jg(3,3)

      integer i,j

      forall(i=1:3,j=1:3)
         Jg(i,j)=kdirac(i,j)/dist+rn(i)*rn(j)/(dist**3.0_8)
      end forall
      end SUBROUTINE Init_Jg


      SUBROUTINE Init_D_Jg(rn,dist,D_Jg)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:kdirac
      !IMPLICIT NONE
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
      !use tensors, only:kdirac
      !IMPLICIT NONE
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
      !use tensors, only:kdirac
      !IMPLICIT NONE
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
      !use tensors, only:kdirac
      !IMPLICIT NONE
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
      !use tensors, only:kdirac
      !IMPLICIT NONE
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
      !use tensors, only:EPS
      !IMPLICIT NONE
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
      !use tensors, only:EPS
      !IMPLICIT NONE
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
      !use tensors, only:EPS
      !IMPLICIT NONE
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
      !IMPLICIT NONE
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
      !IMPLICIT NONE
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
      !IMPLICIT NONE
      real*8,intent(in):: DD2_Jg(3,3,3)
      REAL*8,intent(out):: D2_Kg(3,3,3)

      !integer i,j,k

         D2_Kg=DD2_Jg

      end SUBROUTINE Init_D2_Kg

      SUBROUTINE Init_DD2_Kg(DDD2_Jg,DD2_Kg)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:EPS
      !IMPLICIT NONE
      real*8,intent(in):: DDD2_Jg(3,3,3,3)
      REAL*8,intent(out):: DD2_Kg(3,3,3,3)

      !integer i,j,k

         DD2_Kg=DDD2_Jg

      end SUBROUTINE Init_DD2_Kg

   end Module stokesian_Green