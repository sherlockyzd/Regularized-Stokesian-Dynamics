!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*********************************


   Module Regularization
   use tensors,only:kdirac,EPS,pip5,IIiso,pai
   use LATTICE_BASE,only:phi
      
   IMPLICIT NONE
   private

   public:: Regular_ROTNE_PRAGER_IJ

   contains

!****************************************************************

      SUBROUTINE Regular_ROTNE_PRAGER_IJ(TTP,RRP,TRP,QQP_TR,GPQ_TR,HPQ_TR,RN,aaI,aaJ)
      !use Ewald_summation
      !USE TENSORS
      !USE LATTICE_BASE
      IMPLICIT NONE
      REAL*8,intent(out):: TTP(3,3),RRP(3,3),TRP(3,3)
      real*8,intent(out):: QQP_TR(5,5),GPQ_TR(3,5),HPQ_TR(3,5)
      REAL*8,intent(in):: aaI,aaJ,RN(3)      ! radius
      
      REAL*8 TTPA(3,3),TTPAA(3,3),TTPAB(3,3),TTPB(3,3),TTPBA(3,3),TTPBB(3,3)
      REAL*8 RRPA(3,3),RRPB(3,3),RRPAA(3,3),RRPBA(3,3)
      REAL*8 TRPA(3,3),TRPB(3,3),TRPAA(3,3),TRPBA(3,3),GPQ(3,3,3),HPQ(3,3,3)
      REAL*8 QQPA(3,3,3,3),QQPB(3,3,3,3),QQPAA(3,3,3,3),QQPAB(3,3,3,3)
      REAL*8 QQPBA(3,3,3,3),QQPBB(3,3,3,3),QQP(3,3,3,3)
      REAL*8 GPQA(3,3,3),GPQB(3,3,3),GPQAA(3,3,3),GPQBA(3,3,3)
      real*8 HPQA(3,3,3),HPQB(3,3,3),HPQAA(3,3,3),HPQBA(3,3,3)
      REAL*8 E_def,erf1,erf2,erf3,erf4,erf5
      real*8 Ror(3,3),D2_Ror(3,3),D_Ror(3,3,3),Kr(3,3,3),D_Kr(3,3,3,3),D2_Kr(3,3,3),DD2_Kr(3,3,3,3)
      real*8 Jr(3,3),D_Jr(3,3,3),DD_Jr(3,3,3,3),D2_Jr(3,3),DD2_Jr(3,3,3),DDD2_Jr(3,3,3,3)
      REAL*8 k_ijl(3,3,3),k_ijlm(3,3,3,3),r_ijl(3,3,3),r_ijlm(3,3,3,3)
      real*8 DD_phir(3,3),DDD_phir(3,3,3),DDDD_phir(3,3,3,3)
      real*8 D2_phir,DD2_phir(3),DDD2_phir(3,3),DDDD2_phir(3,3,3),DDDDD2_phir(3,3,3,3)
      real*8 D4_phir,DD4_phir(3),DDD4_phir(3,3)


      INTEGER i,j,k,l,m
      !REAL*8 RN(3),NL(3)
      REAL*8 DIST,nondim_a3,nondim_a,kapa


      nondim_a=1.0_8/(6.0_8*aaI)
      nondim_a3=1.0_8/(6.0_8*aaI*aaI*aaI)

      TTPA=0.0_8
      TTPB=0.0_8
      RRPA=0.D0
      RRPB=0.d0
      QQPA=0.d0
      QQPB=0.D0
      GPQA=0.0_8
      HPQA=0.0_8
      GPQB=0.0_8
      HPQB=0.0_8
      TRPA =0.D0
      TRPB =0.d0
     
       
       !YY=DOT_PRODUCT(RN,RN)/SIG**2/2.D0
       !IF(YY.LE.LOGERR) THEN
        DIST=SQRT( SUM(RN**2))
        !kapa=sqrt(pai)/3.0_8/(1.5_8*(aaI+aaJ)-DIST)
        kapa=sqrt(pai)/3.0_8/(abs((aaI+aaJ)-DIST))
        !kapa=sqrt(pai)/3.0_8*2.0_8**((pai/(6.0_8*phi))**1.0_8/3.0_8)

        call Init_E_def(kapa,dist,E_def)
        call Init_erf1(kapa,dist,E_def,erf1)
        call Init_erf2(kapa,dist,E_def,erf2)
        call Init_erf3(kapa,dist,E_def,erf3)
        call Init_erf4(kapa,dist,E_def,erf4)
        call Init_erf5(kapa,dist,E_def,erf5)
       ! call Init_erfc6(dist,E_def,erfc6)
        call Init_rijl(RN,r_ijl)
        call Init_kijl(RN,k_ijl)
        call Init_rijlm(RN,r_ijlm)
        call Init_kijlm(RN,k_ijlm)
        
        call Init_DD_phir(kapa,RN,dist,erf1,DD_phir)
        call Init_D2_phir(kapa,dist,erf1,D2_phir)
        call Init_D4_phir(kapa,dist,erf2,erf3,D4_phir)
        call Init_DDD_phir(kapa,k_ijl,r_ijl,dist,erf1,erf2,DDD_phir)
        call Init_DD2_phir(kapa,rn,dist,erf1,erf2,DD2_phir)
        call Init_DD4_phir(kapa,rn,dist,erf2,erf3,erf4,DD4_phir)
        call Init_DDDD_phir(kapa,k_ijlm,r_ijlm,Dist,erf1,erf2,erf3,DDDD_phir)
        call Init_DDD2_phir(kapa,RN,dist,erf1,erf2,erf3,DDD2_phir)
        call Init_DDD4_phir(kapa,rn,dist,erf2,erf3,erf4,erf5,DDD4_phir)
        call Init_DDDD2_phir(kapa,k_ijl,r_ijl,dist,erf1,erf2,erf3,erf4,DDDD2_phir) 
        call Init_DDDDD2_phir(kapa,k_ijlm,r_ijlm,dist,erf1,erf2,erf3,erf4,erf5,DDDDD2_phir)      

        call Init_Jr(D2_phir,DD_phir,Jr)
        call Init_D_Jr(DD2_phir,DDD_phir,D_Jr)
        call Init_DD_Jr(DDD2_phir,DDDD_phir,DD_Jr)
        call Init_D2_Jr(D4_phir,DDD2_phir,D2_Jr)
        call Init_DD2_Jr(DD4_phir,DDDD2_phir,DD2_Jr)
        call Init_DDD2_Jr(DDD4_phir,DDDDD2_phir,DDD2_Jr)

        call Init_Ror(D_Jr,Ror)
        call init_D_Ror(DD_Jr,D_Ror)
        call Init_D2_Ror(DD2_Jr,D2_Ror)

        call Init_Kr(D_Jr,Kr)
        call Init_D_Kr(DD_Jr,D_Kr)
        call Init_D2_Kr(DD2_Jr,D2_Kr)
        call Init_DD2_Kr(DDD2_Jr,DD2_Kr)

        TTPAA= 0.75_8*aaI*Jr
        TTPAB= 0.75_8*aaI*(aaI*aaI+aaJ*aaJ)/6.0_8*D2_Jr
        TTPA = TTPAA + TTPAB

        RRPAA=0.0_8
        do k=1,3
         do l=1,3
          forall(i=1:3,j=1:3)
            RRPAA(i,j)=RRPAA(i,j)+EPS(i,k,l)*D_Ror(l,j,k)
          end forall
         enddo
        enddo
        RRPA = 3.0_8*aai*aai*aai/8.0_8*RRPAA

        QQPAA=0.0_8
        QQPAB=0.0_8
        forall(i=1:3,j=1:3,k=1:3,l=1:3)
        QQPAA(i,j,k,l)=0.5_8*(D_Kr(i,k,l,j)+D_Kr(j,k,l,i))
        QQPAB(i,j,k,l)=0.5_8*(DD2_Kr(i,k,l,j)+DD2_Kr(j,k,l,i))
        end forall

        QQPA=-3.0_8/4.0_8*aai*aai*aai*QQPAA &
           -3.0_8/40.0_8*aai*aai*aai*(aai*aai+aaj*aaj)*QQPAB
        GPQAA=Kr+(aai*aai/6.0_8+0.1_8*aaj*aaj)*D2_Kr
        GPQA=-0.75*aai*GPQAA

        HPQAA=0.0_8
        do l=1,3
         do m=1,3
          forall(i=1:3,j=1:3,k=1:3)
           HPQAA(i,j,k)=HPQAA(i,j,k)+EPS(i,l,m)*(D_Kr(m,j,k,l)+0.1_8*aaj*aaj*DD2_Kr(m,j,k,l))
          end forall
         enddo
        enddo
        HPQA=-3.0_8/8.0_8*aai*aai*aai*HPQAA

        TRPAA = aai*(3.0_8/4.0_8)*Ror+aai*aai*aai/8.0_8*D2_Ror
        TRPA =  TRPAA
       !ENDIF


      TTP = (TTPA )*nondim_a
      RRP = (RRPA )*nondim_a3
      QQP = (QQPA )*nondim_a3
      GPQ = (GPQA )*nondim_a
      HPQ = (HPQA )*nondim_a3
      TRP = (TRPA )*nondim_a
      call GH_TR(GPQ_TR,GPQ)
      call GH_TR(HPQ_TR,HPQ)
      call QQ_TR(QQP_TR,QQP)
      END SUBROUTINE Regular_ROTNE_PRAGER_IJ




!************************************************************

      SUBROUTINE Init_Jr(D2_phir,DD_phir,Jr)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:kdirac
      !IMPLICIT NONE
      real*8,intent(in):: D2_phir,DD_phir(3,3)
      REAL*8,intent(out):: Jr(3,3)

      integer i,j

      forall(i=1:3,j=1:3)
         Jr(i,j)=kdirac(i,j)*D2_phir-DD_phir(i,j)
      end forall

      end SUBROUTINE Init_Jr

      SUBROUTINE Init_D_Jr(DD2_phir,DDD_phir,D_Jr)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:kdirac
      !IMPLICIT NONE
      real*8,intent(in):: DD2_phir(3),DDD_phir(3,3,3)
      REAL*8,intent(out):: D_Jr(3,3,3)

      integer i,j,l

      forall(i=1:3,j=1:3,l=1:3)
         D_Jr(i,j,l)=kdirac(i,j)*DD2_phir(l)-DDD_phir(l,i,j)
      end forall

      end SUBROUTINE Init_D_Jr

      SUBROUTINE Init_DD_Jr(DDD2_phir,DDDD_phir,DD_Jr)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:kdirac
      !IMPLICIT NONE
      real*8,intent(in):: DDD2_phir(3,3),DDDD_phir(3,3,3,3)
      REAL*8,intent(out):: DD_Jr(3,3,3,3)

      integer i,j,l,m

      forall(i=1:3,j=1:3,l=1:3,m=1:3)
         DD_Jr(i,j,m,l)=kdirac(i,j)*DDD2_phir(m,l)-DDDD_phir(m,l,i,j)
      end forall

      end SUBROUTINE Init_DD_Jr


      SUBROUTINE Init_D2_Jr(D4_phir,DDD2_phir,D2_Jr)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:kdirac
      !IMPLICIT NONE
      real*8,intent(in):: D4_phir,DDD2_phir(3,3)
      REAL*8,intent(out):: D2_Jr(3,3)

      integer i,j

      forall(i=1:3,j=1:3)
         D2_Jr(i,j)=kdirac(i,j)*D4_phir-DDD2_phir(i,j)
      end forall

      end SUBROUTINE Init_D2_Jr

      SUBROUTINE Init_DD2_Jr(DD4_phir,DDDD2_phir,DD2_Jr)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:kdirac
      !IMPLICIT NONE
      real*8,intent(in):: DD4_phir(3),DDDD2_phir(3,3,3)
      REAL*8,intent(out):: DD2_Jr(3,3,3)

      integer i,j,k

      forall(i=1:3,j=1:3,k=1:3)
         DD2_Jr(i,j,k)=kdirac(i,j)*DD4_phir(k)-DDDD2_phir(k,i,j)
      end forall

      end SUBROUTINE Init_DD2_Jr

      SUBROUTINE Init_DDD2_Jr(DDD4_phir,DDDDD2_phir,DDD2_Jr)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:kdirac
      !IMPLICIT NONE
      real*8,intent(in):: DDD4_phir(3,3),DDDDD2_phir(3,3,3,3)
      REAL*8,intent(out):: DDD2_Jr(3,3,3,3)

      integer i,j,k,l

      forall(i=1:3,j=1:3,k=1:3,l=1:3)
         DDD2_Jr(i,j,l,k)=kdirac(i,j)*DDD4_phir(l,k)-DDDDD2_phir(l,k,i,j)
      end forall

      end SUBROUTINE Init_DDD2_Jr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*********************************

      SUBROUTINE Init_Ror(D_Jr,Ror)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:EPS
      !IMPLICIT NONE
      real*8,intent(in):: D_Jr(3,3,3)
      REAL*8,intent(out):: Ror(3,3)

      integer i,j,k,l

      Ror=0.0_8
      do k=1,3
       do l=1,3
        forall(i=1:3,j=1:3)
         Ror(i,j)=Ror(i,j)-0.5_8*EPS(j,k,l)*D_Jr(i,l,k)
        end forall
       enddo
      enddo

      end SUBROUTINE Init_Ror

      SUBROUTINE Init_Kr(D_Jr,Kr)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:EPS
      !IMPLICIT NONE
      real*8,intent(in):: D_Jr(3,3,3)
      REAL*8,intent(out):: Kr(3,3,3)

      integer i,j,k

      forall(i=1:3,j=1:3,k=1:3)

         Kr(i,j,k)=0.5_8*(D_Jr(i,j,k)+D_Jr(i,k,j))

      end forall

      end SUBROUTINE Init_Kr

      SUBROUTINE Init_D_Ror(DD_Jr,D_Ror)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:EPS
      !IMPLICIT NONE
      real*8,intent(in):: DD_Jr(3,3,3,3)
      REAL*8,intent(out):: D_Ror(3,3,3)

      integer i,j,l,m,n

      D_Ror=0.0_8
      do m=1,3
       do n=1,3
        forall(i=1:3,j=1:3,l=1:3)
         D_Ror(i,j,l)=D_Ror(i,j,l)-0.5_8*EPS(j,m,n)*DD_Jr(i,n,l,m)
        end forall
       enddo
      enddo

      end SUBROUTINE Init_D_Ror


      SUBROUTINE Init_D_Kr(DD_Jr,D_Kr)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:EPS
      !IMPLICIT NONE
      real*8,intent(in):: DD_Jr(3,3,3,3)
      REAL*8,intent(out):: D_Kr(3,3,3,3)

      integer i,j,k,l

      forall(i=1:3,j=1:3,k=1:3,l=1:3)

         D_Kr(i,j,k,l)=0.5_8*(DD_Jr(i,j,l,k)+DD_Jr(i,k,l,j))

      end forall

      end SUBROUTINE Init_D_Kr


      SUBROUTINE Init_D2_Ror(DD2_Jr,D2_Ror)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:EPS
      !IMPLICIT NONE
      real*8,intent(in):: DD2_Jr(3,3,3)
      REAL*8,intent(out):: D2_Ror(3,3)

      integer i,j,k,l

      D2_Ror=0.0_8

      do k=1,3
       do l=1,3
        forall(i=1:3,j=1:3)
        D2_Ror(i,j)=D2_Ror(i,j)-0.5_8*EPS(j,k,l)*DD2_Jr(i,l,k)
        end forall
       enddo
      enddo

      end SUBROUTINE Init_D2_Ror

      SUBROUTINE Init_D2_Kr(DD2_Jr,D2_Kr)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:EPS
     ! IMPLICIT NONE
      real*8,intent(in):: DD2_Jr(3,3,3)
      REAL*8,intent(out):: D2_Kr(3,3,3)

      integer i,j,k

      forall(i=1:3,j=1:3,k=1:3)

         D2_Kr(i,j,k)=0.5_8*(DD2_Jr(i,j,k)+DD2_Jr(i,k,j))

      end forall

      end SUBROUTINE Init_D2_Kr


      SUBROUTINE Init_DD2_Kr(DDD2_Jr,DD2_Kr)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:EPS
      !IMPLICIT NONE
      real*8,intent(in):: DDD2_Jr(3,3,3,3)
      REAL*8,intent(out):: DD2_Kr(3,3,3,3)

      integer i,j,k,l

      forall(i=1:3,j=1:3,k=1:3,l=1:3)

         DD2_Kr(i,j,k,l)=0.5_8*(DDD2_Jr(i,j,l,k)+DDD2_Jr(i,k,l,j))

      end forall

      end SUBROUTINE Init_DD2_Kr




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!******************************

      SUBROUTINE Init_E_def(kapa,r,E_def)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:pip5
      !IMPLICIT NONE
      real*8,intent(in):: r,kapa
      REAL*8,intent(out):: E_def

      real*8 rlam2

      rlam2=r*r*kapa*kapa
      E_def=2.0_8*pip5*kapa*exp(-rlam2)

      end SUBROUTINE Init_E_def

      SUBROUTINE Init_erf1(kapa,r,E_def,erf1)
     ! use LATTICE_BASE,only:kapa
      !IMPLICIT NONE
      real*8,intent(in):: r,E_def,kapa
      REAL*8,intent(out):: erf1

      erf1=E_def
      end SUBROUTINE Init_erf1

      SUBROUTINE Init_erf2(kapa,r,E_def,erf2)
      !use LATTICE_BASE,only:kapa
      !IMPLICIT NONE
      real*8,intent(in):: r,E_def,kapa
      REAL*8,intent(out):: erf2

      erf2=-2.0_8*kapa*kapa*r*E_def

      end SUBROUTINE Init_erf2

      SUBROUTINE Init_erf3(kapa,r,E_def,erf3)
      !use LATTICE_BASE,only:kapa
      !IMPLICIT NONE
      real*8,intent(in):: r,E_def,kapa
      REAL*8,intent(out):: erf3

      real*8 rlam2
      
      rlam2=r*r*kapa*kapa
      erf3=2.0_8*kapa*kapa*E_def*(2.0_8*rlam2-1.0_8)

      end SUBROUTINE Init_erf3

      SUBROUTINE Init_erf4(kapa,r,E_def,erf4)
      !use LATTICE_BASE,only:kapa
      !IMPLICIT NONE
      real*8,intent(in):: r,E_def,kapa
      REAL*8,intent(out):: erf4

      real*8 rlam2
      
      rlam2=r*r*kapa*kapa

      erf4=-4.0_8*kapa*kapa*kapa*kapa*r*E_def*(2.0_8*rlam2-3.0_8)

      end SUBROUTINE Init_erf4

      SUBROUTINE Init_erf5(kapa,r,E_def,erf5)
      !use LATTICE_BASE,only:kapa
      !IMPLICIT NONE
      real*8,intent(in):: r,E_def,kapa
      REAL*8,intent(out):: erf5

      real*8 rlam2,rlam4
      
      rlam2=r*r*kapa*kapa
      rlam4=rlam2*rlam2
      erf5=4.0_8*kapa*kapa*kapa*kapa*E_def*(4.0_8*rlam4-12.0_8*rlam2+3.0_8)

      end SUBROUTINE Init_erf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*********************************


      SUBROUTINE Init_kijl(rn,k_ijl)
      !use LATTICE_BASE,only:kapa
     ! use tensors, only:KDirac
      !IMPLICIT NONE
      real*8,intent(in):: rn(3)
      REAL*8,intent(out):: k_ijl(3,3,3)

      integer i,j,l
      
      forall(i=1:3,j=1:3,l=1:3)
         k_ijl(i,j,l)=kdirac(i,j)*rn(l)+kdirac(i,l)*rn(j)+kdirac(j,l)*rn(i)
      end forall
      end SUBROUTINE Init_kijl

      SUBROUTINE Init_kijlm(rn,k_ijlm)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:KDirac
      !IMPLICIT NONE
      real*8,intent(in):: rn(3)
      REAL*8,intent(out):: k_ijlm(3,3,3,3)

      integer i,j,l,m
      
      forall(i=1:3,j=1:3,l=1:3,m=1:3)
        k_ijlm(i,j,l,m)=kdirac(i,j)*rn(l)*rn(m)+kdirac(i,l)*rn(j)*rn(m)  &
           +kdirac(j,l)*rn(i)*rn(m)+kdirac(i,m)*rn(j)*rn(l)+kdirac(j,m)*rn(i)*rn(l)  &
           +kdirac(l,m)*rn(i)*rn(j)
      end forall
      end SUBROUTINE Init_kijlm

      SUBROUTINE Init_rijl(rn,r_ijl)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:KDirac
      !IMPLICIT NONE
      real*8,intent(in):: rn(3)
      REAL*8,intent(out):: r_ijl(3,3,3)

      integer i,j,l
      
      forall(i=1:3,j=1:3,l=1:3)
         r_ijl(i,j,l)=rn(l)*rn(j)*rn(i)
      end forall
      end SUBROUTINE Init_rijl

      SUBROUTINE Init_rijlm(rn,r_ijlm)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:KDirac
      !IMPLICIT NONE
      real*8,intent(in):: rn(3)
      REAL*8,intent(out):: r_ijlm(3,3,3,3)

      integer i,j,l,m
      
      forall(i=1:3,j=1:3,l=1:3,m=1:3)
         r_ijlm(i,j,l,m)=rn(l)*rn(j)*rn(i)*rn(m)
      end forall
      end SUBROUTINE Init_rijlm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*********************************

      SUBROUTINE Init_D_phir(kapa,RN,r,D_phir)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:KDirac
      !IMPLICIT NONE
      real*8,intent(in):: rn(3),r,kapa
      REAL*8,intent(out):: D_phir(3)

      integer j
      
      forall(j=1:3)
         D_phir(j)=rn(j)/r*erf(kapa*r)
      end forall
      end SUBROUTINE Init_D_phir

      SUBROUTINE Init_DD_phir(kapa,RN,r,erf1,DD_phir)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:KDirac
      !IMPLICIT NONE
      real*8,intent(in):: rn(3),r,erf1,kapa
      REAL*8,intent(out):: DD_phir(3,3)

      integer i,j
      real*8 r2,r3
      r2=r*r
      r3=r*r*r
      
      forall(i=1:3,j=1:3)
         DD_phir(i,j)=(kdirac(i,j)/r-rn(i)*rn(j)/r3)*erf(kapa*r)  &
           +rn(i)*rn(j)/r2*erf1
      end forall

      end SUBROUTINE Init_DD_phir

      SUBROUTINE Init_DDD_phir(kapa,k_ijl,r_ijl,r,erf1,erf2,DDD_phir)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:KDirac
      !IMPLICIT NONE
      real*8,intent(in):: k_ijl(3,3,3),r_ijl(3,3,3),r,erf1,erf2,kapa
      REAL*8,intent(out):: DDD_phir(3,3,3)

      integer i,j,l
      real*8 r2,r3,r4,r5
      r2=r*r
      r3=r*r*r
      r4=r2*r2
      r5=r3*r2
      
      forall(i=1:3,j=1:3,l=1:3)
         DDD_phir(i,j,l)=(-k_ijl(i,j,l)/r3+3.0_8*r_ijl(i,j,l)/r5)*erf(kapa*r)  &
           +(k_ijl(i,j,l)/r2-3.0_8*r_ijl(i,j,l)/r4)*erf1+r_ijl(i,j,l)/r3*erf2
      end forall
      end SUBROUTINE Init_DDD_phir


      SUBROUTINE Init_DDDD_phir(kapa,k_ijlm,r_ijlm,r,erf1,erf2,erf3,DDDD_phir)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:IIiso
      !IMPLICIT NONE
      real*8,intent(in):: k_ijlm(3,3,3,3),r_ijlm(3,3,3,3),r,erf1,erf2,erf3,kapa
      REAL*8,intent(out):: DDDD_phir(3,3,3,3)

      integer i,j,l,m
      real*8 r2,r3,r4,r5,r6,r7
      r2=r*r
      r3=r*r*r
      r4=r2*r2
      r5=r3*r2
      r6=r5*r
      r7=r6*r

      forall(i=1:3,j=1:3,l=1:3,m=1:3)
         DDDD_phir(i,j,l,m)=(-IIiso(i,j,l,m)/r3+3.0_8*k_ijlm(i,j,l,m)/r5 &
            -15.0_8*r_ijlm(i,j,l,m)/r7)*erf(kapa*r)  &
           +(IIiso(i,j,l,m)/r2-3.0_8*k_ijlm(i,j,l,m)/r4+15.0_8*r_ijlm(i,j,l,m)/r6)*erf1 &
           +(k_ijlm(i,j,l,m)/r3-6.0_8*r_ijlm(i,j,l,m)/r5)*erf2 &
           +r_ijlm(i,j,l,m)/r4*erf3
      end forall
      end SUBROUTINE Init_DDDD_phir

      SUBROUTINE Init_D2_phir(kapa,r,erf1,D2_phir)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:KDirac
      !IMPLICIT NONE
      real*8,intent(in):: r,erf1,kapa
      REAL*8,intent(out):: D2_phir

      D2_phir=2.0_8/r*erf(kapa*r)+erf1
      end SUBROUTINE Init_D2_phir




      SUBROUTINE Init_DD2_phir(kapa,rn,r,erf1,erf2,DD2_phir)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:KDirac
      !IMPLICIT NONE
      real*8,intent(in):: rn(3),r,erf1,erf2,kapa
      REAL*8,intent(out):: DD2_phir(3)

      integer l
      real*8 r2,r3
      r2=r*r
      r3=r*r*r
      
      forall(l=1:3)
         DD2_phir(l)=-2.0_8*rn(l)/r3*erf(kapa*r)+ 2.0_8*rn(l)/r2*erf1 &
             +rn(l)/r*erf2
      end forall
      end SUBROUTINE Init_DD2_phir


      SUBROUTINE Init_DDD2_phir(kapa,rn,r,erf1,erf2,erf3,DDD2_phir)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:kdirac
      !IMPLICIT NONE
      real*8,intent(in):: rn(3),r,erf1,erf2,erf3,kapa
      REAL*8,intent(out):: DDD2_phir(3,3)

      integer l,m
      real*8 r2,r3,r4,r5
      r2=r*r
      r3=r*r*r
      r4=r2*r2
      r5=r3*r2

      forall(l=1:3,m=1:3)
         DDD2_phir(l,m)=(-2.0_8*kdirac(l,m)/r3+6.0_8*rn(l)*rn(m)/r5)*erf(kapa*r)  &
           +(2.0_8*kdirac(l,m)/r2-6.0_8*rn(l)*rn(m)/r4)*erf1 &
           +(kdirac(l,m)/r+rn(l)*rn(m)/r3)*erf2 &
           +rn(l)*rn(m)/r2*erf3
      end forall
      end SUBROUTINE Init_DDD2_phir

      SUBROUTINE Init_DDDD2_phir(kapa,k_ijl,r_ijl,r,erf1,erf2,erf3,erf4,DDDD2_phir)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:IIiso
      !IMPLICIT NONE
      real*8,intent(in):: k_ijl(3,3,3),r_ijl(3,3,3),r,erf1,erf2,erf3,erf4,kapa
      REAL*8,intent(out):: DDDD2_phir(3,3,3)

      integer i,j,k
      real*8 r2,r3,r4,r5,r6,r7
      r2=r*r
      r3=r*r*r
      r4=r2*r2
      r5=r3*r2
      r6=r5*r
      r7=r6*r

      forall(i=1:3,j=1:3,k=1:3)
         DDDD2_phir(i,j,k)=(6.0_8*k_ijl(i,j,k)/r5-30.0_8*r_ijl(i,j,k)/r7)*erf(kapa*r)  &
           +(-6.0_8*k_ijl(i,j,k)/r4+30.0_8*r_ijl(i,j,k)/r6)*erf1 &
           +(k_ijl(i,j,k)/r3-9.0_8*r_ijl(i,j,k)/r5)*erf2 &
           +(k_ijl(i,j,k)/r2-r_ijl(i,j,k)/r4)*erf3 &
           +r_ijl(i,j,k)/r3*erf4
      end forall
      end SUBROUTINE Init_DDDD2_phir


      SUBROUTINE Init_DDDDD2_phir(kapa,k_ijlm,r_ijlm,r,erf1,erf2,erf3,erf4,erf5,DDDDD2_phir)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:IIiso
      !IMPLICIT NONE
      real*8,intent(in):: k_ijlm(3,3,3,3),r_ijlm(3,3,3,3),r,erf1,erf2,erf3,erf4,erf5,kapa
      REAL*8,intent(out):: DDDDD2_phir(3,3,3,3)

      integer i,j,l,m
      real*8 r2,r3,r4,r5,r6,r7,r8,r9
      r2=r*r
      r3=r*r*r
      r4=r2*r2
      r5=r3*r2
      r6=r5*r
      r7=r6*r
      r8=r7*r
      r9=r8*r

      forall(i=1:3,j=1:3,l=1:3,m=1:3)
         DDDDD2_phir(i,j,l,m)=(6.0_8*IIiso(i,j,l,m)/r5-30.0_8*k_ijlm(i,j,l,m)/r7 &
            +210.0_8*r_ijlm(i,j,l,m)/r9)*erf(kapa*r)  &
           +(-6.0_8*IIiso(i,j,l,m)/r4+30.0_8*k_ijlm(i,j,l,m)/r6-210.0_8*r_ijlm(i,j,l,m)/r8)*erf1 &
           +(IIiso(i,j,l,m)/r3-9.0_8*k_ijlm(i,j,l,m)/r5+75.0_8*r_ijlm(i,j,l,m)/r7)*erf2 &
           +(IIiso(i,j,l,m)/r2-k_ijlm(i,j,l,m)/r4-5.0_8*r_ijlm(i,j,l,m)/r6)*erf3 &
           +(k_ijlm(i,j,l,m)/r3-4.0_8*r_ijlm(i,j,l,m)/r5)*erf4 &
           +r_ijlm(i,j,l,m)/r4*erf5
      end forall
      end SUBROUTINE Init_DDDDD2_phir

      SUBROUTINE Init_D4_phir(kapa,r,erf2,erf3,D4_phir)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:KDirac
      !IMPLICIT NONE
      real*8,intent(in):: r,erf2,erf3,kapa
      REAL*8,intent(out):: D4_phir


         D4_phir=4.0_8/r*erf2+erf3

      end SUBROUTINE Init_D4_phir

      SUBROUTINE Init_DD4_phir(kapa,rn,r,erf2,erf3,erf4,DD4_phir)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:KDirac
      !IMPLICIT NONE
      real*8,intent(in):: rn(3),r,erf2,erf3,erf4,kapa
      REAL*8,intent(out):: DD4_phir(3)

      integer l
      real*8 r2,r3
      r2=r*r
      r3=r*r*r
      
      forall(l=1:3)
         DD4_phir(l)=-4.0_8*rn(l)/r3*erf2+ 4.0_8*rn(l)/r2*erf3 &
             +rn(l)/r*erf4
      end forall
      end SUBROUTINE Init_DD4_phir


      SUBROUTINE Init_DDD4_phir(kapa,rn,r,erf2,erf3,erf4,erf5,DDD4_phir)
      !use LATTICE_BASE,only:kapa
      !use tensors, only:kdirac
      !IMPLICIT NONE
      real*8,intent(in):: rn(3),r,erf2,erf3,erf4,erf5,kapa
      REAL*8,intent(out):: DDD4_phir(3,3)

      integer l,m
      real*8 r2,r3,r4,r5
      r2=r*r
      r3=r*r*r
      r4=r2*r2
      r5=r3*r2

      forall(l=1:3,m=1:3)
         DDD4_phir(l,m)=(-4.0_8*kdirac(l,m)/r3+12.0_8*rn(l)*rn(m)/r5)*erf2  &
           +(4.0_8*kdirac(l,m)/r2-12.0_8*rn(l)*rn(m)/r4)*erf3 &
           +(kdirac(l,m)/r+3.0_8*rn(l)*rn(m)/r3)*erf4 &
           +rn(l)*rn(m)/r2*erf5
      end forall
      end SUBROUTINE Init_DDD4_phir

   end Module Regularization