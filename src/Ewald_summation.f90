!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*********************************


   Module Ewald_summation
   use tensors,only:kdirac,EPS,pip5,IIiso,pai,II
   use LATTICE_BASE,only:lammda,phi,LI,LR,LV,NR,NI,SIG
      
   IMPLICIT NONE
   private

   public:: PER_ROTNE_PRAGER_SELF,PER_ROTNE_PRAGER_IJ

   contains


!************************************************************

      SUBROUTINE PER_ROTNE_PRAGER_SELF(TTS,RRS,QQS_TR,aaI,LOGERR)
      IMPLICIT NONE
      REAL*8,intent(out):: TTS(3,3),RRS(3,3),QQS_TR(5,5)
      REAL*8,intent(in):: aaI,LOGERR       ! radius
      
      REAL*8 TTSA(3,3),TTSAA(3,3),TTSAB(3,3),TTSB(3,3),TTSR(3,3)
      REAL*8 RRSA(3,3),RRSB(3,3),RRSR(3,3),RRSAA(3,3),RRSBA(3,3)
      REAL*8 QQSA(3,3,3,3),QQSB(3,3,3,3),QQSAA(3,3,3,3),QQSAB(3,3,3,3)
      REAL*8 QQSBA(3,3,3,3),QQSBB(3,3,3,3),QQSR(3,3,3,3),QQS(3,3,3,3)

      INTEGER NX,NY,NZ,i,j,k,l
      REAL*8 RN(3),KN(3),NL(3),nondim_a3,nondim_a
      REAL*8 YY,phix,DIST,aaJ
      REAL*8 E_def,erfc1,erfc2,erfc3,erfc4,erfc5,erfc6
      real*8 Jr(3,3),D2_Jr(3,3),DD_Jr(3,3,3,3),DDD2_Jr(3,3,3,3)
      real*8 D_Ror(3,3,3)
      real*8 D_Kr(3,3,3,3),DD2_Kr(3,3,3,3)

      real*8 F_Jk(3,3),F_D2Jk(3,3),F_DDJk(3,3,3,3)
      real*8 F_DRk(3,3,3)
      real*8 F_DKk(3,3,3,3),F_DD2Kk(3,3,3,3)

      real*8 DD_rerfc(3,3),D2_rerfc,DDD2_rerfc(3,3),D4_rerfc,DDDD_rerfc(3,3,3,3),DDDDD2_rerfc(3,3,3,3),DDD4_rerfc(3,3)
      
      real*8 k_ijlm(3,3,3,3),r_ijlm(3,3,3,3)

      phix=1-0.2_8*phi
      aaJ=aaI
      nondim_a=1.0_8/(6.0_8*aaI)
      nondim_a3=1.0_8/(6.0_8*aaI*aaI*aaI)

      TTSA = 0.D0
      TTSB =0.d0
      RRSA=0.D0
      RRSB=0.d0
      QQSA=0.d0
      QQSB=0.D0
          
! REAL LATTICE
      DO NX=-NR,NR
       DO NY=-NR,NR
        DO NZ=-NR,NR
          IF(NX.EQ.0.AND.NY.EQ.0.AND.NZ.EQ.0) CYCLE
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ

          RN=MATMUL(LR,NL)
          !YY=DOT_PRODUCT(RN,RN)/SIG**2/2.D0
          !IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(RN**2))
           !RW=RN/DIST
           !CALL CALC_RR(RR,RW)
           call Init_E_def(dist,E_def)
           call Init_erfc1(dist,E_def,erfc1)
           call Init_erfc2(dist,E_def,erfc2)
           call Init_erfc3(dist,E_def,erfc3)
           call Init_erfc4(dist,E_def,erfc4)
           call Init_erfc5(dist,E_def,erfc5)
           call Init_erfc6(dist,E_def,erfc6)

           call Init_rijlm(RN,r_ijlm)
           call Init_kijlm(RN,k_ijlm)

           call Init_D2_rerfc(dist,erfc1,erfc2,D2_rerfc)
           call Init_D4_rerfc(dist,erfc2,erfc3,erfc4,D4_rerfc)
           call Init_DD_rerfc(RN,dist,erfc1,erfc2,DD_rerfc)
           call Init_DDDD_rerfc(k_ijlm,r_ijlm,Dist,erfc1,erfc2,erfc3,erfc4,DDDD_rerfc)
           call Init_DDD2_rerfc(RN,Dist,erfc1,erfc2,erfc3,erfc4,DDD2_rerfc)
           call Init_DDD4_rerfc(rn,dist,erfc2,erfc3,erfc4,erfc5,erfc6,DDD4_rerfc)
           call Init_DDDDD2_rerfc(k_ijlm,r_ijlm,dist,erfc1,erfc2,erfc3,erfc4,erfc5,erfc6,DDDDD2_rerfc)

           call Init_Jr(D2_rerfc,DD_rerfc,Jr)
           call Init_DD_Jr(DDD2_rerfc,DDDD_rerfc,DD_Jr)
           call Init_D2_Jr(D4_rerfc,DDD2_rerfc,D2_Jr)
           call Init_DDD2_Jr(DDD4_rerfc,DDDDD2_rerfc,DDD2_Jr)

           call init_D_Ror(DD_Jr,D_Ror)

           call Init_D_Kr(DD_Jr,D_Kr)
           call Init_DD2_Kr(DDD2_Jr,DD2_Kr)

           TTSAA= 0.75_8*aaI*Jr    
           TTSAB= 0.75_8*aaI*(aaI*aaI+aaJ*aaJ)/6.0_8*D2_Jr*phix
           TTSA = TTSA + (TTSAA+TTSAB)

           RRSAA=0.0_8
           do k=1,3
            do l=1,3
             forall(i=1:3,j=1:3)
               RRSAA(i,j)=RRSAA(i,j)+EPS(i,k,l)*D_Ror(l,j,k)
             end forall
            enddo
           enddo

           RRSA = RRSA + 3.0_8*aai*aai*aai/8.0_8*RRSAA

           QQSAA=0.0_8
           QQSAB=0.0_8
           forall(i=1:3,j=1:3,k=1:3,l=1:3)
           QQSAA(i,j,k,l)=0.5_8*(D_Kr(i,k,l,j)+D_Kr(j,k,l,i))
           QQSAB(i,j,k,l)=0.5_8*(DD2_Kr(i,k,l,j)+DD2_Kr(j,k,l,i))
           end forall

          QQSA=QQSA-3.0_8/4.0_8*aai*aai*aai*QQSAA-3.0_8/40.0_8 &
               *aai*aai*aai*(aai*aai+aaj*aaj)*QQSAB
          !ENDIF
        ENDDO
       ENDDO
      ENDDO
      !A1 = A1 - (SIG**2)*P/(2.D0*LV)*KDirac

! INVERSE LATTICE
      DO NX=-NI,NI
       DO NY=-NI,NI
        DO NZ=-NI,NI
         IF(NX.EQ.0.AND.NY.EQ.0.AND.NZ.EQ.0) CYCLE
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          KN=2*PAI*MATMUL(NL,LI)
          !YY=DOT_PRODUCT(KN,KN)*SIG**2/2.D0
          !IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(KN**2) )

           call Init_F_Jk(kn,dist,F_Jk)
           call Init_F_D2Jk(dist,F_Jk,F_D2Jk)
           call Init_F_DDJk(kn,F_Jk,F_DDJk)
           call Init_F_DRk(F_DDJk,F_DRk)
           call Init_F_DKk(F_DDJk,F_DKk)
           call Init_F_DD2Kk(kn,dist,F_Jk,F_DD2Kk)

           TTSB =TTSB +0.75_8*aai*F_Jk+1.0_8/8.0_8*aai* &
                 (aaI*aaI+aaJ*aaJ)*phix*F_D2Jk 

           RRSBA=0.0_8
           do k=1,3
            do l=1,3
             forall(i=1:3,j=1:3)
               RRSBA(i,j)=RRSBA(i,j)+EPS(i,k,l)*F_DRk(l,j,k)
             end forall
            enddo
           enddo

           RRSB = RRSB + 3.0_8/8.0_8*aai*aai*aai*RRSBA   

           forall(i=1:3,j=1:3,k=1:3,l=1:3)
           QQSBA(i,j,k,l)=0.5_8*(F_DKk(i,k,l,j)+F_DKk(j,k,l,i))
           QQSBB(i,j,k,l)=0.5_8*(F_DD2Kk(i,k,l,j)+F_DD2Kk(j,k,l,i))
           end forall

           QQSB=QQSB-3.0_8/4.0_8*aai*aai*aai*QQSBA-3.0_8/40.0_8 &
           *aai*aai*aai*(aai*aai+aaj*aaj)*QQSBB


         ! ENDIF
        ENDDO
       ENDDO
      ENDDO

      TTSR=(1.0_8-6.0_8*Lammda/sqrt(PAI)*aaI+40.0_8/3.0_8/sqrt(PAI)* &
        Lammda*Lammda*Lammda*aaI*aaI*aai*phix)*KDirac
      TTS = (TTSA + TTSB/LV+TTSR)*nondim_a

      RRSR=(0.75_8*(1.0_8-0.0_8*phi)-10.0_8/sqrt(pai)*aai*aai*aai*lammda*lammda*lammda)*KDirac
      RRS = (RRSR+RRSA+RRSB/LV)*nondim_a3

      QQSR=(0.9_8-lammda*lammda*lammda*aaI*aaI*aaI/25.0_8/sqrt(Pai) &
           *(50.0_8-126.0_8*lammda*lammda*aaI*aai)*6.0_8)*II
      QQS =(QQSA+QQSB/LV+QQSR)*nondim_a3
      call QQ_TR(QQS_TR,QQS)

      END SUBROUTINE PER_ROTNE_PRAGER_SELF




      SUBROUTINE PER_ROTNE_PRAGER_IJ(TTP,RRP,TRP,QQP_TR,GPQ_TR,HPQ_TR,R,aaI,aaJ,LOGERR)
      IMPLICIT NONE
      REAL*8,intent(out):: TTP(3,3),RRP(3,3),TRP(3,3)
      real*8,intent(out):: QQP_TR(5,5),GPQ_TR(3,5),HPQ_TR(3,5)
      REAL*8,intent(in):: aaI,aaJ,R(3),LOGERR     ! radius
      
      REAL*8 TTPA(3,3),TTPAA(3,3),TTPAB(3,3),TTPB(3,3),TTPBA(3,3),TTPBB(3,3)
      REAL*8 RRPA(3,3),RRPB(3,3),RRPAA(3,3),RRPBA(3,3)
      REAL*8 TRPA(3,3),TRPB(3,3),TRPAA(3,3),TRPBA(3,3),GPQ(3,3,3),HPQ(3,3,3)
      REAL*8 QQPA(3,3,3,3),QQPB(3,3,3,3),QQPAA(3,3,3,3),QQPAB(3,3,3,3)
      REAL*8 QQPBA(3,3,3,3),QQPBB(3,3,3,3),QQP(3,3,3,3)
      REAL*8 GPQA(3,3,3),GPQB(3,3,3),GPQAA(3,3,3),GPQBA(3,3,3)
      real*8 HPQA(3,3,3),HPQB(3,3,3),HPQAA(3,3,3),HPQBA(3,3,3)
      REAL*8 E_def,erfc1,erfc2,erfc3,erfc4,erfc5,erfc6
      real*8 Ror(3,3),D2_Ror(3,3),D_Ror(3,3,3),Kr(3,3,3),D_Kr(3,3,3,3),D2_Kr(3,3,3),DD2_Kr(3,3,3,3)
      real*8 Jr(3,3),D_Jr(3,3,3),DD_Jr(3,3,3,3),D2_Jr(3,3),DD2_Jr(3,3,3),DDD2_Jr(3,3,3,3)
      REAL*8 k_ijl(3,3,3),k_ijlm(3,3,3,3),r_ijl(3,3,3),r_ijlm(3,3,3,3)
      real*8 DD_rerfc(3,3),DDD_rerfc(3,3,3),DDDD_rerfc(3,3,3,3)
      real*8 D2_rerfc,DD2_rerfc(3),DDD2_rerfc(3,3),DDDD2_rerfc(3,3,3),DDDDD2_rerfc(3,3,3,3)
      real*8 D4_rerfc,DD4_rerfc(3),DDD4_rerfc(3,3)
      real*8 F_Jk(3,3),F_DJk(3,3,3),F_DDJk(3,3,3,3),F_D2Jk(3,3),F_DD2Jk(3,3,3)
      real*8 F_Rk(3,3),F_D2Rk(3,3),F_DRk(3,3,3),F_Kk(3,3,3),F_DKk(3,3,3,3),F_D2Kk(3,3,3),F_DD2Kk(3,3,3,3)

      INTEGER NX,NY,NZ,i,j,k,l,m
      REAL*8 RN(3),KN(3),NL(3)
      REAL*8 YY,phix,DIST,nondim_a3,nondim_a

      phix=1-0.2_8*phi
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
     
! REAL LATTICE
      DO NX=-NR,NR
       DO NY=-NR,NR
        DO NZ=-NR,NR
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          RN=R+MATMUL(LR,NL)

          !YY=DOT_PRODUCT(RN,RN)/SIG**2/2.D0
          !IF(YY.LE.LOGERR) THEN

           DIST=SQRT( SUM(RN**2))
           !RW=RN/DIST
           !CALL CALC_RR(RR,RW)
           call Init_E_def(dist,E_def)
           call Init_erfc1(dist,E_def,erfc1)
           call Init_erfc2(dist,E_def,erfc2)
           call Init_erfc3(dist,E_def,erfc3)
           call Init_erfc4(dist,E_def,erfc4)
           call Init_erfc5(dist,E_def,erfc5)
           call Init_erfc6(dist,E_def,erfc6)
           call Init_rijl(RN,r_ijl)
           call Init_kijl(RN,k_ijl)
           call Init_rijlm(RN,r_ijlm)
           call Init_kijlm(RN,k_ijlm)
           
           call Init_DD_rerfc(RN,dist,erfc1,erfc2,DD_rerfc)
           call Init_D2_rerfc(dist,erfc1,erfc2,D2_rerfc)
           call Init_D4_rerfc(dist,erfc2,erfc3,erfc4,D4_rerfc)
           call Init_DDD_rerfc(k_ijl,r_ijl,dist,erfc1,erfc2,erfc3,DDD_rerfc)
           call Init_DD2_rerfc(rn,dist,erfc1,erfc2,erfc3,DD2_rerfc)
           call Init_DD4_rerfc(rn,dist,erfc2,erfc3,erfc4,erfc5,DD4_rerfc)
           call Init_DDDD_rerfc(k_ijlm,r_ijlm,Dist,erfc1,erfc2,erfc3,erfc4,DDDD_rerfc)
           call Init_DDD2_rerfc(RN,dist,erfc1,erfc2,erfc3,erfc4,DDD2_rerfc)
           call Init_DDD4_rerfc(rn,dist,erfc2,erfc3,erfc4,erfc5,erfc6,DDD4_rerfc)
           call Init_DDDD2_rerfc(k_ijl,r_ijl,dist,erfc1,erfc2,erfc3,erfc4,erfc5,DDDD2_rerfc) 
           call Init_DDDDD2_rerfc(k_ijlm,r_ijlm,dist,erfc1,erfc2,erfc3,erfc4,erfc5,erfc6,DDDDD2_rerfc)      

           call Init_Jr(D2_rerfc,DD_rerfc,Jr)
           call Init_D_Jr(DD2_rerfc,DDD_rerfc,D_Jr)
           call Init_DD_Jr(DDD2_rerfc,DDDD_rerfc,DD_Jr)
           call Init_D2_Jr(D4_rerfc,DDD2_rerfc,D2_Jr)
           call Init_DD2_Jr(DD4_rerfc,DDDD2_rerfc,DD2_Jr)
           call Init_DDD2_Jr(DDD4_rerfc,DDDDD2_rerfc,DDD2_Jr)

           call Init_Ror(D_Jr,Ror)
           call init_D_Ror(DD_Jr,D_Ror)
           call Init_D2_Ror(DD2_Jr,D2_Ror)
    
           call Init_Kr(D_Jr,Kr)
           call Init_D_Kr(DD_Jr,D_Kr)
           call Init_D2_Kr(DD2_Jr,D2_Kr)
           call Init_DD2_Kr(DDD2_Jr,DD2_Kr)

           TTPAA= 0.75_8*aaI*Jr
           TTPAB= 0.75_8*aaI*(aaI*aaI+aaJ*aaJ)/6.0_8*D2_Jr*phix
           TTPA = TTPA + (TTPAA + TTPAB)

           RRPAA=0.0_8
           do k=1,3
            do l=1,3
             forall(i=1:3,j=1:3)
               RRPAA(i,j)=RRPAA(i,j)+EPS(i,k,l)*D_Ror(l,j,k)
             end forall
            enddo
           enddo
           RRPA = RRPA + 3.0_8*aai*aai*aai/8.0_8*RRPAA

           QQPAA=0.0_8
           QQPAB=0.0_8
           forall(i=1:3,j=1:3,k=1:3,l=1:3)
           QQPAA(i,j,k,l)=0.5_8*(D_Kr(i,k,l,j)+D_Kr(j,k,l,i))
           QQPAB(i,j,k,l)=0.5_8*(DD2_Kr(i,k,l,j)+DD2_Kr(j,k,l,i))
           end forall

           QQPA=QQPA-3.0_8/4.0_8*aai*aai*aai*QQPAA &
              -3.0_8/40.0_8*aai*aai*aai*(aai*aai+aaj*aaj)*QQPAB
           GPQAA=Kr+(aai*aai/6.0_8+0.1_8*aaj*aaj)*D2_Kr
           GPQA=GPQA-0.75*aai*GPQAA

           HPQAA=0.0_8
           do l=1,3
            do m=1,3
             forall(i=1:3,j=1:3,k=1:3)
              HPQAA(i,j,k)=HPQAA(i,j,k)+EPS(i,l,m)*(D_Kr(m,j,k,l)+0.1_8*aaj*aaj*DD2_Kr(m,j,k,l))
             end forall
            enddo
           enddo
           HPQA=HPQA-3.0_8/8.0_8*aai*aai*aai*HPQAA

           TRPAA = aai*(3.0_8/4.0_8)*Ror+aai*aai*aai/8.0_8*D2_Ror
           TRPA = TRPA + TRPAA
          !ENDIF
        ENDDO
       ENDDO
      ENDDO


! INVERSE LATTICE
      DO NX=-NI,NI
       DO NY=-NI,NI
        DO NZ=-NI,NI
         IF(NX.EQ.0.AND.NY.EQ.0.AND.NZ.EQ.0) CYCLE
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          KN=2*PAI*MATMUL(NL,LI)
         ! YY=DOT_PRODUCT(KN,KN)*SIG**2/2.D0
          !IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(KN**2) )

           call Init_F_Jk(kn,dist,F_Jk)
           call Init_F_DJk(kn,F_Jk,F_DJk)
           call Init_F_DDJk(kn,F_Jk,F_DDJk)
           call Init_F_D2Jk(dist,F_Jk,F_D2Jk)
           call Init_F_DD2Jk(kn,F_D2Jk,F_DD2Jk)

           call Init_F_Rk(F_DJk,F_Rk)
           call Init_F_DRk(F_DDJk,F_DRk)
           call Init_F_D2Rk(F_DD2Jk,F_D2Rk)

           call Init_F_Kk(F_DJk,F_Kk)
           call Init_F_DKk(F_DDJk,F_DKk)
           call Init_F_D2Kk(F_DD2Jk,F_D2Kk)
           call Init_F_DD2Kk(kn,dist,F_Jk,F_DD2Kk)
           
           !call Init_F_DRk(F_DDJk,F_DRk)
           TTPBA= 0.75_8*aai*F_Jk
           TTPBB= 0.75_8*aaI*(aaI*aaI+aaJ*aaJ)/6.0_8*F_D2Jk*phix
           TTPB = TTPB + (TTPBA + TTPBB)*cos(sum(KN*R))


           RRPBA=0.0_8
           do k=1,3
            do l=1,3
             forall(i=1:3,j=1:3)
               RRPBA(i,j)=RRPBA(i,j)+EPS(i,k,l)*F_DRk(l,j,k)
             end forall
            enddo
           enddo
           RRPB = RRPB + 3.0_8/8.0_8*aai*aai*aai*RRPBA*cos(sum(KN*R))   


           QQPBA=0.0_8
           QQPBB=0.0_8
           forall(i=1:3,j=1:3,k=1:3,l=1:3)
           QQPBA(i,j,k,l)=0.5_8*(F_DKk(i,k,l,j)+F_DKk(j,k,l,i))
           QQPBB(i,j,k,l)=0.5_8*(F_DD2Kk(i,k,l,j)+F_DD2Kk(j,k,l,i))
           end forall
           QQPB=QQPB+(-3.0_8/4.0_8*aai*aai*aai*QQPBA-3.0_8/40.0_8*aai*aai*aai &
            *(aai*aai+aaj*aaj)*QQPBB)*cos(sum(KN*R))

           GPQBA=F_Kk+(aai*aai/6.0_8+0.1_8*aaj*aaj)*F_D2Kk
           GPQB=GPQB+0.75*aai*GPQBA*sin(sum(KN*R))

           HPQBA=0.0_8
           do l=1,3
            do m=1,3
             forall(i=1:3,j=1:3,k=1:3)
              HPQBA(i,j,k)=HPQBA(i,j,k)+EPS(i,l,m)*(F_DKk(m,j,k,l)+0.1_8*aaj*aaj*F_DD2Kk(m,j,k,l))
             end forall
            enddo
           enddo
           HPQB=HPQB-3.0_8/8.0_8*aai*aai*aai*HPQBA*cos(sum(KN*R))

           TRPBA= aai*(3.0_8/4.0_8)*F_Rk+aai*aai*aai/8.0_8*F_D2Rk
           TRPB = TRPB - TRPBA*sin(sum(KN*R))

         ! ENDIF
        ENDDO
       ENDDO
      ENDDO

      TTP = (TTPA + TTPB/LV)*nondim_a
      RRP = (RRPA + RRPB/LV)*nondim_a3
      QQP = (QQPA + QQPB/LV)*nondim_a3
      GPQ = (GPQA + GPQB/LV)*nondim_a
      HPQ = (HPQA + HPQB/LV)*nondim_a3
      TRP = (TRPA + TRPB/LV)*nondim_a
      call GH_TR(GPQ_TR,GPQ)
      call GH_TR(HPQ_TR,HPQ)
      call QQ_TR(QQP_TR,QQP)
      END SUBROUTINE PER_ROTNE_PRAGER_IJ



      SUBROUTINE Init_Jr(D2_rerfc,DD_rerfc,Jr)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:kdirac
      !IMPLICIT NONE
      real*8,intent(in):: D2_rerfc,DD_rerfc(3,3)
      REAL*8,intent(out):: Jr(3,3)

      integer i,j

      forall(i=1:3,j=1:3)
         Jr(i,j)=kdirac(i,j)*D2_rerfc-DD_rerfc(i,j)
      end forall

      end SUBROUTINE Init_Jr

      SUBROUTINE Init_D_Jr(DD2_rerfc,DDD_rerfc,D_Jr)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:kdirac
      !IMPLICIT NONE
      real*8,intent(in):: DD2_rerfc(3),DDD_rerfc(3,3,3)
      REAL*8,intent(out):: D_Jr(3,3,3)

      integer i,j,l

      forall(i=1:3,j=1:3,l=1:3)
         D_Jr(i,j,l)=kdirac(i,j)*DD2_rerfc(l)-DDD_rerfc(l,i,j)
      end forall

      end SUBROUTINE Init_D_Jr

      SUBROUTINE Init_DD_Jr(DDD2_rerfc,DDDD_rerfc,DD_Jr)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:kdirac
      !IMPLICIT NONE
      real*8,intent(in):: DDD2_rerfc(3,3),DDDD_rerfc(3,3,3,3)
      REAL*8,intent(out):: DD_Jr(3,3,3,3)

      integer i,j,l,m

      forall(i=1:3,j=1:3,l=1:3,m=1:3)
         DD_Jr(i,j,m,l)=kdirac(i,j)*DDD2_rerfc(m,l)-DDDD_rerfc(m,l,i,j)
      end forall

      end SUBROUTINE Init_DD_Jr


      SUBROUTINE Init_D2_Jr(D4_rerfc,DDD2_rerfc,D2_Jr)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:kdirac
      !IMPLICIT NONE
      real*8,intent(in):: D4_rerfc,DDD2_rerfc(3,3)
      REAL*8,intent(out):: D2_Jr(3,3)

      integer i,j

      forall(i=1:3,j=1:3)
         D2_Jr(i,j)=kdirac(i,j)*D4_rerfc-DDD2_rerfc(i,j)
      end forall

      end SUBROUTINE Init_D2_Jr

      SUBROUTINE Init_DD2_Jr(DD4_rerfc,DDDD2_rerfc,DD2_Jr)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:kdirac
      !IMPLICIT NONE
      real*8,intent(in):: DD4_rerfc(3),DDDD2_rerfc(3,3,3)
      REAL*8,intent(out):: DD2_Jr(3,3,3)

      integer i,j,k

      forall(i=1:3,j=1:3,k=1:3)
         DD2_Jr(i,j,k)=kdirac(i,j)*DD4_rerfc(k)-DDDD2_rerfc(k,i,j)
      end forall

      end SUBROUTINE Init_DD2_Jr

      SUBROUTINE Init_DDD2_Jr(DDD4_rerfc,DDDDD2_rerfc,DDD2_Jr)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:kdirac
      !IMPLICIT NONE
      real*8,intent(in):: DDD4_rerfc(3,3),DDDDD2_rerfc(3,3,3,3)
      REAL*8,intent(out):: DDD2_Jr(3,3,3,3)

      integer i,j,k,l

      forall(i=1:3,j=1:3,k=1:3,l=1:3)
         DDD2_Jr(i,j,l,k)=kdirac(i,j)*DDD4_rerfc(l,k)-DDDDD2_rerfc(l,k,i,j)
      end forall

      end SUBROUTINE Init_DDD2_Jr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*********************************

      SUBROUTINE Init_Ror(D_Jr,Ror)
      !use LATTICE_BASE,only:lammda
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
      !use LATTICE_BASE,only:lammda
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
      !use LATTICE_BASE,only:lammda
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
      !use LATTICE_BASE,only:lammda
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
      !use LATTICE_BASE,only:lammda
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
      !use LATTICE_BASE,only:lammda
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
      !use LATTICE_BASE,only:lammda
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

      SUBROUTINE Init_E_def(r,E_def)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:pip5
      !IMPLICIT NONE
      real*8,intent(in):: r
      REAL*8,intent(out):: E_def

      real*8 rlam2

      rlam2=r*r*lammda*lammda
      E_def=2.0_8*pip5*lammda*exp(-rlam2)

      end SUBROUTINE Init_E_def

      SUBROUTINE Init_erfc1(r,E_def,erfc1)
     ! use LATTICE_BASE,only:lammda
      !IMPLICIT NONE
      real*8,intent(in):: r,E_def
      REAL*8,intent(out):: erfc1

      erfc1=-E_def
      end SUBROUTINE Init_erfc1

      SUBROUTINE Init_erfc2(r,E_def,erfc2)
      !use LATTICE_BASE,only:lammda
      !IMPLICIT NONE
      real*8,intent(in):: r,E_def
      REAL*8,intent(out):: erfc2

      erfc2=2.0_8*lammda*lammda*r*E_def

      end SUBROUTINE Init_erfc2

      SUBROUTINE Init_erfc3(r,E_def,erfc3)
      !use LATTICE_BASE,only:lammda
      !IMPLICIT NONE
      real*8,intent(in):: r,E_def
      REAL*8,intent(out):: erfc3

      real*8 rlam2
      
      rlam2=r*r*lammda*lammda
      erfc3=-2.0_8*lammda*lammda*E_def*(2.0_8*rlam2-1.0_8)

      end SUBROUTINE Init_erfc3

      SUBROUTINE Init_erfc4(r,E_def,erfc4)
      !use LATTICE_BASE,only:lammda
      !IMPLICIT NONE
      real*8,intent(in):: r,E_def
      REAL*8,intent(out):: erfc4

      real*8 rlam2
      
      rlam2=r*r*lammda*lammda

      erfc4=4.0_8*lammda*lammda*lammda*lammda*r*E_def*(2.0_8*rlam2-3.0_8)

      end SUBROUTINE Init_erfc4

      SUBROUTINE Init_erfc5(r,E_def,erfc5)
      !use LATTICE_BASE,only:lammda
      !IMPLICIT NONE
      real*8,intent(in):: r,E_def
      REAL*8,intent(out):: erfc5

      real*8 rlam2,rlam4
      
      rlam2=r*r*lammda*lammda
      rlam4=rlam2*rlam2
      erfc5=-4.0_8*lammda*lammda*lammda*lammda*E_def*(4.0_8*rlam4-12.0_8*rlam2+3.0_8)

      end SUBROUTINE Init_erfc5

      SUBROUTINE Init_erfc6(r,E_def,erfc6)
      !use LATTICE_BASE,only:lammda
      !IMPLICIT NONE
      real*8,intent(in):: r,E_def
      REAL*8,intent(out):: erfc6

      real*8 rlam2,rlam4
      
      rlam2=r*r*lammda*lammda
      rlam4=rlam2*rlam2
      erfc6=8.0_8*lammda*lammda*lammda*lammda*lammda*lammda*r*E_def*(4.0_8*rlam4-20.0_8*rlam2+15.0_8)
      end SUBROUTINE Init_erfc6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*********************************


      SUBROUTINE Init_kijl(rn,k_ijl)
      !use LATTICE_BASE,only:lammda
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
      !use LATTICE_BASE,only:lammda
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
      !use LATTICE_BASE,only:lammda
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
      !use LATTICE_BASE,only:lammda
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

      SUBROUTINE Init_D_rerfc(RN,r,erfc1,D_rerfc)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:KDirac
      !IMPLICIT NONE
      real*8,intent(in):: rn(3),r,erfc1
      REAL*8,intent(out):: D_rerfc(3)

      integer j
      
      forall(j=1:3)
         D_rerfc(j)=rn(j)/r*erfc(lammda*r)+rn(j)*erfc1
      end forall
      end SUBROUTINE Init_D_rerfc

      SUBROUTINE Init_DD_rerfc(RN,r,erfc1,erfc2,DD_rerfc)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:KDirac
      !IMPLICIT NONE
      real*8,intent(in):: rn(3),r,erfc1,erfc2
      REAL*8,intent(out):: DD_rerfc(3,3)

      integer i,j
      real*8 r2,r3
      r2=r*r
      r3=r*r*r
      
      forall(i=1:3,j=1:3)
         DD_rerfc(i,j)=(kdirac(i,j)/r-rn(i)*rn(j)/r3)*erfc(lammda*r)  &
           +(kdirac(i,j)+rn(i)*rn(j)/r2)*erfc1+rn(i)*rn(j)/r*erfc2
      end forall

      end SUBROUTINE Init_DD_rerfc

      SUBROUTINE Init_DDD_rerfc(k_ijl,r_ijl,r,erfc1,erfc2,erfc3,DDD_rerfc)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:KDirac
      !IMPLICIT NONE
      real*8,intent(in):: k_ijl(3,3,3),r_ijl(3,3,3),r,erfc1,erfc2,erfc3
      REAL*8,intent(out):: DDD_rerfc(3,3,3)

      integer i,j,l
      real*8 r2,r3,r4,r5
      r2=r*r
      r3=r*r*r
      r4=r2*r2
      r5=r3*r2
      
      forall(i=1:3,j=1:3,l=1:3)
         DDD_rerfc(i,j,l)=(-k_ijl(i,j,l)/r3+3.0_8*r_ijl(i,j,l)/r5)*erfc(lammda*r)  &
           -(-k_ijl(i,j,l)/r2+3.0_8*r_ijl(i,j,l)/r4)*erfc1+k_ijl(i,j,l)/r*erfc2+r_ijl(i,j,l)/r2*erfc3
      end forall
      end SUBROUTINE Init_DDD_rerfc


      SUBROUTINE Init_DDDD_rerfc(k_ijlm,r_ijlm,r,erfc1,erfc2,erfc3,erfc4,DDDD_rerfc)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:IIiso
      !IMPLICIT NONE
      real*8,intent(in):: k_ijlm(3,3,3,3),r_ijlm(3,3,3,3),r,erfc1,erfc2,erfc3,erfc4
      REAL*8,intent(out):: DDDD_rerfc(3,3,3,3)

      integer i,j,l,m
      real*8 r2,r3,r4,r5,r6,r7
      r2=r*r
      r3=r*r*r
      r4=r2*r2
      r5=r3*r2
      r6=r5*r
      r7=r6*r

      forall(i=1:3,j=1:3,l=1:3,m=1:3)
         DDDD_rerfc(i,j,l,m)=(-IIiso(i,j,l,m)/r3+3.0_8*k_ijlm(i,j,l,m)/r5 &
            -15.0_8*r_ijlm(i,j,l,m)/r7)*erfc(lammda*r)  &
           -(-IIiso(i,j,l,m)/r2+3.0_8*k_ijlm(i,j,l,m)/r4-15.0_8*r_ijlm(i,j,l,m)/r6)*erfc1 &
           +(IIiso(i,j,l,m)/r-3.0_8*r_ijlm(i,j,l,m)/r5)*erfc2 &
           +(k_ijlm(i,j,l,m)/r2-2.0_8*r_ijlm(i,j,l,m)/r4)*erfc3 &
           +r_ijlm(i,j,l,m)/r3*erfc4
      end forall
      end SUBROUTINE Init_DDDD_rerfc

      SUBROUTINE Init_D2_rerfc(r,erfc1,erfc2,D2_rerfc)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:KDirac
      !IMPLICIT NONE
      real*8,intent(in):: r,erfc1,erfc2
      REAL*8,intent(out):: D2_rerfc


         D2_rerfc=2.0_8/r*erfc(lammda*r)+4.0_8*erfc1+r*erfc2

      end SUBROUTINE Init_D2_rerfc




      SUBROUTINE Init_DD2_rerfc(rn,r,erfc1,erfc2,erfc3,DD2_rerfc)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:KDirac
      !IMPLICIT NONE
      real*8,intent(in):: rn(3),r,erfc1,erfc2,erfc3
      REAL*8,intent(out):: DD2_rerfc(3)

      integer l
      real*8 r2,r3
      r2=r*r
      r3=r*r*r
      
      forall(l=1:3)
         DD2_rerfc(l)=-2.0_8*rn(l)/r3*erfc(lammda*r)+ 2.0_8*rn(l)/r2*erfc1 &
             +5.0_8*rn(l)/r*erfc2+rn(l)*erfc3
      end forall
      end SUBROUTINE Init_DD2_rerfc


      SUBROUTINE Init_DDD2_rerfc(rn,r,erfc1,erfc2,erfc3,erfc4,DDD2_rerfc)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:kdirac
      !IMPLICIT NONE
      real*8,intent(in):: rn(3),r,erfc1,erfc2,erfc3,erfc4
      REAL*8,intent(out):: DDD2_rerfc(3,3)

      integer l,m
      real*8 r2,r3,r4,r5
      r2=r*r
      r3=r*r*r
      r4=r2*r2
      r5=r3*r2

      forall(l=1:3,m=1:3)
         DDD2_rerfc(l,m)=(-2.0_8*kdirac(l,m)/r3+6.0_8*rn(l)*rn(m)/r5)*erfc(lammda*r)  &
           -(-2.0_8*kdirac(l,m)/r2+6.0_8*rn(l)*rn(m)/r4)*erfc1 &
           +(5.0_8*kdirac(l,m)/r-3.0_8*rn(l)*rn(m)/r3)*erfc2 &
           +(kdirac(l,m)+5.0_8*rn(l)*rn(m)/r2)*erfc3 &
           +rn(l)*rn(m)/r*erfc4
      end forall
      end SUBROUTINE Init_DDD2_rerfc

      SUBROUTINE Init_DDDD2_rerfc(k_ijl,r_ijl,r,erfc1,erfc2,erfc3,erfc4,erfc5,DDDD2_rerfc)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:IIiso
      !IMPLICIT NONE
      real*8,intent(in):: k_ijl(3,3,3),r_ijl(3,3,3),r,erfc1,erfc2,erfc3,erfc4,erfc5
      REAL*8,intent(out):: DDDD2_rerfc(3,3,3)

      integer i,j,k
      real*8 r2,r3,r4,r5,r6,r7
      r2=r*r
      r3=r*r*r
      r4=r2*r2
      r5=r3*r2
      r6=r5*r
      r7=r6*r

      forall(i=1:3,j=1:3,k=1:3)
         DDDD2_rerfc(i,j,k)=(6.0_8*k_ijl(i,j,k)/r5-30.0_8*r_ijl(i,j,k)/r7)*erfc(lammda*r)  &
           -(6.0_8*k_ijl(i,j,k)/r4-30.0_8*r_ijl(i,j,k)/r6)*erfc1 &
           +(-3.0_8*k_ijl(i,j,k)/r3+3.0_8*r_ijl(i,j,k)/r5)*erfc2 &
           +(5.0_8*k_ijl(i,j,k)/r2-13.0_8*r_ijl(i,j,k)/r4)*erfc3 &
           +(k_ijl(i,j,k)/r+4.0_8*r_ijl(i,j,k)/r3)*erfc4 +r_ijl(i,j,k)/r2*erfc5
      end forall
      end SUBROUTINE Init_DDDD2_rerfc


      SUBROUTINE Init_DDDDD2_rerfc(k_ijlm,r_ijlm,r,erfc1,erfc2,erfc3,erfc4,erfc5,erfc6,DDDDD2_rerfc)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:IIiso
      !IMPLICIT NONE
      real*8,intent(in):: k_ijlm(3,3,3,3),r_ijlm(3,3,3,3),r,erfc1,erfc2,erfc3,erfc4,erfc5,erfc6
      REAL*8,intent(out):: DDDDD2_rerfc(3,3,3,3)

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
         DDDDD2_rerfc(i,j,l,m)=(6.0_8*IIiso(i,j,l,m)/r5-30.0_8*k_ijlm(i,j,l,m)/r7 &
            +210.0_8*r_ijlm(i,j,l,m)/r9)*erfc(lammda*r)  &
           -(6.0_8*IIiso(i,j,l,m)/r4-30.0_8*k_ijlm(i,j,l,m)/r6+210.0_8*r_ijlm(i,j,l,m)/r8)*erfc1 &
           +(-3.0_8*IIiso(i,j,l,m)/r3+3.0_8*k_ijlm(i,j,l,m)/r5+15.0_8*r_ijlm(i,j,l,m)/r7)*erfc2 &
           +(5.0_8*IIiso(i,j,l,m)/r2-13.0_8*k_ijlm(i,j,l,m)/r4+55.0_8*r_ijlm(i,j,l,m)/r6)*erfc3 &
           +(IIiso(i,j,l,m)/r+4.0_8*k_ijlm(i,j,l,m)/r3-25.0_8*r_ijlm(i,j,l,m)/r5)*erfc4 &
           +(k_ijlm(i,j,l,m)/r2+2.0_8*r_ijlm(i,j,l,m)/r4)*erfc5 &
           +r_ijlm(i,j,l,m)/r3*erfc6
      end forall
      end SUBROUTINE Init_DDDDD2_rerfc

      SUBROUTINE Init_D4_rerfc(r,erfc2,erfc3,erfc4,D4_rerfc)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:KDirac
      !IMPLICIT NONE
      real*8,intent(in):: r,erfc2,erfc3,erfc4
      REAL*8,intent(out):: D4_rerfc


         D4_rerfc=12.0_8/r*erfc2+8.0_8*erfc3+r*erfc4

      end SUBROUTINE Init_D4_rerfc

      SUBROUTINE Init_DD4_rerfc(rn,r,erfc2,erfc3,erfc4,erfc5,DD4_rerfc)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:KDirac
      !IMPLICIT NONE
      real*8,intent(in):: rn(3),r,erfc2,erfc3,erfc4,erfc5
      REAL*8,intent(out):: DD4_rerfc(3)

      integer l
      real*8 r2,r3
      r2=r*r
      r3=r*r*r
      
      forall(l=1:3)
         DD4_rerfc(l)=-12.0_8*rn(l)/r3*erfc2+ 12.0_8*rn(l)/r2*erfc3 &
             +9.0_8*rn(l)/r*erfc4+rn(l)*erfc5
      end forall
      end SUBROUTINE Init_DD4_rerfc


      SUBROUTINE Init_DDD4_rerfc(rn,r,erfc2,erfc3,erfc4,erfc5,erfc6,DDD4_rerfc)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:kdirac
      !IMPLICIT NONE
      real*8,intent(in):: rn(3),r,erfc2,erfc3,erfc4,erfc5,erfc6
      REAL*8,intent(out):: DDD4_rerfc(3,3)

      integer l,m
      real*8 r2,r3,r4,r5
      r2=r*r
      r3=r*r*r
      r4=r2*r2
      r5=r3*r2

      forall(l=1:3,m=1:3)
         DDD4_rerfc(l,m)=(-12.0_8*kdirac(l,m)/r3+36.0_8*rn(l)*rn(m)/r5)*erfc2  &
           -(-12.0_8*kdirac(l,m)/r2+36.0_8*rn(l)*rn(m)/r4)*erfc3 &
           +(9.0_8*kdirac(l,m)/r+3.0_8*rn(l)*rn(m)/r3)*erfc4 &
           +(kdirac(l,m)+9.0_8*rn(l)*rn(m)/r2)*erfc5 &
           +rn(l)*rn(m)/r*erfc6
      end forall
      end SUBROUTINE Init_DDD4_rerfc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*********************************

      SUBROUTINE Init_F_Jk(kn,dist,F_Jk)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:kdirac,pai
      !IMPLICIT NONE
      real*8,intent(in):: kn(3),dist
      REAL*8,intent(out):: F_Jk(3,3)

      integer i,j
      real*8 k2,k4
      k2=dist*dist
      k4=k2*k2

      forall(i=1:3,j=1:3)
         F_Jk(i,j)=(kdirac(i,j)*k2-kn(i)*kn(j))/k4*8.0_8*pai* &
            (1+k2/(4.0_8*lammda*lammda)+k4/(8.0_8*lammda*lammda*lammda*lammda)) &
            *exp(-k2/(4.0_8*lammda*lammda))
      end forall
      end


      SUBROUTINE Init_F_DDJk(kn,F_Jk,F_DDJk)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:kdirac,pai
      !IMPLICIT NONE
      real*8,intent(in):: kn(3),F_Jk(3,3)
      REAL*8,intent(out):: F_DDJk(3,3,3,3)

      integer i,j,m,l

      forall(i=1:3,j=1:3,m=1:3,l=1:3)
         F_DDJk(i,j,m,l)=-kn(m)*kn(l)*F_Jk(i,j)
      end forall
      end

      SUBROUTINE Init_F_DRk(F_DDJk,F_DRk)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:EPS
      !IMPLICIT NONE
      real*8,intent(in):: F_DDJk(3,3,3,3)
      REAL*8,intent(out):: F_DRk(3,3,3)

      integer i,j,m,n,l
         
      F_DRk=0.0_8
      do m=1,3
       do n=1,3
        forall(i=1:3,j=1:3,l=1:3)
         F_DRk(i,j,l)=F_DRk(i,j,l)-0.5_8*EPS(j,m,n)*F_DDJk(i,n,l,m)
        end forall
       enddo
      enddo
      end



      SUBROUTINE Init_F_DKk(F_DDJk,F_DKk)
      !use LATTICE_BASE,only:lammda
     ! use tensors, only:EPS
     ! IMPLICIT NONE
      real*8,intent(in):: F_DDJk(3,3,3,3)
      REAL*8,intent(out):: F_DKk(3,3,3,3)

      integer i,j,k,l

        forall(i=1:3,j=1:3,k=1:3,l=1:3)
         F_DKk(i,j,k,l)=0.5_8*(F_DDJk(i,j,l,k)+F_DDJk(i,k,l,j))
        end forall
  
      end

      SUBROUTINE Init_F_D2Jk(dist,F_Jk,F_D2Jk)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:EPS
      !IMPLICIT NONE
      real*8,intent(in):: dist,F_Jk(3,3)
      REAL*8,intent(out):: F_D2Jk(3,3)

      integer i,j

        forall(i=1:3,j=1:3)
         F_D2Jk(i,j)=-dist*dist*F_Jk(i,j)
        end forall
  
      end

      SUBROUTINE Init_F_DD2Kk(kn,dist,F_Jk,F_DD2Kk)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:EPS
      !IMPLICIT NONE
      real*8,intent(in):: kn(3),dist,F_Jk(3,3)
      REAL*8,intent(out):: F_DD2Kk(3,3,3,3)

      integer i,j,k,l

        forall(i=1:3,j=1:3,k=1:3,l=1:3)
         F_DD2Kk(i,j,k,l)=kn(k)*kn(l)*dist*dist*F_Jk(i,j)
        end forall
  
      end



      SUBROUTINE Init_F_DJk(kn,F_Jk,F_DJk)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:kdirac,pai
      !IMPLICIT NONE
      real*8,intent(in):: kn(3),F_Jk(3,3)
      REAL*8,intent(out):: F_DJk(3,3,3)

      integer i,j,l

      forall(i=1:3,j=1:3,l=1:3)
         F_DJk(i,j,l)=kn(l)*F_Jk(i,j)
      end forall
      end

      SUBROUTINE Init_F_DD2Jk(kn,F_D2Jk,F_DD2Jk)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:kdirac,pai
      !IMPLICIT NONE
      real*8,intent(in):: kn(3),F_D2Jk(3,3)
      REAL*8,intent(out):: F_DD2Jk(3,3,3)

      integer i,j,l

      forall(i=1:3,j=1:3,l=1:3)
         F_DD2Jk(i,j,l)=kn(l)*F_D2Jk(i,j)
      end forall
      end


      SUBROUTINE Init_F_Rk(F_DJk,F_Rk)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:EPS
      !IMPLICIT NONE
      real*8,intent(in):: F_DJk(3,3,3)
      REAL*8,intent(out):: F_Rk(3,3)

      integer i,j,m,n

      F_Rk=0.0_8
      do m=1,3
       do n=1,3
        forall(i=1:3,j=1:3)
         F_Rk(i,j)=F_Rk(i,j)-0.5_8*EPS(j,m,n)*F_DJk(i,n,m)
        end forall
       enddo
      enddo
      end


      SUBROUTINE Init_F_D2Rk(F_DD2Jk,F_D2Rk)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:EPS
      !IMPLICIT NONE
      real*8,intent(in):: F_DD2Jk(3,3,3)
      REAL*8,intent(out):: F_D2Rk(3,3)

      integer i,j,m,n

      F_D2Rk=0.0_8
      do m=1,3
       do n=1,3
        forall(i=1:3,j=1:3)
         F_D2Rk(i,j)=F_D2Rk(i,j)-0.5_8*EPS(j,m,n)*F_DD2Jk(i,n,m)
        end forall
       enddo
      enddo
      end




      SUBROUTINE Init_F_Kk(F_DJk,F_Kk)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:EPS
      !IMPLICIT NONE
      real*8,intent(in):: F_DJk(3,3,3)
      REAL*8,intent(out):: F_Kk(3,3,3)

      integer i,j,k


        forall(i=1:3,j=1:3,k=1:3)
         F_Kk(i,j,k)=0.5_8*(F_DJk(i,j,k)+F_DJk(i,k,j))
        end forall

      end

      SUBROUTINE Init_F_D2Kk(F_DD2Jk,F_D2Kk)
      !use LATTICE_BASE,only:lammda
      !use tensors, only:EPS
      !IMPLICIT NONE
      real*8,intent(in):: F_DD2Jk(3,3,3)
      REAL*8,intent(out):: F_D2Kk(3,3,3)
      integer i,j,k

        forall(i=1:3,j=1:3,k=1:3)
         F_D2Kk(i,j,k)=F_DD2Jk(i,j,k)
        end forall
      end

   end Module Ewald_summation