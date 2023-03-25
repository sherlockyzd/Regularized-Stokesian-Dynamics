      SUBROUTINE FORCE_PER(FEXT,CONFL)
      USE SIZE        ! NN
      USE CONFIG      ! CONFL
      USE FORCE_PAR   ! k_PAR, L0_PAR, A_PAR
      use period_bdy_tools,only:PER_SKEW_CORRECTION
      IMPLICIT NONE
      REAL*8,intent(out):: FEXT(3*NN)
      REAL*8,intent(in)::CONFL(3,NN)
      
      REAL*8 FLJ(3*NN)    ! STRECHING and BENDING FORCES
      REAL*8 L(1:NN+1),LINV(1:NN+1),T(3,NN+2)
      REAL*8 L0,T0(3),FIJ(6)
      INTEGER I,J

      T=0.D0
      L   =0.D0
      LINV=0.D0
      DO I=1,NN-1
       T(:,I)=CONFL(:,I+1)-CONFL(:,I)
       CALL PER_SKEW_CORRECTION(T(:,I),L(I))      
        IF (L(I).GT.0.D0) THEN
         T(:,I)=T(:,I)/L(I)
        ENDIF
         LINV(I)=1.D0/L(I)
      ENDDO

      FLJ=0.D0
 
      DO I=1,NN-1
       DO J=I+1,NN
        T0=CONFL(:,J)-CONFL(:,I)
        CALL PER_SKEW_CORRECTION(T0,L0)
        CALL FORCE_LJ(T0,FIJ)
        FLJ(3*(I-1)+1:3*I) = FLJ(3*(I-1)+1:3*I) + FIJ(1:3)
        FLJ(3*(J-1)+1:3*J) = FLJ(3*(J-1)+1:3*J) + FIJ(4:6)
       ENDDO
      ENDDO

!C TOTAL FORCE -------------------------------------------------

      FEXT=FLJ

      RETURN
      END


!***********************************************************
!***********************************************************
!***********************************************************

      SUBROUTINE FORCE_LJ(R,F)
      USE FORCE_PAR   ! LJ_SIGMA, LJ_EPS 
      IMPLICIT NONE
      REAL*8,intent(in)::R(3)
      REAL*8,intent(out)::F(6)
      
      REAL*8 T0(3),L0

      T0=R
      L0=SQRT(SUM(T0**2))
      IF(L0.LE.LJ_CUT) THEN
       T0=T0/L0
       F(1:3) = &
         +4.D0*LJ_EPS/LJ_SIGMA*( &
           12.D0*(LJ_SIGMA/L0)**13 &
          -6.D0*(LJ_SIGMA/L0)**7 &
          )*T0
       F(4:6)= &
         -4.D0*LJ_EPS/LJ_SIGMA*( &
           12.D0*(LJ_SIGMA/L0)**13 &
          -6.D0*(LJ_SIGMA/L0)**7 &
          )*T0
      ELSE
       F=0.D0
      ENDIF

      RETURN
      END

      subroutine  repulsiveforce(ntotal,p_pos,radii,Frp)
      use prutil,only:cp
      implicit none
      integer,intent(in) :: ntotal
      real(kind=cp),dimension(3,ntotal),intent(in) :: p_pos
      real(kind=cp),dimension(ntotal),intent(in) :: radii
      real(kind=cp),dimension(3,ntotal),intent(out)  :: Frp
      integer :: i,j
      real(kind=cp) :: r,r1,r2,r3,tau,epsila,Fr,F0,r_ave
      !real(kind=cp),dimension(ntotal,ntotal):: Fr
      Frp=0.0_cp
      Fr=0.0_cp
      tau=20.0_cp
      F0=0.05_cp

        !write(*,*) 'r_ave=',r_ave
      do i=1,ntotal
          do j=i+1,ntotal
              r1=p_pos(1,j)-p_pos(1,i)
              r2=p_pos(2,j)-p_pos(2,i)
              r3=p_pos(3,j)-p_pos(3,i)
              r=sqrt(r1*r1+r2*r2+r3*r3)
              r_ave=0.5_8*(radii(i)+radii(j))
              epsila=(r-2.0_8*r_ave)/r_ave
              !write(*,*) 'epsila=',epsila
              if (epsila.lt. 2.5_8) then
                 Fr=F0*tau*exp(-tau*epsila)/(1-exp(-tau*epsila))
                 Frp(1,i)=Frp(1,i)-Fr*r1/r
                 Frp(2,i)=Frp(2,i)-Fr*r2/r
                 Frp(3,i)=Frp(3,i)-Fr*r3/r
                 Frp(1,j)=Frp(1,j)+Fr*r1/r
                 Frp(2,j)=Frp(2,j)+Fr*r2/r
                 Frp(3,j)=Frp(3,j)+Fr*r3/r
              end if
            end do
        end do
        end subroutine repulsiveforce


      SUBROUTINE ppiclf_collision(conf,RADII,U_pos,ppiclf_F)
      USE SIZE,only: NN,Nb     ! NN
      USE method,only: Isperiod
      use SYS_property,only:ksp,erest,pho_par
      use tensors,only:pai2,pai
      use CONFIG,only:Floc_index
      IMPLICIT NONE
      real*8,intent(in)::RADII(NN),CONF(3,NN),U_pos(6*NN)
      real*8,intent(out)::ppiclf_F(3,NN)

      INTEGER I,J
      real*8 rthresh, rx(3),rxdiff, rydiff, rzdiff, rdiff, rm1, rm2
      real*8 rmult, eta, rbot, rn_12x, rn_12y, rn_12z, rdelta12
      real*8 rv12_mag, rv12_mage, rksp_max, rnmag


      ppiclf_F=0.0_8

      do i=1,NN

          if(i.gt.NN-Nb) then
         ! write(*,*) 'boundayr pass-------------------------'
           cycle
          endif
          
        do j=i+1,NN
         
          if(Floc_index(i).eq.0.and.Floc_index(j).eq.0) then
           cycle
          endif
       !if (j .ne. 0) then
         rthresh  = RADII(i) + RADII(j)
         rx=conf(:,j)-conf(:,i)
 
        !if(IsPeriod) then
         ! CALL PER_SKEW_CORRECTION(rx,rdiff)
        !else
          rdiff=sqrt(sum(rx*rx))
        !endif
 
        rxdiff = rx(1)
        rydiff = rx(2)
        rzdiff = rx(3)


         if (rdiff .lt. rthresh) then
          write(*,*) 'collision Contact overlapping, delta='
          write(*,*) rdiff-rthresh,i,j         
         rm1 = pho_par*4.0_8/3.0_8*pai*RADII(i)**3
         rm2 = pho_par*4.0_8/3.0_8*pai*RADII(j)**3
         
         rmult = 1.0d0/sqrt(1.0d0/rm1+1.0d0/rm2)
         eta   = 2.0d0*sqrt(ksp)*log(erest)/sqrt(log(erest)**2+pai2)*rmult
         
         rbot = 1.0d0/rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         rn_12z = rzdiff*rbot
         
         rdelta12 = rthresh - rdiff
         
         rv12_mag=(U_pos(3*(i-1)+1)-U_pos(3*(j-1)+1))*rn_12x &
                 +(U_pos(3*(i-1)+2)-U_pos(3*(j-1)+2))*rn_12y &
                 +(U_pos(3*i)-U_pos(3*j))*rn_12z

         rv12_mage = rv12_mag*eta
         rksp_max  = ksp*rdelta12
         rnmag     = -rksp_max - rv12_mage
         
         ppiclf_F(1,i) = ppiclf_F(1,i)+ rnmag*rn_12x
         ppiclf_F(2,i) = ppiclf_F(2,i)+ rnmag*rn_12y
         ppiclf_F(3,i) = ppiclf_F(3,i)+ rnmag*rn_12z
         ppiclf_F(1,j) = ppiclf_F(1,j)- rnmag*rn_12x
         ppiclf_F(2,j) = ppiclf_F(2,j)- rnmag*rn_12y
         ppiclf_F(3,j) = ppiclf_F(3,j)- rnmag*rn_12z

         endif

         enddo
       enddo

      end subroutine ppiclf_collision


      SUBROUTINE GAUSS_nb(G)
      USE TENSORS,ONLY:PAI
      IMPLICIT NONE

      REAL*8 G
      REAL*8 A,B

      CALL RANDOM_NUMBER(A)
      CALL RANDOM_NUMBER(B)

      G=SQRT(-2*LOG(A)) * COS(2*Pai*B)

      RETURN
      END

      SUBROUTINE GAUSS_VB(APP,VB)
      use size,only:NN
      !USE TENSORS,ONLY:PAI
      IMPLICIT NONE
      real*8,intent(in)::APP(6*NN,6*NN)
      real*8,intent(out)::VB(3*NN)

      integer i
      real*8 A(3*NN,3*NN),G(3*NN)

      A=APP(1:3*NN,1:3*NN)
      CALL CHOLESKY(A,3*NN)
      DO I=1,3*NN
       CALL GAUSS_nb(G(I))
      ENDDO
      VB=MATMUL(A,G)
      end




      SUBROUTINE source_inter_F(conf,RADII,U_pos,Fe)
      !use SYS_property,only:gravity
      use size,only:NN,Nb
      use method,only:usecollision,wall_method,useDLVO,usebond
      use filament,only:F_rb,solve_implicit
      use filament_math,only:InternalForcesAndTorques
      IMPLICIT NONE
      real*8,intent(in)::RADII(NN),conf(3,NN),U_pos(6*NN)
      real*8,intent(out)::Fe(6*NN)!,SijN_FP(5*NN)

      real*8 collision_F(3,NN),ppiclf_F_wall(3,NN),F_DLVO(3,NN),bond_F(3,NN)
      real*8 F(3*NN),F_rep(3*NN)!,g
      real*8 Filament_internal_force_torque(6*NN)
      integer i,j,ii


      Fe=0.0_8
      F_rep=0.0_8
      F_DLVO=0.0_8
      bond_F=0.0_8

 

      F=0.0_8

      !SijN_FP=0.0_8
      if(usecollision)then
          ppiclf_F_wall=0.0_8
          collision_F=0.0_8
          call ppiclf_collision(conf,RADII,U_pos,collision_F)
          !call JKR_collision(conf,RADII,collision_F)
          if(usebond)then  
            call FENE_bond(conf,RADII,bond_F)
          endif
          if(wall_method) then
            call ppiclf_collision_wall(conf,RADII,U_pos,ppiclf_F_wall)
          endif

          forall(i=1:NN-Nb,j=1:3)
             F(3*(i-1)+j)=ppiclf_F_wall(j,i)+collision_F(j,i)+bond_F(j,i)
          end forall
      endif
      
      F_DLVO=0.0_8
      if(useDLVO)then
        call DLVO(conf,RADII,F_DLVO)
        forall(i=1:NN-Nb,j=1:3)
           F(3*(i-1)+j)=F(3*(i-1)+j)+F_DLVO(j,i)
        end forall
        !write(*,*) "DLVO-------------------------"
      endif

      !call repulsiveforce(NN,conf,radii,F_rep)

      Fe(1:3*NN)=Fe(1:3*NN)+F !+F_rep!
      do ii=1,NN
          write(*,*) 'i,internal_Fe_torque===',ii,Fe(3*NN+3*(ii-1)+1:3*NN+3*ii)
              !write(*,*) 'i, torue===',ii,Ftotal(3*NN+3*(ii-1)+1:3*NN+3*ii)
      enddo  

      if (F_rb.ne.0.and.(.not.solve_implicit))then
        call InternalForcesAndTorques(NN,conf,Filament_internal_force_torque)            
        Fe=Fe+Filament_internal_force_torque
      endif


      end  SUBROUTINE source_inter_F

      SUBROUTINE source_body_F(conf,RADII,Fe)
      use SYS_property,only:gravity
      use size,only:NN,Nb
      !use method,only:usecollision,wall_method,useDLVO
      IMPLICIT NONE
      real*8,intent(in)::RADII(NN),conf(3,NN)!,U_pos(3*NN)
      real*8,intent(out)::Fe(6*NN)!,SijN_FP(5*NN)

      !real*8 collision_F(3,NN),ppiclf_F_wall(3,NN),F_DLVO(3,NN),bond_F(3,NN)
      !real*8 F(3*NN),F_rep(3*NN)!,g
      integer I


      Fe=0.0_8


      
      !g=2.62e-3
      !g=2.62_8
      do I=1,NN-Nb
         Fe(3*(i-1)+3)=0.0_8-1.0_8*gravity
      enddo

      END SUBROUTINE source_body_F




      SUBROUTINE ppiclf_collision_wall(conf,RADII,U_pos,ppiclf_F)
      USE SIZE,only: NN      ! NN
      use SYS_property,only:ksp,erest,pho_par
      use tensors,only:pai2,pai
      IMPLICIT NONE
      real*8,intent(in)::RADII(NN),CONF(3,NN),U_pos(6*NN)
      real*8,intent(out)::ppiclf_F(3,NN)

      INTEGER I
      real*8 rthresh, rx(3),rxdiff, rydiff, rzdiff, rdiff, rm1
      real*8 rmult, eta, rbot, rn_12x, rn_12y, rn_12z, rdelta12
      real*8 rv12_mag, rv12_mage, rksp_max, rnmag


      ppiclf_F=0.0_8

      do i=1,NN
      !do j=i+1,NN
      !if (j .ne. 0) then
         rthresh  = RADII(i) !+ RADII(j)
         rx=0.0_8
         rx(3)=-conf(3,i)
        rdiff=sqrt(sum(rx*rx))
         if (rdiff .lt. rthresh) then
         rxdiff = rx(1)
         rydiff = rx(2)
         rzdiff = rx(3)
          write(*,*) 'collision Wall Contact overlapping, delta='
          write(*,*) rdiff-rthresh,i        
         rm1 = pho_par*4.0_8/3.0_8*pai*RADII(i)**3
         !rm2 = pho_par*4.0_8/3.0_8*pai*RADII(j)**3
         !rmult = 1.0d0/sqrt(1.0d0/rm1+1.0d0/rm2)
         rmult = sqrt(rm1)
         eta   = 2.0d0*sqrt(ksp)*log(erest)/sqrt(log(erest)**2+pai2)*rmult
         
         rbot = 1.0d0/rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         rn_12z = rzdiff*rbot
         
         rdelta12 = rthresh - rdiff
         
         rv12_mag=(U_pos(3*(i-1)+1))*rn_12x &
                 +(U_pos(3*(i-1)+2))*rn_12y &
                 +(U_pos(3*i))*rn_12z

         rv12_mage = rv12_mag*eta
         rksp_max  = ksp*rdelta12
         rnmag     = -rksp_max - rv12_mage
         
         ppiclf_F(1,i) = ppiclf_F(1,i)+ rnmag*rn_12x
         ppiclf_F(2,i) = ppiclf_F(2,i)+ rnmag*rn_12y
         ppiclf_F(3,i) = ppiclf_F(3,i)+ rnmag*rn_12z
         endif
       enddo

      end subroutine ppiclf_collision_wall




      SUBROUTINE DLVO(conf,RADII,F_DLVO)
      USE SIZE,only: NN,Nb    ! NN
      USE method,only: simplePeriod
      use DLVO_property
      use tensors,only:PAI
      use CONFIG,only:Floc_index
      use period_bdy_tools,only:PER_SKEW_CORRECTION
      IMPLICIT NONE
      real*8,intent(in)::RADII(NN),CONF(3,NN)
      real*8,intent(out)::F_DLVO(3,NN)
      !real*8,intent(inout)::SijN_FP(5*NN)

      INTEGER i,j
      real*8 p0,p,f_hVdW,F_EDL,F_VDW,h
      real*8 rthresh, rx(3),rxdiff, rydiff, rzdiff, rdiff
      real*8 rbot,rnmag,rn_12x,rn_12y,rn_12z,DF1,DF2,DF3
      !real*8 Sij(5,5),Sij_FP_part(5)

      F_DLVO=0.0_8
      p0=0.5709_8

      do i=1,NN

          if(i.gt.NN-Nb) then
         ! write(*,*) 'boundayr pass-------------------------'
           cycle
          endif
          
        do j=i+1,NN
         
          if(Floc_index(i).eq.0.and.Floc_index(j).eq.0) then
          !if(Floc_index(i).eq.0.or.Floc_index(j).eq.0) then
           cycle
          endif
       !if (j .ne. 0) then
         rthresh  = RADII(i) + RADII(j)
         rx=conf(:,j)-conf(:,i)
         rdiff=sqrt(sum(rx*rx))
 
        if(simplePeriod) then
          CALL PER_SKEW_CORRECTION(rx,rdiff)   
        endif
 
        rxdiff = rx(1)
        rydiff = rx(2)
        rzdiff = rx(3)

        h= rdiff-rthresh

         if (h .gt. 0.0_8 .and. h .lt. 5.5_8*rthresh) then
          !write(*,*) 'DLVO using, delta='
          !write(*,*) rdiff-rthresh,i,j         
         F_EDL=-2.0_8*PAI*kapa_EDL*rthresh*0.5_8*Permittivity &
              & *Potential*Potential*exp(-kapa_EDL*h)

         p=2.0_8*PAI*h/Lambda_l
         if(p.le.p0)then
            f_hVdW=(1.0_8+3.54_8*p)/(1.0_8+1.77_8*p)**2
         else
            f_hVdW=0.98_8/p-0.434_8/p**2+0.067429_8/p**3
         endif

         F_VdW=AH*rthresh*0.5_8/12.0_8/((h+z0)*(h+z0))*f_hVdW
         
         rbot = 1.0d0/rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         rn_12z = rzdiff*rbot
         rnmag  = (F_EDL+F_VDW)

         DF1=rnmag*rn_12x
         DF2=rnmag*rn_12y
         DF3=rnmag*rn_12z
         F_DLVO(1,i) = F_DLVO(1,i)+ DF1
         F_DLVO(2,i) = F_DLVO(2,i)+ DF2
         F_DLVO(3,i) = F_DLVO(3,i)+ DF3
         F_DLVO(1,j) = F_DLVO(1,j)- DF1
         F_DLVO(2,j) = F_DLVO(2,j)- DF2
         F_DLVO(3,j) = F_DLVO(3,j)- DF3

#ifdef SijN_FP
         Sij=0.0_8
         Sij(1,2)=-0.25*(DF1*rydiff+DF2*rxdiff)
         Sij(1,3)=-0.25*(DF1*rzdiff+DF3*rxdiff)
         Sij(2,3)=-0.25*(DF2*rzdiff+DF3*rydiff)
         Sij(2,1)=Sij(1,2)
         Sij(3,1)=Sij(1,3)
         Sij(3,2)=Sij(2,3)
         call EI_transY2(Sij_FP_part,Sij)

          !Do i=1,NN 
            SijN_FP(5*(i-1)+1:5*i) = SijN_FP(5*(i-1)+1:5*i)+Sij_FP_part(1:5)
            SijN_FP(5*(j-1)+1:5*j) = SijN_FP(5*(j-1)+1:5*j)+Sij_FP_part(1:5)
          !enddo
#endif




         endif

         enddo
       enddo

      end subroutine DLVO


      SUBROUTINE JKR_collision(conf,RADII,collision_F)
      USE SIZE,only: NN,Nb     ! NN
      USE method,only: Isperiod
      use SYS_property,only:ksp,erest,pho_par,F_adh,b0,hp0
      use tensors,only:pai2,pai
      use CONFIG,only:Floc_index
      !use DLVO
      IMPLICIT NONE
      real*8,intent(in)::RADII(NN),CONF(3,NN)!,U_pos(3*NN)
      real*8,intent(out)::collision_F(3,NN)
      !real*8,intent(inout)::SijN_FP(5*NN)

      INTEGER I,J
      real*8 rthresh, rx(3),rxdiff, rydiff, rzdiff, rdiff, rm1, rm2
      real*8  rbot, rn_12x, rn_12y, rn_12z, h_delta12, rnmag,DF1,DF2,DF3
     ! real*8 Sij(5,5),Sij_FP_part(5)
      !real*8 rv12_mag, rv12_mage, rksp_max


      collision_F=0.0_8

      do i=1,NN

          if(i.gt.NN-Nb) then
         ! write(*,*) 'boundayr pass-------------------------'
           cycle
          endif
          
        do j=i+1,NN
         
          if(Floc_index(i).eq.0.and.Floc_index(j).eq.0) then
           cycle
          endif
       !if (j .ne. 0) then
         rthresh  = RADII(i) + RADII(j)
         rx=conf(:,j)-conf(:,i)
 
        !if(IsPeriod) then
         ! CALL PER_SKEW_CORRECTION(rx,rdiff)
        !else
          rdiff=sqrt(sum(rx*rx))
        !endif
 
        rxdiff = rx(1)
        rydiff = rx(2)
        rzdiff = rx(3)


         if (rdiff .lt. (rthresh+hp0) )then
          write(*,*) 'JKR_collision method using, delta====='
          write(*,*) rdiff-rthresh,i,j         
         !rm1 = pho_par*4.0_8/3.0_8*pai*RADII(i)**3
         !rm2 = pho_par*4.0_8/3.0_8*pai*RADII(j)**3
         
         !rmult = 1.0d0/sqrt(1.0d0/rm1+1.0d0/rm2)
         !eta   = 2.0d0*sqrt(ksp)*log(erest)/sqrt(log(erest)**2+pai2)*rmult
         
         rbot = 1.0d0/rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         rn_12z = rzdiff*rbot

         
         h_delta12 = rdiff-rthresh
         rnmag     = -(-1.0_8+1.10_8*(rthresh*0.5_8/(b0*b0) &
                     & *(hp0-h_delta12))**(5.0_8/3.0_8))*abs(F_adh)


         !rv12_mag=(U_pos(3*(i-1)+1)-U_pos(3*(j-1)+1))*rn_12x &
          !       +(U_pos(3*(i-1)+2)-U_pos(3*(j-1)+2))*rn_12y &
          !       +(U_pos(3*i)-U_pos(3*j))*rn_12z

         !rv12_mage = rv12_mag*eta
         !rksp_max  = ksp*rdelta12
         
         DF1=rnmag*rn_12x
         DF2=rnmag*rn_12y
         DF3=rnmag*rn_12z
         collision_F(1,i) = collision_F(1,i)+ DF1
         collision_F(2,i) = collision_F(2,i)+ DF2
         collision_F(3,i) = collision_F(3,i)+ DF3
         collision_F(1,j) = collision_F(1,j)- DF1
         collision_F(2,j) = collision_F(2,j)- DF2
         collision_F(3,j) = collision_F(3,j)- DF3

#ifdef SijN_FP
         Sij=0.0_8
         Sij(1,2)=-0.25*(DF1*rydiff+DF2*rxdiff)
         Sij(1,3)=-0.25*(DF1*rzdiff+DF3*rxdiff)
         Sij(2,3)=-0.25*(DF2*rzdiff+DF3*rydiff)
         Sij(2,1)=Sij(1,2)
         Sij(3,1)=Sij(1,3)
         Sij(3,2)=Sij(2,3)
         call EI_transY2(Sij_FP_part,Sij)

          !Do i=1,NN 
            SijN_FP(5*(i-1)+1:5*i) = SijN_FP(5*(i-1)+1:5*i)+Sij_FP_part(1:5)
            SijN_FP(5*(j-1)+1:5*j) = SijN_FP(5*(j-1)+1:5*j)+Sij_FP_part(1:5)
#endif         !enddo
         


         endif

         enddo
       enddo

      end subroutine JKR_collision


      SUBROUTINE FENE_bond(conf,RADII,bond_F)
      USE SIZE,only: NN,Nb     ! NN
      !USE method,only: Isperiod
      !use SYS_property,only:ksp,erest,pho_par,F_adh,b0,hp0
      use SYS_property,only:hp0
      !8use tensors,only:pai2,pai
      use CONFIG,only:Floc_index
      use BROWN_proerty,only:k_LJ
      !use DLVO
      IMPLICIT NONE
      real*8,intent(in)::RADII(NN),CONF(3,NN)!,U_pos(3*NN)
      real*8,intent(out)::bond_F(3,NN)

      INTEGER I,J
      real*8 rthresh, rx(3),rxdiff, rydiff, rzdiff, rdiff, rm1, rm2
      real*8  rbot, rn_12x, rn_12y, rn_12z, h_delta12, rnmag,sigma_LJ,R0_LJ
      !real*8 rv12_mag, rv12_mage, rksp_max


      bond_F=0.0_8

      do i=1,NN-Nb
          !if(i.gt.NN-Nb) then
         ! write(*,*) 'boundayr pass-------------------------'
           !cycle
          !endif
        do j=i+1,NN
         
          if(Floc_index(i).eq.0.and.Floc_index(j).eq.0) then
           cycle
          endif
          if(Floc_index(i).ne.0.and.Floc_index(j).ne.0) then
           cycle
          endif
       !if (j .ne. 0) then
         rthresh  = RADII(i) + RADII(j)
         sigma_LJ=rthresh/2.0_8
         R0_LJ=0.5_8*sigma_LJ;
         !k_LJ=30.0_8*kBT/sigma_LJ/sigma_LJ;
         rx=conf(:,j)-conf(:,i)
 
        !if(IsPeriod) then
         ! CALL PER_SKEW_CORRECTION(rx,rdiff)
        !else
          rdiff=sqrt(sum(rx*rx))
        !endif
        rxdiff = rx(1)
        rydiff = rx(2)
        rzdiff = rx(3)


         if (rdiff .lt. (rthresh+R0_LJ) )then
          write(*,*) 'LJ_bond method using, delta='
          write(*,*) rdiff-rthresh,i,j         
         !rm1 = pho_par*4.0_8/3.0_8*pai*RADII(i)**3
         !!rm2 = pho_par*4.0_8/3.0_8*pai*RADII(j)**3
         
         !rmult = 1.0d0/sqrt(1.0d0/rm1+1.0d0/rm2)
         !eta   = 2.0d0*sqrt(ksp)*log(erest)/sqrt(log(erest)**2+pai2)*rmult
         
         rbot = 1.0d0/rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         rn_12z = rzdiff*rbot

         
         h_delta12 = rdiff-rthresh-hp0
         rnmag     = k_LJ*h_delta12/(1-(h_delta12/R0_LJ)**2)

         !rv12_mag=(U_pos(3*(i-1)+1)-U_pos(3*(j-1)+1))*rn_12x &
          !       +(U_pos(3*(i-1)+2)-U_pos(3*(j-1)+2))*rn_12y &
          !       +(U_pos(3*i)-U_pos(3*j))*rn_12z

         !rv12_mage = rv12_mag*eta
         !rksp_max  = ksp*rdelta12
         
         
         bond_F(1,i) = bond_F(1,i)+ rnmag*rn_12x
         bond_F(2,i) = bond_F(2,i)+ rnmag*rn_12y
         bond_F(3,i) = bond_F(3,i)+ rnmag*rn_12z
         bond_F(1,j) = bond_F(1,j)- rnmag*rn_12x
         bond_F(2,j) = bond_F(2,j)- rnmag*rn_12y
         bond_F(3,j) = bond_F(3,j)- rnmag*rn_12z

         endif

         enddo
       enddo

      end subroutine FENE_bond








!****************************************************

#ifdef collision
         rksp_wall = ksp

         ! give a bit larger collision threshold for walls
         rextra   = 0.5d0
         rthresh  = (0.5d0+rextra)*rpropi(PPICLF_R_JDP)
         
         rxdiff = yj(PPICLF_JX) - yi(PPICLF_JX)
         rydiff = yj(PPICLF_JY) - yi(PPICLF_JY)
         rzdiff = yj(PPICLF_JZ) - yi(PPICLF_JZ)
         
         rdiff = sqrt(rxdiff**2 + rydiff**2 + rzdiff**2)
         
         if (rdiff .gt. rthresh) return
         
         rm1 = rpropi(PPICLF_R_JRHOP)*rpropi(PPICLF_R_JVOLP)
         
         rmult = sqrt(rm1)
         eta   = 2.0d0*sqrt(rksp_wall)*log(erest)/sqrt(log(erest)**2+rpi2)*rmult
         
         rbot = 1.0d0/rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         rn_12z = rzdiff*rbot
         
         rdelta12 = rthresh - rdiff
         
         rv12_mag = -1.0d0*(yi(PPICLF_JVX)*rn_12x +yi(PPICLF_JVY)*rn_12y
          &        +yi(PPICLF_JVZ)*rn_12z)

         rv12_mage = rv12_mag*eta
         rksp_max  = rksp_wall*rdelta12
         rnmag     = -rksp_max - rv12_mage
         
         ppiclf_ydotc(PPICLF_JVX,i) = ppiclf_ydotc(PPICLF_JVX,i)
          &                          + rnmag*rn_12x
         ppiclf_ydotc(PPICLF_JVY,i) = ppiclf_ydotc(PPICLF_JVY,i)
          &                          + rnmag*rn_12y
         ppiclf_ydotc(PPICLF_JVZ,i) = ppiclf_ydotc(PPICLF_JVZ,i)
          &                          + rnmag*rn_12z
#endif