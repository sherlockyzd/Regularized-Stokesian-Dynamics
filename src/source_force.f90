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
      USE SIZE,only: NN      ! NN
      USE method,only: Isperiod
      use SYS_property,only:ksp,erest,pho_par
      use tensors,only:pai2,pai
      IMPLICIT NONE
      real*8,intent(in)::RADII(NN),CONF(3,NN),U_pos(3*NN)
      real*8,intent(out)::ppiclf_F(3,NN)

      INTEGER I,J
      real*8 rthresh, rx(3),rxdiff, rydiff, rzdiff, rdiff, rm1, rm2
      real*8 rmult, eta, rbot, rn_12x, rn_12y, rn_12z, rdelta12
      real*8 rv12_mag, rv12_mage, rksp_max, rnmag


      ppiclf_F=0.0_8

      do i=1,NN
        do j=i+1,NN
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




      SUBROUTINE source_F(conf,RADII,U_pos,Fe)
      use size,only:NN
      use method,only:usecollision
      IMPLICIT NONE
      real*8,intent(in)::RADII(NN),CONF(3,NN),U_pos(3*NN)
      real*8,intent(out)::Fe(6*NN)

      real*8 ppiclf_F(3,NN),F(3*NN),F_rep(3*NN)
      integer i,j


      Fe=0.0_8
      F_rep=0.0_8
      !do I=1,NN
       !  Fe(3*(i-1)+2)=-1.0_8
      !enddo
      !call FORCE_PER(Fe,CONF)
      F=0.0_8

      if(usecollision)then
        call ppiclf_collision(conf,RADII,U_pos,ppiclf_F)

        forall(i=1:NN,j=1:3)
           F(3*(i-1)+j)=ppiclf_F(j,i)
        end forall
      endif
      !call repulsiveforce(NN,conf,radii,F_rep)

      Fe(1:3*NN)=Fe(1:3*NN)+F !+F_rep!


      end

