! WRITES THE CURRENT CONFIGURATION ON FILE

      PROGRAM MAIN
      USE DEM_gz
      USE SIZE           ! NN,NX,NY,NZ
      USE CONFIG         
      USE LATTICE_BASE
      USE FORCE_PAR      ! GAMMA,SIGMA
      USE TENSORS,ONLY:EPS
      implicit none
      INTEGER I,J,K,K_WRITE,SEEDSIZE, INP,iLb,NLB 
      INTEGER,ALLOCATABLE :: SEED(:) 
      REAL*8 T,DT,END_TIME,BASE_DT,T0
      REAL*8 EI_bg(6)
      LOGICAL BROWN, THERE
      real*8 yeta_mu(7)

      CALL INIT_KDirac
      CALL INIT_EPS
      CALL INIT_Y2()
      !call INIT_Y21()
      CALL INIT_II


      
      INQUIRE( FILE='control_file.dat', EXIST=THERE ) 
       IF ( THERE ) THEN
       OPEN(10,FILE='control_file.dat')
       read(10,*) IsCaseLattice
       read(10,*) NLB
       read(10,*) LB0
       READ(10,*) GAMMA
       READ(10,*) T0 
       read(10,*) u_bg
       read(10,*) omega_bg
       read(10,*) EI_bg
       READ(10,*) LINERR
       READ(10,*) BROWN
       READ(10,*) muCalculate
       IF ( BROWN ) THEN
        CALL RANDOM_SEED(SIZE=SEEDSIZE)
        ALLOCATE(SEED(SEEDSIZE))
        READ(10,*) INP
        SEED = INP
        CALL RANDOM_SEED(PUT=SEED)
       ENDIF
       READ(10,*) END_TIME
       READ(10,*) BASE_DT
       READ(10,*) K_WRITE
       write(*,*) 'K_write=',K_WRITE
       READ(10,*) LJ_EPS
       READ(10,*) LJ_SIGMA
       READ(10,*) LJ_CUT
       CLOSE(10)
      ELSE
       WRITE(0,*) "error: no control_file.dat file "  &
        // "present in this folder"
       WRITE(0,*) "        Please provide one"
       CALL EXIT()
      ENDIF

      ! background Einf11,Einf22,Einf33,Einf12,Einf13,Einf23

      EI_bg(5)=GAMMA*0.5_8
      Eij=0.0_8
      Eij(1,1)=EI_bg(1)
      Eij(2,2)=EI_bg(2)
      Eij(3,3)=EI_bg(3)
      Eij(1,2)=EI_bg(4)
      Eij(1,3)=EI_bg(5)
      Eij(2,3)=EI_bg(6)
      Eij(2,1)=EI_bg(4)
      Eij(3,1)=EI_bg(5)
      Eij(3,2)=EI_bg(6)
      
      omega_bg(2)=GAMMA*0.5_8
      omegaT=0.0_8
      do I=1,3
        do J=1,3
        omegaT(I,J)=EPS(J,I,1)*omega_bg(1)+ &
          EPS(J,I,2)*omega_bg(2)+EPS(J,I,3)*omega_bg(3)
        enddo
      enddo
      write(*,*)'Eij=',Eij
      write(*,*)'omegaT=',omegaT


      DT = BASE_DT
      
      CALL INIT_VMD_WRITE()

      !klb=0
      !LB0=(/2.01:0.05:3,3.1:0.5:8/)!1.9_8
      LB=0.0_8
      !LB0=10.0_8
      do iLb =1,NLB !100

      T=T0
      K=0

      if(LB0.lt.2.1_8)then
      LB0=LB0+0.002_8
      else if(LB0.lt.3.5_8)then
      LB0=LB0+0.05_8
      else if(LB0.lt.10.0_8)then
      LB0=LB0+0.5_8
      else
      LB0=LB0+2.0_8
      endif

      CALL INITIAL_CONF()


      !u_bg(1)=-GAMMA*0.5_8*LB(3)
      
      WRITE(*,*) 'N: ',NN,' dt: ',BASE_DT,' shear: ',GAMMA,  &
      ' X: ',LB(1),' Y: ',LB(2), &
      ' Z: ',LB(3)

      !CALL INIT_VMD_WRITE(RADII,DT)


      DO

       DO I = 1, NN
        CALL PERIOD_CORRECTION(CONF(:,I),LB,GAMMA,T)
       ENDDO
       
       CALL STEPPER_NOBROWN(CONF,RADII,DT,T,yeta_mu)
       
       !====================== DEM ========================!       
       call DEM_calc_initial_position(CONF,LB,RADII*1.0_8,NN,10000)
       !====================== DEM ========================!
              
      IF (mod(K,K_WRITE) .eq. 0) then
        CALL VMD_WRITE(CONF,T,RADII)
        !Call xyz_write(CONF,T)
        write(20,*)T,yeta_mu,NI,NR

       !ENDIF
     ENDIF


      K=K+1
      T=T+DT

       IF (T.GT.END_TIME) THEN
        EXIT
       ENDIF
          

      write(*,*)'T=',T


      ENDDO

      DEALLOCATE( CONF )
      DEALLOCATE( RADII)

      enddo
     
      CALL END_VMD_WRITE()

END


!***********************************************************
!***********************************************************
!***********************************************************

! VTF file for vmd reader header

      SUBROUTINE INIT_VMD_WRITE()
      USE SIZE
      USE LATTICE_BASE
      USE FORCE_PAR
      IMPLICIT NONE
      !REAL*8,intent(in)::  RADII(NN),DT
      !INTEGER I

      OPEN(31,FILE='output.VTF') ! xyz file for vmd
      open(20,file='yeta_mu.txt',position='append')

      !WRITE(31,*) '# ',NN*1,' ',L0_PAR,' ',K_PAR,' ',LJ_SIGMA, &
      !           ' ',LJ_EPS
      !WRITE(31,*) '# ',LB(1),' ',LB(2),' ',LB(3),' ',DT,' ',GAMMA
      !DO I=1,NN
      ! WRITE(31,*) 'atom ',I-1,'radius',RADII(I),'name H'
      !ENDDO

      !WRITE(31,'(A)') ' '
      !WRITE(31,*) 'pbc',LB(1),LB(2),LB(3)
      !WRITE(31,*)

      RETURN
      END

!***********************************************************
!***********************************************************
!***********************************************************

! VTF file for vmd reader timestep write

      SUBROUTINE VMD_WRITE(CONF,T,RADII)
      USE SIZE
      IMPLICIT NONE
      REAL(8),intent(in):: CONF(3,NN),T,RADII(NN)
      !INTEGER,intent(out):: K
      INTEGER I

      !WRITE(31,'(A)') " "
      !WRITE(31,*) "timestep"
      !WRITE(31,*) "#",T,K
      DO I = 1, NN
       WRITE(31,*) T,CONF(:,I),RADII(I)
      ENDDO

      RETURN
      END

!***********************************************************
!***********************************************************
!***********************************************************

! close VTF file for vmd reader

      SUBROUTINE END_VMD_WRITE()
      IMPLICIT NONE

      CLOSE(31) ! xyz file for vmd
      close(20)
      RETURN
      END

      
