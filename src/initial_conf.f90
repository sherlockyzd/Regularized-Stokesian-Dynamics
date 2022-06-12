

  module Init_IO_MOD 
  USE SIZE
  use method
  use control
  USE CONFIG         
  USE LATTICE_BASE
  USE FORCE_PAR      ! LJ
  use SYS_property
  use tensors,only:pai

  implicit none
  private



  public::INIT_READ_CONTROL,INIT_VMD_WRITE,VMD_WRITE,yetamu_WRITE,END_VMD_WRITE, &
      & INITIAL_SIZE,INITIAL_LB0,INITIAL_CONF

  contains
      subroutine INIT_READ_CONTROL()

      LOGICAL::THERE

      INQUIRE( FILE='control_file.dat', EXIST=THERE ) 
       IF ( THERE ) THEN
       OPEN(10,FILE='control_file.dat')
       read(10,*) correction_method
       read(10,*) FTS_method       
       read(10,*) uselub
       read(10,*) useDEM 
       read(10,*) usecollision
       read(10,*) useGMRESmethod
       read(10,*) Xshear       
       read(10,*) IsPeriod
       read(10,*) Oscillation_shear
       read(10,*) fix
       read(10,*) IsCaseLattice
       read(10,*) beginstep
       read(10,*) interval
       read(10,*) NUM_THREADS
       read(10,*) timestep
       !read(10,*) NLB
       !read(10,*) LB0
       read(10,*) radii_test
       read(10,*) u_bg
       read(10,*) omega_bg
       read(10,*) EI_bg
       read(10,*) mu_f
       read(10,*) pho_f
       READ(10,*) Re_number
       READ(10,*) Gamma_alpha
       read(10,*) ksp
       read(10,*) erest
       read(10,*) pho_par 
       READ(10,*) LINERR
       READ(10,*) BROWN
       READ(10,*) collision_condition
       READ(10,*) muCalculate
       IF ( BROWN ) THEN
        CALL RANDOM_SEED(SIZE=SEEDSIZE)
        ALLOCATE(SEED(SEEDSIZE))
        READ(10,*) INP
        SEED = INP
        CALL RANDOM_SEED(PUT=SEED)
       ENDIF
       READ(10,*) T0 
       READ(10,*) END_TIME
       READ(10,*) BASE_DT
       READ(10,*) K_WRITE
      ! write(*,*) 'K_write=',K_WRITE
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

     end subroutine INIT_READ_CONTROL

! VTF file for vmd reader header

      SUBROUTINE INIT_VMD_WRITE()
      IMPLICIT NONE
      !REAL*8,intent(in)::  RADII(NN),DT
      !INTEGER I

      open(20,file='yeta_mu.txt')
      OPEN(31,FILE='output.VTF') ! xyz file for vmd
      open(50,file='SYS_property.txt')

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

      SUBROUTINE VMD_WRITE(CONF,T,RADII,U_pos)
      IMPLICIT NONE
      REAL(8),intent(in):: CONF(3,NN),T,RADII(NN),U_pos(3*NN)
      INTEGER I
      DO I = 1, NN
       write(31,"(8ES24.15)")T,CONF(:,I),RADII(I),U_pos(3*(I-1)+1),U_pos(3*(I-1)+2),U_pos(3*(I-1)+3)
      ENDDO

      RETURN
      END


      SUBROUTINE yetamu_WRITE(T,yeta_mu,LB,frequency,GAMMA,GAMMAangle)
      USE SIZE
      IMPLICIT NONE
      REAL(8),intent(in):: T,yeta_mu(5),LB,frequency,GAMMA,GAMMAangle
      
       write(20,"(11ES24.15)")T,yeta_mu(:),LB,frequency,GAMMA,GAMMAangle

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
      close(50)

      END SUBROUTINE END_VMD_WRITE

!***********************************************************

      subroutine INITIAL_SIZE(T)
      IMPLICIT NONE
      REAL*8,intent(out):: T

      
      INTEGER:: I,istat
      LOGICAL:: THERE



      INQUIRE( FILE='initial_config.dat', EXIST=THERE ) 
      IF ( THERE ) THEN
       OPEN(11,FILE='initial_config.dat')
       READ(11,*) NX,NY,NZ
         if(IsCaseLattice==0) then
           READ(11,*) NN
           READ(11,*) LB
           ALLOCATE( CONF(3,NN),STAT=istat )
           ALLOCATE( RADII(NN),STAT=istat )
           DO I=1,NN
             READ(11,*) T,CONF(1,I),CONF(2,I),CONF(3,I),RADII(I)
           ENDDO
         elseif(IsCaseLattice==1) then      !SC
           NN=NX*NY*NZ
           ALLOCATE( CONF(3,NN) ,STAT=istat)
           ALLOCATE( RADII(NN) ,STAT=istat)

         elseif(IsCaseLattice==2) then        !BCC
           NN=NX*NY*NZ*2
           ALLOCATE( CONF(3,NN) ,STAT=istat)
           ALLOCATE( RADII(NN) ,STAT=istat)

         elseif(IsCaseLattice==3) then        !FCC
           NN=NX*NY*NZ*4
           ALLOCATE( CONF(3,NN) ,STAT=istat)
           ALLOCATE( RADII(NN) ,STAT=istat)
         endif
       CLOSE(11)
      ELSE
       WRITE(0,*) "error: no initial_config.dat file " &
         // "present in this folder"
       WRITE(0,*) "        Please provide one"
       CALL EXIT()
      ENDIF

      END subroutine INITIAL_SIZE


      SUBROUTINE INITIAL_CONF(T)
      IMPLICIT NONE
      REAL*8,intent(out):: T
      
      
      INTEGER I,J,K,IJK,nl,mm,num
      
      real*8 pos(3,3),Ireal1

      T=T0

       nl=NX+0
      if(IsCaseLattice==0) then
         return
       elseif(IsCaseLattice==1) then      !SC
         LB(1)=NX*LB0*radii_test
         LB(2)=NY*LB0*radii_test
         LB(3)=NZ*LB0*radii_test
         conf=0.0_8
         radii=0.0_8
         IJK=0
         DO I=1,NX
           DO J=1,NY
             DO K=1,NZ
               IJK=IJK+1
               CONF(1,IJK)=LB(1)/NX*(real(I-1,8)+0.5_8)
               CONF(2,IJK)=LB(2)/NY*(real(J-1,8)+0.5_8)
               CONF(3,IJK)=LB(3)/NZ*(real(K-1,8)+0.5_8)
               RADII(IJK)=radii_test 
            enddo
          enddo          
         ENDDO

       elseif(IsCaseLattice==2) then        !BCC
         LB(1)=NX*LB0*radii_test
         LB(2)=NY*LB0*radii_test
         LB(3)=NZ*LB0*radii_test
         conf=0.0_8
         radii=0.0_8
         IJK=0
         DO I=1,NX
           DO J=1,NY
             DO K=1,NZ
               IJK=IJK+1
               CONF(1,IJK)=LB(1)/NX*(real(I-1,8))
               CONF(2,IJK)=LB(2)/NY*(real(J-1,8))
               CONF(3,IJK)=LB(3)/NZ*(real(K-1,8))
               RADII(IJK)=1.0_8*radii_test
            enddo
          enddo          
         ENDDO

         DO I=1,NX
           DO J=1,NY
             DO K=1,NZ
              IJK=IJK+1
              if(i+j+k==3) then
              CONF(1,IJK)=0.5_8*LB(1)/NX
              CONF(2,IJK)=0.5_8*LB(2)/NY
              CONF(3,IJK)=0.5_8*LB(3)/NZ
              else
               CONF(1,IJK)=0.5_8*LB(1)/NX+LB(1)/NX*(real(I-1,8))!*0.5_8*nl
               CONF(2,IJK)=0.5_8*LB(2)/NY+LB(2)/NY*(real(J-1,8))!*0.5_8*nl
               CONF(3,IJK)=0.5_8*LB(3)/NZ+LB(3)/NZ*(real(K-1,8))!*0.5_8*nl
              endif
               RADII(IJK)=1.0_8*radii_test
            enddo
          enddo          
         ENDDO
         conf(1,:)=conf(1,:)+LB(1)/NX*0.2_8
         conf(2,:)=conf(2,:)+LB(2)/NY*0.2_8
         conf(3,:)=conf(3,:)+LB(3)/NZ*0.2_8

       elseif(IsCaseLattice==3) then        !FCC
         LB(1)=NX*LB0*radii_test
         LB(2)=NY*LB0*radii_test
         LB(3)=NZ*LB0*radii_test
         conf=0.0_8
         radii=0.0_8

         pos=0.5_8
         pos(1,1)=0
         pos(2,2)=0
         pos(3,3)=0
         IJK=0
         DO I=1,NX
           DO J=1,NY
             DO K=1,NZ
               IJK=IJK+1
               CONF(1,IJK)=LB(1)/NX*(real(I-1,8))
               CONF(2,IJK)=LB(2)/NY*(real(J-1,8))
               CONF(3,IJK)=LB(3)/NZ*(real(K-1,8))
               RADII(IJK)=1.0_8*radii_test
            enddo
          enddo          
         ENDDO

         DO I=1,NX
           DO J=1,NY
             DO K=1,NZ
               num=nx*ny*nz+3*((i-1)*nz*ny+(j-1)*nz+k-1)
               DO mm=1,3
               IJK=num+mm
               call Logtoreal(Ireal1,(i==1))
               CONF(1,IJK)=LB(1)/NX*(real(I-1,8)+Ireal1*pos(mm,1) &
               +(1.0_8-Ireal1)*pos(mm,1))
               call Logtoreal(Ireal1,(j==1))
               CONF(2,IJK)=LB(2)/NY*(real(J-1,8)+Ireal1*pos(mm,2) &
               +(1.0_8-Ireal1)*pos(mm,2))
               call Logtoreal(Ireal1,(k==1))
               CONF(3,IJK)=LB(3)/NZ*(real(K-1,8)+Ireal1*pos(mm,3) &
               +(1.0_8-Ireal1)*pos(mm,3))
               RADII(IJK)=1.0_8*radii_test
               !write(*,*) CONF(:,IJK)
               enddo  
            enddo
          enddo          
         ENDDO
         conf(1,:)=conf(1,:)+LB(1)/NX*0.2_8
         conf(2,:)=conf(2,:)+LB(2)/NY*0.2_8
         conf(3,:)=conf(3,:)+LB(3)/NZ*0.2_8

       endif

      END     

!********************************************************


      SUBROUTINE INITIAL_LB0
      IMPLICIT NONE

      integer::i,nt,istatus

      if(IsCaseLattice==1) then
        length_Lb0=31
        phi_test=0.52_8
      elseif(IsCaseLattice==2) then
        length_Lb0=38
        phi_test=0.66_8
      elseif(IsCaseLattice==3) then
        length_Lb0=41
        phi_test=0.72_8
      else
        length_Lb0=1
      endif

      ALLOCATE(Lb0_list(length_Lb0),stat=istatus)

      if(IsCaseLattice.ne.0) then
        nt=length_Lb0-6
        Lb0_list(1:5)=(/(i,i=0,4)/)*1.0d0/4*(0.018_8-0.002_8)+0.002_8
        Lb0_list(6:length_Lb0)=(/(i,i=0,nt)/)*1.0d0/nt*(phi_test-0.02_8)+0.02_8
        if(IsCaseLattice==1) then
          Lb0_list=(1.d0/Lb0_list*4.0_8/3.0_8*pai)**(1/3.d0)
        elseif(IsCaseLattice==2) then
          Lb0_list=(2.d0/Lb0_list*4.0_8/3.0_8*pai)**(1/3.d0)
        elseif(IsCaseLattice==3) then
          Lb0_list=(4.d0/Lb0_list*4.0_8/3.0_8*pai)**(1/3.d0)
        endif
        !do i=1,length_Lb0
        !write(*,*) Lb0_list(i)
        !enddo
      endif

      end SUBROUTINE INITIAL_LB0

  end module Init_IO_MOD 




