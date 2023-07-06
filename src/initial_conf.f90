

  module Init_IO_MOD 
  USE SIZE
  use method
  use control
  USE CONFIG
  USE LATTICE_BASE
  USE FORCE_PAR      ! LJ
  use SYS_property
  use DLVO_property
  use BROWN_proerty
  use lambda_Dim
  use tensors,only:pai
  use rb_conglomerate
  use filament
  use filament_implicit2
  use filament_explicit_or_implicit1

  implicit none
  private



  public::INIT_READ_CONTROL,INIT_VMD_WRITE,VMD_WRITE,yetamu_WRITE,END_VMD_WRITE, &
      & INITIAL_SIZE,INITIAL_LB0,INITIAL_CONF,boundary_condition,VMD_WRITE_BDY,&
      & INIT_DIM_WRITE,boundary_correct_condition,INITIAL_LB,INITIAL_CONF_ALLOCATION, &
      & END_CONF_DEALLOCATION

  contains



    subroutine INITIAL_CONF_ALLOCATION()
    IMPLICIT NONE
    INTEGER status

      if(K_rb.ne.0)then
        ALLOCATE (conf_rb(3,K_rb) ,STAT=status)
        ALLOCATE (conf_rb_vector(3,Nswimer) ,STAT=status)
        ALLOCATE (U_par_rb(6*K_rb) ,STAT=status)
        ALLOCATE (q_rb(K_rb) ,STAT=status)
        ALLOCATE (rbmconn_Inertial_body(3*K_rb,3*K_rb) ,STAT=status)
        ALLOCATE (rbmconn_Inertial_body_inverse(3*K_rb,3*K_rb) ,STAT=status)
      endif

      if(F_rb.ne.0)then

        !ALLOCATE (Filament_obj(F_rb) ,STAT=status)
        !if(filament_implicit_method.ne.2) then
        if(.not.filament_solve_implicit.or.filament_implicit_method.eq.1) then
        ALLOCATE (index1(F_rb) ,STAT=status)
        ALLOCATE (Filament_conf_past(3,Nfilament) ,STAT=status)
        ALLOCATE (Filament_conf_now(3,Nfilament) ,STAT=status)
        ALLOCATE (Filament_tau_now(3,Nfilament) ,STAT=status)
        !ALLOCATE (Filament_tau_next(3,Nfilament) ,STAT=status)
        ALLOCATE (Filament_U1_now(3*F_rb) ,STAT=status)
        !ALLOCATE (Filament_X1_next(3*F_rb) ,STAT=status)
        ALLOCATE (Filament_X1_now(3*F_rb) ,STAT=status)
        ALLOCATE (Filament_X1_past(3*F_rb) ,STAT=status)
        !ALLOCATE (Filament_Inertial_lambda(3*Nfilament),STAT=status)
        !ALLOCATE (Filament_Interal_force(6*Nfilament),STAT=status)
        ALLOCATE (Filament_Lie_algebra_now(3*Nfilament),STAT=status)
        !ALLOCATE (Filament_Lie_algebra_next(3*Nfilament),STAT=status)
        ALLOCATE (Filament_q(Nfilament) ,STAT=status)
        endif
      endif

    end subroutine INITIAL_CONF_ALLOCATION



  subroutine END_CONF_DEALLOCATION()
  IMPLICIT NONE
  if(K_rb.ne.0)then 
      if(allocated(KK_rbmconn)) DEALLOCATE(KK_rbmconn)
      if(allocated(conf_rb_vector))  DEALLOCATE( conf_rb_vector)
      if(allocated(U_par_rb)) DEALLOCATE(U_par_rb)
      if(allocated(q_rb)) DEALLOCATE( q_rb)  
      if(allocated(rbmconn_Inertial_body)) DEALLOCATE( rbmconn_Inertial_body)
      if(allocated(rbmconn_Inertial_body_inverse)) DEALLOCATE( rbmconn_Inertial_body_inverse)
      if(allocated(conf_rb))  DEALLOCATE( conf_rb)
  endif
  if(F_rb.ne.0)then
      if(allocated(Filament_obj)) DEALLOCATE(Filament_obj)
      if(allocated(Filament_num)) DEALLOCATE(Filament_num)
      if(allocated(index1)) DEALLOCATE(index1)
      if(allocated(Filament_conf_past))  DEALLOCATE( Filament_conf_past)
      if(allocated(Filament_conf_now))  DEALLOCATE( Filament_conf_now)
      if(allocated(Filament_tau_base)) DEALLOCATE( Filament_tau_base)
      if(allocated(Filament_tau_now)) DEALLOCATE( Filament_tau_now)
      !if(allocated(Filament_tau_next)) DEALLOCATE( Filament_tau_next)
      if(allocated(Filament_tau_now)) DEALLOCATE( Filament_tau_now)
      if(allocated(Filament_U1_now))  DEALLOCATE( Filament_U1_now)
      !if(allocated(Filament_X1_next))  DEALLOCATE( Filament_X1_next)
      if(allocated(Filament_X1_now))  DEALLOCATE( Filament_X1_now)
      if(allocated(Filament_X1_past))  DEALLOCATE( Filament_X1_past)
      !if(allocated(Filament_Interal_force))  DEALLOCATE( Filament_Interal_force)
      !if(allocated(Filament_Inertial_torque))  DEALLOCATE( Filament_Inertial_torque)
      if(allocated(Filament_Lie_algebra_now)) DEALLOCATE(Filament_Lie_algebra_now)
      !if(allocated(Filament_Lie_algebra_next)) DEALLOCATE(Filament_Lie_algebra_next)
      if(allocated(Filament_q))  DEALLOCATE( Filament_q)
  endif

  end subroutine END_CONF_DEALLOCATION



      subroutine INIT_READ_CONTROL()

      LOGICAL::THERE

      INQUIRE( FILE='control_file.dat', EXIST=THERE ) 
       IF ( THERE ) THEN
       OPEN(10,FILE='control_file.dat')
       read(10,*) correction_method
       read(10,*) FTS_method 
       read(10,*) filament_solve_implicit      
       read(10,*) uselub
       read(10,*) useDEM 
       read(10,*) usecollision
       read(10,*) useDEMstepper
       read(10,*) useGMRESmethod
       read(10,*) Xshear       
       read(10,*) IsPeriod
       read(10,*) simplePeriod
       read(10,*) Oscillation_shear
       !read(10,*) useCongl
       read(10,*) fix
       read(10,*) IsCaseLattice
       read(10,*) beginstep
       read(10,*) interval
       read(10,*) NUM_THREADS
       read(10,*) timestep
       !read(10,*) NLB
       !read(10,*) LB0
       read(10,*) radii_test
       read(10,*) radii_true
       read(10,*) g_gravity
       read(10,*) u_bg
       read(10,*) omega_bg
       read(10,*) EI_bg
       read(10,*) mu_f
       read(10,*) pho_f
       READ(10,*) Re_number
       READ(10,*) Gamma_alpha
       READ(10,*) Youngs_module
       READ(10,*) Poisson_ratio
       read(10,*) ksp
       read(10,*) erest
       read(10,*) pho_par 
       READ(10,*) LINERR
       READ(10,*) collision_condition
       READ(10,*) muCalculate
       READ(10,*) T0 
       READ(10,*) END_TIME
       READ(10,*) BASE_DT
       READ(10,*) K_WRITE
      ! write(*,*) 'K_write=',K_WRITE
       READ(10,*) LJ_EPS
       READ(10,*) LJ_SIGMA
       READ(10,*) LJ_CUT
       READ(10,*) Boundary
       READ(10,*) wall_method
       READ(10,*) useDLVO
       READ(10,*) DebyeLength
       READ(10,*) z0
       READ(10,*) Lambda_l
       READ(10,*) Permittivity
       READ(10,*) Potential
       READ(10,*) AH
       READ(10,*) Surfaceenergy
       READ(10,*) BROWN
        IF ( BROWN ) THEN
         READ(10,*) kBT
         CALL RANDOM_SEED(SIZE=SEEDSIZE)
         ALLOCATE(SEED(SEEDSIZE))
         SEED = INP
         CALL RANDOM_SEED(PUT=SEED)
        ENDIF
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
      !open(50,file='SYS_property.txt')
      open(51,file='DIM_property.txt')

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

      SUBROUTINE INIT_DIM_WRITE()
      IMPLICIT NONE
      !REAL*8,intent(in)::  RADII(NN),DT
      !INTEGER I

      WRITE(51,*) 'lambda_pho,lambda_Length,lambda_mass,lambda_Viscosity,lambda_Permittivity,lambda_kB'
      WRITE(51,"(6ES24.15)") lambda_pho,lambda_Length,lambda_mass,lambda_Viscosity,lambda_Permittivity,lambda_kB
      WRITE(51,*) 'lambda_velocity,lambda_acceleration,lambda_time,lambda_Shearrate,lambda_Force'
      WRITE(51,"(5ES24.15)") lambda_velocity,lambda_acceleration,lambda_time,lambda_Shearrate,lambda_Force
      WRITE(51,*) 'lambda_potential,lambda_DebyeLength,lambda_AH,lambda_Surfaceenergy'
      WRITE(51,"(4ES24.15)") lambda_potential,lambda_DebyeLength,lambda_AH,lambda_Surfaceenergy
      WRITE(51,*) 'lambda_temperature,lambda_Pressure,lambda_Elasticmodulu'
      WRITE(51,"(3ES24.15)") lambda_temperature,lambda_Pressure,lambda_Elasticmodulu
      WRITE(51,*) 'lambda_ksp,lambda_erest,lambda_Renold'
      WRITE(51,"(3ES24.15)") lambda_ksp,lambda_erest,lambda_Renold
      WRITE(51,*) 'ksp,pho_par,pho_f,mass_par,gravity'
      WRITE(51,"(5ES24.15)") ksp,pho_par,pho_f,mass_par,gravity
      WRITE(51,*) 'mu_f,END_TIME,BASE_DT,T0,dt_dem,NtDEM,kBT=kBT'
      WRITE(51,"(7ES24.15)") mu_f,END_TIME,BASE_DT,T0,dt_dem0,NtDEM,kBT
      WRITE(51,*) 'F_adh,b0,hp0,heq'
      WRITE(51,"(4ES24.15)") F_adh,b0,hp0,heq
      

      close(51)
      RETURN
      END

!***********************************************************
!***********************************************************
!***********************************************************

! VTF file for vmd reader timestep write

      SUBROUTINE VMD_WRITE(CONF,T,RADII,U_pos)
      IMPLICIT NONE
      REAL(8),intent(in):: CONF(3,NN),T,RADII(NN),U_pos(6*NN)!,SijN(5*NN)
      INTEGER I
      DO I = 1, NN-Nb
       write(31,"(13e24.15e3)")T,CONF(:,I),RADII(I),U_pos(3*(I-1)+1), &
       & U_pos(3*(I-1)+2),U_pos(3*(I-1)+3)
       !,SijN(5*(I-1)+1),SijN(5*(I-1)+2), &
       !& SijN(5*(I-1)+3),SijN(5*(I-1)+4),SijN(5*(I-1)+5)
      ENDDO

      RETURN
      END


      SUBROUTINE VMD_WRITE_BDY(CONF,T,RADII)
      IMPLICIT NONE
      REAL(8),intent(in):: CONF(3,NN),T,RADII(NN)
      INTEGER I
      OPEN(32,FILE='output_bdy.VTF') ! xyz file for vmd
      DO I = NN-Nb+1,NN
       write(32,"(13ES24.15)")T,CONF(:,I),RADII(I)
       !,SijN(5*(I-1)+1),SijN(5*(I-1)+2), &
       !& SijN(5*(I-1)+3),SijN(5*(I-1)+4),SijN(5*(I-1)+5)
      ENDDO
      close(32)
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
      

      END SUBROUTINE END_VMD_WRITE

!***********************************************************

      subroutine boundary_condition()
      INTEGER:: I,istat
      LOGICAL:: THERE
      INQUIRE( FILE='initial_bdy.dat', EXIST=THERE ) 
      IF ( THERE ) THEN
       OPEN(11,FILE='initial_bdy.dat')
       READ(11,*) Nb
       ALLOCATE( CONFbdy(3,Nb),STAT=istat )
       ALLOCATE( RADIIBDY(Nb),STAT=istat )
       DO I=1,Nb
          READ(11,*) CONFBDY(1,I),CONFBDY(2,I),CONFBDY(3,I),RADIIBDY(I)
       ENDDO
       CLOSE(11)
      ELSE
       WRITE(0,*) "error: no initial_bdy.dat file " &
         // "present in this folder"
       WRITE(0,*) "        Please provide one"
       CALL EXIT()
      ENDIF
      end subroutine boundary_condition


      subroutine boundary_correct_condition()
      INTEGER:: I,J,istat,Nbsqrtx,Nbsqrty,K

       Nbsqrtx=ceiling(LB(1)/2.2_8)
       Nbsqrty=ceiling(LB(2)/2.2_8)
       Nb=Nbsqrtx*Nbsqrty
       ALLOCATE( CONFbdy(3,Nb),STAT=istat )
       ALLOCATE( RADIIBDY(Nb),STAT=istat )
       K=0
       DO I=1,Nbsqrtx
         DO J=1,Nbsqrty
           K=K+1
           CONFBDY(1,K)=1.0+(I-1)*2.2
           CONFBDY(2,K)=1.0+(J-1)*2.2
           CONFBDY(3,K)=1.0_8
           RADIIBDY(K)=1.0_8
         ENDDO
       ENDDO

       write(*,*) "CONFBDY(1,K),CONFBDY(2,K),CONFBDY(3,K),RADIIBDY(K)",K
      DO K=1,Nb
        write(*,*) CONFBDY(1,K),CONFBDY(2,K),CONFBDY(3,K),RADIIBDY(K)
      enddo

     !endif
      end subroutine boundary_correct_condition



      subroutine INITIAL_SIZE()
      IMPLICIT NONE
      !REAL*8,intent(out):: T
      
      INTEGER:: I,istat
      LOGICAL:: THERE 
      !T=T0
      INQUIRE( FILE='initial_config.dat', EXIST=THERE ) 
      IF ( THERE ) THEN
       OPEN(11,FILE='initial_config.dat')
       
       READ(11,*) NX,NY,NZ
       call INITIAL_LB()
        if(Boundary) then 
          !call boundary_condition()
          call boundary_correct_condition()
        else
          Nb=0
        endif

        if(IsCaseLattice==0) then
           READ(11,*) Np
           READ(11,*) K_rb,Nswimer
           write(*,*) 'K_rb,Nswimer====',K_rb,Nswimer
           if(K_rb.ne.0)then
              ALLOCATE (KK_rbmconn(K_rb) ,STAT=istat)
              READ(11,*) KK_rbmconn(:)
              if(sum(KK_rbmconn(:)).eq.Nswimer) then
               write(*,*) 'KK_rbmconn',KK_rbmconn(:) 
              else
               write(*,*) 'error--------KK_rbmconn',KK_rbmconn(:) 
               CALL EXIT()
              endif
            else
              READ(11,*) 
            endif

           READ(11,*) F_rb,Nfilament
           write(*,*) 'F_rb,Nfilament====',F_rb,Nfilament
           if(F_rb.ne.0)then
             ALLOCATE (Filament_obj(F_rb) ,STAT=istat)
             ALLOCATE (Filament_num(F_rb) ,STAT=istat)
             READ(11,*) (Filament_num(i),i=1,F_rb)
             if(sum(Filament_num(:)).eq.Nfilament) then
               write(*,*) 'Filament_num',Filament_num(:) 
              else
               write(*,*) 'error--------Filament_num',Filament_num(:) 
               CALL EXIT()
             endif
           else
             READ(11,*)
           endif

           if(Nfilament+Nswimer.gt.Np) then
              write(*,*) 'error:Nswimer,Nfilament,Np====',Nswimer,Nfilament,Np
              CALL EXIT()
           endif
           READ(11,*) alphaX,betaY,gamaZ
           alphaX=alphaX*pai/180.0_8
           betaY=betaY*pai/180.0_8
           gamaZ=gamaZ*pai/180.0_8
           READ(11,*) LB
           NN=Np+Nb
        elseif(IsCaseLattice==1) then      !SC
           Np=NX*NY*NZ
           NN=Np+Nb
        elseif(IsCaseLattice==2) then        !BCC
           Np=NX*NY*NZ*2
           NN=Np+Nb
        elseif(IsCaseLattice==3) then        !FCC
           Np=NX*NY*NZ*4
           NN=Np+Nb
        endif

        CLOSE(11)
      ELSE
         WRITE(0,*) "error: no initial_config.dat file " &
           // "present in this folder"
         WRITE(0,*) "        Please provide one"
         CALL EXIT()
      ENDIF

      END subroutine INITIAL_SIZE



      SUBROUTINE INITIAL_CONF(T,U_pos)
      IMPLICIT NONE
      REAL*8,intent(out):: T
      REAL*8,intent(inout):: U_pos(6*NN)    
      INTEGER I,J,K,IJK,mm,num,istat   
      real*8 pos(3,3),Ireal1
      LOGICAL hot_start
      
      T=T0

      if(IsCaseLattice==0) then
         
        ALLOCATE( CONF(3,NN),STAT=istat )
        ALLOCATE( RADII(NN),STAT=istat )
        OPEN(11,FILE='initial_config.dat')
        READ(11,*)
        READ(11,*)
        READ(11,*)
        READ(11,*)
        READ(11,*)
        READ(11,*)
        READ(11,*)
        READ(11,*)
        READ(11,*) hot_start
        
        if(hot_start) then
            DO I=1,Np
              READ(11,*) T,CONF(1,I),CONF(2,I),CONF(3,I),RADII(I),U_pos(3*(I-1)+1), &
                  & U_pos(3*(I-1)+2),U_pos(3*(I-1)+3)
            ENDDO
        else
            DO I=1,Np
              READ(11,*) T,CONF(1,I),CONF(2,I),CONF(3,I),RADII(I)
            ENDDO
        endif
        if(K_rb.ne.0) then
          call swim_rotate_rb(Nswimer,CONF(:,1:Nswimer),alphaX,betaY,gamaZ)
        endif

        write(*,*) 'hot_start===',hot_start,'conf===================='
        do i=1,Np
            write(*,*) i,CONF(:,i),U_pos(3*(I-1)+1), &
              & U_pos(3*(I-1)+2),U_pos(3*(I-1)+3)
        enddo
        CLOSE(11)

      elseif(IsCaseLattice==1) then      !SC
         ALLOCATE( CONF(3,NN) ,STAT=istat)
         ALLOCATE( RADII(NN) ,STAT=istat)
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
         
         do i=1,Np
           conf(1,i)=conf(1,i)+LB(1)/NX*0.2_8
           conf(2,i)=conf(2,i)+LB(2)/NY*0.2_8
           conf(3,i)=conf(3,i)+LB(3)/NZ*0.2_8
         enddo

      elseif(IsCaseLattice==2) then        !BCC
         ALLOCATE( CONF(3,NN) ,STAT=istat)
         ALLOCATE( RADII(NN) ,STAT=istat)
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
         do i=1,Np
           conf(1,i)=conf(1,i)+LB(1)/NX*0.2_8
           conf(2,i)=conf(2,i)+LB(2)/NY*0.2_8
           conf(3,i)=conf(3,i)+LB(3)/NZ*0.2_8
         enddo

      elseif(IsCaseLattice==3) then        !FCC
         ALLOCATE( CONF(3,NN) ,STAT=istat)
         ALLOCATE( RADII(NN) ,STAT=istat)
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

         do i=1,Np
           conf(1,i)=conf(1,i)+LB(1)/NX*0.2_8
           conf(2,i)=conf(2,i)+LB(2)/NY*0.2_8
           conf(3,i)=conf(3,i)+LB(3)/NZ*0.2_8
         enddo

      endif


        if(useDEM)then
          ALLOCATE( CONF_DEM(3,Np),STAT=istat)
          CONF_DEM=conf(:,1:Np)
        endif
        
   !  write(*,*) "CONF(1,K),CONF(2,K),CONF(3,K),RADII(K)",IJK
    !  DO I=1,Np
     !   write(*,*) CONF(1,I),CONF(2,I),CONF(3,I),RADII(I)
     ! enddo
     ALLOCATE( Floc_index(NN),STAT=istat )
     Floc_index=1
     if(Boundary) then 
        do I=1,Nb
          CONF(1,Np+I)=CONFBDY(1,I)
          CONF(2,Np+I)=CONFBDY(2,I)
          CONF(3,Np+I)=CONFBDY(3,I)
          RADII(Np+I)=RADIIBDY(I)
          Floc_index(Np+I)=0
        enddo
        

        do I=1,Np
        conf(3,I)=conf(3,I)+LB(3)*2.0_8+1.0_8
        enddo
        LB(3)=LB(3)*5.0_8+1.0_8
      endif
      
     ! write(*,*) "CONF_corr(1,K),CONF(2,K),CONF(3,K),RADII(K)",IJK
      !DO I=1,Np
       ! write(*,*) CONF(1,I),CONF(2,I),CONF(3,I),RADII(I)
      !enddo

     ! write(*,*) "CONF_BDY(1,K),CONF_BDY(2,K),CONF_BDY(3,K),RADII(K)",IJK
      !DO I=1,Nb
      !  write(*,*) CONF(1,Np+I),CONF(2,Np+I),CONF(3,Np+I),RADII(I)
      !enddo
      !stop

      
      END SUBROUTINE INITIAL_CONF  

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



      SUBROUTINE INITIAL_LB
      IMPLICIT NONE

      if(IsCaseLattice==0) then
         return
       else     !SC
         LB(1)=NX*LB0*radii_test
         LB(2)=NY*LB0*radii_test
         LB(3)=NZ*LB0*radii_test
         LB_DEM=LB
      endif
       write(*,*) 'LB',LB(2),LB(2),LB(3)
      END SUBROUTINE INITIAL_LB

  end module Init_IO_MOD 




