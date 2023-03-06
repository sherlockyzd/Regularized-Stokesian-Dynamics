! WRITES THE CURRENT CONFIGURATION ON FILE

  PROGRAM MAIN
  use method
  use control
  USE DEM_gz
  USE SIZE           ! NN,NX,NY,NZ
  USE CONFIG         
  USE LATTICE_BASE
  USE SYS_property      ! GAMMA,SIGMA
  use master_time_m, only: tic, toc, time_print
  use Init_IO_MOD 
  use STEPPER_MOD
  use hydro_tools,only:Init_frequency,INIT_u_bg
  USE BROWN_proerty
  use dimentionless
  use lambda_Dim
  use Mobilitymatrix 
  use tensors,only:pai2
  USE OMP_LIB
  use rb_conglomerate
  use conglomerate,only:rbmconn_Init


  implicit none
  INTEGER K_time,iLb,NLB,thread_id,status,I,J
  REAL*8 T,DT,yeta_mu(5),GAMMAangle,ompstart,ompend
  !REAL*8 SijN(5*NN)
  REAL*8,ALLOCATABLE :: U_pos(:),SijN(:)

  write(*,*) 'start-------------'
  call tic()
  call time_print() 
  ompstart= omp_get_wtime()   
  CALL INIT_KDirac()
  CALL INIT_EPS()
  CALL INIT_Y2()
  !call INIT_Y20()
  CALL INIT_II()
  call INIT_READ_CONTROL()
  CALL INIT_VMD_WRITE()
  call dimentionless_calc()
  call dimentionless_translate()
  CALL INIT_DIM_WRITE()


  usebond=.False.
  thread_id = OMP_GET_MAX_THREADS()
  PRINT *, "max processos: ", thread_id

  call OMP_SET_NUM_THREADS(NUM_THREADS)
  PRINT *, "we use processos: ", NUM_THREADS
  DT = BASE_DT
  Nb=0

  call INITIAL_LB0()


  yeta_mu=0.0_8
  NLB=length_Lb0

  !NLB=1
  LB_size: do iLb =beginstep,NLB,interval !100

  !iLb=20;
  if(IsCaseLattice.ne.0) then
      LB0=LB0_list(length_Lb0-iLb+1)
  endif
  
  CALL INITIAL_SIZE(T)
  ALLOCATE(U_pos(6*NN),STAT=status) 
  ALLOCATE(SijN(5*NN),STAT=status)
  U_pos=0.0_8
  CALL INITIAL_CONF()

  if(K_rb.ne.0)then
    ALLOCATE (conf_rb(3,K_rb) ,STAT=status)
    ALLOCATE (conf_rb_vector(3,Np) ,STAT=status)
    ALLOCATE (U_par_rb(6*K_rb) ,STAT=status)
    ALLOCATE (q_rb(K_rb) ,STAT=status)
    ALLOCATE (rbmconn_Inertial_body(3*K_rb,3*K_rb) ,STAT=status)
    ALLOCATE (rbmconn_Inertial_body_inverse(3*K_rb,3*K_rb) ,STAT=status)
    CALL rbmconn_Init(NN,CONF,RADII,U_pos,conf_rb,U_par_rb,q_rb)
  endif
  call VMD_WRITE_BDY(CONF,T,RADII)
   
    K_time=0

    ALLOCATE(grmobmxsave(6*NN,6*NN),STAT=status) 
    Hydrosave=.false.
    Brownsave=.false.

     !====================== DEM ========================!  
     if(useDEM) then 
       !conf0=conf(:,1:Np);
       call DEM_calc_initial_position(CONF_DEM,LB_DEM,LB0-2.0_8,RADII*1.0_8,Np,timestep)

        do I=1,Np
        conf(1,I)=CONF_DEM(1,I)
        conf(2,I)=CONF_DEM(2,I)
        conf(3,I)=CONF_DEM(3,I)!+LB_DEM(3)*2.0_8+1.0_8
        enddo
       print*,'LB_DEM=',LB_DEM
       print*,'============= DEM ============'
#ifdef blockdem0
       block
       integer::pid,ia,ib
       real(8)::Rnew(3),DMIN
       do pid=1,Np
         print*,pid,CONF_DEM(1:3,pid)
       enddo
       do ia=1,Np-1
         do ib=ia+1,Np
           Rnew=conf_DEM(:,ib)-conf_DEM(:,ia)
           call PER_SKEW_CORRECTION(Rnew,DMIN)
           print*, ia,ib,Dmin
         enddo
       enddo
       end block
#endif
       print*,'============= DEM ============'
     endif
     !====================== DEM ========================!

    if(Oscillation_shear)then
      Call Init_frequency(LB,radii,frequency)
    else
      frequency=1.d0
    endif

    GAMMA=GAMMA_alpha*frequency
    WRITE(*,*) 'N: ',NN,' dt: ',BASE_DT,' shear: ',GAMMA,  &
    ' X: ',LB(1),' Y: ',LB(2), &
    ' Z: ',LB(3)
    write(*,*) 'LB0:',LB0,'Re_number=',Re_number, &
               'pho_f=',pho_f,'frequency=',frequency,'gamma_alpha=',GAMMA_alpha

    WRITE(50,*) 'BASE_DT,GAMMA,LB(1),LB0,Re_number,pho_f, & 
              & frequency,GAMMA_alpha'

    WRITE(50,"(15ES24.15)") BASE_DT,GAMMA,LB(1),LB0,Re_number,pho_f, & 
              frequency,GAMMA_alpha
    close(50)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!important!!!!!
    

    Time: DO



      if(Oscillation_shear) then
        GAMMA=GAMMA_alpha*frequency*cos(frequency*T)
        GAMMAangle=GAMMA_alpha*sin(frequency*T)
      endif

      !if(Xshear)then   !!!!error
      !u_bg(1)=-GAMMA*0.5_8*LB(3)
      !endif

      IF (mod(K_time,K_WRITE) .eq. 0) then
        CALL VMD_WRITE(CONF,T/lambda_time,RADII,U_pos)
        if(Nfloc_sum.eq.NN) then
          write(*,*)'Floc all complete'
          stop
        endif
      ENDIF

      CALL INIT_u_bg()
      !if(BROWN) then
        CALL STEPPER_Stokesian(CONF,RADII,DT,T,yeta_mu,U_pos,SijN)
      !else

      write(*,*) "rigid5==================================="

      IF (mod(K_time,K_WRITE) .eq. 0) then
        write(20,"(12ES24.15)") T/lambda_time,yeta_mu(:),LB(1),frequency,GAMMA,GAMMAangle
      ENDIF



      K_time=K_time+1
      T=T+DT

       IF (T.GT.END_TIME) THEN
        EXIT
       ENDIF

       !do I=1,NN-Nb
        !  do J=1,3
         ! if(abs(CONF(J,I)).gt.2E2) then
          !  write(*,*) 'CONF error'
           ! stop
         ! endif
         ! enddo
       ! enddo

      write(*,*)'T=========================',T

      !write(*,*)'LR=',LR
      !write(*,*)'LI=',LI
    ENDDO Time

  if(allocated(grmobmxsave))  DEALLOCATE( grmobmxsave)
  if(allocated(CONF))  DEALLOCATE( CONF)
  if(allocated(CONFBDY))  DEALLOCATE( CONFBDY)
  if(allocated(CONF_DEM)) DEALLOCATE( CONF_DEM)
  if(allocated(RADII)) DEALLOCATE( RADII)
  if(allocated(RADIIBDY)) DEALLOCATE( RADIIBDY)
  if(allocated(U_pos)) DEALLOCATE( U_pos)
  if(allocated(Floc_index)) DEALLOCATE(Floc_index) 
  if(allocated(KK_rbmconn)) DEALLOCATE(KK_rbmconn)
  if(allocated(conf_rb_vector))  DEALLOCATE( conf_rb_vector)
  if(allocated(U_par_rb)) DEALLOCATE(U_par_rb)
  if(allocated(q_rb)) DEALLOCATE( q_rb)  
  if(allocated(rbmconn_Inertial_body)) DEALLOCATE( rbmconn_Inertial_body)
  if(allocated(rbmconn_Inertial_body_inverse)) DEALLOCATE( rbmconn_Inertial_body_inverse)
  if(allocated(conf_rb))  DEALLOCATE( conf_rb)

  enddo LB_size
 
  CALL END_VMD_WRITE()

  if(allocated(LB0_list)) DEALLOCATE(LB0_list)
  call toc()
  ompend=omp_get_wtime()
  print*,' >>> omp Total time (s)',ompend-ompstart

  END