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
  use STEPPER_NOBROWN_MOD,only:STEPPER_NOBROWN
  use hydro_tools,only:Init_frequency,INIT_u_bg
  USE OMP_LIB


  implicit none
  INTEGER K_time,iLb,NLB,thread_id,status
  REAL*8 T,DT,yeta_mu(5),GAMMAangle,ompstart,ompend
  REAL*8,ALLOCATABLE :: U_pos(:),conf0(:,:)

  write(*,*) 'start-------------'
  call tic()
  call time_print() 
  ompstart= omp_get_wtime()   
  CALL INIT_KDirac
  CALL INIT_EPS
  CALL INIT_Y2()
  !call INIT_Y20()
  CALL INIT_II
  call INIT_READ_CONTROL


  thread_id = OMP_GET_MAX_THREADS()
  PRINT *, "max processos: ", thread_id

  call OMP_SET_NUM_THREADS(NUM_THREADS)
  PRINT *, "we use processos: ", NUM_THREADS
  DT = BASE_DT

  CALL INIT_VMD_WRITE()
  call INITIAL_LB0
  CALL INITIAL_SIZE(T)
  ALLOCATE( CONF0(3,NN),STAT=status )
  ALLOCATE(U_pos(3*NN),STAT=status)  

  yeta_mu=0.0_8
  NLB=length_Lb0

  !NLB=1
  LB_size: do iLb =beginstep,NLB,interval !100


    if(IsCaseLattice.ne.0) then
      LB0=LB0_list(length_Lb0-iLb+1)
    endif

    CALL INITIAL_CONF(T) 
    !conf0=conf
    K_time=0


     !====================== DEM ========================!  
     if(useDEM) then 
       !conf=conf0    
       call DEM_calc_initial_position(CONF,LB,LB0-2.0_8,RADII*1.0_8,NN,timestep)
       print*,'LB=',LB
       print*,'============= DEM ============'
#ifdef blockdem0
       block
       integer::pid,ia,ib
       real(8)::Rnew(3),DMIN
       do pid=1,NN
         print*,pid,CONF(1:3,pid)
       enddo
       do ia=1,NN-1
         do ib=ia+1,NN
           Rnew=conf(:,ib)-conf(:,ia)
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

    WRITE(50,"(15ES24.15)") BASE_DT,GAMMA,LB(1),LB0,Re_number,pho_f, & 
              frequency,GAMMA_alpha
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!important!!!!!
    U_pos=0.0_8

    Time: DO



      if(Oscillation_shear) then
        GAMMA=GAMMA_alpha*frequency*cos(frequency*T)
        GAMMAangle=GAMMA_alpha*sin(frequency*T)
      endif

      !if(Xshear)then   !!!!error
      !u_bg(1)=-GAMMA*0.5_8*LB(3)
      !endif
      IF (mod(K_time,K_WRITE) .eq. 0) then
        CALL VMD_WRITE(CONF,T,RADII,U_pos)
      ENDIF

      CALL INIT_u_bg()
      CALL STEPPER_NOBROWN(CONF,RADII,DT,T,yeta_mu,U_pos)

      IF (mod(K_time,K_WRITE) .eq. 0) then
        write(20,"(12ES24.15)") T,yeta_mu(:),LB(1),frequency,GAMMA,GAMMAangle
      ENDIF

      K_time=K_time+1
      T=T+DT

       IF (T.GT.END_TIME) THEN
        EXIT
       ENDIF

          

      write(*,*)'T=========================',T

      !write(*,*)'LR=',LR
      !write(*,*)'LI=',LI
    ENDDO Time

  enddo LB_size



  if(allocated(CONF))  DEALLOCATE( CONF)
  if(allocated(CONF0)) DEALLOCATE( CONF0)
  if(allocated(RADII)) DEALLOCATE( RADII)
  if(allocated(U_pos)) DEALLOCATE( U_pos)
  if(allocated(LB0_list)) DEALLOCATE(LB0_list)
  CALL END_VMD_WRITE()

  call toc()
  ompend=omp_get_wtime()
  print*,' >>> omp Total time (s)',ompend-ompstart

  END