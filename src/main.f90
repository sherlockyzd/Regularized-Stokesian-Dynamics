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
  use filament
  use filament_math


  implicit none
  INTEGER K_time,iLb,NLB,thread_id,status,I,J
  REAL*8 T,DT,yeta_mu(5),GAMMAangle,ompstart,ompend
  !REAL*8 SijN(5*NN)
  REAL*8,ALLOCATABLE :: U_pos(:),SijN(:),U_pos_swimer(:),U_pos_filament(:)


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
  CALL INITIAL_CONF_ALLOCATION()

  if(K_rb.ne.0)then
    ALLOCATE(U_pos_swimer(6*Nswimer),STAT=status) 
    CALL rbmconn_Init(Nswimer,CONF(:,1:Nswimer),RADII(1:Nswimer),U_pos_swimer,conf_rb,U_par_rb,q_rb)
  endif

  if(F_rb.ne.0)then
    ALLOCATE(U_pos_filament(6*Nfilament),STAT=status) 
    CALL filament_Init(Nfilament,CONF(:,Nswimer+1:Nswimer+Nfilament), &
      & RADII(Nswimer+1:Nswimer+Nfilament),U_pos_filament)
  endif

  if(Nswimer.ne.0.and.Nfilament.ne.0) then
    call subToAll_particle(NN,Nswimer,Nfilament,U_pos_swimer,U_pos_filament,U_pos)
  elseif(Nswimer.ne.0.and.Nfilament.eq.0) then
    U_pos(1:6*Nswimer)=U_pos_swimer
  elseif(Nswimer.eq.0.and.Nfilament.ne.0) then
    U_pos(1:6*Nfilament)=U_pos_filament
  endif


  if(boundary) then
    call VMD_WRITE_BDY(CONF,T,RADII)
  endif
   write(*,*)'1---------------------check'
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

    WRITE(*,*) 'BASE_DT,GAMMA,LB(1),LB0,Re_number,pho_f, & 
              & frequency,GAMMA_alpha'

    WRITE(*,"(15ES24.15)") BASE_DT,GAMMA,LB(1),LB0,Re_number,pho_f, & 
              frequency,GAMMA_alpha
    !close(50)
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
      
      CALL STEPPER_Stokesian(CONF,RADII,DT,T,yeta_mu,U_pos,SijN)
    
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

      if(isnan(U_pos(2)).or.abs(U_pos(2)).gt.1.0_8) then 
       write(*,*) "U_pos has error",U_pos(1:3)
       call exit()
      endif
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
    if(allocated(U_pos_swimer)) DEALLOCATE( U_pos_swimer)
    if(allocated(U_pos_filament)) DEALLOCATE( U_pos_filament)
    if(allocated(Floc_index)) DEALLOCATE(Floc_index)

    CALL END_CONF_DEALLOCATION()

  enddo LB_size
 
  CALL END_VMD_WRITE()

  if(allocated(LB0_list)) DEALLOCATE(LB0_list)
  if(allocated(Filament_num)) DEALLOCATE(Filament_num)
  call toc()
  ompend=omp_get_wtime()
  print*,' >>> omp Total time (s)',ompend-ompstart

  END