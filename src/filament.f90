  Module filament_math
   use quaternions
   use method,only:FTS_method,simplePeriod
  ! USE CONFIG,only:u_bg,omega_bg,omegaT,Eij   ! CONF,POLY_LEN
  ! use hydro_tools,only:calc_uo_bg
   use period_bdy_tools,only:PERIOD_CORRECTION
   use SYS_property,only:mass_par,DT_DEM
   USE LATTICE_BASE,only:LB
   use filament!,only:F_rb,fila_rbmconn,conf_rb_vector,rbmconn_Inertial_body,rbmconn_Inertial_body_inverse

   IMPLICIT NONE
   private

   public:: filament_Init,InternalForcesAndTorques,new_U1_filament,new_conf_filament
   !,new_U_par_filament,new_conf_filament,new_q_filament

   contains

!****************************************************************


    subroutine filament_Init(NN,CONF,RADII,U_pos)
    IMPLICIT NONE
    INTEGER,intent(in)::NN
    REAL*8,intent(in):: CONF(3,NN),RADII(NN)
    real*8,intent(out):: U_pos(6*NN)

    integer :: ii, jj, KK_fila
    integer ::num_fila_sum,num_fila
    real*8 :: Ixx, Iyy,Izz,RADII2

      U_pos=0.0_8
      Filament_tau_now(:,:)=0.0_8
      num_fila_sum=0
      do  KK_fila=1,F_rb       
          num_fila=Filament_num(KK_fila)       
          do ii = 1, num_fila
              jj=ii+num_fila_sum
              if(ii.eq.1)then
                Filament_tau_now( :, jj) = conf(:,jj+1)-conf(:,jj)
              elseif(ii.eq.num_fila)then
                Filament_tau_now( :, jj) = conf(:,jj)-conf(:,jj-1)
              else
                Filament_tau_now( :, jj) = conf(:,jj+1)-conf(:,jj)
              endif
              Filament_q(jj)=quat(1.0_8,0.0_8,0.0_8,0.0_8)
          end do
          index1(KK_fila)=num_fila_sum+1
         ! Filament_X1(3*(KK_fila-1):3*KK_fila)=conf(:,num_fila_sum+1)
         ! Filament_U1(3*(KK_fila-1):3*KK_fila)=U_pos(3*num_fila_sum+1:3*num_fila_sum+3)
          num_fila_sum=num_fila_sum+num_fila
      end do
      !Filament_tau_confold=Filament_tau_conf
      Filament_tau_base=Filament_tau_now
      Filament_Inertial_body_inverse=0.0_8
      Filament_Inertial_body=0.0_8
      RADII2=0.4_8*RADII(1)*RADII(1)
      Ixx=RADII2
      Iyy=RADII2
      Izz=RADII2
      Filament_Inertial_body( 1, 1 ) =Ixx 
      Filament_Inertial_body( 2, 2 ) =Iyy
      Filament_Inertial_body( 3, 3 ) =Izz 
      Filament_Inertial_body_inverse=Filament_Inertial_body
      CALL MATREV(Filament_Inertial_body_inverse,3)

    end subroutine filament_Init


    subroutine InternalForcesAndTorques(NN,conf,Filament_internal_force_torque)
    ! Filament.INTERNALFORCESANDTORQUES()  Place gravity, elastic and 
    !                                      constraint forces and torques 
    !                                      on the filament. They are stored
    !                                      in Filament.F
    IMPLICIT NONE
    INTEGER,intent(in)::NN
    real*8,intent(in)::conf(3,NN)!,Filament_tau_base(3,NN)
    real*8,intent(out)::Filament_internal_force_torque(6*NN)
    integer :: ii, jj, KK_fila
    integer ::num_fila_sum,num_fila
    real*8 :: dr(3),ds,tauj(3),tauj1(3),M(3),F0(3),F1(3),tauF(3)
    type(quaternion)::q_mid,dqds,qb
    real*8 :: internal_force(3),internal_torque(3),angle(3)


      Filament_internal_force_torque=0.0_8
      num_fila_sum=0
      do  KK_fila=1,F_rb       
          num_fila=Filament_num(KK_fila)       
          do ii = 1, num_fila-1
              jj=ii+num_fila_sum
              if(ii.eq.num_fila)then
                Filament_internal_force_torque( 3*(jj-1)+1:3*jj) = 0.0_8
                Filament_internal_force_torque( 3*NN+3*(jj-1)+1:3*NN+3*jj) = 0.0_8
              else
                tauj=0.5_8*Filament_tau_now(:,jj)
                tauj1=0.5_8*Filament_tau_now(:,jj+1)
                dr=conf(:,jj+1)-conf(:,jj)
                !ds=SQRT(tauj1**2)+sqrt()
                ds=SQRT(SUM((Filament_tau_base(:,jj))**2))
                q_mid = qMidpoint(Filament_q(jj),Filament_q(jj+1))
                ! Elastic torques
                !dqds =dqsub(Filament_q(jj+1),Filament_q(jj))/ds
                !qb = 2.0_8*qmul(qconj(q_mid),dqds)
                angle=qangle(Filament_q(jj+1),Filament_q(jj))/ds
                !M(1)=GI*qb%i;M(2)=EI*qb%j;M(3)=EI*qb%k
                M(1)=GI*angle(1);M(2)=EI*angle(2);M(3)=EI*angle(3)               
                internal_torque = qrot(q_mid,M)


                !internal Force
                
                tauF=dr/ds!(tauj1+tauj)/dL
                F0=qrot_inverse(q_mid,tauF)
                F0=F0-Filament_tau_base(:,jj)/ds
                F1=0.0_8
                F1(1)=F0(1)*EA;
                !F1(2)=F0(2)*GA;F1(3)=F0(3)*GA
                internal_force=qrot(q_mid,F1)
                write(*,*) 'MMMMMMMMMMMMMMM===========',M(:)!,GI,EI
                write(*,*) 'internal_torque===========',internal_torque(:)
                write(*,*) 'FFFFFFFFFFFFFFF===========',F1(:)
                write(*,*) 'internal_Force ===========',internal_force(:)

                Filament_internal_force_torque( 3*(jj-1)+1:3*jj) = Filament_internal_force_torque( 3*(jj-1) &
                  & +1:3*jj) +internal_force(3)
                Filament_internal_force_torque( 3*jj+1:3*jj+3) = Filament_internal_force_torque( 3*jj+1:3*jj+3)&
                & -internal_force(3)
                Filament_internal_force_torque( 3*NN+3*(jj-1)+1:3*NN+3*jj) = Filament_internal_force_torque& 
                  & ( 3*NN+3*(jj-1)+1:3*NN+3*jj) +internal_torque(:)+CrossProduct3D(0.5*dr,internal_force)
                Filament_internal_force_torque( 3*NN+3*jj+1:3*NN+3*jj+3) = Filament_internal_force_torque &
                  & ( 3*NN+3*jj+1:3*NN+3*jj+3) -internal_torque(:)+CrossProduct3D(0.5*dr,internal_force)
              endif
          end do
          num_fila_sum=num_fila_sum+num_fila
      end do
    do ii=1,NN
        write(*,*) 'i,internal_filament_torque===',ii,Filament_internal_force_torque(3*NN+3*(ii-1)+1:3*NN+3*ii)
            !write(*,*) 'i, torue===',ii,Ftotal(3*NN+3*(ii-1)+1:3*NN+3*ii)
    enddo  


      end subroutine InternalForcesAndTorques



    subroutine new_U1_filament(NN,RADII,Ftotal,U_pos)
    IMPLICIT NONE
    INTEGER,intent(in)::NN
    REAL*8,intent(in):: RADII(NN),Ftotal(6*NN)
    real*8,intent(inout):: U_pos(6*NN)

    integer :: ii,KK_fila

    real*8 :: filament_torque(3*NN)


    do ii=1,NN
      filament_torque(3*(ii-1)+1:3*ii)=matmul(Filament_Inertial_body_inverse,Ftotal(3*NN+3*(ii-1)+1:3*NN+3*ii))
    enddo

    do ii=1,NN
        write(*,*) 'i,filament_torque/Internal_moment===',ii,filament_torque(3*(ii-1)+1:3*ii)
            !write(*,*) 'i, torue===',ii,Ftotal(3*NN+3*(ii-1)+1:3*NN+3*ii)
    enddo
   !CALL MATREV(rbmconn_Inertial,3*K_rb)
    U_pos(3*NN+1:6*NN)=U_pos(3*NN+1:6*NN)+filament_torque*DT_DEM/mass_par
    U_pos(1:3*NN)=U_pos(1:3*NN)+Ftotal(1:3*NN)*DT_DEM/mass_par


    !do  KK_fila=1,F_rb
    !do  ii=1,NN
        !Filament_U1(3*(KK_fila-1)+1:3*KK_fila)=Filament_U1(3*(KK_fila-1)+1:3*KK_fila)+ &
        !    & Ftotal(3*(index1(KK_fila)-1):3*index1(KK_fila))*DT_DEM/mass_par  
        !U_pos(3*(ii-1)+1:3*ii)=Filament_U1(3*(KK_fila-1)+1:3*KK_fila)
    !end do            

    end subroutine new_U1_filament



    subroutine new_conf_filament(NN,conf,T,U_pos)
    IMPLICIT NONE
    INTEGER,intent(in)::NN
    REAL*8,intent(in):: T,U_pos(6*NN)!,U_pos(6*NN)!,RADII(NN)
    !real*8,intent(in):: U_par_rb(6*K_rb)
    real*8,intent(inout):: conf(3,NN)
   ! type(quaternion),intent(inout)::q_rb(K_rb)

    integer :: KK_fila,ii, jj!, jj2, kk, ll, mm,num_rb
    integer ::num_fila_sum,num_fila
    real*8 :: tauj(3),tauj1(3),taujold(3),tauj1old(3)
    !integer ::n6,n5!, idostep
    !real*8, dimension( 3 ) :: swimcm
    real*8 ::DCONF1(3),w(3)
    type(quaternion) ::dq,q_norm

    !Filament_tau_confold=Filament_tau_now

    ! write(*,*) "q==000000000000000"
     ! do ii=1,NN
     !   CALL qprint(Filament_q(ii))
      !  Filament_tau_conf(:,ii)=qrot(Filament_q(ii),Filament_tau_base(:,ii))
      !enddo

    do ii=1,NN
      jj=3*NN+3*(ii-1)
      w=U_pos(jj+1:jj+3)
      dq=qscal(DT_DEM,qderivat(Filament_q(ii),w))
      q_norm=qadd(Filament_q(ii),dq)
      Filament_q(ii)=qscal(1.0_8/qnorm(q_norm) , q_norm)
      !write(*,*) "ii==",ii
      !CALL qprint(Filament_q(ii))
      Filament_tau_now(:,ii)=qrot(Filament_q(ii),Filament_tau_base(:,ii))
    enddo

      forall(ii=1:NN,jj=1:3)
        CONF(jj,ii)=CONF(jj,ii)+DT_DEM*U_pos(3*(ii-1)+jj)
      end forall


#ifdef U1
    num_fila_sum=0
    do KK_fila=1,F_rb
      num_fila=Filament_num(KK_fila)
      conf(:,num_fila_sum+1)=Filament_X1(3*(KK_fila-1)+1:3*KK_fila)
      U_pos(3*num_fila_sum+1:3*num_fila_sum+3)=Filament_U1(3*(KK_fila-1)+1:3*KK_fila)
      do ii = 2, num_fila
          jj=ii+num_fila_sum
          tauj=0.5_8*Filament_tau_conf(:,jj-1)
          tauj1=0.5_8*Filament_tau_conf(:,jj)
          conf(:,jj)=conf(:,jj-1)+tauj+tauj1
          taujold=0.5_8*Filament_tau_confold(:,jj-1)
          tauj1old=0.5_8*Filament_tau_confold(:,jj)
          U_pos(3*(jj-1)+1:3*jj)=U_pos(3*(jj-2)+1:3*(jj-1))+(tauj1+tauj-tauj1old-taujold)/DT_DEM
      end do
      num_fila_sum=num_fila_sum+num_fila
    enddo
#endif

      if(simplePeriod) then
         DO ii = 1, NN
          CALL PERIOD_CORRECTION(conf(:,ii),LB,T)
         ENDDO
      endif


    end subroutine new_conf_filament

    end  Module filament_math


  Module filament_solve_Jacobian
  use quaternions
  use method,only:FTS_method,simplePeriod
  ! USE CONFIG,only:u_bg,omega_bg,omegaT,Eij   ! CONF,POLY_LEN
  ! use hydro_tools,only:calc_uo_bg
  use period_bdy_tools,only:PERIOD_CORRECTION
  use SYS_property,only:mass_par,DT_DEM
  USE LATTICE_BASE,only:LB
  use filament!,only:F_rb,fila_rbmconn,conf_rb_vector,rbmconn_Inertial_body,rbmconn_Inertial_body_inverse

  IMPLICIT NONE
  private

 ! public:: 
 !,new_U_par_filament,new_conf_filament,new_q_filament

  contains
  !filament_function(X1_NP1,lamda,cita,X1_NM1,X1_N,torque,Fe,Te,)






  end Module filament_solve_Jacobian
