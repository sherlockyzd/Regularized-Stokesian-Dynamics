  module STEPPER_MOD
  USE SIZE           ! NN
  USE CONFIG,only:Eij,omegaT,u_bg,omega_bg
  use method
  USE TENSORS,ONLY:Y21,PAI,front
  USE LATTICE_BASE
  use SYS_property
  use gmres_mod_gmres
  use R1d_mod, only: norm_2
  use hydro_tools
  use period_bdy_tools,only:PERIOD_CORRECTION,SET_LATTICE_SKEW
  use hydro_lub_mod
  use Wall_interaction
  use BROWN_proerty
  use Mobilitymatrix ,only:Brownsave
  use conglomerate
  use rb_conglomerate,only:K_rb,conf_rb,U_par_rb,q_rb
  use filament
  use filament_explicit_or_implicit1
  use filament_implicit2
  use filament_solve_Explicit
  use filament_solve_Implicit_method1
  use filament_solve_Implicit_method2
  implicit none
  private


    public::STEPPER_Stokesian

  contains

      SUBROUTINE STEPPER_Stokesian(CONF,RADII,DT,T,yeta_mu,U_pos,SijN)    ! RETURNS NEW T=T+DT
      IMPLICIT NONE
      REAL*8,intent(inout):: CONF(3,NN),U_pos(6*NN)
      REAL*8,intent(in):: RADII(NN),DT,T
      real(8),intent(out):: yeta_mu(5),SijN(5*NN)

      REAL*8 EIN(5*NN),TPP(3*NN,3*NN)
      REAL*8 APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8 RPP(6*NN,6*NN),SijN_hyd(5*NN),SijN_FP(5*NN),SijN_B(5*NN)!,,SijN(5*NN)
      REAL*8 LATTICE(3,3),EWS_ERR,Fe_inter(6*NN),Fe_body(6*NN),Fe(6*NN)
      REAL*8 uo_bg(6*NN),U_relative(6*NN),u_brown(3*NN)
      real*8 rfu(6*NN,6*NN),rfe(6*NN,5*NN),rse(5*NN,5*NN)
      real*8 A(3*NN,3*NN),Z(3*NN),VB(3*NN),DMDX(3*NN),V_DM(3*NN)
      logical pos_collision
      INTEGER I,J,Nt

      !character(128)::myformat

      call EIN_transY2(EIN,Eij,NN)
      call pos_collision_judge(conf,radii,pos_collision)
      call ppiclf_collision_timestep(RADII,DT,pos_collision)
      if(simplePeriod) then
        CALL SET_LATTICE_SKEW(T,LATTICE)
        !write(*,*) T,LATTICE
      endif      
      if(.not.IsPeriod) then
       CALL GRPERY_INV_FRI(APP,APQ,AQQ,CONF,RADII,NN)
       write(*,*) 'notperiod'
      else
       EWS_ERR = LINERR
       CALL PER_GRPERY_INV_FRI(APP,APQ,AQQ,CONF,RADII,NN,LATTICE,EWS_ERR)
       write(*,*) 'period'
      endif
      
      if(Wall_method) then
        call Wall_GRPERY_INV_FRI(TPP,CONF,RADII,NN)
        APP(1:3*NN,1:3*NN)=APP(1:3*NN,1:3*NN)+TPP
      endif

      if(BROWN) then      
        A= APP(1:3*NN,1:3*NN)
        CALL CHOLESKY(A,3*NN)
        Z=0.0_8
        DO I=1,3*NN
        CALL GAUSS(Z(I))
        ENDDO
        VB=sqrt(2.0_8*kBT)*MATMUL(A,Z)/DT 
        !V_DM=0.0_8
        call DMDX_CAL(DMDX,CONF,RADII,Z,T)
        V_DM=kBT*DMDX
        u_brown=VB(1:3*NN)+V_DM(1:3*NN)
      else
        u_brown=0.0_8
      endif


        write(*,*) 'u_brown==',BROWN
       ! do I=1,NN
        !  write(*,*) 'u_brown',u_brown(3*(I-1)+1),u_brown(3*(I-1)+2),u_brown(3*(I-1)+3)
        !enddo


      NT_DEMdo:do Nt=1,NT_DEM
        call source_inter_F(conf,RADII,U_pos,Fe_inter)
        call source_body_F(conf,RADII,Fe_body)
        !do I=1,NN
        !  write(*,*) Fe_body(3*(I-1)+1:3*I)
        !enddo
        Fe=Fe_body + Fe_inter
        rfu=0.0_8
        rfe=0.0_8
        rse=0.0_8
        if(uselub)then
         call lubmxcalc(NN,conf,rfu,rfe,rse)
        endif

        do i=1,NN-Nb
            write(*,*) 'i,FeFtotal===',i,Fe(3*(i-1)+1:3*i)
            write(*,*) 'i,Fetorue====',i,Fe(3*NN+3*(i-1)+1:3*NN+3*i)
        enddo

        !if(T.gt.0.0_8.and.useDEMstepper) then
        if(useDEMstepper) then
            call INIT_uo_bg(conf,uo_bg)
            U_relative=0.0_8 !angular velocity
            U_relative=U_pos-uo_bg!-u_brown
            
            call STEPPER_UDEM(radii,APP,APQ,AQQ,rfu,rfe,rse,Fe,EIN,RPP,U_relative,U_pos,u_brown,SijN_hyd &
               & ,Fe_inter,SijN_FP,CONF,SijN_B)
            write(*,*) "Using STEPPER_UDEM success"
             
        else
          if(useGMRESmethod) then
            call STEPPER_Gmres(APP,APQ,AQQ,rfu,rfe,rse,Fe,EIN,U_relative,SijN_hyd)
            write(*,*) "Using STEPPER_Gmres"
          else
            call STEPPER_NOGmres(APP,APQ,AQQ,rfu,rfe,rse,Fe,EIN,RPP,U_relative,u_brown,SijN_hyd &
               & ,Fe_inter,SijN_FP,CONF,SijN_B)
            write(*,*) "Using STEPPER_NOGmres"
          endif
        endif

          
        CALL STEPPER_conf(conf,radii,U_relative,U_pos,T)
         


      enddo NT_DEMdo

      
      !SijN=SijN_FP+SijN_hyd
      !SijN_FP=0.0_8
      call yetamu_solve(APP,RPP,Fe,U_relative,RADII,SijN_hyd,SijN_FP,SijN_B,EIN,yeta_mu)

      write(*,*) 'check________STEPPER_Stokesian__________Success!'
      END SUBROUTINE STEPPER_Stokesian

!***********************************************************

      SUBROUTINE STEPPER_conf(conf,radii,U_par,U_pos,T)
      REAL*8,intent(inout):: CONF(3,NN),U_pos(6*NN)
      REAL*8,intent(in)::U_par(6*NN),radii(NN),T
      !REAL*8,intent(inout):: 
   
      REAL*8 DCONF_1(3*NN),uo_bg(6*NN),U_pos_filament(6*Nfilament)
      INTEGER I,J


      if(K_rb.ne.0) then
        call new_conf_swim(Nswimer,conf(:,1:Nswimer),T,U_par_rb,conf_rb,q_rb)
        write(*,*) 'new_conf_rb_rb_rb____okok'
      endif


      if(F_rb.ne.0) then
        if(.not.filament_solve_implicit)then
          U_pos_filament(1:3*Nfilament)=U_pos(3*Nswimer+1:3*Nswimer+3*Nfilament)
          U_pos_filament(3*Nfilament+1:6*Nfilament)=U_pos(3*NN+3*Nswimer+1:3*NN+3*Nswimer+3*Nfilament)
          call new_conf_explicit_filament(Nfilament,conf(:,Nswimer+1:Nswimer+Nfilament),T,U_pos_filament)
        else
          if(filament_implicit_method.ne.2)then
              CONF=Filament_conf_now
              if(simplePeriod) then
                 DO i = 1, Nfilament
                  CALL PERIOD_CORRECTION(conf(:,i),LB,T)
                 ENDDO
              endif
          else
            call new_conf_Implicit2_filament(Nfilament,conf,T)
          endif
        endif

        write(*,*) 'new_conf_filament_____filament_____filament_____okok'
      endif

      if(K_rb.eq.0.and.F_rb.eq.0) then
        call INIT_uo_bg(conf,uo_bg)
        U_pos=U_par+uo_bg
        !call par_index_floc(CONF,radii)
        !call U_BDY_CORR(U_pos)
        if(.not.fix) then
          DCONF_1=DT_DEM*U_pos(1:3*NN)

          forall(i=1:NN-Nb,j=1:3)
            CONF(j,i)=CONF(j,i)+DCONF_1(3*(i-1)+j)
          end forall
          if(simplePeriod) then
             DO I = 1, Np
              CALL PERIOD_CORRECTION(CONF(:,I),LB,T)
             ENDDO
             call INIT_uo_bg(conf,uo_bg)
             !if(K_rb.ne.0)then
             !   call  swim_uo_bg(NN,CONF,uo_bg)
             !endif
             U_pos=U_par+uo_bg
          endif
        endif
       ! call par_index_floc(CONF,radii)
       ! call U_BDY_CORR(U_pos)
       write(*,*) 'check________STEPPER_conf__________Success!'
      endif

      END SUBROUTINE STEPPER_conf

      SUBROUTINE STEPPER_UDEM(radii,APP,APQ,AQQ,rfu,rfe,rse,Fe,EIN,RPP,U_par,U_pos,u_brown,SijN_hyd,Fe_inter,SijN_FP,CONF,SijN_B)
      IMPLICIT NONE
      !REAL*8,intent(in):: CONF(3,NN)
      REAL*8,intent(in):: APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8,intent(in):: rfu(6*NN,6*NN),rfe(6*NN,5*NN),rse(5*NN,5*NN)
      real*8,intent(in):: EIN(5*NN),Fe(6*NN),radii(NN),u_brown(3*NN),Fe_inter(6*NN),CONF(3,NN)
      real(8),intent(inout):: U_par(6*NN),U_pos(6*NN)
      real(8),intent(out):: SijN_hyd(5*NN),SijN_FP(5*NN),SijN_B(5*NN),RPP(6*NN,6*NN)

      REAL*8 RPQ(6*NN,5*NN),RQQ(5*NN,5*NN),UB(6*NN),Fh(6*NN)
      REAL*8 RPP0(6*NN,6*NN),RPQ0(6*NN,5*NN),RQQ0(5*NN,5*NN)
      real*8 mu(5*NN,5*NN),APPR(6*NN,6*NN),U_add(6*NN)
      REAL*8 DU_1(6*NN),uo_bg(6*NN),Ftotal(6*NN)
      REAL*8 Ftotal_filament(6*Nfilament),U_pos_filament(6*Nfilament)
      REAL*8 U_pos_swimer(6*Nswimer),Ftotal_swimer(6*Nswimer),Filament_Interal_force(6*Nfilament)

      integer I,J

        RPQ=0.0_8
        RQQ=0.0_8
        UB=0.0_8
        do I=1,3*NN
          UB(I)=u_brown(I)
        enddo
        if(FTS_method)then
         call MobToResist_FTS(APP,APQ,AQQ,NN,RPP0,RPQ0,RQQ0)
        else
         call MobToResist_FT(APP,NN,RPP0)
        endif
        call arrangResist(NN,rfu,rfe,rse,RPP0,RPQ0,RQQ0,RPP,RPQ,RQQ)

        U_add=U_par
        U_add(1:3*NN)=U_add(1:3*NN)-u_brown
        if(FTS_method)then
         Fh = matmul(RPQ,EIN)-matmul(RPP,U_add)
        else
         Fh = -matmul(RPP,U_add)
        endif
        do i=1,NN-Nb
            write(*,*) 'i,Fhforce===',i,Fh(3*(i-1)+1:3*i)
            write(*,*) 'i,Fhtorue===',i,Fh(3*NN+3*(i-1)+1:3*NN+3*i)
        enddo
        Ftotal = Fe+Fh

        if(Nswimer.ne.0.and.Nfilament.ne.0) then
          call AllTosub_particle(NN,Nswimer,Nfilament,U_pos_swimer,U_pos_filament,U_pos)
          call AllTosub_particle(NN,Nswimer,Nfilament,Ftotal_swimer,Ftotal_filament,Ftotal)
        elseif(Nswimer.ne.0.and.Nfilament.eq.0) then
          U_pos_swimer=U_pos(1:6*Nswimer)
          Ftotal_swimer=Ftotal(1:6*Nswimer)
        elseif(Nswimer.eq.0.and.Nfilament.ne.0) then
          U_pos_filament=U_pos(1:6*Nfilament)
          Ftotal_filament=Ftotal(1:6*Nfilament)
        endif
        
        if(K_rb.ne.0)then
           call new_U_par_rb(Nswimer,CONF(:,1:Nswimer),RADII(1:Nswimer),Ftotal_swimer,U_pos_swimer,U_par_rb,q_rb)
           write(*,*) 'U_par_swimmer_____swimmer_____swimmer_____swimmer_____okok'
        endif
           
        if(F_rb.ne.0) then
            if(.not.filament_solve_implicit) then
              write(*,*) 'filament_____using EXPLICIT'
              call filament_explicit_solve(Nfilament,RADII(Nswimer+1:Nswimer+Nfilament),Ftotal_filament,U_pos_filament)
            else
              if(filament_implicit_method.ne.2) then
                write(*,*) 'filament_____using IMPLICIT filament_solve_Implicit_method1'
                call filament_implicit_solve(Nfilament,Filament_X1_now,Filament_X1_past,Filament_U1_now, &
                      & Filament_tau_now,Filament_Lie_algebra_now,Filament_Interal_force,Ftotal_filament, &
                      & Filament_conf_now,Filament_conf_past,U_pos_filament,Filament_q)
              else
                write(*,*) 'filament_____using IMPLICIT filament_solve_Implicit_method2'
                call filament_implicit_solve2(Nfilament,Filament_Interal_force,Ftotal_filament,U_pos_filament)
              endif
            endif
            write(*,*) 'U_par_filament_____filament_____filament_____okok'
        endif

        if(K_rb.eq.0.and.F_rb.eq.0) then
            U_par=U_par+Ftotal*DT_DEM/mass_par
        endif

        if(Nswimer.ne.0.and.Nfilament.ne.0) then
          call subToAll_particle(NN,Nswimer,Nfilament,U_pos_swimer,U_pos_filament,U_pos)
        elseif(Nswimer.ne.0.and.Nfilament.eq.0) then
          U_pos(1:6*Nswimer)=U_pos_swimer
        elseif(Nswimer.eq.0.and.Nfilament.ne.0) then
          U_pos(1:6*Nfilament)=U_pos_filament
        endif
         

        call arrangement(RPP,RPQ,RQQ,mu,NN)
        
        SijN_hyd=matmul(mu,EIN)
        
        call Rheology(CONF,RPQ,RPP,Fe_inter,SijN_FP)

        SijN_B=matmul(transpose(RPQ),UB)
       !do i=1,NN-Nb
         ! mass = pho_par*4.0_8/3.0_8*pai*RADII(i)**3
         ! do j=1,3 
            !DU_1=Ftotal*DT_DEM/mass_par
           ! DU_1(3*NN+3*(i-1)+j)=Ftotal(3*NN+3*(i-1)+j)*DT_DEM/mass
          !enddo
        !enddo
        !do i=1,NN-Nb
         ! write(*,*) "U_pos====",i,U_pos(3*(i-1)+1:3*i)
        !enddo
        write(*,*) 'check________STEPPER_UDEM__________Success!'
      end SUBROUTINE STEPPER_UDEM

      SUBROUTINE STEPPER_NOGmres(APP,APQ,AQQ,rfu,rfe,rse,Fe,EIN,RPP,U_par,u_brown,SijN_hyd,Fe_inter,SijN_FP,CONF,SijN_B)
      IMPLICIT NONE
      REAL*8,intent(in):: APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8,intent(in):: rfu(6*NN,6*NN),rfe(6*NN,5*NN),rse(5*NN,5*NN)
      real*8,intent(in):: EIN(5*NN),Fe(6*NN),Fe_inter(6*NN),CONF(3,NN),u_brown(3*NN)
      real(8),intent(out):: SijN_hyd(5*NN),SijN_FP(5*NN),SijN_B(5*NN),RPP(6*NN,6*NN),U_par(6*NN)

      REAL*8 RPQ(6*NN,5*NN),RQQ(5*NN,5*NN),UB(6*NN)
      REAL*8 RPP0(6*NN,6*NN),RPQ0(6*NN,5*NN),RQQ0(5*NN,5*NN)
      real*8 mu(5*NN,5*NN),APPR(6*NN,6*NN),uo_bg_rb(6*K_rb)!,rbmconn(6*NN, 6*K_rb),rfu_swim( 6*K_rb, 6*K_rb )
      INTEGER I,J

        RPQ=0.0_8
        RQQ=0.0_8
        UB=0.0_8
        do I=1,3*NN
          UB(I)=u_brown(I)
        enddo
        if(FTS_method)then
         call MobToResist_FTS(APP,APQ,AQQ,NN,RPP0,RPQ0,RQQ0)
        else
         call MobToResist_FT(APP,NN,RPP0)
        endif

        call arrangResist(NN,rfu,rfe,rse,RPP0,RPQ0,RQQ0,RPP,RPQ,RQQ)

        if(K_rb.ne.0)then
           call new_vel(NN,CONF,RPP,RPQ,Fe,EIN,U_par,U_par_rb,q_rb)
           call calc_uo_bg(K_rb,conf_rb,uo_bg_rb)
           U_par_rb=U_par_rb+uo_bg_rb
           !write(*,*) 'swimmer_____swimmer_____swimmer_____swimmer_____okok'
           !write(*,*) 'K_rb===',K_rb  
        else

          APPR=RPP  
          !write(*,*) '1111111111111111111okok'
          CALL MATREV(APPR,6*NN)
          !write(*,*) '222222222222222222222222okok'
          if(FTS_method)then
           U_par = matmul(APPR,Fe+matmul(RPQ,EIN))
          else
           U_par = matmul(APPR,Fe)
          endif

        endif

        call arrangement(RPP,RPQ,RQQ,mu,NN)
        SijN_hyd=matmul(mu,EIN)
        
        call Rheology(CONF,RPQ,RPP,Fe_inter,SijN_FP)

        SijN_B=matmul(transpose(RPQ),UB)

      end SUBROUTINE STEPPER_NOGmres

      SUBROUTINE Rheology(CONF,RPQ,RPP,Fe_inter,SijN_FP)
      IMPLICIT NONE
      REAL*8,intent(in):: CONF(3,NN),RPP(6*NN,6*NN),RPQ(6*NN,5*NN),Fe_inter(6*NN)
      real(8),intent(out):: SijN_FP(5*NN)

      real*8 Sij(3,3),Sij_FP_part(5)
      INTEGER I,J,KK
      !real*8 RPPINV(6*NN,6*NN)!,RQQTemp(5*NN,5*NN),muTemp(5,5)

      CALL MATREV(RPP,6*NN)

      SijN_FP=matmul(matmul(transpose(RPQ),RPP),Fe_inter)

      Do KK=1,NN 
         Sij=0.0_8
         DO I=1,3
           Do J=1,3
            Sij(I,J)=0.5_8*(CONF(I,KK)*Fe_inter(3*(KK-1)+J)+CONF(J,KK)*Fe_inter(3*(KK-1)+I))
           enddo
         enddo

        ! Sij(1,2)=-0.25*(DF1*rydiff+DF2*rxdiff)
        ! Sij(1,3)=-0.25*(DF1*rzdiff+DF3*rxdiff)
         !Sij(2,3)=-0.25*(DF2*rzdiff+DF3*rydiff)
        ! Sij(2,1)=Sij(1,2)
        ! Sij(3,1)=Sij(1,3)
         !Sij(3,2)=Sij(2,3)
         call EI_transY2(Sij_FP_part,Sij)

          
          SijN_FP(5*(KK-1)+1:5*KK) = SijN_FP(5*(KK-1)+1:5*KK)+Sij_FP_part(1:5)
          !SijN_FP(5*(j-1)+1:5*j) = SijN_FP(5*(j-1)+1:5*j)+Sij_FP_part(1:5)
      ENDDO


      END SUBROUTINE Rheology


      SUBROUTINE STEPPER_Gmres(APP,APQ,AQQ,rfu,rfe,rse,Fe,EIN,U_par,SijN)
      IMPLICIT NONE
      REAL*8,intent(in):: APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8,intent(in):: rfu(6*NN,6*NN),rfe(6*NN,5*NN),rse(5*NN,5*NN)
      real*8,intent(in):: EIN(5*NN),Fe(6*NN)
      real(8),intent(out):: SijN(5*NN),U_par(6*NN)
      if(FTS_method)then
       call hydro_SD_FTS(APP,APQ,AQQ,rfu,rfe,rse,Fe,EIN,U_Par,SijN) 
      else
       call hydro_SD_FT(APP,rfu,Fe,U_Par)  
      endif

      end SUBROUTINE STEPPER_Gmres

      SUBROUTINE hydro_SD_FTS(APP,APQ,AQQ,rfu,rfe,rse,Fe,EIN,U_Par,SijN)     ! RETURNS NEW T=T+DT
      IMPLICIT NONE
      REAL*8,intent(in):: APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8,intent(in):: rfu(6*NN,6*NN),rfe(6*NN,5*NN),rse(5*NN,5*NN)
      REAL*8,intent(in):: EIN(5*NN),Fe(6*NN)
      real(8),intent(out):: U_Par(6*NN),SijN(5*NN)

      REAL*8 eval_RHS(11*NN),eval_LHS(11*NN)
      REAL*8 LPP(6*NN,6*NN),LPQ(6*NN,5*NN),LQQ(5*NN,5*NN)
      REAL*8 A_grmb(11*NN,11*NN),Rlub(11*NN,11*NN),grmb_LHS(11*NN,11*NN),grmb_RHS(11*NN,11*NN)
      real*8 RlubM(11*NN,11*NN),RM_LHS(11*NN,11*NN),RM_RHS(11*NN,11*NN)
      real*8 F1_RHS(11*NN),F2_RHS(11*NN),U1LHS(11*NN),U2LHS(11*NN),Fff(6*NN),Sff(5*NN)
      INTEGER I

      integer:: iter, nb
      real(8)  :: res,err

      call Part_to_Whole(APP,APQ,AQQ,A_grmb)
      call M_Tr_front(APP,APQ,AQQ,grmb_LHS,grmb_RHS)
      call Part_to_Whole(rfu,rfe,rse,Rlub)
      RlubM=matmul(Rlub,A_grmb)
      forall(i=1:11*NN)
        RlubM(i,i)=RlubM(i,i)+1.0_8
      end forall
      call Whole_to_Part(LPP,LPQ,LQQ,RlubM)
      call M_Tr_front(LPP,LPQ,LQQ,RM_LHS,RM_RHS)

      Fff=Fe
      F1_RHS(6*NN+1:11*NN)=-EIN
      F2_RHS(1:6*NN)=Fe

      !if(FTS_method)then

      Gmres: do i=1,50
        F1_RHS(1:6*NN)=Fff
        call simplegmres(U1LHS,11*NN,matmul(grmb_RHS,F1_RHS),grmb_LHS,iter,nb,res)
        F2_RHS(6*NN+1:11*NN)=U1LHS(6*NN+1:11*NN)
        call simplegmres(U2LHS,11*NN,matmul(RM_LHS,F2_RHS),RM_RHS,iter,nb,res)
        Fff=U2LHS(1:6*NN)
        U_Par=U1LHS(1:6*NN)
        Sff=U1LHS(6*NN+1:11*NN)


        eval_RHS(1:6*NN)=Fff+matmul(rfu,U_Par)+matmul(rfe,-EIN)
        eval_RHS(6*NN+1:11*NN)=matmul(transpose(APQ),Fff)+matmul(AQQ,Sff)
        eval_LHS(1:6*NN)=Fe
        eval_LHS(6*NN+1:11*NN)=-EIN
        err=norm_2(eval_LHS-eval_RHS)
       if ( err.lt.LINERR*1e-2 ) then
          !  do j=1,NN
          !    write(*,*) 'Fff,U0',Fff(3*(j-1)+1:3*j), U_Par(3*(j-1)+1:3*j)
          !  enddo
            U_Par=matmul(APP,Fff)+matmul(APQ,-EIN)
            SijN=U2LHS(6*NN+1:11*NN)
          !  do j=1,NN
          !    write(*,*) 'Fff,U1',Fff(3*(j-1)+1:3*j), U_Par(3*(j-1)+1:3*j)
          !  enddo
           write(*,*) 'err=',err
           print*, "gmres_iter success! iter=",i

           return

        endif
      enddo Gmres

         write(*,*) 'err=',err
         print*, "gmres_iter failed! iter=",i

      END



      SUBROUTINE hydro_SD_FT(APP,rfu,Fe,U_Par)      ! RETURNS NEW T=T+DT
      IMPLICIT NONE
      REAL*8,intent(in):: APP(6*NN,6*NN),rfu(6*NN,6*NN),Fe(6*NN)
      real(8),intent(out):: U_Par(6*NN)

      REAL*8 Fe_eval(6*NN)
      real*8 RlubM(6*NN,6*NN)
      real*8 Fff(6*NN)

      integer:: iter, nb,i
      real(8)  :: res,err

      RlubM=matmul(rfu,APP)
      forall(i=1:6*NN)
        RlubM(i,i)=RlubM(i,i)+1.0_8
      end forall

      call simplegmres(Fff,6*NN,Fe,RlubM,iter,nb,res)

      U_Par=matmul(APP,Fff)
      Fe_eval=Fff+matmul(rfu,U_Par)
      err=norm_2(Fe_eval-Fe)
      if ( err<LINERR*1e-6  ) then
        print*, "gmres_mod_gmres_test = ok"
      endif
      END

      SUBROUTINE GAUSS(G)
      IMPLICIT NONE
      REAL*8 G,PI
      PARAMETER(PI=3.141592653589793D0)
      REAL*8 A,B
      ! CALL RANDOM_SEED()
      CALL RANDOM_NUMBER(A)
      ! CALL RANDOM_SEED()
      CALL RANDOM_NUMBER(B)
      G=SQRT(-2*LOG(A)) * COS(2*PI*B)
      END

      SUBROUTINE DMDX_CAL(DMDX,CONF,RADII,Z,T)
      IMPLICIT NONE
      real*8,intent(in)::CONF(3,NN),RADII(NN),Z(3*NN),T
      real*8,intent(out)::DMDX(3*NN)

      REAL*8 CONF1(3,NN),CONF2(3,NN)
      REAL*8 APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      real*8 A1(3*NN,3*NN),A2(3*NN,3*NN)
      REAL*8 LATTICE(3,3),deta
      INTEGER::I,J
      PARAMETER(deta=0.001)


       forall(I=1:3,J=1:NN)
       CONF1(i,j)=CONF(i,j)+0.5_8*deta*Z(3*(j-1)+i)
       CONF2(i,j)=CONF(i,j)-0.5_8*deta*Z(3*(j-1)+i)
       end forall
      
      Brownsave=.true.
      if(.not.IsPeriod) then
       CALL GRPERY_INV_FRI(APP,APQ,AQQ,CONF1,RADII,NN)
       A1=APP(1:3*NN,1:3*NN)

       CALL GRPERY_INV_FRI(APP,APQ,AQQ,CONF2,RADII,NN)
       A2=APP(1:3*NN,1:3*NN)

       write(*,*) 'BROWN notperiod'
      else
       CALL SET_LATTICE_SKEW(T,LATTICE)
       CALL PER_GRPERY_INV_FRI(APP,APQ,AQQ,CONF1,RADII,NN,LATTICE,LINERR)
       A1=APP(1:3*NN,1:3*NN)
       CALL PER_GRPERY_INV_FRI(APP,APQ,AQQ,CONF2,RADII,NN,LATTICE,LINERR)
       A2=APP(1:3*NN,1:3*NN)
      endif
      Brownsave=.false.
      DMDX=matmul((A1-A2),Z)
      DMDX=DMDX/deta

      END SUBROUTINE DMDX_CAL
  end module STEPPER_MOD




#ifdef stokesian
      SUBROUTINE STEPPER_NOBROWN(CONF,RADII,DT,T,yeta_mu,U_pos,SijN)    ! RETURNS NEW T=T+DT
      IMPLICIT NONE
      REAL*8,intent(inout):: CONF(3,NN),U_pos(3*NN)
      REAL*8,intent(in):: RADII(NN),DT,T
      real(8),intent(out):: yeta_mu(5),SijN(5*NN)

      REAL*8 EIN(5*NN),TPP(3*NN,3*NN)
      REAL*8 APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8 RPP(6*NN,6*NN),Fe(6*NN)
      REAL*8 LATTICE(3,3),EWS_ERR
      REAL*8 uo_bg(6*NN),U_par(6*NN),u_add(3*NN)
      real*8 rfu(6*NN,6*NN),rfe(6*NN,5*NN),rse(5*NN,5*NN)
      logical pos_collision

      INTEGER I,J,Nt

      !character(128)::myformat

      call EIN_transY2(EIN,Eij,NN)
      call pos_collision_judge(conf,radii,pos_collision)
      call ppiclf_collision_timestep(RADII,DT,pos_collision)
      if(simplePeriod) then
        CALL SET_LATTICE_SKEW(T,LATTICE)
      endif
      if(.not.IsPeriod) then
       CALL GRPERY_INV_FRI(APP,APQ,AQQ,CONF,RADII,NN)
       write(*,*) 'notperiod'
      else
       EWS_ERR = LINERR
       CALL PER_GRPERY_INV_FRI(APP,APQ,AQQ,CONF,RADII,NN,LATTICE,EWS_ERR)
      endif

      if(Wall_method) then
        call Wall_GRPERY_INV_FRI(TPP,CONF,RADII,NN)
        APP(1:3*NN,1:3*NN)=APP(1:3*NN,1:3*NN)+TPP
      endif

      u_add=0.0_8

      NT_DEMdo:do Nt=1,NT_DEM
        call source_F(conf,RADII,U_pos,Fe)
        rfu=0.0_8
        rfe=0.0_8
        rse=0.0_8
        if(uselub)then
         call lubmxcalc(NN,conf,rfu,rfe,rse)
        endif

        if(T.gt.0.0_8.and.useDEMstepper) then
            call INIT_uo_bg(conf,uo_bg)
            U_par=0.0_8
            U_par(1:3*NN)=U_pos-uo_bg(1:3*NN)-u_add
            call STEPPER_UDEM(conf,radii,APP,APQ,AQQ,rfu,rfe,rse,Fe,EIN,RPP,U_par,SijN)
        else
          if(useGMRESmethod) then
            call STEPPER_Gmres(APP,APQ,AQQ,rfu,rfe,rse,Fe,EIN,U_par,SijN)
          else
            call STEPPER_NOGmres(APP,APQ,AQQ,rfu,rfe,rse,Fe,EIN,RPP,U_par,SijN)
          endif
        endif
   
        CALL STEPPER_conf(conf,radii,U_par,U_pos,u_add,T)

      enddo NT_DEMdo

      call yetamu_solve(APP,RPP,Fe,U_Par,RADII,SijN,EIN,yeta_mu)

      END SUBROUTINE STEPPER_NOBROWN

#endif