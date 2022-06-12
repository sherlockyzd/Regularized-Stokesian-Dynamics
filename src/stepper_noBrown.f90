  module STEPPER_NOBROWN_MOD
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

  implicit none
  private


    public::STEPPER_NOBROWN

  contains


      SUBROUTINE STEPPER_NOBROWN(CONF,RADII,DT,T,yeta_mu,U_pos)    ! RETURNS NEW T=T+DT
      IMPLICIT NONE
      REAL*8,intent(inout):: CONF(3,NN),U_pos(3*NN)
      REAL*8,intent(in):: RADII(NN),DT,T
      real(8),intent(out):: yeta_mu(5)

      REAL*8 EIN(5*NN),SijN(5*NN)
      REAL*8 APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8 RPP(6*NN,6*NN),Fe(6*NN)
      REAL*8 LATTICE(3,3),EWS_ERR
      REAL*8 uo_bg(6*NN),U_par(6*NN)
      REAL*8 DCONF_1(3*NN)
      real*8 rfu(6*NN,6*NN),rfe(6*NN,5*NN),rse(5*NN,5*NN)

      INTEGER I,J,Nt

      !character(128)::myformat

      call EIN_transY2(EIN,Eij,NN)
      call ppiclf_collision_timestep(RADII,DT)
      
      if(.not.IsPeriod) then
       CALL GRPERY_INV_FRI(APP,APQ,AQQ,CONF,RADII,NN)
       write(*,*) 'notperiod'
      else
       CALL SET_LATTICE_SKEW(T,LATTICE)
       EWS_ERR = LINERR
       CALL PER_GRPERY_INV_FRI(APP,APQ,AQQ,CONF,RADII,NN,LATTICE,EWS_ERR)
      endif
      call INIT_uo_bg(conf,uo_bg)

      NT_DEMdo:do Nt=1,NT_DEM
        call source_F(conf,RADII,U_pos,Fe)
        rfu=0.0_8
        rfe=0.0_8
        rse=0.0_8
        if(uselub)then
         call lubmxcalc(NN,conf,rfu,rfe,rse)
        endif

        if(useGMRESmethod) then
          call STEPPER_Gmres(APP,APQ,AQQ,rfu,rfe,rse,Fe,EIN,U_par,SijN)
        else
          call STEPPER_NOGmres(APP,APQ,AQQ,rfu,rfe,rse,Fe,EIN,RPP,U_par,SijN)
        endif
        
        U_pos=U_par(1:3*NN)+uo_bg(1:3*NN)
        if(.not.fix) then
          DCONF_1=DT_DEM*U_pos
          forall(i=1:NN,j=1:3)
            CONF(j,i)=CONF(j,i)+DCONF_1(3*(i-1)+j)
          end forall
          if(IsPeriod) then
           DO I = 1, NN
            CALL PERIOD_CORRECTION(CONF(:,I),LB,T)
           ENDDO
          endif
          call INIT_uo_bg(conf,uo_bg)
          U_pos=U_par(1:3*NN)+uo_bg(1:3*NN)
        endif
      enddo NT_DEMdo

      call yetamu_solve(APP,RPP,Fe,U_Par,RADII,SijN,EIN,yeta_mu)

      END SUBROUTINE STEPPER_NOBROWN

!***********************************************************

      SUBROUTINE STEPPER_NOGmres(APP,APQ,AQQ,rfu,rfe,rse,Fe,EIN,RPP,U_par,SijN)
      IMPLICIT NONE
      REAL*8,intent(in):: APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8,intent(in):: rfu(6*NN,6*NN),rfe(6*NN,5*NN),rse(5*NN,5*NN)
      real*8,intent(in):: EIN(5*NN),Fe(6*NN)
      real(8),intent(out):: SijN(5*NN),RPP(6*NN,6*NN),U_par(6*NN)

      REAL*8 RPQ(6*NN,5*NN),RQQ(5*NN,5*NN)
      REAL*8 RPP0(6*NN,6*NN),RPQ0(6*NN,5*NN),RQQ0(5*NN,5*NN)
      real*8 mu(5*NN,5*NN),APPR(6*NN,6*NN)

        RPQ=0.0_8
        RQQ=0.0_8
        if(FTS_method)then
         call MobToResist_FTS(APP,APQ,AQQ,NN,RPP0,RPQ0,RQQ0)
        else
         call MobToResist_FT(APP,NN,RPP0)
        endif

        call arrangResist(NN,rfu,rfe,rse,RPP0,RPQ0,RQQ0,RPP,RPQ,RQQ)
        APPR=RPP  
        CALL MATREV(APPR,6*NN)
        if(FTS_method)then
         U_par = matmul(APPR,Fe+matmul(RPQ,EIN))
        else
         U_par = matmul(APPR,Fe)
        endif
        call arrangement(RPP,RPQ,RQQ,mu,NN)
        SijN=matmul(mu,EIN)

      end SUBROUTINE STEPPER_NOGmres


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


  end module STEPPER_NOBROWN_MOD