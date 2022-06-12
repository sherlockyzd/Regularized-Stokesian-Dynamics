!***********************************************************
!***********************************************************
!***********************************************************

    module period_bdy_tools
    use size
    use CONFIG, only:RADII 
    USE LATTICE_BASE
    use method,only:Oscillation_shear
    USE SYS_property     ! GAMMA0,GAMMA,frequency
    use method,only:Oscillation_shear
    USE FORCE_PAR      ! GAMMA
    use tensors, only:pai
      
       
    implicit none
    private


    public::PER_SKEW_CORRECTION,PERIOD_CORRECTION,SET_LATTICE_SKEW,DISTANCE_PER_3D,CALC_LATTICE_INVERSE

  contains

      SUBROUTINE PER_SKEW_CORRECTION(R,DMIN)
      IMPLICIT NONE
      REAL*8,intent(inout)::R(3),DMIN
      REAL*8 RC(3),RL(3),RMIN(3),X(3),M(3),D
      INTEGER I,M1,M2,M3
      
      X=MATMUL(LI,R)

      X=MOD(X,1.D0)

      DO I=1,3
      IF(X(I).LT.0.D0) X(I)=X(I)+1.D0
      ENDDO

      RC=MATMUL(LR,X)
      DMIN=SUM(RC**2)
      RMIN=RC

      DO M1=0,1
      DO M2=0,1
      DO M3=0,1
      IF(M1==0.AND.M2==0.AND.M3==0) CYCLE
       M(1)=M1
       M(2)=M2
       M(3)=M3
       RL=MATMUL(LR,M)
       D=SUM( (RC-RL)**2 )
       IF(D.LT.DMIN) THEN
        DMIN=D
        RMIN=RC-RL
       ENDIF
      ENDDO
      ENDDO
      ENDDO

      R=RMIN
      DMIN=SQRT(DMIN)

      END SUBROUTINE PER_SKEW_CORRECTION

      SUBROUTINE PERIOD_CORRECTION(R,LB,T)
      IMPLICIT NONE
      REAL*8,intent(out)::R(3)
      REAL*8,intent(in):: T,LB(3)
      
      real*8 Lshere
!   CORRECTIONS FOR PERIOD

      if(Oscillation_shear) then
        Lshere=GAMMA_alpha*sin(frequency*T)*LB(3)
      else
        Lshere=LB(3)*GAMMA*T
      end if


      IF(R(3).GT.LB(3)) THEN
       R(3) = R(3) - LB(3)
       R(1) = MOD(R(1) - Lshere,LB(1))
      ENDIF
      IF(R(3).LT.0) THEN
       R(3) = R(3) + LB(3)
       R(1) = MOD(R(1) + Lshere,LB(1))
      ENDIF

      IF(R(1).GT.LB(1)) THEN
       R(1)=R(1)-LB(1)
      ELSE IF(R(1).LT.0) THEN
       R(1)=R(1)+LB(1)
      ENDIF

      IF(R(2).GT.LB(2)) THEN
       R(2)=R(2)-LB(2)
      ELSE IF(R(2).LT.0) THEN
       R(2)=R(2)+LB(2)
      ENDIF


      END SUBROUTINE PERIOD_CORRECTION


      SUBROUTINE SET_LATTICE_SKEW(T,LATTICE)
      IMPLICIT NONE
      REAL*8,intent(in):: T
      REAL*8,intent(out):: LATTICE(3,3)  ! MM = METRIC MATRIX, LE = EIGENVALUS
      
      
      REAL*8 TSHEAR,TMOD,angle_shere,Lshere

      LATTICE=0.D0
      LATTICE(1,1)=LB(1)
      LATTICE(2,2)=LB(2)
      LATTICE(3,3)=LB(3)

      if(Oscillation_shear) then
        angle_shere=GAMMA_alpha*sin(frequency*T)
        Lshere=angle_shere*LB(3)
        LATTICE(1,3)=MOD(Lshere,LB(1))
      else
        IF (GAMMA .NE. 0.0) THEN
          TSHEAR=LB(1)/(LB(3)*GAMMA)
          TMOD=MOD(T,TSHEAR)
          LATTICE(1,3)=LB(3)*TMOD*GAMMA
        ELSE
          LATTICE(1,3)=0.0
        ENDIF 
      endif

      END SUBROUTINE SET_LATTICE_SKEW

      SUBROUTINE DISTANCE_PER_3D(CONF)
      IMPLICIT NONE
      REAL*8,intent(in):: CONF(3,NN)
      REAL*8 R(3),DIST
      INTEGER I,J

      DO I=1,NN-1
      DO J=I+1,NN
      R=CONF(:,I)-CONF(:,J)      
      CALL PER_SKEW_CORRECTION(R,DIST)      
      IF( DIST.LE.1.D0 ) THEN
       WRITE(*,*) 'WARNING! Periodic distance between particles'
       WRITE(*,*) I,' and ',J,' is less or equal 1.0'
      ENDIF
      ENDDO
      ENDDO

      END SUBROUTINE DISTANCE_PER_3D


      SUBROUTINE CALC_LATTICE_INVERSE(LATTICE,EWS_ERR)
      IMPLICIT NONE
      REAL*8,intent(in)::LATTICE(3,3),EWS_ERR
      
      REAL*8 MM(3,3),LE(3)  ! MM = METRIC MATRIX, LE = EIGENVALUS
      REAL*8 RNR,RNI,LMIN

      LR = LATTICE
      CALL MAT_INV3(LI,LR,LV)
      LV=ABS(LV)
      SIG=(LV)**(1.D0/3.D0)/SQRT(2.D0*pai)
      !Lammda=1.0_8/(sqrt(2.0_8)*SIG)
      Lammda=SQRT(pai)/(LV)**(1.D0/3.D0)
      LOGERR=-LOG(EWS_ERR)
      MM=MATMUL(TRANSPOSE(LR),LR)      
      CALL EIGEN(MM,LE,3)
      LMIN=MINVAL(LE)
      RNR=2.D0*LOGERR*SIG**2/LMIN
      NR=int(SQRT(RNR)+1)

      LMIN=1.D0/MAXVAL(LE)
      RNI=LOGERR/(2.D0*pai**2*SIG**2*LMIN)
      NI=int(SQRT(RNI)+1)

      phi=sum(RADII(1:NN)**3)*4.0_8/3.0_8*Pai/LV

      END
  end module period_bdy_tools

!********************************************************

      module master_time_m

      use, intrinsic :: iso_fortran_env, only: stdout => output_unit
      implicit none
      private

      public :: tic, toc, time_print

      real(8), save :: time_save !! save the time

      contains

        subroutine tic()
        !! start the timer
            call cpu_time(time_save)
        end subroutine tic

        subroutine toc(t)
        !! stop the timer and return the time in seconds
            real(8), intent(out), optional :: t
            real(8) :: time_now

            call cpu_time(time_now)
            time_now = time_now - time_save

            if (present(t)) then
                t = time_now
            else
                write (stdout, "(A, ES20.10, A)") 'Time elapsed: ', time_now, " s"
            end if
        end subroutine toc

        subroutine time_print()
        !! print the current time and date
            character(len=8) :: datstr
            character(len=10) :: timstr

            ! Get the current date and time.
            call date_and_time(datstr, timstr)

            ! Write out the date and time.
            write (stdout, "(A)") "Date = "//datstr(1:4)//"/"// &
                datstr(5:6)//"/"// &
                datstr(7:8)
            write (stdout, "(A)") "Time = "//timstr(1:2)//":"// &
                timstr(3:4)//":"// &
                timstr(5:10)

        end subroutine time_print

      end module master_time_m
!*********************************************************************


  module hydro_tools
  use size
  use tensors
  use method
  use LATTICE_BASE
  use SYS_property
  use tensors,only:pai,pai2,EPS
  USE CONFIG    ! CONF,POLY_LEN
  use period_bdy_tools,only:PER_SKEW_CORRECTION
  implicit none
  private


    public::conf_overlap_correction,yetamu_solve,Init_frequency, &
      & ppiclf_collision_timestep,INIT_u_bg,INIT_uo_bg

  contains



     subroutine conf_overlap_correction(CONF,U_Par,NN,radii)
      implicit none
      INTEGER,intent(in):: NN
      REAL*8,intent(in):: conf(3,NN),radii(NN)
      real(8),intent(inout)::U_Par(3*NN)

      real*8 rx(3),valpha_per(3),vbeta_per(3),valpha(3),vbeta(3)
      real*8 r,r_ave,epsila
      integer alpha,beta

      do alpha=1,NN
        do beta=alpha+1,NN
        rx=conf(:,beta)-conf(:,alpha)

        if(IsPeriod) then
         CALL PER_SKEW_CORRECTION(rx,r)
        else
         r=sqrt(sum(rx*rx))
        endif
        r_ave=0.5_8*(radii(alpha)+radii(beta))
        epsila=(r-2.0_8*r_ave)/r_ave

        if(epsila.lt.0.0_8) then
                write(*,*) 'Contact Error, r=',r,alpha,beta
                write(*,*) 'make correction'
          valpha=U_Par(3*(alpha-1)+1:3*alpha)
          vbeta=U_Par(3*(beta-1)+1:3*beta)

          valpha_per=sum(valpha*rx)/r*rx
          vbeta_per=sum(vbeta*rx)/r*rx
          if(collision_condition==0)  then         !perfectly elastic collision
           U_Par(3*(alpha-1)+1:3*alpha)=U_Par(3*(alpha-1)+1:3*alpha)-2*valpha_per
           U_Par(3*(beta-1)+1:3*beta)=U_Par(3*(beta-1)+1:3*beta)-2*vbeta_per
          elseif(collision_condition==1) then      !Viscoelastic collision
           U_Par(3*(alpha-1)+1:3*alpha)=U_Par(3*(alpha-1)+1:3*alpha)-valpha_per
           U_Par(3*(beta-1)+1:3*beta)=U_Par(3*(beta-1)+1:3*beta)-vbeta_per
          else
           U_Par(3*(alpha-1)+1:3*alpha)=U_Par(3*(alpha-1)+1:3*alpha)-valpha_per
           U_Par(3*(beta-1)+1:3*beta)=U_Par(3*(beta-1)+1:3*beta)-vbeta_per
          endif
        endif
        enddo
      enddo
     end subroutine conf_overlap_correction


      subroutine yetamu_solve(APP,RPP,FE,U_Par,RADII,SijN,EIN,yeta_mu)
      IMPLICIT none
      REAL*8,intent(in):: APP(6*NN,6*NN),RPP(6*NN,6*NN),radii(NN)
      REAL*8,intent(in):: Fe(6*NN)
      real(8),intent(in):: U_Par(6*NN),SijN(5*NN),EIN(5*NN)
      real*8,intent(out)::yeta_mu(5)

      real*8 U_ave,KperTT,KperRR
      integer i

      yeta_mu=0.0_8

      do I=1,NN
        yeta_mu(5)=yeta_mu(5)+SijN(4+5*(I-1))
      enddo
      yeta_mu(5)=yeta_mu(5)/mu_f
      do I=1,NN
          yeta_mu(1)=yeta_mu(1)+SijN(4+5*(I-1))/(-2.0_8*EiN(4+5*(I-1)))
      enddo
      yeta_mu(1)=yeta_mu(1)/LV/mu_f
      yeta_mu(2)= phi

      call averageU(APP,U_ave,NN) 
      yeta_mu(3)=1./U_ave

      if(useGmresmethod) then
        U_ave=0.0_8
        do I=1,NN
          yeta_mu(4)=yeta_mu(4)+Fe(3*(i-1)+1)/NN
          U_ave=U_ave+(U_Par(3*(i-1)+1)*RADII(i))/NN
        enddo

        yeta_mu(4)=yeta_mu(4)/(U_ave*6*PAI*mu_f)
      else
        !call Mualphabeta(mu,alpha,beta,NN)
        !yeta_mu(4)=alpha
        !yeta_mu(5)=beta
        call Kpermeability(KperTT,KperRR,RPP,NN)
        yeta_mu(4)=KperTT
      endif

      end subroutine yetamu_solve


      SUBROUTINE Init_frequency(LB,radii,frequency)
      IMPLICIT NONE
      real*8,intent(in)::LB(3),radii(NN)
      real*8,intent(out)::frequency
      real*8 Lm,aa,Gamma_amp
      
        Lm=minval(LB)
        aa=minval(radii)
        Gamma_amp=Re_number*mu_f/(pho_f*aa*Lm)
        frequency=Gamma_amp/Gamma_alpha
      end SUBROUTINE Init_frequency


      SUBROUTINE ppiclf_collision_timestep(RADII,DT)
      IMPLICIT NONE
      real*8,intent(in)::RADII(NN),DT
      !real*8,intent(out)::DT_DEM
      !integer,intent(out):: NT_DEM

      real*8 rm1, rm2,rmult,dt_dem0
        if(usecollision)then
          rm1 = pho_par*4.0_8/3.0_8*pai*minval(RADII)**3
          rm2 = pho_par*4.0_8/3.0_8*pai*(0.5_8*(minval(RADII)+maxval(radii)))**3
           
          rmult = 1.0d0/sqrt(1.0d0/rm1+1.0d0/rm2)
          dt_dem0 = sqrt(rmult/ksp*(log(erest)**2+pai2))/10.0_8
          if(DT.lt.dt_dem0)then
            NT_DEM=10
          elseif(DT.gt.10.0_8*dt_dem0)then
            NT_DEM=10
          else
            NT_DEM=max(floor(DT/dt_dem0),5)
          endif
          WRITE(*,*) 'Nt_DEM=',Nt_DEM
          WRITE(*,*) 'dt_dem0=',dt_dem0
        else
          NT_DEM=1
        endif
        DT_DEM=DT/real(NT_DEM)
      end SUBROUTINE ppiclf_collision_timestep


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE INIT_u_bg()
      !USE SIZE      ! NN
      !USE LATTICE_BASE ! LB LB0
      !use method
      IMPLICIT NONE
      INTEGER I,J

      ! background Einf11,Einf22,Einf33,Einf12,Einf13,Einf23
        ! fill in the -Einf Vector
        ! E vector conversion
        ! EV1=E11-E33, EV2=2E12, EV3=2E13, EV4=2E23, EV5=E22-E33
        ! or EIN
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
      !write(*,*)'T=, GAMMA=', T,GAMMA   
     ! write(*,*)'Eij=',Eij
      !write(*,*)'omegaT=',omegaT

     end



      SUBROUTINE INIT_uo_bg(conf,uo_bg)
      IMPLICIT NONE
      real*8,intent(in)::CONF(3,NN)
      real*8,intent(out)::uo_bg(6*NN)

      real*8 u_bg_local(3)
      INTEGER I,J

      do i=1,NN
        do j=1,3
        ! U_bg_local=U_bg + Omega_bg * x + Einf_bg . x
        u_bg_local(j)=u_bg(j)+(OmegaT(j,1)*CONF(1,i)          &
          +OmegaT(j,2)*CONF(2,i)+OmegaT(j,3)*CONF(3,i)) &
          +(Eij(j,1)*CONF(1,i)+Eij(j,2)*CONF(2,i)       &
          +Eij(j,3)*CONF(3,i))
        ! fill in the U-Uinf, Omega-Omegainf vector
        uo_bg(3*(i-1)+j)=u_bg_local(j)
        uo_bg(3*NN+3*(i-1)+j)=omega_bg(j)
        end do
      end do

     end


 

  end module hydro_tools











#ifdef test
write(myformat,'(a,I6,a)')'(a, ',6*NN,'E24.15)'
  open(unit=25,file='APP.txt',position='append',form='formatted') 
   do i=1,6*NN
    write(25,myformat)'APP0',APP(i,:)
  enddo
   do i=1,6*NN
    write(25,myformat)'APP1',APP1(i,:)
  enddo
  close(25)

write(myformat,'(a,I6,a)')'(a, ',5*NN,'E24.15)'
  open(unit=25,file='APQ.txt',position='append',form='formatted') 
   do i=1,6*NN
    write(25,myformat)'APQ0',APQ(i,:)
  enddo
   do i=1,6*NN
    write(25,myformat)'APQ1',APQ1(i,:)
  enddo
  close(25)

write(myformat,'(a,I6,a)')'(a, ',5*NN,'E24.15)'
  open(unit=25,file='AQQ.txt',position='append',form='formatted') 
   do i=1,5*NN
    write(25,myformat)'AQQ0',AQQ(i,:)
  enddo
   do i=1,5*NN
    write(25,myformat)'AQQ1',AQQ1(i,:)
  enddo
  close(25)
#endif

!  open(unit=25,file='grmobmx.txt',position='append',form='formatted') 
!  do i=1,11*ntotal
!    write(25,myformat)filetime,grmobmx(i,:)
!  enddo
!  close(25)
      !write(*,*) 'NR,NI=',NR,NI

#ifdef test
      write(*,*) 'AQQ=',AQQ(1,1:5)
      write(*,*) AQQ(1,6:10)
      write(*,*) AQQ(2,1:5)
      !write(*,*) AQQ(6,6:10)
      write(*,*) 'RQQ=',RQQ(1,1:5)
      write(*,*) RQQ(1,6:10)
      write(*,*) RQQ(2,1:5)
#endif

#ifdef test
      write(*,*) 'mu=',mu(1,1:5)
      write(*,*) mu(2,1:5)
      write(*,*) mu(3,1:5)
      write(*,*) mu(4,1:5)
      write(*,*) mu(5,1:5)
      write(*,*) 'mu=',mu(6,1:5)
      write(*,*) mu(7,1:5)
      write(*,*) mu(8,1:5)
      write(*,*) mu(9,1:5)
      write(*,*) mu(10,1:5)
#endif
