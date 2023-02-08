!****************************************
      MODULE method   ! CONF, W
      implicit none
      save
      INTEGER:: IsCaseLattice,collision_condition,NUM_THREADS,correction_method
      Logical:: IsPeriod,simplePeriod,useDEM,uselub,FTS_method,Boundary,Wall_method
      Logical:: Oscillation_shear,fix,useCollision,useGMRESmethod,useDEMstepper,useDLVO,usebond
      END MODULE method

!****************************************
      MODULE control   ! CONF, W
      implicit none
      save
      INTEGER:: K_WRITE,SEEDSIZE,INP,timestep,beginstep,interval
      INTEGER,ALLOCATABLE :: SEED(:) 
      LOGICAL THERE,Xshear
      REAL*8 END_TIME,BASE_DT,T0,dt_dem0,NtDEM
      END MODULE control

!****************************************
      MODULE CONFIG   ! CONF, W
      implicit none
      !save
      REAL*8, ALLOCATABLE :: CONF(:,:),RADII(:),CONFBDY(:,:),RADIIBDY(:),CONF_DEM(:,:)
      INTEGER, ALLOCATABLE :: Floc_index(:)
      REAL*8 W,u_bg(3),omega_bg(3),omegaT(3,3),Eij(3,3),EI_bg(6),radii_test,radii_true
      INTEGER POLY_LEN
      END MODULE CONFIG

!*************************************************
      MODULE LATTICE_BASE   ! LB(3)
      implicit none
      real*8,allocatable::Lb0_list(:)
      REAL*8 LB(3),LB0,LB_DEM(3)
      REAL*8 LR(3,3),LI(3,3),LV
      REAL*8 SIG,LINERR,LOGERR,Lammda,phi,phi_test
      INTEGER NR,NI,length_Lb0
      logical:: muCalculate
      END MODULE LATTICE_BASE

!****************************************
      MODULE TENSORS  
      implicit none
      real*8,parameter :: PAI=3.141592653589793238462643383279502884_8
      real*8,parameter :: sqrtPI=1.7724538509055160272981674833411_8
      real*8,parameter :: pai2=9.8696044010893586188344909998762_8
      real*8,parameter :: PIP5=0.564189583547756286948079451560772585844_8
      REAL*8 :: KDirac(3,3),EPS(3,3,3),Y2(5,3,3),Y21(5,3,3),Y22(5,3,3)
      REAL*8 :: II(3,3,3,3),Iunit(3,3,3,3),IIiso(3,3,3,3)
      logical:: front=.false.
      END MODULE TENSORS

!****************************************
      MODULE FORCE_PAR   ! k_PAR, L0_PAR, A_PAR, GAMMA
      implicit none
      REAL*8 :: K_PAR  = 10.D0           
      REAL*8 :: L0_PAR = 1.20D0   
      REAL*8 :: A_PAR  = 0.5D0
      REAL*8 :: A0_PAR  = 0.5D0
      REAL*8 :: LJ_EPS,LAMBDA,LJ_SIGMA,LJ_CUT
      INTEGER :: EQUILIBRATE = 0
      END MODULE FORCE_PAR

!****************************************
      MODULE SIZE       
      implicit none
      INTEGER NN,Nb,Np,NX,NY,NZ,Nfloc,Nfloc_sum,Nsuspen_sum
      END MODULE SIZE

      MODULE rb_conglomerate     
      implicit none
      INTEGER K_rb
      INTEGER,ALLOCATABLE :: KK_rbmconn(:) 
      !LOGICAL useCongl
      END MODULE rb_conglomerate    
!****************************************
      module prutil     !precision and utils
      implicit none
      integer, parameter :: SP=8
      integer, parameter :: cp=8
      real*8 ::  kd(3,3),per(3,3,3)
      END MODULE prutil
!**************************************************


!****************************************
      MODULE SYS_property   ! CONF, W
      implicit none
      !save
      REAL*8 ksp,erest,pho_par,t_collision,DT_DEM,F_adh,b0,hp0,heq
      REAL*8 pho_f,Re_number,mu_f,mass_par,gravity,Poisson_ratio,Youngs_module,g_gravity
      real*8 :: GAMMA,GAMMA_alpha,frequency  
      integer::Nt_DEM
      END MODULE SYS_property

      MODULE DLVO_property    ! CONF, W
      implicit none
      !save
      REAL*8 DebyeLength,kapa_EDL,Permittivity,z0,Lambda_l
      REAL*8 Potential,AH,Surfaceenergy
      END MODULE DLVO_property

      MODULE lambda_Dim   ! dimentionless
      implicit none
      !save
      REAL*8 lambda_pho,lambda_Length,lambda_mass,lambda_Viscosity,lambda_Permittivity,lambda_kB
      REAL*8 lambda_velocity,lambda_acceleration,lambda_time,lambda_Shearrate,lambda_Force
      REAL*8 lambda_potential,lambda_DebyeLength,lambda_AH,lambda_Surfaceenergy
      REAL*8 lambda_temperature,lambda_Pressure,lambda_Elasticmodulu
      REAL*8 lambda_ksp,lambda_erest,lambda_Renold
      END MODULE lambda_Dim

!****************************************
      MODULE BROWN_proerty   ! CONF, W
      implicit none
      !save
      Logical BROWN
      REAL*8 kBT,temperature,kB,k_LJ
      END MODULE BROWN_proerty
!****************************************

      MODULE Mobilitymatrix   ! CONF, W
      implicit none
      !save
      Logical Brownsave,Hydrosave
      REAL*8, ALLOCATABLE :: grmobmxsave(:,:)
      END MODULE Mobilitymatrix
!****************************************
  module dimentionless
  use SYS_property
  use DLVO_property
  use BROWN_proerty
  use lambda_Dim
  use CONFIG,only:radii_true,radii_test
  use TENSORS,only:PAI,pai2
  use control,only: END_TIME,BASE_DT,T0,dt_dem0,NtDEM
  implicit none
  private

  public::dimentionless_calc,dimentionless_translate
  contains

      SUBROUTINE dimentionless_calc()    ! dimentionless
      !IMPLICIT NONE
      !real*8 ::
      !integer:: i,j
      
      lambda_Renold=1.0_8
      lambda_Permittivity=1.0_8
      lambda_kB=1.0_8
      lambda_Viscosity=1.0_8/PAI*1.0e3
      lambda_velocity=1.0_8
      lambda_Length=radii_test/radii_true

      lambda_pho=lambda_Renold*lambda_Viscosity/(lambda_velocity*lambda_Length)
      lambda_mass=lambda_pho*lambda_Length*lambda_Length*lambda_Length
      lambda_time=lambda_Length/lambda_velocity
      lambda_acceleration=lambda_Length/(lambda_time*lambda_time)
      lambda_Force=lambda_mass*lambda_acceleration
      lambda_potential=lambda_Viscosity/(sqrt(lambda_pho))
      lambda_Surfaceenergy=lambda_Viscosity*lambda_Viscosity/(lambda_pho*lambda_Length)
      lambda_temperature=lambda_Length*lambda_Viscosity*lambda_Viscosity/(lambda_pho*lambda_kB)
      lambda_Pressure=lambda_Force/(lambda_Length*lambda_Length)
      lambda_DebyeLength=lambda_Length
      lambda_AH=lambda_Length*lambda_Viscosity*lambda_Viscosity/lambda_pho
      lambda_Shearrate=1.0_8/lambda_time
      lambda_ksp=lambda_Pressure*lambda_Length
      lambda_Elasticmodulu=lambda_Pressure

      end SUBROUTINE dimentionless_calc

      SUBROUTINE dimentionless_translate()    ! dimentionless
      !real*8 ::
      !integer:: i,j
      
      ksp=ksp*lambda_ksp
      pho_par=pho_par*lambda_pho
      pho_f=pho_f*lambda_pho
      mass_par=pho_par*4*PAI/3.0_8*1.0_8**3
      gravity=mass_par*g_gravity*lambda_acceleration
      mu_f=mu_f*lambda_Viscosity
      END_TIME=END_TIME*lambda_time
      BASE_DT=BASE_DT*lambda_time
      T0=T0*lambda_time
      

      dt_dem0 = sqrt(mass_par/ksp*(log(erest)**2+pai2))/10.0_8
      NtDEM=floor(BASE_DT/dt_dem0)

      DebyeLength=DebyeLength*lambda_Length
      kapa_EDL=1.0_8/DebyeLength
      z0=z0*lambda_Length
      Lambda_l=Lambda_l*lambda_Length
      Permittivity=Permittivity*lambda_Permittivity
      Potential=Potential*lambda_potential
      AH=AH*lambda_AH

      Youngs_module=Youngs_module*lambda_Elasticmodulu
      Surfaceenergy=Surfaceenergy*lambda_Surfaceenergy
      F_adh=1.50_8*PAI*Surfaceenergy*1.0_8
      b0=((9.0_8*PAI*Surfaceenergy*1.0_8*1.0_8*(1.0_8-Poisson_ratio*Poisson_ratio)) &
          & /(0.5_8*Youngs_module))**(1.0_8/3.0_8)
      hp0=b0*b0/(2**(1.0_8/3.0_8)*3.0_8*1.0_8)
      heq=hp0-(1.1_8)**(-3.0_8/5.0_8)*b0*b0/1.0_8
      kBT=kBT*lambda_kB*lambda_temperature
      k_LJ=30.0_8*kBT/1.0_8/1.0_8
      !kB=kB*lambda_kB


      end SUBROUTINE dimentionless_translate

  end module dimentionless