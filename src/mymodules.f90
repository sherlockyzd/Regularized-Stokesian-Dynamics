!****************************************
      MODULE method   ! CONF, W
      implicit none
      save
      INTEGER:: IsCaseLattice,collision_condition,NUM_THREADS,correction_method
      Logical:: IsPeriod,useDEM,uselub,FTS_method
      Logical:: Oscillation_shear,fix,useCollision,useGMRESmethod
      END MODULE method

!****************************************
      MODULE control   ! CONF, W
      implicit none
      save
      INTEGER:: K_WRITE,SEEDSIZE,INP,timestep,beginstep,interval
      INTEGER,ALLOCATABLE :: SEED(:) 
      LOGICAL BROWN, THERE,Xshear
      REAL*8 END_TIME,BASE_DT,T0
      END MODULE control

!****************************************
      MODULE CONFIG   ! CONF, W
      implicit none
      !save
      REAL*8, ALLOCATABLE :: CONF(:,:),RADII(:)
      REAL*8 W,u_bg(3),omega_bg(3),omegaT(3,3),Eij(3,3),EI_bg(6),radii_test
      INTEGER POLY_LEN
      END MODULE CONFIG

!****************************************
      MODULE SYS_property   ! CONF, W
      implicit none
      !save
      REAL*8 ksp,erest,pho_par,t_collision,DT_DEM
      REAL*8 pho_f,Re_number,mu_f 
      real*8 :: GAMMA,GAMMA_alpha,frequency  
      integer::Nt_DEM
      END MODULE SYS_property

!****************************************
      MODULE LATTICE_BASE   ! LB(3)
      implicit none
      real*8,allocatable::Lb0_list(:)
      REAL*8 LB(3),LB0
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
      INTEGER NN,NX,NY,NZ
      END MODULE SIZE

!****************************************
      module prutil     !precision and utils
      implicit none
      integer, parameter :: SP=8
      integer, parameter :: cp=8
      real*8 ::  kd(3,3),per(3,3,3)
      END MODULE prutil
!**************************************************