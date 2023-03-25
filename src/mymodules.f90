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
      REAL*8, ALLOCATABLE :: CONF(:,:),RADII(:)!,U_pos(:)
      REAL*8, ALLOCATABLE :: CONFBDY(:,:),RADIIBDY(:),CONF_DEM(:,:)
      INTEGER, ALLOCATABLE :: Floc_index(:)
      REAL*8 W,u_bg(3),omega_bg(3),omegaT(3,3),Eij(3,3),EI_bg(6),radii_test,radii_true
      INTEGER POLY_LEN
      real*8 alphaX,betaY,gamaZ 
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
      INTEGER NN,Nb,Np,NX,NY,NZ,Nfilament,Nswimer
      INTEGER Nfloc,Nfloc_sum,Nsuspen_sum
      END MODULE SIZE 
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
      REAL*8 pho_f,Re_number,mu_f,mass_par,gravity,Poisson_ratio,Youngs_module,Shear_module,g_gravity
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
    module quaternions

      implicit none

      type quaternion
        double precision :: r
        double precision :: i
        double precision :: j
        double precision :: k
      end type quaternion


    interface operator(+)
      module procedure qadd
    end interface 

    interface operator(-)
      module procedure qsub
    end interface 

    interface operator(*)
      module procedure qmul
    end interface 

    interface operator(*)
      module procedure qscal
    end interface 

    interface operator(/)
      module procedure qdivision
    end interface 

    contains



    function quat(r, i, j, k) result(q)

      implicit none
      double precision, intent(in) :: r, i, j, k
      type(quaternion) :: q

      q%r = r
      q%i = i
      q%j = j
      q%k = k

    end function quat

    function qconst(c) result(q)

      implicit none
      double precision, intent(in) :: c
      type(quaternion) :: q

      q%r = c
      q%i = c
      q%j = c
      q%k = c

    end function qconst

    function qadd(x, y) result(q)

      implicit none
      type(quaternion), intent(in) :: x, y
      type(quaternion) :: q

      q%r = x%r + y%r
      q%i = x%i + y%i
      q%j = x%j + y%j
      q%k = x%k + y%k

    end function qadd


    function qsub(x, y) result(q)

      implicit none
      type(quaternion), intent(in) :: x, y
      type(quaternion) :: q

      q%r = x%r - y%r
      q%i = x%i - y%i
      q%j = x%j - y%j
      q%k = x%k - y%k

    end function qsub


    function qscal(a, x) result(q)

      implicit none
      type(quaternion), intent(in) :: x
      double precision, intent(in) :: a
      type(quaternion) :: q

      q%r = a * x%r
      q%i = a * x%i
      q%j = a * x%j
      q%k = a * x%k

    end function qscal


    function qconj(x) result(q)

      implicit none
      type(quaternion), intent(in) :: x
      type(quaternion) :: q

      q%r =  x%r
      q%i = -x%i
      q%j = -x%j
      q%k = -x%k

    end function qconj


    function qnorm(x) result(res)

      implicit none
      type(quaternion), intent(in) :: x
      double precision :: res

      res = sqrt(x%r * x%r + x%i * x%i + x%j * x%j + x%k * x%k)

    end function qnorm

    subroutine qprint(q)

      implicit none
      type(quaternion) :: q

      write(*,*) "quaternion and norm:",q%r, q%i, q%j, q%k,qnorm(q)
        
    end subroutine qprint

    function qnorm2(x) result(res)

      implicit none
      type(quaternion), intent(in) :: x
      double precision :: res

      res = x%r * x%r + x%i * x%i + x%j * x%j + x%k * x%k

    end function qnorm2


    function qinv(x) result(q)

      implicit none
      type(quaternion), intent(in) :: x
      type(quaternion) :: q

      q = qscal(1.0d0/qnorm2(x) , qconj(x))

    end function qinv


    function qmul(x, y) result(q)

      implicit none
      type(quaternion), intent(in) :: x, y
      type(quaternion) :: q

      q%r = x%r*y%r - x%i*y%i - x%j*y%j - x%k*y%k
      q%i = x%r*y%i + x%i*y%r + x%j*y%k - x%k*y%j
      q%j = x%r*y%j - x%i*y%k + x%j*y%r + x%k*y%i
      q%k = x%r*y%k + x%i*y%j - x%j*y%i + x%k*y%r

    end function qmul

    function dqsub(q1, q0) result(dq)

      implicit none
      type(quaternion), intent(in) :: q1, q0
      type(quaternion) :: dq

      dq=qmul(q1,qconj(q0)) 

    end function dqsub

    function qequal(q1,q0) result(YN)
    implicit none
    type(quaternion), intent(in) :: q1, q0
    type(quaternion) :: dq
    real*8:: q_norm
    logical:: YN
    dq=q1-q0
    q_norm=qnorm(dq)
    YN=q_norm.lt.1e-6
    end function qequal


    function qangle(q1,q0) result(angle)
    implicit none
    type(quaternion), intent(in) :: q1, q0
    type(quaternion) :: dq0,dq
    real*8:: angle(3),coshalfsita,sinhalfsita(3),nn(3),halfsita

    dq0=dqsub(q1, q0)
    !call qprint(dq0)
    if(qnorm(dq0).lt.1e-6)then
      angle=0.0_8
    else
      dq=qscal(1.0_8/qnorm(dq0) , dq0)
      coshalfsita=dq%r
      halfsita=acos(coshalfsita)
      sinhalfsita(1)=dq%i;sinhalfsita(2)=dq%j;sinhalfsita(3)=dq%k;
      nn=sinhalfsita/sqrt(sum(sinhalfsita**2))
      angle=2.0*halfsita*nn
      write(*,*) 'angle=',angle(:)
    endif
    end function qangle


    function qdivision(x, a) result(q)

      implicit none
      type(quaternion), intent(in) :: x
      double precision, intent(in) :: a
      type(quaternion) :: q

      q%r = x%r / a
      q%i = x%i / a
      q%j = x%j / a
      q%k = x%k / a

    end function qdivision




    function qderivat(q, w) result(dq)

      implicit none
      type(quaternion), intent(in) :: q
      double precision, intent(in) :: w(3)
      type(quaternion) :: dq,qw!,q_norm

      qw=quat(0.0_8, w(1), w(2), w(3))
      dq=0.5_8*qmul(qw, q)
      !dq=0.5_8*qmul(q,qw)

    end function qderivat


    function qMidpoint(q1, q2) result(qmid)

      implicit none
      type(quaternion), intent(in) :: q1,q2
      !double precision, intent(in) :: w(3)
      type(quaternion) :: q_conj,q_rot,q_sqrt,qmid!,q_norm
      q_conj=qconj(q1)
      q_rot=qmul(q2,q_conj)
      q_sqrt=qsqrt(q_rot)
      qmid=qmul(q_sqrt,q1)

    end function qMidpoint

    function qsqrt(q) result(q_sqrt)
      implicit none
      type(quaternion), intent(in) :: q
      double precision :: s,w(3),s1,w1(3)
      type(quaternion) :: q_sqrt
      s=q%r;w(1)=q%i;w(2)= q%j; w(3) = q%k;
      if(s.lt.-0.99999) then
         q_sqrt=quat(0.0_8,0.0_8,0.0_8,1.0_8)
      else
         s1=sqrt((s+1.0_8)*0.5_8)
         w1=w/(s1*2.0_8)
         q_sqrt=quat(s1,w1(1),w1(2),w1(3))
      endif
    end function qsqrt



    function qrot(q,v0)   result(v)
    implicit none
    type(quaternion), intent(in) :: q
    real*8, intent(in)::v0(3)
    !real*8, intent(out)::
    real*8 :: R(3,3),s,vx,vy,vz,v(3)
    type(quaternion) :: q_norm
    q_norm=qscal(1.0_8/qnorm(q) , q)
    R=quat2rotmat(q_norm)
    v=matmul(R,v0)

    end function qrot

    function qrot_inverse(q,v0)   result(v)
    implicit none
    type(quaternion), intent(in) :: q
    real*8, intent(in)::v0(3)
    !real*8, intent(out)::
    real*8 :: R(3,3),s,vx,vy,vz,v(3)
    type(quaternion) :: q_norm
    q_norm=qscal(1.0_8/qnorm(q) , q)
    R=transpose(quat2rotmat(q_norm))
    v=matmul(R,v0)

    end function qrot_inverse

    function quat2rotmat(q_norm)   result(R)
    implicit none
    type(quaternion), intent(in) :: q_norm
    
    !real*8, intent(out)::
    real*8 :: R(3,3),s,vx,vy,vz!,v(3)
    !type(quaternion) :: q_norm
    !q_norm=qscal(1.0_8/qnorm(q) , q)
    s = q_norm%r; vx = q_norm%i; vy = q_norm%j; vz = q_norm%k;
    R(1,1) = 1.0_8 - 2.0_8*vy**2 - 2.0_8*vz**2
    R(1,2) = 2.0_8*vx*vy - 2.0_8*s*vz
    R(1,3) = 2.0_8*vx*vz + 2.0_8*s*vy
    R(2,1) = 2.0_8*vx*vy + 2.0_8*s*vz
    R(2,2) = 1.0_8 - 2.0_8*vx**2 - 2.0_8*vz**2
    R(2,3) = 2.0_8 *vy*vz - 2.0_8 *s*vx
    R(3,1) = 2.0_8 *vx*vz - 2.0_8 *s*vy
    R(3,2) = 2.0_8 *vy*vz + 2.0_8 *s*vx
    R(3,3) = 1.0_8  - 2.0_8 *vx**2 - 2.0_8 *vy**2
    end function quat2rotmat

    function quat2rotmatinverse(q_norm)   result(R)
    implicit none
    type(quaternion), intent(in) :: q_norm
    
    !real*8, intent(out)::
    real*8 :: R(3,3),s,vx,vy,vz!,v(3)
    !type(quaternion) :: q_norm
    !q_norm=qscal(1.0_8/qnorm(q) , q)
    s = q_norm%r; vx = q_norm%i; vy = q_norm%j; vz = q_norm%k;
    R(1,1) = 1.0_8 - 2.0_8*vy**2 - 2.0_8*vz**2
    R(1,2) = 2.0_8*vx*vy + 2.0_8*s*vz
    R(1,3) = 2.0_8*vx*vz - 2.0_8*s*vy
    R(2,1) = 2.0_8*vx*vy - 2.0_8*s*vz
    R(2,2) = 1.0_8 - 2.0_8*vx**2 - 2.0_8*vz**2
    R(2,3) = 2.0_8 *vy*vz + 2.0_8 *s*vx
    R(3,1) = 2.0_8 *vx*vz + 2.0_8 *s*vy
    R(3,2) = 2.0_8 *vy*vz - 2.0_8 *s*vx
    R(3,3) = 1.0_8  - 2.0_8 *vx**2 - 2.0_8 *vy**2
    end function quat2rotmatinverse



    subroutine qsetzero(q)

      implicit none
      type(quaternion) :: q

      q%r = 0.0d0
      q%i = 0.0d0
      q%j = 0.0d0
      q%k = 0.0d0
        
    end subroutine qsetzero


    subroutine qset(r, i, j, k, q)

      implicit none
      double precision, intent(in) :: r, i, j, k
      type(quaternion) :: q

      q%r = r
      q%i = i
      q%j = j
      q%k = k
        
    end subroutine qset


    function CrossProduct3D (X,Y) result(Z)
    IMPLICIT NONE
      REAL(8),INTENT(IN)  :: X(3), Y(3)
      REAL(8) :: Z(3)
      Z=0.0_8
      Z(1) = X(2)*Y(3)-X(3)*Y(2)
      Z(2) = X(3)*Y(1)-X(1)*Y(3)
      Z(3) = X(1)*Y(2)-X(2)*Y(1)
    END function CrossProduct3D

    function Mat_Cross(P) result(P_m)
    IMPLICIT NONE
      REAL(8),INTENT(IN)  :: P(3)! Y(3)
      REAL(8) :: P_m(3,3)
      P_m=0.0_8
      P_m(1,2) = -P(3)
      P_m(1,3) =  P(2)
      P_m(2,1) =  P(3)
      P_m(2,3) =  -P(1)
      P_m(3,1) =  -P(2)
      P_m(3,2) =  P(1)
    END function Mat_Cross

    function Mat_inv_Cross(P) result(P_m)
    IMPLICIT NONE
      REAL(8),INTENT(IN)  :: P(3)! Y(3)
      REAL(8) :: P_m(3,3)
      P_m=transpose(Mat_Cross(P))

    END function Mat_inv_Cross



    function Lie_algebra_update(u,q0)  result(q_norm)
    type(quaternion), intent(in) :: q0
    real*8,intent(in)::u(3)
    type(quaternion) :: q1,qu,q_norm
    !real*8, intent(out)::
    real*8 :: ucos,usin,unorm

    unorm=sqrt(sum(u**2))+1e-16
    ucos=cos(unorm*0.5_8)
    !if(unorm.lt.1e-5) then
     ! usin=0.0_8
    !else
    usin=sin(unorm*0.5_8)/unorm
    !endif
    qu%r=ucos
    qu%i=usin*u(1);qu%j=usin*u(2);qu%k=usin*u(3)
    q1=qu*q0
    q_norm=qscal(1.0_8/qnorm(q1) , q1)

    end function Lie_algebra_update


    function dexpu_inv(u,w) result(dexpu)
    real*8, intent(in) :: u(3),w(3)
    real*8::dexpu(3)

    !type(quaternion) :: q1,qu
    !real*8, intent(out)::
    real*8 :: ucot,uscalar,unorm,uw(3)

    unorm=sqrt(sum(u**2))

    if(unorm.lt.1e-16) then
      uscalar=-1.0_8/12.0_8
    else
      ucot=unorm/tan(unorm*0.5_8)-2.0_8
      uscalar=0.5_8*ucot/(unorm*unorm)
    endif
    uw(:)= CrossProduct3D(u,CrossProduct3D(u,w))
    dexpu(:)=w-0.5*CrossProduct3D(u,w)-uscalar*uw
    end function dexpu_inv

    end module quaternions  




!****************************************************************

      MODULE rb_conglomerate
      use quaternions     
      implicit none
      INTEGER K_rb
      INTEGER,ALLOCATABLE :: KK_rbmconn(:)
      REAL*8, ALLOCATABLE :: conf_rb(:,:),conf_rb_vector(:,:),U_par_rb(:)
      REAL*8, ALLOCATABLE ::rbmconn_Inertial_body(:,:),rbmconn_Inertial_body_inverse(:,:)
      type(quaternion), ALLOCATABLE:: q_rb(:)
      !LOGICAL useCongl
      END MODULE rb_conglomerate   


!****************************************************************
      MODULE filament
      use quaternions     
      implicit none
      INTEGER F_rb
      INTEGER,ALLOCATABLE :: Filament_num(:),index1(:)
      REAL*8, ALLOCATABLE :: Filament_conf_past(:,:),Filament_conf_now(:,:)
      REAL*8, ALLOCATABLE :: Filament_tau_base(:,:),Filament_tau_now(:,:)
      REAL*8, ALLOCATABLE :: Filament_U1_now(:),Filament_X1_now(:),Filament_X1_past(:)
      REAL*8, ALLOCATABLE :: Filament_Interal_force(:)!,Filament_Inertial_torque(:)
      REAL*8, ALLOCATABLE :: Filament_Lie_algebra_now(:)
      type(quaternion), ALLOCATABLE:: Filament_q(:)
      logical::solve_implicit
      REAL*8 ::Filament_Inertial_body_inverse(3,3),Filament_Inertial_body(3,3)
      REAL*8 ::GI,EI,GA,EA
      
      !real*8 alphaX,betaY,gamaZ 
      !LOGICAL useCongl
      END MODULE filament   
!****************************************************************
  module dimentionless
  use SYS_property
  use DLVO_property
  use BROWN_proerty
  use lambda_Dim
  use CONFIG,only:radii_true,radii_test
  use TENSORS,only:PAI,pai2
  use control,only: END_TIME,BASE_DT,T0,dt_dem0,NtDEM
  use filament,only:GI,EI,GA,EA
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
      lambda_Viscosity=1.0_8/PAI/mu_f
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
      Shear_module=Youngs_module/(2.0_8*(1.0_8+Poisson_ratio))
      Surfaceenergy=Surfaceenergy*lambda_Surfaceenergy
      F_adh=1.50_8*PAI*Surfaceenergy*1.0_8
      b0=((9.0_8*PAI*Surfaceenergy*1.0_8*1.0_8*(1.0_8-Poisson_ratio*Poisson_ratio)) &
          & /(0.5_8*Youngs_module))**(1.0_8/3.0_8)
      hp0=b0*b0/(2**(1.0_8/3.0_8)*3.0_8*1.0_8)
      heq=hp0-(1.1_8)**(-3.0_8/5.0_8)*b0*b0/1.0_8
      kBT=kBT*lambda_kB*lambda_temperature
      k_LJ=30.0_8*kBT/1.0_8/1.0_8

      GA=pai*Shear_module
      EA=pai*Youngs_module
      GI=0.25*Shear_module
      EI=0.25*Youngs_module
      !kB=kB*lambda_kB


      end SUBROUTINE dimentionless_translate

  end module dimentionless