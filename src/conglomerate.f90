!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*********************************

   Module conglomerate
   !use tensors,only:kdirac,EPS,pip5,IIiso,pai
   !use LATTICE_BASE,only:kapa
   use quaternions
   use method,only:FTS_method,simplePeriod
   USE CONFIG,only:u_bg,omega_bg,omegaT,Eij   ! CONF,POLY_LEN
   use hydro_tools,only:calc_uo_bg
   use period_bdy_tools,only:PERIOD_CORRECTION
   use SYS_property,only:mass_par,DT_DEM
   USE LATTICE_BASE,only:LB
   use rb_conglomerate,only:K_rb,KK_rbmconn,conf_rb_vector,rbmconn_Inertial_body,rbmconn_Inertial_body_inverse

   IMPLICIT NONE
   private

   public:: new_vel,new_U_par_rb,new_conf_swim,rbmconn_Init,new_q_rb_swim!,swim_rbmconn,swim_uo_bg,

   contains

!****************************************************************


    subroutine new_conf_swim(NN,conf,T,U_par_rb,conf_rb,q_rb)
    IMPLICIT NONE
    INTEGER,intent(in)::NN
    REAL*8,intent(in):: T!,U_pos(6*NN)!,RADII(NN)
    real*8,intent(in):: U_par_rb(6*K_rb)
    real*8,intent(inout):: conf_rb(3,K_rb)
    type(quaternion),intent(inout)::q_rb(K_rb)
    real*8,intent(out):: CONF(3,NN)


    integer :: ii, jj!, jj2, kk, ll, mm,num_rb
    !integer ::n6,n5!, idostep
    !real*8, dimension( 3 ) :: swimcm
    real*8 ::DCONF(3*K_rb),w(3)
    type(quaternion) ::dq,q_norm


    !conf_rb(:,:),conf_rb_vector(:,:),U_par_rb(:)

      DCONF=DT_DEM*U_par_rb(1:3*K_rb)

      forall(ii=1:K_rb,jj=1:3)
        conf_rb(jj,ii)=conf_rb(jj,ii)+DCONF(3*(ii-1)+jj)
      end forall

      if(simplePeriod) then
         DO ii = 1, K_rb
          CALL PERIOD_CORRECTION(conf_rb(:,ii),LB,T)
         ENDDO
      endif

      call conf_rb2conf(NN,conf_rb,conf,q_rb)

      do ii=1,K_rb
        jj=3*K_rb+3*(ii-1)
        w=U_par_rb(jj+1:jj+3)
        !dq=qderivat(q_rb(ii),w)*DT_DEM
        dq=qscal(DT_DEM,qderivat(q_rb(ii),w))
        q_norm=qadd(q_rb(ii),dq)
        q_rb(ii)=qscal(1.0_8/qnorm(q_norm) , q_norm)
        !write(*,*) "ii==",ii
        CALL qprint(q_rb(ii))
      enddo
      
        write(*,*) "rigid3==================================="
    end subroutine new_conf_swim


    subroutine conf_rb2conf(NN,conf_rb,conf,q_rb)
    IMPLICIT NONE
    INTEGER,intent(in)::NN
    REAL*8,intent(in):: conf_rb(3,K_rb)
    type(quaternion),intent(in)::q_rb(K_rb)
    real*8,intent(out):: CONF(3,NN)

    integer :: ii, jj, kk, KK_rb
    integer ::num_rb_sum,num_rb
    real*8 :: swimcm(3),swimcm_rotation(3,3)
    real*8 :: dx1, dx2, dx3

    num_rb_sum=0
    do  KK_rb=1,K_rb
        num_rb=KK_rbmconn(KK_rb)
        swimcm_rotation=quat2mat(q_rb(KK_rb))
       ! do jj=1,3
       ! write(*,*) "KK_rb and swimcm_rotation(jj,:)",KK_rb,jj,swimcm_rotation(jj,:)
       ! enddo
        do ii = num_rb_sum+1, num_rb_sum+num_rb     
            swimcm( : ) = matmul(swimcm_rotation,conf_rb_vector(:,ii))
       !     write(*,*) "conf_rb_vector(:,ii)=",ii,conf_rb_vector(:,ii)
         !   write(*,*) "swimcm(:,ii)=",ii,swimcm(:)
            conf(:,ii)=conf_rb(:,KK_rb)+swimcm(:)
        end do  
        num_rb_sum=num_rb_sum+num_rb    
    end do
    end subroutine conf_rb2conf


    subroutine new_q_rb_swim(NN,conf,U_pos,T,U_par_rb,conf_rb,q_rb)
    IMPLICIT NONE
    INTEGER,intent(in)::NN
    REAL*8,intent(in):: T,U_pos(6*NN)!,RADII(NN)
    real*8,intent(in):: U_par_rb(6*K_rb)
    real*8,intent(inout):: conf_rb(3,K_rb)
    type(quaternion),intent(inout)::q_rb(K_rb)
    real*8,intent(inout):: CONF(3,NN)


    integer :: ii, jj!, jj2, kk, ll, mm,num_rb
    !integer ::n6,n5!, idostep
    !real*8, dimension( 3 ) :: swimcm
    real*8 ::DCONF(3*K_rb),w(3),DCONF_1(3*NN)
    type(quaternion) ::dq


    !conf_rb(:,:),conf_rb_vector(:,:),U_par_rb(:)
        
      DCONF_1=DT_DEM*U_pos(1:3*NN)

      forall(ii=1:NN,jj=1:3)
        CONF(jj,ii)=CONF(jj,ii)+DCONF_1(3*(ii-1)+jj)
      end forall

      DCONF=DT_DEM*U_par_rb(1:3*K_rb)

      forall(ii=1:K_rb,jj=1:3)
        conf_rb(jj,ii)=conf_rb(jj,ii)+DCONF(3*(ii-1)+jj)
      end forall

      if(simplePeriod) then
         DO ii = 1, K_rb
          CALL PERIOD_CORRECTION(conf_rb(:,ii),LB,T)
         ENDDO
      endif

      !call conf_rb2conf(NN,conf_rb,conf,q_rb)

      do ii=1,K_rb
        jj=3*K_rb+3*(ii-1)
        w=U_par_rb(jj+1:jj+3)
        !dq=qderivat(q_rb(ii),w)*DT_DEM
        dq=qscal(DT_DEM,qderivat(q_rb(ii),w))
        q_rb(ii)=q_rb(ii)+dq
        !write(*,*) "ii==",ii
        CALL qprint(q_rb(ii))
      enddo
      
        write(*,*) "rigid3==================================="
    end subroutine new_q_rb_swim






    subroutine new_U_par_rb(NN,CONF,RADII,Ftotal,U_pos,U_par_rb,q_rb)
    IMPLICIT NONE
    INTEGER,intent(in)::NN
    REAL*8,intent(in):: CONF(3,NN),RADII(NN),Ftotal(6*NN)!
    type(quaternion),intent(in)::q_rb(K_rb)
    real*8,intent(inout):: U_par_rb(6*K_rb)
    real*8,intent(out):: U_pos(6*NN)

    integer :: ii, jj, jj2, kk, ll, mm,num_rb
    !integer ::n6,n5!, idostep
    !real*8, dimension( 3 ) :: swimcm
    real*8, dimension( 6*NN, 6*K_rb ) :: rbmconn
    real*8, dimension( 3*K_rb, 3*K_rb ) :: rbmconn_Inertial,rbmconn_Inertial_inverse
    real*8, dimension( 6*K_rb ) :: swimhold_Ftotal,Ftotal_rb!,swimhold2
    real*8, dimension( 6*K_rb,6*K_rb ) :: mass_inverse_matrix
    real*8, dimension( 3*K_rb ) :: rbmconn_torque


    mass_inverse_matrix=0.0_8
    do ii=1,K_rb  
        num_rb=KK_rbmconn(ii)
        jj=3*(ii-1)
        mass_inverse_matrix(jj+1,jj+1)=1.0_8/num_rb
        mass_inverse_matrix(jj+2,jj+2)=1.0_8/num_rb
        mass_inverse_matrix(jj+3,jj+3)=1.0_8/num_rb
    enddo
    CALL swim_rbmconn_Inertial(q_rb,rbmconn_Inertial,rbmconn_Inertial_inverse)
    CALL swim_rbmconn_torque(rbmconn_Inertial,U_par_rb,rbmconn_torque)

    CALL swim_rbmconn(NN,CONF,rbmconn,q_rb)
    swimhold_Ftotal=matmul(transpose(rbmconn),Ftotal)
    swimhold_Ftotal(3*K_rb+1:6*K_rb)=swimhold_Ftotal(3*K_rb+1:6*K_rb)-rbmconn_torque

   !CALL MATREV(rbmconn_Inertial,3*K_rb)
    mass_inverse_matrix(3*K_rb+1:6*K_rb,3*K_rb+1:6*K_rb)=rbmconn_Inertial_inverse
    Ftotal_rb=matmul(mass_inverse_matrix,swimhold_Ftotal)

    U_par_rb=U_par_rb+Ftotal_rb*DT_DEM/mass_par
    U_pos=matmul(rbmconn,U_par_rb)


    end subroutine new_U_par_rb




    subroutine new_vel(NN,CONF,rfu,rfe,fext,einf,u,U_par_rb,q_rb)
    ! ------------------------------------------------------------ !
    ! This routine computes the contributions to the particle
    ! velocities due to a swimming gait and an external force/torque.
    !
    ! Inputs: none
    ! Outputs: none
    ! Changes: u - the velocities due to the swimming gait 
    !          sh - the hydrodynamic contribution to the stresslet
    ! ------------------------------------------------------------ !
      IMPLICIT NONE
      INTEGER,intent(in)::NN
      REAL*8,intent(in):: CONF(3,NN),rfu(6*NN,6*NN),rfe(6*NN,5*NN),fext(6*NN),einf(5*NN)
      type(quaternion),intent(in)::q_rb(K_rb)
      real*8,intent(out):: u(6*NN),U_par_rb(6*K_rb)

        integer :: ii, jj, jj2, kk, ll, mm
        !integer ::n6,n5!, idostep
        !real*8, dimension( 3 ) :: swimcm
        real*8, dimension( 6*NN, 6*K_rb ) :: rbmconn
        real*8, dimension( 6*K_rb, 6*K_rb ) :: rfu_swim
        real*8, dimension( 6*NN, 6*K_rb ) :: swimhold
        real*8, dimension( 6*K_rb ) :: fswim
        !real*8, dimension( 5*NN ) :: einf
        real*8, dimension( 6*NN ) :: fconst,uext,udiff,ugate, fgate
        !real*8, dimension( 3, 3 ) :: ee, xf
        !real*8 :: omega_gate

        CALL swim_rbmconn(NN,CONF,rbmconn,q_rb)
        !do ii=1,6*NN
         !   write(*,*) rbmconn(ii,:)
        !enddo
        !CALL swim_rbmconn(NN,CONF,rbmconn)

        swimhold=matmul(rfu,rbmconn)
        !call dgemm( 't', 'n', 6, 6, n6, 1.d0, &
        !    rbmconn, n6, swimhold, n6, 0.d0, rfu_swim, 6 )
        rfu_swim=matmul(transpose(rbmconn),swimhold)

        !rfu_swim=matmul(transpose(rbmconn),matmul(rfu,rbmconn))
        !write(*,*) 'swim_rbmconnswim_rbmconnswim_rbmconn'
        CALL MATREV(rfu_swim,6*K_rb)
        !write(*,*) 'swim_rbmconnswim_rbmconnswim_rbmconn'
        !idostep = 3
        !call cholesky( rfu_swim, t, 6, istop, idostep ) 
        !if ( istop .eq. 0 ) then
            !write( *, * ) 'Stopped after chol. 1 in vel.f90'
           ! stop
        !end if
    !!!!
    ! Beyond this point, all the resistance tensors and their inverses have been constructed
    ! one simply need compute the appropriate matrix-vector products (this is done with dgemv,
    ! from the blas package though one can hard-code this as well) to study the motion of
    ! colloids/rigid bodies/swimmers.  Note, the rate of strain (einf) has been set to zero.
    ! to study swimmers, this must be changed as suggested in the tutorial.
    ! At this point:
    ! rfu_swim = ( sigma * rfu * sigma^T ) ^ -1
    ! The implicit gate
    ! einf( : ) = 0.d0
    ! The explicit gate
        ugate( : ) = 0.d0 !tuoluo bianmao velocity
        fgate( : ) = 0.d0       
#ifdef ugate        
   ! omega_gate=1.0_8
    do ii=1,NN
       ugate(3*(ii-1)+1)=-omega_gate*conf(2,ii)
       ugate(3*(ii-1)+2)=omega_gate*conf(1,ii)
       ugate(3*NN+3*(ii-1)+3)=omega_gate

    enddo
   fgate=-matmul(rfu,ugate)
#endif
        ! Begin calculating the velocities of the particles...
        
        ! Calculate the propulsive thrust due to an explicit gate
        !call dgemv( 'n', n6, n6, -1.d0, rfu, n6, ugate, 1, 0.d0, fgate, 1 )
        ! Calculate the propulsive thrust due to an implicit gate
        !call dgemv( 'n', n6, n5, 1.d0, rfe, n6, einf, 1, 0.d0, uext, 1 )
        ! Add in any external forces(6*NN,1)
        if(FTS_method)then
           uext=matmul(rfe,einf)
        else
           uext = 0.d0
        endif

        uext( : ) =  uext( : )+ fext( : ) +fgate( : )
        ! Compute Sigma * u and store in fswim(6,1)
        !call dgemv( 't', n6, 6, 1.d0, rbmconn, n6, uext, 1, 0.d0, fswim, 1 )
        !fswim( : ) = 0.d0
        fswim=matmul(transpose(rbmconn),uext)
       ! uswim( : ) = 0.d0
        ! Compute rfu_swim(6,6) * fswim(6,1)
        !call dgemv( 'n', 6, 6, 1.d0, rfu_swim, 6, fswim, 1, 0.d0, uswim, 1 )
        U_par_rb=matmul(rfu_swim,fswim)
        ! Compute Sigma^T * uswim and store in udiff 
        !call dgemv( 't', n6, 6, 1.d0, rbmconn, n6, uswim, 1, 0.d0, udiff, 1 )
        !udiff(:) = 0.d0
        udiff=matmul(rbmconn,U_par_rb)
        u(:)=udiff( : ) + ugate( : )
    !******************************************************************************!
    end subroutine new_vel




    subroutine swim_rbmconn(NN,CONF,rbmconn,q_rb)
    ! ------------------------------------------------------------ !
    ! This routine computes the contributions to the particle
    ! velocities due to a swimming gait and an external force/torque.
    !
    ! Inputs: none
    ! Outputs: none
    ! Changes: u - the velocities due to the swimming gait 
    !          sh - the hydrodynamic contribution to the stresslet
    ! ------------------------------------------------------------ !
      IMPLICIT NONE
      INTEGER,intent(in)::NN
      REAL*8,intent(in):: CONF(3,NN)
      type(quaternion),intent(in)::q_rb(K_rb)
      real*8,intent(out):: rbmconn(6*NN, 6*K_rb)

      integer :: ii, jj, kk, KK_rb
      integer ::num_rb_sum,num_rb
      real*8 :: swimcm(3),swimcm_rotation(3,3)
      real*8 :: dx1, dx2, dx3

        ! Find the center of mass of each swimmer
        rbmconn( :, : ) = 0.d0
        num_rb_sum=0
        do  KK_rb=1,K_rb
            num_rb=KK_rbmconn(KK_rb)
            swimcm_rotation=quat2mat(q_rb(KK_rb))
            do ii = num_rb_sum+1, num_rb_sum+num_rb     
                swimcm( : ) = matmul(swimcm_rotation,conf_rb_vector(:,ii))               
                dx1 = swimcm( 1 )
                dx2 = swimcm( 2 )
                dx3 = swimcm( 3 )
                jj = 3 * ( ii - 1 ) 
                do kk = 1, 3
                    rbmconn( jj+kk,      3*(KK_rb-1)+kk ) = 1.d0
                    rbmconn( 3*NN+jj+kk, 3*(K_rb+KK_rb-1)+kk ) = 1.d0
                end do
                rbmconn( jj +1, 3*(K_rb+KK_rb-1) + 2 ) = -dx3
                rbmconn( jj +1, 3*(K_rb+KK_rb-1) + 3 ) = dx2
                rbmconn( jj +2, 3*(K_rb+KK_rb-1) + 1 ) = dx3
                rbmconn( jj +2, 3*(K_rb+KK_rb-1) + 3 ) = -dx1
                rbmconn( jj +3, 3*(K_rb+KK_rb-1) + 1 ) = -dx2
                rbmconn( jj +3, 3*(K_rb+KK_rb-1) + 2 ) = dx1
            end do  
            num_rb_sum=num_rb_sum+num_rb    
        end do

        ! Sigma^T = rbmconn(6*NN,6)
        ! Calculate rfu_swim ( Sigma * RFU * Sigma^T)-->(6,6)

        !call dgemm( 'n', 'n', n6, 6, n6, 1.d0, rfu, n6, rbmconn, &
        !    n6, 0.d0, swimhold, n6 )
        !swimhold=matmul(rfu,rbmconn)
        !call dgemm( 't', 'n', 6, 6, n6, 1.d0, &
        !    rbmconn, n6, swimhold, n6, 0.d0, rfu_swim, 6 )
        !rfu_swim=matmul(transpose(rbmconn),matmul(rfu,rbmconn))

        ! Invert rfu_swim(6,6)
        

    end subroutine swim_rbmconn

    subroutine swim_rbmconn_Inertial(q_rb,rbmconn_Inertial,rbmconn_Inertial_inverse)
    ! ------------------------------------------------------------ !
    ! This routine computes the contributions to the particle
    ! velocities due to a swimming gait and an external force/torque.
    !
    ! Inputs: none
    ! Outputs: none
    ! Changes: u - the velocities due to the swimming gait 
    !          sh - the hydrodynamic contribution to the stresslet
    ! ------------------------------------------------------------ !
      IMPLICIT NONE
      !INTEGER,intent(in)::NN
      !REAL*8,intent(in):: CONF(3,NN),RADII(NN)
      type(quaternion),intent(in)::q_rb(K_rb)
      real*8,intent(out):: rbmconn_Inertial(3*K_rb, 3*K_rb),rbmconn_Inertial_inverse(3*K_rb, 3*K_rb)

      integer :: ii, jj, KK_rb
      integer ::num_rb_sum,num_rb
      real*8, dimension( 3 ) :: dr!swimcm
      real*8, dimension( 3,3 ) :: swimcm_rotation,Inertial,Inertial_component!,Inertial_component_inverse
      real*8 :: dx, dy, dz, Ixx, Ixy,Ixz,Iyy,Iyz,Izz,RADII2

        ! Find the center of mass of each swimmer
        rbmconn_Inertial( :, : ) = 0.d0
        rbmconn_Inertial_inverse( :, : ) = 0.d0
        do  KK_rb=1,K_rb
            swimcm_rotation=quat2mat(q_rb(KK_rb))
            jj=3*(KK_rb-1)
            Inertial_component=rbmconn_Inertial_body(jj +1:jj+3, jj +1:jj+3)
            Inertial=matmul(swimcm_rotation,Inertial_component)
            rbmconn_Inertial( jj +1:jj+3, jj +1:jj+3 ) =matmul(Inertial,transpose(swimcm_rotation))
            Inertial_component=rbmconn_Inertial_body_inverse(jj +1:jj+3, jj +1:jj+3)
            Inertial=matmul(swimcm_rotation,Inertial_component)
            rbmconn_Inertial_inverse( jj +1:jj+3, jj +1:jj+3 ) =matmul(Inertial,transpose(swimcm_rotation))
        end do
    end subroutine swim_rbmconn_Inertial
 


    subroutine swim_rbmconn_torque(rbmconn_Inertial,U_par_rb,rbmconn_torque)
      IMPLICIT NONE
      !INTEGER,intent(in)::NN
      REAL*8,intent(in):: rbmconn_Inertial(3*K_rb,3*K_rb),U_par_rb(6*K_rb)!,RADII(NN)
      real*8,intent(out):: rbmconn_torque(3*K_rb)

      integer :: ii, jj, kk,KK_rb
      integer ::num_rb_sum,num_rb
      real*8, dimension( 3 ) :: swimcm,angular_velocity,angular_Inertial
      real*8, dimension( 3 ,3) :: swimcm_Inertial
      !real*8 ::angular_velocity(3)
      !real*8 :: dx, dy, dz, Ixx, Ixy,Ixz,Iyy,Iyz,Izz,RADII2

        ! Find the center of mass of each swimmer



        rbmconn_torque= 0.0_8
        angular_velocity=0.0_8
        num_rb_sum=0
        do  KK_rb=1,K_rb
            jj=3*(KK_rb-1)
            kk=3*K_rb+3*(KK_rb-1)
            angular_velocity=U_par_rb(kk+1:kk+3)
             
            swimcm_Inertial=rbmconn_Inertial(jj+1:jj+3,jj+1:jj+3)*mass_par
            angular_Inertial=matmul(swimcm_Inertial,angular_velocity)
            CALL CrossProduct3D (angular_velocity, angular_Inertial,swimcm)
            rbmconn_torque(jj+1:jj+3)=swimcm
        end do

    end subroutine swim_rbmconn_torque
 







    subroutine lub_drag_cor(dR,q)
    ! ------------------------------------------------------------ !
    ! This routine computes the contributions to the particle
    ! velocities due to a swimming gait and an external force/torque.
    !
    ! Inputs: none
    ! Outputs: none
    ! Changes: u - the velocities due to the swimming gait 
    !          sh - the hydrodynamic contribution to the stresslet
    ! ------------------------------------------------------------ !
    !use global

      !use Ewald_summation
      !USE TENSORS
      !USE LATTICE_BASE
      IMPLICIT NONE
      !INTEGER,intent(in)::NN
      REAL*8,intent(in):: dR(3)
      real*8,intent(out):: q(12,12)

        integer :: ii!, jj,jj2, kk, ll, mm
        real*8:: rbmconn(3,3),ONEI(3,3)
        real*8 :: dx1, dx2, dx3



        ! Find the center of mass of each swimmer
        ONEI( :, : ) = 0.d0
        do ii = 1, 3
            ONEI( ii,ii ) = 1.0_8
        end do
        ! Assemble the rigid body motion connectivity tensor (Sigma)
        rbmconn( :, : ) = 0.d0
        dx1 = dR(1)
        dx2 = dR(2)
        dx3 = dR(3)
        rbmconn( 1, 2 ) = -dx3
        rbmconn( 1, 3 ) = dx2
        rbmconn( 2, 1 ) = dx3
        rbmconn( 2, 3 ) = -dx1
        rbmconn( 3, 1 ) = -dx2
        rbmconn( 3, 2 ) = dx1

        q(:,:)=0.0_8
        q(1:3,1:3)=ONEI
        q(1:3,4:6)=-ONEI
        q(4:6,1:3)=-ONEI
        q(4:6,4:6)=ONEI
        q(7:9,7:9)=ONEI
        q(7:9,10:12)=-ONEI
        q(10:12,7:9)=-ONEI
        q(10:12,10:12)=ONEI
        q(1:3,7:9)=-rbmconn
        q(1:3,10:12)=-rbmconn
        q(4:6,7:9)=rbmconn
        q(4:6,10:12)=rbmconn
    end subroutine lub_drag_cor


    subroutine rbmconn_Init(NN,CONF,RADII,U_pos,conf_rb,U_par_rb,q_rb)
    IMPLICIT NONE
    INTEGER,intent(in)::NN
    REAL*8,intent(in):: CONF(3,NN),RADII(NN)
    real*8,intent(out):: U_pos(6*NN),conf_rb(3,K_rb),U_par_rb(6*K_rb)
    type(quaternion),intent(out):: q_rb(K_rb)

    integer :: ii, jj, kk, KK_rb
    integer ::num_rb_sum,num_rb
    real*8, dimension( 3 ) :: swimcm
    real*8 :: dx, dy, dz, Ixx, Ixy,Ixz,Iyy,Iyz,Izz,RADII2
    real*8 ::rbmconn_Inertial(3,3),rbmconn(6*NN, 6*K_rb)

      rbmconn_Inertial_body=0.0_8
      rbmconn_Inertial_body_inverse=0.0_8
      num_rb_sum=0
      do  KK_rb=1,K_rb       
          num_rb=KK_rbmconn(KK_rb)       
          swimcm( : ) = 0.d0
          do ii = num_rb_sum+1, num_rb_sum+num_rb
              swimcm( 1:3 ) = swimcm( 1:3 ) + conf(1:3,ii)
          end do
          swimcm = swimcm / num_rb
          conf_rb(:,KK_rb)=swimcm(:)

          Ixx=0.0_8
          Ixy=0.0_8
          Ixz=0.0_8
          Iyy=0.0_8
          Iyz=0.0_8
          Izz=0.0_8 
          do ii = num_rb_sum+1, num_rb_sum+num_rb 
              conf_rb_vector(:,ii)= conf( :,ii )- swimcm(:)
              dx = conf_rb_vector(1,ii)
              dy = conf_rb_vector(2,ii)
              dz = conf_rb_vector(3,ii)
              RADII2=0.4_8*RADII(ii)*RADII(ii)
              Ixx=Ixx+dy*dy+dz*dz+RADII2
              Iyy=Iyy+dx*dx+dz*dz+RADII2
              Izz=Izz+dx*dx+dy*dy+RADII2
              Ixy=Ixy-dx*dy
              Ixz=Ixz-dx*dz
              Iyz=Iyz-dy*dz 
          end do  
          rbmconn_Inertial=0.0_8   
          rbmconn_Inertial( 1, 1 ) =Ixx 
          rbmconn_Inertial( 2, 2 ) =Iyy
          rbmconn_Inertial( 3, 3 ) =Izz 
          rbmconn_Inertial( 1, 2 ) =Ixy
          rbmconn_Inertial( 2, 1 ) =Ixy  
          rbmconn_Inertial( 1, 3 ) =Ixz
          rbmconn_Inertial( 3, 1 ) =Ixz
          rbmconn_Inertial( 2, 3 ) =Iyz
          rbmconn_Inertial( 3, 2 ) =Iyz
          jj=3*(KK_rb-1)
          rbmconn_Inertial_body(jj+1:jj+3,jj+1:jj+3)=rbmconn_Inertial
          CALL MATREV(rbmconn_Inertial,3)
          rbmconn_Inertial_body_inverse(jj+1:jj+3,jj+1:jj+3)=rbmconn_Inertial
          num_rb_sum=num_rb_sum+num_rb
          q_rb(KK_rb)=quat(1.0_8,0.0_8,0.0_8,0.0_8)
      end do

      do ii=1,K_rb
        write(*,*) "CONF_rb(:,ii)=",ii,conf_rb(:,ii)
      enddo
      do ii=1,NN
        write(*,*) "conf_rb_vector(:,ii)=",ii,conf_rb_vector(:,ii)
      enddo


        CALL swim_rbmconn(NN,CONF,rbmconn,q_rb)
        U_par_rb=0.0_8
        U_pos=matmul(rbmconn,U_par_rb)
    end subroutine rbmconn_Init







end  Module conglomerate


















#ifdef dissipation
        ! Compute the constraining forces
        call dgemv( 'n', n6, n6, 1.d0, rfu, n6, udiff, 1, 0.d0, fconst, 1 )
        fconst( : ) = fconst( : ) + uext( : ) - 2.d0 * fext( : )
        sh( : ) = 0.d0
        ! Compute the hydrodynamic contribution to the stresslet
        call dgemv( 't', n6, n5, 1.d0, rfe, n6, udiff, 1, 0.d0, sh, 1 )
        call dgemv( 't', n6, n5, 1.d0, rfe, n6, ugate, 1, 1.d0, sh, 1 )
        call dgemv( 'n', n5, n5, -1.d0, rse, n5, einf, 1, 1.d0, sh, 1 )    
        dissipation = 0.d0
        do ii = 1, NN

            jj = 5 * ii - 4
            kk = 6 * ii - 5

            ee( 1, 1 ) = 2.d0 / 3.d0 * ( einf( jj ) - 0.5d0 * einf( jj + 4 ) )
            ee( 2, 2 ) = 2.d0 / 3.d0 * ( einf( jj + 4 ) - 0.5d0 * einf( jj ) )
            ee( 3, 3 ) = -ee( 1, 1 ) - ee( 2, 2 )
            ee( 1, 2 ) = 0.5d0 * einf( jj + 1 )
            ee( 2, 1 ) = 0.5d0 * einf( jj + 1 )
            ee( 1, 3 ) = 0.5d0 * einf( jj + 2 )
            ee( 3, 1 ) = 0.5d0 * einf( jj + 2 )
            ee( 2, 3 ) = 0.5d0 * einf( jj + 3 )
            ee( 3, 2 ) = 0.5d0 * einf( jj + 3 )

            xf( 1, 1 ) = -x( kk ) * fconst( kk ) + sh( jj )
            xf( 2, 2 ) = -x( kk + 1 ) * fconst( kk + 1 ) + sh( jj + 4 )
            xf( 3, 3 ) = -x( kk + 2 ) * fconst( kk + 2 ) - sh( jj ) - sh( jj + 4 )
            xf( 1, 2 ) = -x( kk ) * fconst( kk + 1 ) + sh( jj + 1 )
            xf( 2, 1 ) = -x( kk + 1 ) * fconst( kk ) + sh( jj + 1 )
            xf( 1, 3 ) = -x( kk ) * fconst( kk + 2 ) + sh( jj + 2 )
            xf( 3, 1 ) = -x( kk + 2 ) * fconst( kk ) + sh( jj + 2 )
            xf( 2, 3 ) = -x( kk + 1 ) * fconst( kk + 2 ) + sh( jj + 3 )
            xf( 3, 2 ) = -x( kk + 2 ) * fconst( kk + 1 ) + sh( jj + 3 )

            dissipation = dissipation + dot_product( ee( 1, : ), xf( 1, : ) ) + dot_product( ee( 2, : ), xf( 2, : ) ) &
                + dot_product( ee( 3, : ), xf( 3, : ) )
        end do
        x0( : ) = x( : )
        
     subroutine conf_rb_bg(NN,CONF,uo_bg)
        ! ------------------------------------------------------------ !
        ! This routine computes the contributions to the particle
        ! velocities due to a swimming gait and an external force/torque.
        !
        ! Inputs: none
        ! Outputs: none
        ! Changes: u - the velocities due to the swimming gait 
        !          sh - the hydrodynamic contribution to the stresslet
        ! ------------------------------------------------------------ !

          INTEGER,intent(in)::NN
          REAL*8,intent(in):: CONF(3,NN)
          real*8,intent(inout):: uo_bg(6*NN)

          real*8 rbmconn(6*NN, 6*K_rb)
          real*8 uo_rb_bg(6*K_rb),swimcm(6)
          INTEGER KK_rb,num_rb_sum,num_rb,ii
          !real*8 rbm( 6*NN, 6*K_rb )
          !uo_bgnew=0.0_8
          !write(*,*) 'K_rb===kkkkkkk',K_rb
          !write(*,*) 'NN===kkkkkkk',NN
          !write(*,*) 'CONF',CONF(:,:)

          CALL swim_rbmconn(NN,CONF,rbmconn)
          uo_rb_bg=0.0_8
          num_rb_sum=0
            do  KK_rb=1,K_rb
                
                num_rb=KK_rbmconn(KK_rb)
                
                swimcm= 0.d0

                do ii = num_rb_sum+1, num_rb_sum+num_rb
                    swimcm( 1:3 ) = swimcm( 1:3 ) + uo_bg(3*(ii-1)+1:3*ii)
                    swimcm( 4:6 ) = swimcm( 4:6 ) + uo_bg(3*NN+3*(ii-1)+1:3*NN+3*ii)
                end do

                uo_rb_bg( 3*(KK_rb-1)+1:3*KK_rb) = swimcm( 1:3 ) / num_rb
                uo_rb_bg( 3*K_rb+3*(KK_rb-1)+1:3*K_rb+3*KK_rb ) = swimcm( 4:6 ) / num_rb

                num_rb_sum=num_rb_sum+num_rb  
            enddo
            uo_bg=matmul(rbmconn,uo_rb_bg)


     end subroutine conf_rb_bg

       subroutine swim_uo_bg(NN,CONF,uo_bg)
        ! ------------------------------------------------------------ !
        ! This routine computes the contributions to the particle
        ! velocities due to a swimming gait and an external force/torque.
        !
        ! Inputs: none
        ! Outputs: none
        ! Changes: u - the velocities due to the swimming gait 
        !          sh - the hydrodynamic contribution to the stresslet
        ! ------------------------------------------------------------ !

          INTEGER,intent(in)::NN
          REAL*8,intent(in):: CONF(3,NN)
          real*8,intent(out):: uo_bg(6*NN)

          real*8 rbmconn(6*NN, 6*K_rb)
          real*8 uo_rb_bg(6*K_rb)!,CONF_rb(3,K_rb)
          !INTEGER KK_rb,num_rb_sum,num_rb,ii

          
          !CALL swim_conf_rb(NN,CONF,CONF_rb)

         !uo_rb_bg=0.0_8
         call calc_uo_bg(K_rb,conf_rb,uo_rb_bg)

         CALL swim_rbmconn(NN,CONF,rbmconn)
         
         uo_bg=matmul(rbmconn,uo_rb_bg)


         end subroutine swim_uo_bg

#endif