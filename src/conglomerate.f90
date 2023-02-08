!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*********************************

   Module conglomerate
   !use tensors,only:kdirac,EPS,pip5,IIiso,pai
   !use LATTICE_BASE,only:kapa
   use method,only:FTS_method
   USE CONFIG,only:u_bg,omega_bg,omegaT,Eij   ! CONF,POLY_LEN
   !use hydro_tools,only:calc_uo_bg
   use rb_conglomerate   

   IMPLICIT NONE
   private

   public:: new_vel,swim_uo_bg!,swim_rbmconn,

   contains

!****************************************************************

    subroutine new_vel(NN,CONF,rfu,rfe,fext,einf,u)
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
      real*8,intent(out):: u(6*NN)

        integer :: ii, jj, jj2, kk, ll, mm
        !integer ::n6,n5!, idostep
        !real*8, dimension( 3 ) :: swimcm
        real*8, dimension( 6*NN, 6*K_rb ) :: rbmconn
        real*8, dimension( 6*K_rb, 6*K_rb ) :: rfu_swim
        real*8, dimension( 6*NN, 6*K_rb ) :: swimhold
        real*8, dimension( 6*K_rb ) :: uswim, fswim
        !real*8, dimension( 5*NN ) :: einf
        real*8, dimension( 6*NN ) :: fconst,uext,udiff,ugate, fgate
        !real*8, dimension( 3, 3 ) :: ee, xf
        !real*8 :: omega_gate

        CALL swim_rbmconn(NN,CONF,rbmconn)
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
         !if(FTS_method)then
          ! uext=matmul(rfe,einf)
          !else
           uext = 0.d0
          !endif

        uext( : ) =  fext( : ) +uext( : )+fgate( : )
        ! Compute Sigma * u and store in fswim(6,1)
        !call dgemv( 't', n6, 6, 1.d0, rbmconn, n6, uext, 1, 0.d0, fswim, 1 )
        !fswim( : ) = 0.d0
        fswim=matmul(transpose(rbmconn),uext)
       ! uswim( : ) = 0.d0
        ! Compute rfu_swim(6,6) * fswim(6,1)
        !call dgemv( 'n', 6, 6, 1.d0, rfu_swim, 6, fswim, 1, 0.d0, uswim, 1 )
        uswim=matmul(rfu_swim,fswim)
        ! Compute Sigma^T * uswim and store in udiff 
        !call dgemv( 't', n6, 6, 1.d0, rbmconn, n6, uswim, 1, 0.d0, udiff, 1 )
        !udiff(:) = 0.d0
        udiff=matmul(rbmconn,uswim)
        u(:)=udiff( : ) !+ ugate( : )
    !******************************************************************************!
    end subroutine new_vel

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
          real*8 uo_rb_bg(6*K_rb),CONF_rb(3,K_rb)
          INTEGER KK_rb,num_rb_sum,num_rb,ii

          CALL swim_rbmconn(NN,CONF,rbmconn)
          CONF_rb=0.0_8
          uo_rb_bg=0.0_8
          num_rb_sum=0
            do  KK_rb=1,K_rb
                num_rb=KK_rbmconn(KK_rb)
                do ii = num_rb_sum+1, num_rb_sum+num_rb
                    CONF_rb( 1:3,KK_rb ) = CONF_rb( 1:3,KK_rb ) + CONF(1:3,ii)
                end do
                CONF_rb( 1:3,KK_rb ) = CONF_rb( 1:3,KK_rb )/ num_rb
                num_rb_sum=num_rb_sum+num_rb
            enddo
            call calc_uo_bg(K_rb,CONF_rb,uo_rb_bg)

            uo_bg=matmul(rbmconn,uo_rb_bg)


     end subroutine swim_uo_bg


    subroutine swim_rbmconn(NN,CONF,rbmconn)
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
      real*8,intent(out):: rbmconn(6*NN, 6*K_rb)

      integer :: ii, jj, kk, KK_rb
      integer ::num_rb_sum,num_rb
      real*8, dimension( 3 ) :: swimcm
      real*8 :: dx1, dx2, dx3

        ! Find the center of mass of each swimmer
        rbmconn( :, : ) = 0.d0
        num_rb_sum=0
        do  KK_rb=1,K_rb
            
            num_rb=KK_rbmconn(KK_rb)
            
            swimcm( : ) = 0.d0

            do ii = num_rb_sum+1, num_rb_sum+num_rb
                swimcm( 1:3 ) = swimcm( 1:3 ) + conf(1:3,ii)
            end do

            swimcm = swimcm / num_rb
         

            ! Assemble the rigid body motion connectivity tensor (Sigma)
            
            
            !do kk_rb=1,
                do ii = num_rb_sum+1, num_rb_sum+num_rb                    
                    dx1 = conf( 1,ii ) - swimcm( 1 )
                    dx2 = conf( 2,ii ) - swimcm( 2 )
                    dx3 = conf( 3,ii ) - swimcm( 3 )

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

      SUBROUTINE calc_uo_bg(MM,conf,uo_bg)
      IMPLICIT NONE
      integer,intent(in) ::MM
      real*8, intent(in) ::CONF(3,MM)
      real*8, intent(out)::uo_bg(6*MM)

      real*8 u_bg_local(3)
      INTEGER i,j
   
        do i=1,MM
          u_bg_local=0.0_8
          do j=1,3
          !U_bg_local=U_bg + Omega_bg * x + Einf_bg . x
          u_bg_local(j)=u_bg(j)+(OmegaT(j,1)*CONF(1,i)    &
            +OmegaT(j,2)*CONF(2,i)+OmegaT(j,3)*CONF(3,i)) &
            +(Eij(j,1)*CONF(1,i)+Eij(j,2)*CONF(2,i)       &
            +Eij(j,3)*CONF(3,i)) 
          !fill in the U-Uinf, Omega-Omegainf vector
          uo_bg(3*(i-1)+j)=u_bg_local(j)
          uo_bg(3*MM+3*(i-1)+j)=omega_bg(j)
          end do
        end do

      end SUBROUTINE calc_uo_bg



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

#endif