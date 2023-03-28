  Module filament_solve_Jacobian
  use quaternions
  use TENSORS,only:kDirac
  use method,only:FTS_method,simplePeriod
  ! USE CONFIG,only:u_bg,omega_bg,omegaT,Eij   ! CONF,POLY_LEN
  ! use hydro_tools,only:calc_uo_bg
  use period_bdy_tools,only:PERIOD_CORRECTION
  use SYS_property,only:mass_par,DT_DEM
  USE LATTICE_BASE,only:LB
  use filament,only:F_rb,Filament_num,Filament_tau_base,Filament_Inertial_body_inverse,index1,GI,EI,GA,EA

  IMPLICIT NONE
  private

  public:: filament_implicit_solve
 !,new_U_par_filament,new_conf_filament,new_q_filament

  contains

  subroutine filament_implicit_solve(Nfilament,X1_now,X1_past,U1_now,tau_now,Lie_algebra_now, &
                    & Internal_force,Fe,conf_now,conf_past,U_now,q_now)
  IMPLICIT NONE
  integer,intent(in):: Nfilament
  real*8,intent(in):: Fe(6*Nfilament)
  real*8,intent(inout)::X1_now(3*F_rb),X1_past(3*F_rb),U1_now(3*F_rb),Lie_algebra_now(3*Nfilament)
  real*8,intent(inout)::tau_now(3,Nfilament),conf_now(3,Nfilament),conf_past(3,Nfilament),U_now(6*Nfilament)
  type(quaternion),intent(inout):: q_now(Nfilament)
  real*8,intent(out)::Internal_force(6*Nfilament)

  integer::ii,iter
  real*8:: Internal_torque(3*Nfilament),conf_next(3,Nfilament),omega_next(3*Nfilament),error,Yk_next(6*Nfilament)
  real*8:: Lie_algebra_next(3*Nfilament),X1_next(3*Nfilament),Internal_lamda(3*Nfilament),tau_next(3,Nfilament)
  real*8:: Internal_lamda_test(3*Nfilament)
  type(quaternion):: q_next(Nfilament)

  iter=200;error=1e-6
  call Calc_Face_InternalTorques(Nfilament,conf_now,q_now,Internal_torque,Internal_lamda_test)
  call sub_To_Yk(Nfilament,Yk_next,X1_now,Internal_lamda_test,Lie_algebra_now)
  call BDF_solver(Nfilament,iter,error,Yk_next,X1_now,X1_past,U1_now,tau_now, &
                    & Lie_algebra_now,Internal_torque,Fe,conf_now,conf_past,U_now,q_now, &
                    & X1_next,Internal_lamda,Lie_algebra_next, &
                    & tau_next,conf_next,omega_next,q_next)
  X1_past=X1_now
  X1_now=X1_next
  U1_now=(X1_now-X1_past)/DT_DEM
  conf_past=conf_now
  conf_now=conf_next
  tau_now=tau_next
  q_now=q_next
  Lie_algebra_now=Lie_algebra_next
  do ii=1,Nfilament
    U_now(3*(ii-1)+1:3*ii)=(conf_now(:,ii)-conf_past(:,ii))/DT_DEM
  enddo
  U_now(3*Nfilament+1:6*Nfilament)=omega_next


  Internal_force(1:3*Nfilament)=Internal_lamda
  Internal_force(3*Nfilament+1:6*Nfilament)=Internal_torque

  end subroutine filament_implicit_solve




  subroutine Calc_Face_InternalTorques(Nfilament,conf_now,q_now,Internal_torque,Internal_force)
  !  Place gravity, elastic and constraint forces and torques 
  !  on the filament. They are stored in Filament.F
  IMPLICIT NONE
  INTEGER,intent(in)::Nfilament
  real*8,intent(in)::conf_now(3,Nfilament)!,Filament_tau_base(3,NN)
  type(quaternion),intent(in)::q_now(Nfilament)
  real*8,intent(out)::Internal_torque(3*Nfilament),Internal_force(3*Nfilament)
  integer :: ii, jj, KK_fila
  integer ::num_fila_sum,num_fila
  real*8 :: dr(3),ds,M(3),tauF(3),F0(3),F1(3)
  type(quaternion)::q_mid,dqds,qb
  real*8 :: torque(3),angle(3),force(3)


  Internal_torque=0.0_8
  Internal_force=0.0_8
  num_fila_sum=0
  do  KK_fila=1,F_rb       
      num_fila=Filament_num(KK_fila)       
      do ii = 1, num_fila-1
          jj=ii+num_fila_sum

          dr=conf_now(:,jj+1)-conf_now(:,jj)
          ds=SQRT(SUM((Filament_tau_base(:,jj))**2))
          q_mid = qMidpoint(q_now(jj),q_now(jj+1))
          tauF=dr/ds!(tauj1+tauj)/dL
          F0=qrot_inverse(q_mid,tauF)
          F0=F0-Filament_tau_base(:,jj)/ds
          F1=0.0_8
          F1(1)=F0(1)*EA;
          F1(2)=F0(2)*GA;F1(3)=F0(3)*GA
          force=qrot(q_mid,F1)
          ! Elastic torques
          !dqds =dqsub(Filament_q(jj+1),Filament_q(jj))/ds
          !qb = 2.0_8*qmul(qconj(q_mid),dqds)
          angle=qangle(q_now(jj+1),q_now(jj))/ds
          !M(1)=GI*qb%i;M(2)=EI*qb%j;M(3)=EI*qb%k
          M(1)=GI*angle(1);M(2)=EI*angle(2);M(3)=EI*angle(3)               
          torque = qrot(q_mid,M)


          write(*,*) 'MMMMMMMMMMMMMMM===========',M(:)!,GI,EI
          write(*,*) 'internal_torque===========',torque(:)

          Internal_torque( 3*(jj-1)+1:3*jj) =torque(:)
          Internal_force ( 3*(jj-1)+1:3*jj) =force(:)
      end do
      num_fila_sum=num_fila_sum+num_fila
  end do

    end subroutine Calc_Face_InternalTorques

  subroutine BDF_solver(Nfilament,iter,error,Yk_next,X1_now,X1_past,U1_now,tau_now,Lie_algebra_now, &
                    & Inertial_torque,Fe,conf_now,conf_past,U_now,q_now,X1_next,Inertial_lamda,Lie_algebra_next, &
                    tau_next,conf_next,omega_next,q_next)
  IMPLICIT NONE
  integer,intent(in):: Nfilament,iter
  real*8,intent(in)::X1_now(3*F_rb),X1_past(3*F_rb),U1_now(3*F_rb),Lie_algebra_now(3*Nfilament)
  real*8,intent(in)::Inertial_torque(3*Nfilament),error
  real*8,intent(in):: tau_now(3,Nfilament),Fe(6*Nfilament),conf_now(3,Nfilament),conf_past(3,Nfilament),U_now(6*Nfilament)
  type(quaternion),intent(in):: q_now(Nfilament)
  real*8,intent(out)::X1_next(3*F_rb),Inertial_lamda(3*Nfilament),Lie_algebra_next(3*Nfilament)
  real*8,intent(out)::tau_next(3,Nfilament),conf_next(3,Nfilament),omega_next(3*Nfilament)
  type(quaternion),intent(out):: q_next(Nfilament)
  real*8,intent(inout)::Yk_next(6*Nfilament)

  integer::ii,jj
  real*8:: dYk_norm,fy_norm
  real*8:: dYk(6*Nfilament),filament_Fy(6*Nfilament),Ck(6*Nfilament),Dk(6*Nfilament),jac_mat(6*Nfilament,6*Nfilament)
  real*8:: Jacobian_filament(6*Nfilament,6*Nfilament),CDk_sum(6*Nfilament,6*Nfilament),dk_mat(6*Nfilament,6*Nfilament)
  

  Ck=0.0_8;Dk=0.0_8;CDk_sum=0.0_8
  do ii=1,iter
  write(*,*) "iter==================",ii
    call Yk_To_sub(Nfilament,Yk_next,X1_next,Inertial_lamda,Lie_algebra_next)
    call Lie_algebra_To_q(Nfilament,q_now,Lie_algebra_next,q_next)
    call q_To_tau(Nfilament,q_next,tau_next)
    call filament_X1ToCONF(Nfilament,X1_next,tau_next,conf_next)
    call Rotate_omeganew(Nfilament,tau_now,Inertial_torque,Inertial_lamda,Fe(3*Nfilament+1:6*Nfilament), &
        & U_now(3*Nfilament+1:6*Nfilament),omega_next)

    call filament_function(Nfilament,X1_now,X1_past,U1_now,tau_now,Lie_algebra_now, &
                    & Lie_algebra_next,X1_next,Inertial_lamda, &
                    & Inertial_torque,Fe,conf_now,conf_past,conf_next,U_now,omega_next,filament_Fy)
    fy_norm=sqrt(sum(filament_Fy**2))/Nfilament
    if(fy_norm.lt.error)then
       write(*,*) "fy_norm===",fy_norm
       write(*,*) "iter complete successs!!!!!!!!!!!!!!!!!!!!!"
       return
    endif
    !if(mod(ii,5).eq.0)then
      call Jacobian_filament_function(Nfilament,Lie_algebra_next,tau_next,omega_next, &
              & tau_now,Jacobian_filament)
    !endif

    CDk_sum=CDk_sum+matrix_mul77(Ck,Dk,6*Nfilament)
    jac_mat=Jacobian_filament+CDk_sum

    call inverse_squarematrix(6*Nfilament,jac_mat,dk_mat)
    dYk=-matmul(dk_mat,filament_Fy)
    Yk_next=Yk_next+dYk
    dYk_norm=sqrt(sum(dYk**2))
    Ck=filament_Fy/dYk_norm;Dk=dYk/dYk_norm

  enddo

  write(*,*) "fy_norm===",fy_norm
  write(*,*) "iter failed!!!!!!!!!!!!!!!!!!!!!"

  end subroutine BDF_solver



  function matrix_mul77(a,b,bc)  result(c)   
  implicit none
  integer,intent(in)::bc
  real*8,intent(in)::a(bc),b(bc)
  real*8:: c(bc,bc)
  integer i,j
  do i=1,bc 
      do j=1,bc
          c(i,j)=a(i)*b(j)
      end do 
  end do
  return
  end

  subroutine Yk_To_sub(Nfilament,Yk_next,X1_next,Inertial_lamda,Lie_algebra_next)
  IMPLICIT NONE
  integer:: Nfilament
  real*8,intent(in)::Yk_next(6*Nfilament)
  real*8,intent(out):: X1_next(3*F_rb),Inertial_lamda(3*Nfilament),Lie_algebra_next(3*Nfilament)

  X1_next=Yk_next(1:3*F_rb)
  call exactToInertial_lamda(Nfilament,Inertial_lamda,Yk_next(3*F_rb+1:3*Nfilament))
  Lie_algebra_next=Yk_next(3*Nfilament+1:6*Nfilament)
  end subroutine Yk_To_sub

  subroutine sub_To_Yk(Nfilament,Yk_next,X1_next,Inertial_lamda,Lie_algebra_next)
  IMPLICIT NONE
  integer:: Nfilament
  real*8,intent(out)::Yk_next(6*Nfilament)
  real*8,intent(in):: X1_next(3*F_rb),Inertial_lamda(3*Nfilament),Lie_algebra_next(3*Nfilament)

  real*8:: exactInertial_lamda(3*(Nfilament-F_rb))

  Yk_next(1:3*F_rb)=X1_next
  !call exactToInertial_lamda(Nfilament,Inertial_lamda,Yk_next(3*F_rb+1:3*Nfilament))
  call Find_exactInertial_lamda(Nfilament,Inertial_lamda,exactInertial_lamda)
  Yk_next(3*F_rb+1:3*Nfilament)=exactInertial_lamda
  Yk_next(3*Nfilament+1:6*Nfilament)=Lie_algebra_next
  end subroutine sub_To_Yk


  subroutine Jacobian_filament_function(Nfilament,Lie_algebra_next,tau_next,omega_next, &
                    & tau_now,Jacobian_filament)
  IMPLICIT NONE
  integer:: Nfilament
  real*8,intent(in)::tau_next(3,Nfilament),omega_next(3*Nfilament)
  real*8,intent(in)::tau_now(3,Nfilament),Lie_algebra_next(3*Nfilament)
  real*8,intent(out)::Jacobian_filament(6*Nfilament,6*Nfilament)

  integer::ii,jj,kk,mm,nn,index1_internal,index2_internal
  real*8::Jacobian_11_filament(3*F_rb,3*F_rb),Jacobian_12_filament(3*F_rb,3*(Nfilament-F_rb))
  real*8::Jacobian_21_filament(3*(Nfilament-F_rb),3*F_rb),Jacobian_23_filament(3*(Nfilament-F_rb),3*Nfilament)
  real*8::Jacobian_22_filament(3*(Nfilament-F_rb),3*(Nfilament-F_rb))
  real*8::Jacobian_32_filament(3*Nfilament,3*(Nfilament-F_rb)),Jacobian_33_filament(3*Nfilament,3*Nfilament)
  real*8::Jacobian_sub(3,3)

  Jacobian_filament=0.0_8
  Jacobian_11_filament=0.0_8
  Jacobian_12_filament=0.0_8
  Jacobian_21_filament=0.0_8
  Jacobian_22_filament=0.0_8
  Jacobian_23_filament=0.0_8
  Jacobian_32_filament=0.0_8
  Jacobian_33_filament=0.0_8

  do ii=1,F_rb
    index1_internal=index1(ii)-ii+1
    Jacobian_11_filament(3*(ii-1)+1:3*(ii),3*(ii-1)+1:3*(ii))=kDirac
    jj=index1_internal
    Jacobian_12_filament(3*(ii-1)+1:3*(ii),3*(jj-1)+1:3*(jj))=-2.0_8/3.0_8*DT_DEM*DT_DEM/mass_par*kDirac
  enddo

  do jj=1,F_rb
    index1_internal=index1(jj)-jj+1
    index2_internal=index1_internal+Filament_num(jj)-2
    do ii=index1_internal,index2_internal
      Jacobian_21_filament(3*(ii-1)+1:3*(ii),3*(jj-1)+1:3*(jj))=kDirac
    enddo
  enddo


  do mm=1,F_rb
    index1_internal=index1(mm)-mm+1
    index2_internal=index1_internal+Filament_num(mm)-2
    do ii=index1_internal,index2_internal
      jj=ii
      Jacobian_22_filament(3*(ii-1)+1:3*(ii),3*(jj-1)+1:3*(jj))=2.0_8/3.0_8*DT_DEM*DT_DEM/mass_par*kDirac
      if(ii.lt.index2_internal) then
      jj=ii+1
      Jacobian_22_filament(3*(ii-1)+1:3*(ii),3*(jj-1)+1:3*(jj))=-2.0_8/3.0_8*DT_DEM*DT_DEM/mass_par*kDirac
      endif

      kk=ii+mm
      do jj=1,Nfilament
        if(jj.lt.kk.and.jj.gt.index1_internal) then
          Jacobian_23_filament(3*(ii-1)+1:3*(ii),3*(jj-1)+1:3*(jj))=Mat_inv_Cross(tau_next(:,jj))
        elseif(jj.eq.kk)then
          Jacobian_23_filament(3*(ii-1)+1:3*(ii),3*(jj-1)+1:3*(jj))=0.5_8*Mat_inv_Cross(tau_next(:,jj))
        elseif(jj.lt.kk.and.jj.eq.index1_internal)then
          Jacobian_23_filament(3*(ii-1)+1:3*(ii),3*(jj-1)+1:3*(jj))=0.5_8*Mat_inv_Cross(tau_next(:,jj))
        endif
      enddo

    enddo
  enddo


  do mm=1,F_rb
    index1_internal=index1(mm)-mm+1
    index2_internal=index1_internal+Filament_num(mm)-2
    !do ii=1,Nfilament-F_rb
    do jj=index1_internal,index2_internal
      kk=jj+mm-1
      !if((ii.ge.index1_internal).and.(ii.lt.index2_internal)) then
      ii=kk
      call Jacobian_32_sub(tau_now(:,ii),Lie_algebra_next(3*(ii-1)+1:3*ii),Jacobian_sub)
      Jacobian_32_filament(3*(ii-1)+1:3*(ii),3*(jj-1)+1:3*(jj))=-1.0_8/3.0_8*DT_DEM*DT_DEM*Jacobian_sub
      if(jj.ne.index2_internal) then
        ii=kk+1
        call Jacobian_32_sub(tau_now(:,ii),Lie_algebra_next(3*(ii-1)+1:3*ii),Jacobian_sub)
        Jacobian_32_filament(3*(ii-1)+1:3*(ii),3*(jj-1)+1:3*(jj))=-1.0_8/3.0_8*DT_DEM*DT_DEM*Jacobian_sub
      endif
    enddo
  enddo

  do ii=1,Nfilament
    jj=ii
    call Jacobian_33_sub(omega_next(3*(ii-1)+1:3*ii),Lie_algebra_next(3*(ii-1)+1:3*ii),Jacobian_sub)
    Jacobian_33_filament(3*(ii-1)+1:3*(ii),3*(jj-1)+1:3*(jj))=kDirac-1.0_8/3.0_8*DT_DEM*Jacobian_sub
  enddo



  Jacobian_filament(1:3*F_rb,1:3*F_rb)=Jacobian_11_filament
  Jacobian_filament(1:3*F_rb,3*F_rb+1:3*Nfilament)  =Jacobian_12_filament
  Jacobian_filament(3*F_rb+1:3*Nfilament,1:3*F_rb)=Jacobian_21_filament
  Jacobian_filament(3*F_rb+1:3*Nfilament,3*F_rb+1:3*Nfilament)=Jacobian_22_filament
  Jacobian_filament(3*F_rb+1:3*Nfilament,3*Nfilament+1:6*Nfilament)=Jacobian_23_filament
  Jacobian_filament(3*Nfilament+1:6*Nfilament,3*F_rb+1:3*Nfilament)=Jacobian_32_filament
  Jacobian_filament(3*Nfilament+1:6*Nfilament,3*Nfilament+1:6*Nfilament)=Jacobian_33_filament

  !write(*,*) "det(Jacobian_filament)===============",det(Jacobian_filament)

  end  subroutine Jacobian_filament_function



  subroutine Jacobian_32_sub(tau_now,Lie_algebra_next,Jacobian_sub)
  IMPLICIT NONE
  real*8,intent(in)::tau_now(3),Lie_algebra_next(3)
  real*8,intent(out)::Jacobian_sub(3,3)

  real*8 submat_Lie_algebra(3,3),submat_tau(3,3),submat0(3,3),submat1(3,3)


  submat_Lie_algebra=Mat_Cross(Lie_algebra_next)
  submat_tau=Mat_Cross(tau_now)
  submat0=matmul(Filament_Inertial_body_inverse,submat_tau)
  submat1=kDirac-0.5_8*submat_Lie_algebra+1.0_8/12.0_8*matmul(submat_Lie_algebra,submat_Lie_algebra)
  Jacobian_sub=matmul(submat1,submat0)

  end subroutine Jacobian_32_sub

  subroutine Jacobian_33_sub(omega_next,Lie_algebra_next,Jacobian_sub)
  IMPLICIT NONE
  real*8,intent(in)::omega_next(3),Lie_algebra_next(3)
  real*8,intent(out)::Jacobian_sub(3,3)

  real*8 submat_Lie_algebra(3,3),submat_omega(3,3),submat0(3,3),submat1(3,3)

  submat_Lie_algebra=Mat_Cross(Lie_algebra_next)
  submat_omega=Mat_inv_Cross(omega_next)

  submat0=Mat_inv_Cross(CrossProduct3D(Lie_algebra_next,omega_next))
  submat1=matmul(submat_Lie_algebra,submat_omega)
  Jacobian_sub=-0.5_8*submat_omega+1.0_8/12.0_8*(submat0+submat1)

  end subroutine Jacobian_33_sub

  subroutine filament_function(Nfilament,X1_now,X1_past,U1_now,tau_now,Lie_algebra_now, &
                    & Lie_algebra_next,X1_next,Inertial_lamda,&
                    & Inertial_torque,Fe,conf_now,conf_past,conf_next,U_now,omega_next,filament_Fy)
  IMPLICIT NONE
  integer:: Nfilament
  real*8,intent(in)::Lie_algebra_next(3*Nfilament),X1_next(3*F_rb),Inertial_lamda(3*Nfilament)
  real*8,intent(in)::X1_now(3*F_rb),X1_past(3*F_rb),U1_now(3*F_rb),U_now(6*Nfilament),omega_next(3*Nfilament)
  real*8,intent(in)::tau_now(3,Nfilament),Inertial_torque(3*Nfilament),Fe(6*Nfilament)
  real*8,intent(in)::Lie_algebra_now(3*Nfilament)
  real*8,intent(in)::conf_now(3,Nfilament),conf_past(3,Nfilament),conf_next(3,Nfilament)
  !type(quaternion),intent(in):: q(Nfilament)
  real*8,intent(out)::filament_Fy(6*Nfilament)

  integer::ii,jj,kk,index_Y2,KK_fila,num_fila,num_fila_sum

  !type(quaternion):: q_next(Nfilament)
  real*8 Fe_ext(3*Nfilament)!,Te_ext(3*Nfilament)
  real*8 X1Inertial_lamda(3*F_rb),dexpu(3)!,Lie_algebra_next(3*Nfilament)
  !,tau_next(3,Nfilament)

  !call Yk_To_sub(Nfilament,Yk_next,X1_next,Inertial_lamda,Lie_algebra_next)

  filament_Fy=0.0_8
  Fe_ext=Fe(1:3*Nfilament)
  

  !Te_ext=Fe(3*Nfilament+1:6*Nfilament)
  call Find_X1Inertial_lamda(Nfilament,Inertial_lamda,X1Inertial_lamda)

  !do jj=1,Nfilament
  !    write(*,*) 'Inertial_lamda==',jj,Inertial_lamda(3*(jj-1)+1:3*jj)
  !enddo

  do ii=1,F_rb
      filament_Fy(3*(ii-1)+1:3*ii)=X1_next(3*(ii-1)+1:3*ii)-4.0_8/3.0_8*X1_now(3*(ii-1)+1:3*ii) &
      & +1.0_8/3.0_8*X1_past(3*(ii-1)+1:3*ii) &
      & -2.0_8/3.0_8*DT_DEM*DT_DEM/mass_par*X1Inertial_lamda(3*(ii-1)+1:3*ii) &
      & -2.0_8/3.0_8*DT_DEM*(U1_now(3*(ii-1)+1:3*ii)+DT_DEM/mass_par*Fe_ext(3*(index1(ii)-1)+1:3*index1(ii)))
     ! write(*,*) 'filament_Fy(3*(ii-1)+1:3*ii)==',ii,filament_Fy(3*(ii-1)+1:3*ii)

  enddo


  num_fila_sum=0
  index_Y2=0
  do  KK_fila=1,F_rb
      num_fila=Filament_num(KK_fila)
      do ii = 2, num_fila
        index_Y2=index_Y2+1
        jj=ii+num_fila_sum
        if(ii.ne.num_fila)then
          filament_Fy(3*F_rb+3*(index_Y2-1)+1:3*F_rb+3*index_Y2)=conf_next(:,jj) &
        & -4.0_8/3.0_8*conf_now(:,jj)+1.0_8/3.0_8*conf_past(:,jj) &
        & -2.0_8/3.0_8*DT_DEM*DT_DEM/mass_par*(Inertial_lamda(3*(jj-1)+1:3*jj)-Inertial_lamda(3*(jj-2)+1:3*jj-3))&
        & -2.0_8/3.0_8*DT_DEM*(U_now(3*(jj-1)+1:3*jj)+DT_DEM/mass_par*Fe_ext(3*(jj-1)+1:3*jj))
        else
          filament_Fy(3*F_rb+3*(index_Y2-1)+1:3*F_rb+3*index_Y2)=conf_next(:,jj) &
        & -4.0_8/3.0_8*conf_now(:,jj)+1.0_8/3.0_8*conf_past(:,jj) &
        & +2.0_8/3.0_8*DT_DEM*DT_DEM/mass_par*Inertial_lamda(3*(jj-2)+1:3*jj-3)&
        & -2.0_8/3.0_8*DT_DEM*(U_now(3*(jj-1)+1:3*jj)+DT_DEM/mass_par*Fe_ext(3*(jj-1)+1:3*jj))
        endif
      !write(*,*) 'filament_Fy(3*(ii-1)+1:3*ii)==',index_Y2,filament_Fy(3*F_rb+3*(index_Y2-1)+1:3*F_rb+3*index_Y2)
      enddo
      num_fila_sum=num_fila_sum+num_fila
  enddo

  write(*,*) 'Nfilament==',Nfilament,'index_Y2==',index_Y2



  do ii=1,Nfilament
      dexpu=dexpu_inv(Lie_algebra_next(3*(ii-1)+1:3*ii),omega_next(3*(ii-1)+1:3*ii))
      filament_Fy(3*Nfilament+3*(ii-1)+1:3*Nfilament+3*ii)=Lie_algebra_next(3*(ii-1)+1:3*ii) &
      & -1.0_8/3.0_8*Lie_algebra_now(3*(ii-1)+1:3*ii) &
      & -2.0_8/3.0_8*DT_DEM*dexpu
  enddo

  end  subroutine filament_function



  subroutine Rotate_omeganew(Nfilament,tau_now,Inertial_torque,Inertial_lamda,Te_ext,omega_now,omega_new)
  integer,intent(in)::Nfilament
  real*8,intent(in)::tau_now(3,Nfilament),Inertial_torque(3*Nfilament),Inertial_lamda(3*Nfilament)
  real*8,intent(in)::Te_ext(3*Nfilament),omega_now(3*Nfilament)
  real*8,intent(out)::omega_new(3*Nfilament)

  integer:: ii,jj,KK_fila,num_fila,num_fila_sum
  real*8:: angular_momentum(3),Inertial_lamda_left(3),Inertial_lamda_right(3)
  real*8:: inertial_torque_right(3),Inertial_torque_left(3)

  do  KK_fila=1,F_rb       
      num_fila=Filament_num(KK_fila)
      do ii = 1, num_fila
        jj=ii+num_fila_sum
        if(ii.eq.1) then
          Inertial_lamda_left=0.0_8
          Inertial_lamda_right=Inertial_lamda(3*(jj-1)+1:3*jj)
          Inertial_torque_left=0.0_8
          Inertial_torque_right=Inertial_torque(3*(jj-1)+1:3*jj)
        elseif(ii.eq.num_fila) then
          Inertial_lamda_left=Inertial_lamda(3*(jj-2)+1:3*jj-3)
          Inertial_lamda_right=0.0_8
          Inertial_torque_left=Inertial_torque(3*(jj-2)+1:3*jj-3)
          Inertial_torque_right=0.0_8
        else
          Inertial_lamda_left=Inertial_lamda(3*(jj-2)+1:3*jj-3)
          Inertial_lamda_right=Inertial_lamda(3*(jj-1)+1:3*jj)
          Inertial_torque_left=Inertial_torque(3*(jj-2)+1:3*jj-3)
          Inertial_torque_right=Inertial_torque(3*(jj-1)+1:3*jj)
        endif
        angular_momentum=0.5*CrossProduct3D(tau_now(:,jj),(Inertial_lamda_left+Inertial_lamda_right)) &
        & +Inertial_torque_right-Inertial_torque_left+Te_ext(3*(ii-1)+1:3*ii)
        omega_new(3*(ii-1)+1:3*ii)=omega_now(3*(ii-1)+1:3*ii)+DT_DEM*matmul(Filament_Inertial_body_inverse, &
         & angular_momentum)
      enddo
      num_fila_sum=num_fila_sum+num_fila
  enddo

  end subroutine Rotate_omeganew

  subroutine filament_X1ToCONF(Nfilament,X1_next,tau_next,conf_next)
  integer,intent(in):: Nfilament
  real*8,intent(in) ::X1_next(3*F_rb),tau_next(3,Nfilament)
  real*8,intent(out)::conf_next(3,Nfilament)

  integer::ii,jj,KK_fila,num_fila,num_fila_sum

  num_fila_sum=0
  do  KK_fila=1,F_rb       
      num_fila=Filament_num(KK_fila)
      conf_next(:,num_fila_sum+1)=X1_next(3*(KK_fila-1)+1:3*KK_fila)
      do ii = 1, num_fila-1
          jj=ii+num_fila_sum
          conf_next(:,jj+1)=conf_next(:,jj)+0.5_8*(tau_next(:,jj)+tau_next(:,jj+1))
      enddo
      num_fila_sum=num_fila_sum+num_fila
  enddo

  end  subroutine filament_X1ToCONF


  subroutine Lie_algebra_To_q(Nfilament,q,Lie_algebra_next,q_next)
  IMPLICIT NONE
  integer:: Nfilament
  real*8,intent(in)::Lie_algebra_next(3*Nfilament)
  type(quaternion),intent(in):: q(Nfilament)
  type(quaternion),intent(out):: q_next(Nfilament)

  integer:: ii

  do ii=1,Nfilament
      q_next(ii)=Lie_algebra_update(Lie_algebra_next(3*(ii-1)+1:3*ii),q(ii))
  enddo
  end subroutine Lie_algebra_To_q


  subroutine q_To_tau(Nfilament,q,tau)
  IMPLICIT NONE
  integer:: Nfilament
  type(quaternion),intent(in):: q(Nfilament)
  real*8,intent(out)::tau(3,Nfilament)

  integer:: ii
  do ii=1,Nfilament
      tau(:,ii)=qrot(q(ii),Filament_tau_base(:,ii))
  enddo
  end subroutine q_To_tau

  subroutine Find_X1Inertial_lamda(Nfilament,Inertial_lamda,X1Inertial_lamda)
  integer,intent(in):: Nfilament
  real*8,intent(in) ::Inertial_lamda(3*Nfilament)
  real*8,intent(out)::X1Inertial_lamda(3*F_rb)

  integer::ii

  do ii=1,F_rb
      X1Inertial_lamda(3*(ii-1)+1:3*ii)=Inertial_lamda(3*(index1(ii)-1)+1:3*index1(ii))
  enddo

  end  subroutine Find_X1Inertial_lamda


  subroutine Find_exactInertial_lamda(Nfilament,Inertial_lamda,exactInertial_lamda)
  integer,intent(in):: Nfilament
  real*8,intent(in) ::Inertial_lamda(3*Nfilament)
  real*8,intent(out)::exactInertial_lamda(3*(Nfilament-F_rb))

  integer::ii,indexstart,indexend,length

  exactInertial_lamda=0.0_8
  indexstart=1
  do ii=1,F_rb
    length=Filament_num(ii)-1
    indexend=indexstart+length-1
    exactInertial_lamda(3*(indexstart-1)+1:3*indexend)=Inertial_lamda(3*(index1(ii)-1)+1:3*(index1(ii)-1+length))
    indexstart=indexend+1
  enddo

  end  subroutine Find_exactInertial_lamda

  subroutine exactToInertial_lamda(Nfilament,Inertial_lamda,exactInertial_lamda)
  integer,intent(in):: Nfilament
  real*8,intent(in) ::exactInertial_lamda(3*(Nfilament-F_rb))
  real*8,intent(out)::Inertial_lamda(3*Nfilament)

  integer::ii,indexstart,indexend,length

  Inertial_lamda=0.0_8
  indexstart=1
  do ii=1,F_rb
    length=Filament_num(ii)-1
    indexend=indexstart+length-1
    Inertial_lamda(3*(index1(ii)-1)+1:3*(index1(ii)-1+length))=exactInertial_lamda(3*(indexstart-1)+1:3*indexend)
    indexstart=indexend+1
  enddo

  end  subroutine exactToInertial_lamda

  end Module filament_solve_Jacobian



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


    subroutine filament_Init(Nfilament,CONF,RADII,U_pos)
    IMPLICIT NONE
    INTEGER,intent(in)::Nfilament
    REAL*8,intent(in):: CONF(3,Nfilament),RADII(Nfilament)
    real*8,intent(out):: U_pos(6*Nfilament)

    integer :: ii, jj, KK_fila
    integer ::num_fila_sum,num_fila
    real*8 :: Ixx, Iyy,Izz,RADII2

      U_pos=0.0_8
      Filament_conf_now=CONF
      Filament_conf_past=CONF
      Filament_U1_now=0.0_8
      Filament_tau_now(:,:)=0.0_8
      Filament_Lie_algebra_now=0.0_8

      num_fila_sum=0
      do  KK_fila=1,F_rb       
          num_fila=Filament_num(KK_fila)
          Filament_X1_now(3*(KK_fila-1)+1:3*KK_fila) = conf(:,num_fila_sum+1)
          do ii = 1, num_fila
              jj=ii+num_fila_sum
              if(ii.eq.1)then
                Filament_tau_now( :, jj) = conf(:,jj+1)-conf(:,jj)
              elseif(ii.eq.num_fila)then
                Filament_tau_now( :, jj) = conf(:,jj)-conf(:,jj-1)
              else
                Filament_tau_now( :, jj) = conf(:,jj+1)-conf(:,jj)
              endif
              write(*,*) "Filament_tau_base====",jj,Filament_tau_now( :, jj)
              Filament_q(jj)=quat(1.0_8,0.0_8,0.0_8,0.0_8)
          end do
          index1(KK_fila)=num_fila_sum+1
         ! Filament_X1(3*(KK_fila-1):3*KK_fila)=conf(:,num_fila_sum+1)
         ! Filament_U1(3*(KK_fila-1):3*KK_fila)=U_pos(3*num_fila_sum+1:3*num_fila_sum+3)
          num_fila_sum=num_fila_sum+num_fila
      end do

      Filament_X1_past=Filament_X1_now
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


    subroutine InternalForcesAndTorques(Nfilament,conf,Filament_internal_force_torque)
    ! Filament.INTERNALFORCESANDTORQUES()  Place gravity, elastic and 
    !                                      constraint forces and torques 
    !                                      on the filament. They are stored
    !                                      in Filament.F
    IMPLICIT NONE
    INTEGER,intent(in)::Nfilament
    real*8,intent(in)::conf(3,Nfilament)!,Filament_tau_base(3,NN)
    real*8,intent(out)::Filament_internal_force_torque(6*Nfilament)
    integer :: ii, jj, KK_fila
    integer ::num_fila_sum,num_fila
    real*8 :: dr(3),ds,tauj(3),tauj1(3),M(3),F0(3),F1(3),tauF(3)
    type(quaternion)::q_mid!,dqds,qb
    real*8 :: internal_force(3),internal_torque(3),angle(3)


      Filament_internal_force_torque=0.0_8
      num_fila_sum=0
      do  KK_fila=1,F_rb       
          num_fila=Filament_num(KK_fila)       
          do ii = 1, num_fila-1
              jj=ii+num_fila_sum
              if(ii.eq.num_fila)then
                Filament_internal_force_torque( 3*(jj-1)+1:3*jj) = 0.0_8
                Filament_internal_force_torque( 3*Nfilament+3*(jj-1)+1:3*Nfilament+3*jj) = 0.0_8
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
                write(*,*) "dr,ds,angle===",dr(:),ds,angle(:)
                F0=qrot_inverse(q_mid,tauF)
                F0=F0-Filament_tau_base(:,jj)/ds
                F1=0.0_8
                F1(1)=F0(1)*EA;
                F1(2)=F0(2)*GA;F1(3)=F0(3)*GA
                internal_force=qrot(q_mid,F1)
                write(*,*) 'index_filament,MMMMMMMMMMMMMMM===========',jj,M(:)!,GI,EI
                write(*,*) 'index_filament,internal_torque===========',jj,internal_torque(:)
                write(*,*) 'index_filament,FFFFFFFFFFFFFFF===========',jj,F1(:)
                write(*,*) 'index_filament,internal_Force ===========',jj,internal_force(:)

                Filament_internal_force_torque( 3*(jj-1)+1:3*jj) = Filament_internal_force_torque( 3*(jj-1) &
                  & +1:3*jj) +internal_force(:)

                Filament_internal_force_torque( 3*jj+1:3*jj+3) = Filament_internal_force_torque( 3*jj+1:3*jj+3)&
                & -internal_force(:)

                Filament_internal_force_torque( 3*Nfilament+3*(jj-1)+1:3*Nfilament+3*jj) = &
                & Filament_internal_force_torque( 3*Nfilament+3*(jj-1)+1:3*Nfilament+3*jj)  & 
                & +internal_torque(:)+CrossProduct3D(0.5*dr,internal_force)

                Filament_internal_force_torque( 3*Nfilament+3*jj+1:3*Nfilament+3*jj+3) =  &
                & Filament_internal_force_torque( 3*Nfilament+3*jj+1:3*Nfilament+3*jj+3)  &
                & -internal_torque(:)+CrossProduct3D(0.5*dr,internal_force)
              endif
          end do
          num_fila_sum=num_fila_sum+num_fila
      end do
    do ii=1,Nfilament
        write(*,*) 'i,internal_filament_torque===',ii, &
        & Filament_internal_force_torque(3*Nfilament+3*(ii-1)+1:3*Nfilament+3*ii)
            !write(*,*) 'i, torue===',ii,Ftotal(3*NN+3*(ii-1)+1:3*NN+3*ii)

    enddo  


      end subroutine InternalForcesAndTorques



    subroutine new_U1_filament(Nfilament,RADII,Ftotal,U_pos)
    IMPLICIT NONE
    INTEGER,intent(in)::Nfilament
    REAL*8,intent(in):: RADII(Nfilament),Ftotal(6*Nfilament)
    real*8,intent(inout):: U_pos(6*Nfilament)

    integer :: ii,KK_fila

    real*8 :: filament_torque(3*Nfilament)


    do ii=1,Nfilament
      filament_torque(3*(ii-1)+1:3*ii)=matmul(Filament_Inertial_body_inverse, &
        & Ftotal(3*Nfilament+3*(ii-1)+1:3*Nfilament+3*ii))
    enddo

    do ii=1,Nfilament
        write(*,*) 'i,filament_torque/Internal_moment===',ii,filament_torque(3*(ii-1)+1:3*ii)
            !write(*,*) 'i, torue===',ii,Ftotal(3*NN+3*(ii-1)+1:3*NN+3*ii)
    enddo
   !CALL MATREV(rbmconn_Inertial,3*K_rb)
    U_pos(3*Nfilament+1:6*Nfilament)=U_pos(3*Nfilament+1:6*Nfilament)+filament_torque*DT_DEM/mass_par
    U_pos(1:3*Nfilament)=U_pos(1:3*Nfilament)+Ftotal(1:3*Nfilament)*DT_DEM/mass_par


    !do  KK_fila=1,F_rb
    !do  ii=1,NN
        !Filament_U1(3*(KK_fila-1)+1:3*KK_fila)=Filament_U1(3*(KK_fila-1)+1:3*KK_fila)+ &
        !    & Ftotal(3*(index1(KK_fila)-1):3*index1(KK_fila))*DT_DEM/mass_par  
        !U_pos(3*(ii-1)+1:3*ii)=Filament_U1(3*(KK_fila-1)+1:3*KK_fila)
    !end do            

    end subroutine new_U1_filament



    subroutine new_conf_filament(Nfilament,conf,T,U_pos)
    IMPLICIT NONE
    INTEGER,intent(in)::Nfilament
    REAL*8,intent(in):: T,U_pos(6*Nfilament)!,U_pos(6*NN)!,RADII(NN)
    !real*8,intent(in):: U_par_rb(6*K_rb)
    real*8,intent(inout):: conf(3,Nfilament)
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
    if(.not.solve_implicit) then
        do ii=1,Nfilament
          jj=3*Nfilament+3*(ii-1)
          w=U_pos(jj+1:jj+3)
          dq=qscal(DT_DEM,qderivat(Filament_q(ii),w))
          q_norm=qadd(Filament_q(ii),dq)
          Filament_q(ii)=qscal(1.0_8/qnorm(q_norm) , q_norm)
          !write(*,*) "ii==",ii
          !CALL qprint(Filament_q(ii))
          Filament_tau_now(:,ii)=qrot(Filament_q(ii),Filament_tau_base(:,ii))
        enddo

        forall(ii=1:Nfilament,jj=1:3)
          CONF(jj,ii)=CONF(jj,ii)+DT_DEM*U_pos(3*(ii-1)+jj)
        end forall
    else
        !filament_implicit_solve(Nfilament,Filament_X1_now,Filament_X1_past,Filament_U1_now,Filament_tau_now, &
        !    & Filament_Lie_algebra_now,Filament_Interal_force,Fe,Filament_conf_now,Filament_conf_past,U_pos,Filament_q)

        CONF=Filament_conf_now

    endif


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
         DO ii = 1, Nfilament
          CALL PERIOD_CORRECTION(conf(:,ii),LB,T)
         ENDDO
      endif


    end subroutine new_conf_filament

    end  Module filament_math
