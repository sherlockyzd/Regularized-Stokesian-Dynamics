!******************************************************
!******************************************************
!******************************************************
  module hydro_lub_mod
  use prutil
  use size,only:NB
  use TENSORS,only:Y21,Y22,Y2,front
  USE method
  use CONFIG, only:RADII
  use period_bdy_tools,only:PER_SKEW_CORRECTION

  implicit none
  private


  public::lubmxcalc

  contains


    subroutine lubmxcalc(ntotal,P_pos,rfu,rfe,rse)
    implicit none
    integer, intent(in) :: ntotal
    real(kind=cp), dimension(3,ntotal), intent(in) :: P_pos
    real(kind=cp), dimension(6*ntotal,6*ntotal), intent(out) :: rfu
    real(kind=cp), dimension(6*ntotal,5*ntotal), intent(out) :: rfe
    real(kind=cp), dimension(5*ntotal,5*ntotal), intent(out) :: rse

    real(kind=cp) :: r,rx(3)
    real(kind=cp) :: x11a,x12a,y11a,y12a,y11b,y12b,x11c,x12c,y11c,y12c
    real(kind=cp) :: x11g,x12g,y11g,y12g,y11h,y12h,xm,ym,zm
    real(kind=cp), dimension(3) :: e
    real(kind=cp), dimension(2,2,3,3) :: alub
    real(kind=cp), dimension(2,2,3,3) :: blub
    real(kind=cp), dimension(2,2,3,3) :: club
    real(kind=cp), dimension(2,2,3,3,3) :: glubtensor
    real(kind=cp), dimension(2,2,3,3,3) :: hlubtensor
    real(kind=cp), dimension(2,2,5,3) :: glubtensor_Y2,hlubtensor_Y2
    real(kind=cp), dimension(3,3,3,3) :: mlubtensor
    real(kind=cp), dimension(5,5) :: mmx
    !real(kind=cp), dimension(ntotal,3) :: p_pos
    real(kind=cp), dimension(6*ntotal,6*ntotal) :: APPl
    real(kind=cp), dimension(5*ntotal,6*ntotal) :: AQPl
    real(kind=cp), dimension(5*ntotal,5*ntotal) :: AQQl
    !real(kind=cp), dimension(3*ntotal,3*ntotal) :: APPl_B

    integer :: alpha,beta,NN
    integer :: i,j,k,l,m,n

    real*8 V0(3,3,3),V1(5,3),r_ave,s


    NN=ntotal
    rfu=0.0_8
    rfe=0.0_8
    rse=0.0_8
    APPl=0.0_8
    AQPl=0.0_8
    AQQl=0.0_8
    glubtensor_Y2=0.0_8
    hlubtensor_Y2=0.0_8
    !APPl_B=0.0_8


     do alpha=1,ntotal
       do beta=alpha+1,ntotal
          if(alpha.gt.NN-Nb) then
         ! write(*,*) 'boundayr pass-------------------------'
           cycle
          endif
            rx=p_pos(:,alpha)-p_pos(:,beta)
            r=sqrt(sum(rx*rx))
          if(IsPeriod) then
            CALL PER_SKEW_CORRECTION(rx,r)
          endif
            alub=0.0_8
            blub=0.0_8
            club=0.0_8
            glubtensor=0.0_8
            hlubtensor=0.0_8
            mlubtensor=0.0_8
            r_ave=0.5_8*(radii(alpha)+radii(beta))
            s=r/r_ave
            if(s.le.2.0_8) then
                  write(*,*) 'lubrication Contact Error, r='
                   write(*,*) s,alpha,beta
                    !cycle 
                    s=4.00001_8-s
                    call lubcalc(s,x11a,x12a,y11a,y12a,y11b,y12b,x11c,x12c,y11c,y12c,&
                            & x11g,x12g,y11g,y12g,y11h,y12h,xm,ym,zm)
                    !stop
                    !cycle
            else if(s.gt.4.0_8) then
                    cycle                
            else
            ! initialize the lub2bmx matrix
            ! build the 2body lubrication resistance matrix lub2bmx(22*22) 
            !write(*,*) 's-2',s-2.0_8
            call lubcalc(s,x11a,x12a,y11a,y12a,y11b,y12b,x11c,x12c,y11c,y12c,&
                            & x11g,x12g,y11g,y12g,y11h,y12h,xm,ym,zm)

            endif



      
            !e(1)=r1/r
            !e(2)=r2/r
            !e(3)=r3/r
            e=rx/r
            forall(i=1:3,j=1:3)
            alub(1,1,i,j)=x11a*e(i)*e(j)+y11a*(kd(i,j)-e(i)*e(j))
            alub(1,2,i,j)=x12a*e(i)*e(j)+y12a*(kd(i,j)-e(i)*e(j))
            blub(1,1,i,j)=y11b*(per(i,j,1)*e(1)+per(i,j,2)*e(2)+per(i,j,3)*e(3))
            blub(1,2,i,j)=y12b*(per(i,j,1)*e(1)+per(i,j,2)*e(2)+per(i,j,3)*e(3))
            club(1,1,i,j)=x11c*e(i)*e(j)+y11c*(kd(i,j)-e(i)*e(j))
            club(1,2,i,j)=x12c*e(i)*e(j)+y12c*(kd(i,j)-e(i)*e(j))
            end forall

            forall(i=1:3,j=1:3,k=1:3)
            glubtensor(1,1,i,j,k)=x11g*(e(i)*e(j)-kd(i,j)/3.0_8)*e(k)& 
            &+y11g*(e(i)*kd(j,k)+e(j)*kd(i,k)-2*e(i)*e(j)*e(k))
            glubtensor(1,2,i,j,k)=x12g*(e(i)*e(j)-kd(i,j)/3.0_8)*e(k)& 
            &+y12g*(e(i)*kd(j,k)+e(j)*kd(i,k)-2*e(i)*e(j)*e(k))
            hlubtensor(1,1,i,j,k)=y11h*((per(i,k,1)*e(1)+per(i,k,2)*e(2)+&
            &per(i,k,3)*e(3))*e(j) + (per(j,k,1)*e(1)+per(j,k,2)*e(2)+&
            &per(j,k,3)*e(3))*e(i))
            hlubtensor(1,2,i,j,k)=y12h*((per(i,k,1)*e(1)+per(i,k,2)*e(2)+&
            &per(i,k,3)*e(3))*e(j) + (per(j,k,1)*e(1)+per(j,k,2)*e(2)+&
            &per(j,k,3)*e(3))*e(i))
            end forall 
            
            forall(i=1:3,j=1:3,k=1:3,l=1:3)
            mlubtensor(i,j,k,l)=xm*(3.0_8*(e(i)*e(j)-kd(i,j)/3.0_8)*(e(k)*e(l)&
            & -kd(k,l)/3.0_8)/2.0_8)+ym*(e(i)*kd(j,l)*e(k)+e(j)*kd(i,l)*e(k)&
            & +e(i)*kd(j,k)*e(l)+e(j)*kd(i,k)*e(l)-4.0_8*e(i)*e(j)*e(k)*e(l))&
            & /2.0_8+zm*(kd(i,k)*kd(j,l)+kd(j,k)*kd(i,l)-kd(i,j)*kd(k,l)& 
            & +e(i)*e(j)*kd(k,l)+kd(i,j)*e(k)*e(l)-e(i)*kd(j,l)*e(k)&
            & -e(j)*kd(i,l)*e(k)-e(i)*kd(j,k)*e(l)-e(j)*kd(i,k)*e(l)&
            & +e(i)*e(j)*e(k)*e(l))/2.0_8
            end forall
          
            forall(i=1:3,j=1:3)
            alub(2,2,i,j)=alub(1,1,i,j) !x22a=x11a,y22a=y11a
            alub(2,1,i,j)=alub(1,2,i,j) !x12a=x21a,y12a=y21a
            blub(2,2,i,j)=-blub(1,1,i,j) !y11b=-y22b
            blub(2,1,i,j)=-blub(1,2,i,j) !y12b=-y21b
            club(2,2,i,j)=club(1,1,i,j) !x11c=x22c
            club(2,1,i,j)=club(1,2,i,j) !x21c=x12c
            end forall

            forall(i=1:3,j=1:3,k=1:3)
            glubtensor(2,2,i,j,k)=-glubtensor(1,1,i,j,k) ! x22g=-x11g, y22g=-y11g
            glubtensor(2,1,i,j,k)=-glubtensor(1,2,i,j,k) ! x21g=-x12g, y21g=-y12g
            hlubtensor(2,2,i,j,k)=hlubtensor(1,1,i,j,k) !  y22h=y11h 
            hlubtensor(2,1,i,j,k)=hlubtensor(1,2,i,j,k) ! y21h=y12h
            end forall


            call DD_transY2(mmx,mlubtensor)



           do m=1,2
            do n=1,2
               !CALL mulT3T2(glubtensor(m,n,1:3,1:3,1:3),Y2(I,1:3,1:3),V)
                V0=glubtensor(m,n,1:3,1:3,1:3)
                call RTD_front_transY2(V1,V0)
               glubtensor_Y2(m,n,1:5,1:3)=V1
               !CALL mulT3T2(hlubtensor(m,n,1:3,1:3,1:3),Y2(I,1:3,1:3),V)
                V0=hlubtensor(m,n,1:3,1:3,1:3)
               call RTD_front_transY2(V1,V0)
               hlubtensor_Y2(m,n,1:5,1:3)=V1
            enddo
          enddo
               
          

           !DO I=1,2
           !CALL PER_ROTNE_PRAGER_TT_SELF(C1PP,RADII(I))
           APPl(3*(alpha-1)+1:3*(alpha-1)+3,3*(alpha-1)+1:3*(alpha-1)+3) =&
           APPl(3*(alpha-1)+1:3*(alpha-1)+3,3*(alpha-1)+1:3*(alpha-1)+3)+ alub(1,1,1:3,1:3)
           APPl(3*(beta-1)+1:3*(beta-1)+3,3*(beta-1)+1:3*(beta-1)+3) =&
           APPl(3*(beta-1)+1:3*(beta-1)+3,3*(beta-1)+1:3*(beta-1)+3)+ alub(2,2,1:3,1:3)

           !CALL PER_ROTNE_PRAGER_RR_SELF(C1PP,RADII(I))
           APPl(3*NN+3*(alpha-1)+1:3*NN+3*(alpha-1)+3,3*NN+3*(alpha-1)+1:3*NN+3*(alpha-1)+3)= &
           APPl(3*NN+3*(alpha-1)+1:3*NN+3*(alpha-1)+3,3*NN+3*(alpha-1)+1:3*NN+3*(alpha-1)+3)+club(1,1,1:3,1:3)
           APPl(3*NN+3*(beta-1)+1:3*NN+3*(beta-1)+3,3*NN+3*(beta-1)+1:3*NN+3*(beta-1)+3)= &
           APPl(3*NN+3*(beta-1)+1:3*NN+3*(beta-1)+3,3*NN+3*(beta-1)+1:3*NN+3*(beta-1)+3)+club(2,2,1:3,1:3)

           !CALL PER_ROTNE_PRAGER_RT_SELF(C1PP,RADII(I))
           APPl(3*NN+3*(alpha-1)+1:3*NN+3*alpha,3*(alpha-1)+1:3*alpha)= &
           APPl(3*NN+3*(alpha-1)+1:3*NN+3*alpha,3*(alpha-1)+1:3*alpha)+blub(1,1,1:3,1:3)        
           APPl(3*NN+3*(beta-1)+1:3*NN+3*beta,3*(beta-1)+1:3*beta)= &
           APPl(3*NN+3*(beta-1)+1:3*NN+3*beta,3*(beta-1)+1:3*beta)+blub(2,2,1:3,1:3)  

           !CALL PER_ROTNE_PRAGER_TD_SELF_Y2(C1PQ,RADII(I))

           AQPl(5*(alpha-1)+1:5*(alpha-1)+5,3*(alpha-1)+1:3*(alpha-1)+3) = &
           AQPl(5*(alpha-1)+1:5*(alpha-1)+5,3*(alpha-1)+1:3*(alpha-1)+3) + glubtensor_Y2(1,1,1:5,1:3)
           AQPl(5*(beta-1)+1:5*(beta-1)+5,3*(beta-1)+1:3*(beta-1)+3) = &
           AQPl(5*(beta-1)+1:5*(beta-1)+5,3*(beta-1)+1:3*(beta-1)+3) + glubtensor_Y2(2,2,1:5,1:3)

           AQPl(5*(alpha-1)+1:5*(alpha-1)+5,3*NN+3*(alpha-1)+1:3*NN+3*(alpha-1)+3) = &
           AQPl(5*(alpha-1)+1:5*(alpha-1)+5,3*NN+3*(alpha-1)+1:3*NN+3*(alpha-1)+3) + hlubtensor_Y2(1,1,1:5,1:3)
           AQPl(5*(beta-1)+1:5*(beta-1)+5,3*NN+3*(beta-1)+1:3*NN+3*(beta-1)+3) = &
           AQPl(5*(beta-1)+1:5*(beta-1)+5,3*NN+3*(beta-1)+1:3*NN+3*(beta-1)+3) + hlubtensor_Y2(2,2,1:5,1:3)

#ifdef back
           APQl(3*(alpha-1)+1:3*(alpha-1)+3,5*(alpha-1)+1:5*(alpha-1)+5) = &
           APQl(3*(alpha-1)+1:3*(alpha-1)+3,5*(alpha-1)+1:5*(alpha-1)+5) + glubtensor_Y2(1,1,1:3,1:5)
           APQl(3*(beta-1)+1:3*(beta-1)+3,5*(beta-1)+1:5*(beta-1)+5) = &
           APQl(3*(beta-1)+1:3*(beta-1)+3,5*(beta-1)+1:5*(beta-1)+5) + glubtensor_Y2(2,2,1:3,1:5)

           !CALL PER_ROTNE_PRAGER_RD_SELF_Y2(C1PQ,RADII(I))
           APQl(3*NN+3*(alpha-1)+1:3*NN+3*(alpha-1)+3,5*(alpha-1)+1:5*(alpha-1)+5) =&
           APQl(3*NN+3*(alpha-1)+1:3*NN+3*(alpha-1)+3,5*(alpha-1)+1:5*(alpha-1)+5)+ hlubtensor_Y2(1,1,1:3,1:5)
           APQl(3*NN+3*(beta-1)+1:3*NN+3*(beta-1)+3,5*(beta-1)+1:5*(beta-1)+5) =&
           APQl(3*NN+3*(beta-1)+1:3*NN+3*(beta-1)+3,5*(beta-1)+1:5*(beta-1)+5) + hlubtensor_Y2(2,2,1:3,1:5)
#endif

           !CALL PER_ROTNE_PRAGER_DD_SELF_Y2(C1QQ,RADII(I))         
           AQQl(5*(alpha-1)+1:5*(alpha-1)+5,5*(alpha-1)+1:5*(alpha-1)+5) = &
           AQQl(5*(alpha-1)+1:5*(alpha-1)+5,5*(alpha-1)+1:5*(alpha-1)+5)+mmx
           AQQl(5*(beta-1)+1:5*(beta-1)+5,5*(beta-1)+1:5*(beta-1)+5) = &
           AQQl(5*(beta-1)+1:5*(beta-1)+5,5*(beta-1)+1:5*(beta-1)+5)+mmx
          !ENDDO


          !DO I=1,NN-1           
          !DO J=I+1,NN
          !aaI=RADII(I)
          !aaJ=RADII(J)
          !R=CONF(:,I)-CONF(:,J)
          !CALL PER_SKEW_CORRECTION(R,DMIN)
          !A1PP = 0.D0
            !CALL PER_ROTNE_PRAGER_TT_IJ(C1PP,R,aaI,aaJ)
            APPl(3*(alpha-1)+1:3*alpha,3*(beta-1)+1:3*beta)=&
            APPl(3*(alpha-1)+1:3*alpha,3*(beta-1)+1:3*beta)+alub(1,2,1:3,1:3)
            
            !CALL PER_ROTNE_PRAGER_RR_IJ(C1PP,R,aaI,aaJ)
            APPl(3*NN+3*(alpha-1)+1:3*NN+3*alpha,3*NN+3*(beta-1)+1:3*NN+3*beta)=&
            APPl(3*NN+3*(alpha-1)+1:3*NN+3*alpha,3*NN+3*(beta-1)+1:3*NN+3*beta)+club(1,2,1:3,1:3)



            !CALL PER_ROTNE_PRAGER_RT_IJ(C1PP,R,aaI,aaJ)
            APPl(3*NN+3*(alpha-1)+1:3*NN+3*alpha,3*(beta-1)+1:3*beta)=&
            APPl(3*NN+3*(alpha-1)+1:3*NN+3*alpha,3*(beta-1)+1:3*beta)+blub(1,2,1:3,1:3)

            !CALL PER_ROTNE_PRAGER_RT_IJ(C1PP,-R,aaJ,aaI)
            APPl(3*NN+3*(beta-1)+1:3*NN+3*beta,3*(alpha-1)+1:3*alpha)=&
            APPl(3*NN+3*(beta-1)+1:3*NN+3*beta,3*(alpha-1)+1:3*alpha)+blub(2,1,1:3,1:3)       

           AQPl(5*(alpha-1)+1:5*(alpha-1)+5,3*(beta-1)+1:3*(beta-1)+3) = &
           AQPl(5*(alpha-1)+1:5*(alpha-1)+5,3*(beta-1)+1:3*(beta-1)+3) + glubtensor_Y2(1,2,1:5,1:3)
           AQPl(5*(beta-1)+1:5*(beta-1)+5,3*(alpha-1)+1:3*(alpha-1)+3) = &
           AQPl(5*(beta-1)+1:5*(beta-1)+5,3*(alpha-1)+1:3*(alpha-1)+3) + glubtensor_Y2(2,1,1:5,1:3)

           AQPl(5*(alpha-1)+1:5*(alpha-1)+5,3*NN+3*(beta-1)+1:3*NN+3*(beta-1)+3) = &
           AQPl(5*(alpha-1)+1:5*(alpha-1)+5,3*NN+3*(beta-1)+1:3*NN+3*(beta-1)+3) + hlubtensor_Y2(1,2,1:5,1:3)
           AQPl(5*(beta-1)+1:5*(beta-1)+5,3*NN+3*(alpha-1)+1:3*NN+3*(alpha-1)+3) = &
           AQPl(5*(beta-1)+1:5*(beta-1)+5,3*NN+3*(alpha-1)+1:3*NN+3*(alpha-1)+3) + hlubtensor_Y2(2,1,1:5,1:3)

#ifdef back       
            !CALL PER_ROTNE_PRAGER_TD_IJ_Y2(C1PQ,R,aaI,aaJ)
            APQl(3*(alpha-1)+1:3*(alpha-1)+3,5*(beta-1)+1:5*(beta-1)+5)= &
            APQl(3*(alpha-1)+1:3*(alpha-1)+3,5*(beta-1)+1:5*(beta-1)+5)+glubtensor_Y2(1,2,1:3,1:5)
            !CALL PER_ROTNE_PRAGER_TD_IJ_Y2(C1PQ,-R,aaI,aaJ)
            APQl(3*(beta-1)+1:3*(beta-1)+3,5*(alpha-1)+1:5*(alpha-1)+5) = &
            APQl(3*(beta-1)+1:3*(beta-1)+3,5*(alpha-1)+1:5*(alpha-1)+5)+glubtensor_Y2(2,1,1:3,1:5)              
            !CALL PER_ROTNE_PRAGER_RD_IJ_Y2(C1PQ,R,aaI,aaJ)
            APQl(3*NN+3*(alpha-1)+1:3*NN+3*(alpha-1)+3,5*(beta-1)+1:5*(beta-1)+5)=&
            APQl(3*NN+3*(alpha-1)+1:3*NN+3*(alpha-1)+3,5*(beta-1)+1:5*(beta-1)+5)+hlubtensor_Y2(1,2,1:3,1:5)
            !CALL PER_ROTNE_PRAGER_RD_IJ_Y2(C1PQ,-R,aaI,aaJ)
            APQl(3*NN+3*(beta-1)+1:3*NN+3*(beta-1)+3,5*(alpha-1)+1:5*(alpha-1)+5)= &
            APQl(3*NN+3*(beta-1)+1:3*NN+3*(beta-1)+3,5*(alpha-1)+1:5*(alpha-1)+5)+hlubtensor_Y2(2,1,1:3,1:5)
#endif

            !CALL PER_ROTNE_PRAGER_DD_IJ_Y2(C1QQ,R,aaI,aaJ)              
            AQQl(5*(alpha-1)+1:5*(alpha-1)+5,5*(beta-1)+1:5*(beta-1)+5) = &
            AQQl(5*(alpha-1)+1:5*(alpha-1)+5,5*(beta-1)+1:5*(beta-1)+5) + mmx

            AQQl(5*(beta-1)+1:5*(beta-1)+5,5*(alpha-1)+1:5*(alpha-1)+5) = &
            AQQl(5*(beta-1)+1:5*(beta-1)+5,5*(alpha-1)+1:5*(alpha-1)+5)+mmx 
           
           !ENDDO
          !ENDDO
         enddo
      enddo



          APPl(1:3*NN,3*NN+1:6*NN)=transpose(APPl(3*NN+1:6*NN,1:3*NN))

          DO I=1,6*NN-1       
           DO J=I+1,6*NN
           APPl(J,I)=APPl(I,J)
           enddo
          enddo

      rfu=APPl*6.0_8
      rfe=transpose(AQPl)*6.0_8
      rse=AQQl*6.0_8
    write(*,*) 'ok lub-------------------------'
    end subroutine lubmxcalc

    subroutine lubcalc(r,x11a,x12a,y11a,y12a,y11b,y12b,x11c,x12c,y11c,y12c,&
     & x11g,x12g,y11g,y12g,y11h,y12h,xm,ym,zm)
    implicit none
    real(kind=cp), intent(in) :: r
    real(kind=cp), intent(out) :: x11a
    real(kind=cp), intent(out) :: x12a
    real(kind=cp), intent(out) :: y11a
    real(kind=cp), intent(out) :: y12a
    real(kind=cp), intent(out) :: y11b
    real(kind=cp), intent(out) :: y12b
    real(kind=cp), intent(out) :: x11c
    real(kind=cp), intent(out) :: x12c
    real(kind=cp), intent(out) :: y11c
    real(kind=cp), intent(out) :: y12c
    real(kind=cp), intent(out) :: x11g
    real(kind=cp), intent(out) :: x12g
    real(kind=cp), intent(out) :: y11g
    real(kind=cp), intent(out) :: y12g
    real(kind=cp), intent(out) :: y11h
    real(kind=cp), intent(out) :: y12h
    real(kind=cp), intent(out) :: xm
    real(kind=cp), intent(out) :: ym
    real(kind=cp), intent(out) :: zm

    real(kind=cp) :: xi
    real(kind=cp) :: invxi
    real(kind=cp) :: loginvxi
    real(kind=cp) :: xiloginvxi

    if (r.lt.2.1_cp) then
            xi=r-2.0_cp
            invxi=1.0_cp/xi
            loginvxi=log(invxi)
            xiloginvxi=xi*loginvxi
            
            x11a=-1.23041_cp+0.25_cp*invxi+1.8918_cp*xi+9.0_cp*loginvxi/40.0_cp & 
            &+3.0_cp*xiloginvxi/112.0_cp
            x12a=-x11a+0.00312_cp-0.0011_cp*xi

            y11a=-0.39394_cp+0.95665_cp*xi+loginvxi/6.0_cp
            y12a=-y11a+0.004636_cp-0.007049_cp*xi

            y11b=0.408286_cp-0.84055_cp*xi-loginvxi/6.0_cp-xiloginvxi/12.0_cp

            x11c=0.0479_cp+0.12494_cp*xi-xiloginvxi/6.0_cp
            x12c=-x11c+0.016869_cp+0.049536_cp*xi

            y11c=-0.605434_cp+0.939139_cp*xi+4.0_cp*loginvxi/15.0_cp & 
            &+94.0_cp*xiloginvxi/375.0_cp
            !y12c corrected, compare to the old code. 62/375
            !y12c tabulated data is the same
            y12c=loginvxi/15.0_cp+62.0_cp*xiloginvxi/375.0_cp-0.219032_cp&
            &+0.332843_cp*xi
            y12b=-0.405978_cp+0.833042_cp*xi+loginvxi/6.0_cp+xiloginvxi/12.0_cp

            x11g=-1.16897_cp+1.47882_cp*xi+invxi/4.0_cp & 
            &+9.0_cp*loginvxi/40.0_cp+39.0_cp*xiloginvxi/280.0_cp
            x12g=-x11g+0.01_cp-0.00167_cp*xi

            y11g=-0.2041_cp+0.44226_cp*xi+loginvxi/12.0_cp+xiloginvxi/24.0_cp
            y12g=-y11g+0.012265_cp-0.02757_cp*xi

            y11h=-0.143777_cp+0.264207_cp*xi+loginvxi/30.0_cp & 
            &+137.0_cp*xiloginvxi/1500.0_cp
            y12h=-0.298166_cp+0.534123_cp*xi+2.0_cp*loginvxi/15.0_cp &
            &+113.0_cp*xiloginvxi/1500.0_cp
            
            ! pay attention to m terms. here, in the equation,
            ! xm=(x11m+x12m)/2, ym=(y11m+y12m)/2, zm=(z11m+z12m)/2
            ! however, in the tabulated data, the data are 
            ! (x11m+x12m), (y11m+y12m), (z11m+z12m)
            ! so, in the tabcalc function, 
            ! xm ym zm tab values should be divided by 2
            ! this gives the same results as all the old codes

            xm=invxi/6.0_cp+3.0_cp*loginvxi/20.0_cp+47.0_cp*xiloginvxi/280.0_cp&
            & -0.740815_cp+0.706802_cp*xi
            ! the ym equation contains an extra term (xiloginvxi)
            ! compared to the old code, so the O(1), O(xi) terms are modified
            ym=-0.193518_cp+0.553883_cp*xi+loginvxi/12.0_cp&
            & -29.0_cp*xiloginvxi/500.0_cp
            ! the zm equation & tabulation are not modified
            zm = 0.00645755_cp - 0.021142_cp * xi
    else
            call lubtabcalc(r,x11a,x12a,y11a,y12a,y11b,y12b,x11c,x12c,y11c,y12c,&
     & x11g,x12g,y11g,y12g,y11h,y12h,xm,ym,zm)
    end if

    end subroutine lubcalc

    subroutine lubtabcalc(r,x11a,x12a,y11a,y12a,y11b,y12b,x11c,x12c,y11c,y12c,&
     & x11g,x12g,y11g,y12g,y11h,y12h,xm,ym,zm)
    implicit none
     
    real(kind=cp), intent(in) :: r
    real(kind=cp), intent(out) :: x11a
    real(kind=cp), intent(out) :: x12a
    real(kind=cp), intent(out) :: y11a
    real(kind=cp), intent(out) :: y12a
    real(kind=cp), intent(out) :: y11b
    real(kind=cp), intent(out) :: y12b
    real(kind=cp), intent(out) :: x11c
    real(kind=cp), intent(out) :: x12c
    real(kind=cp), intent(out) :: y11c
    real(kind=cp), intent(out) :: y12c
    real(kind=cp), intent(out) :: x11g
    real(kind=cp), intent(out) :: x12g
    real(kind=cp), intent(out) :: y11g
    real(kind=cp), intent(out) :: y12g
    real(kind=cp), intent(out) :: y11h
    real(kind=cp), intent(out) :: y12h
    real(kind=cp), intent(out) :: xm
    real(kind=cp), intent(out) :: ym
    real(kind=cp), intent(out) :: zm

    integer, parameter :: ntabc=39
    real(kind=cp),dimension(ntabc) , save :: trabc
    real(kind=cp),dimension(ntabc) , save:: tx11a
    real(kind=cp),dimension(ntabc) , save:: tx12a
    real(kind=cp),dimension(ntabc) , save:: ty11a
    real(kind=cp),dimension(ntabc) , save:: ty12a
    real(kind=cp),dimension(ntabc) , save:: ty11b
    real(kind=cp),dimension(ntabc) , save:: ty12b
    real(kind=cp),dimension(ntabc) , save:: tx11c
    real(kind=cp),dimension(ntabc) , save:: tx12c
    real(kind=cp),dimension(ntabc) , save:: ty11c
    real(kind=cp),dimension(ntabc) , save:: ty12c

    integer, parameter :: ntgh=47
    real(kind=cp),dimension(ntgh) , save:: trgh
    real(kind=cp),dimension(ntgh) , save:: tx11g
    real(kind=cp),dimension(ntgh) , save:: tx12g
    real(kind=cp),dimension(ntgh) , save:: ty11g
    real(kind=cp),dimension(ntgh) , save:: ty12g
    real(kind=cp),dimension(ntgh) , save:: ty11h
    real(kind=cp),dimension(ntgh) , save:: ty12h

    integer, parameter :: ntm=47
    real(kind=cp),dimension(ntm) , save:: trm
    real(kind=cp),dimension(ntm) , save:: txm
    real(kind=cp),dimension(ntm) , save:: tym
    real(kind=cp),dimension(ntm) , save:: tzm

    integer :: i

    integer, save :: readflag=0

    if (readflag.eq.0) then

    call readlubabc(ntabc,trabc,tx11a,tx12a,ty11a,ty12a,ty11b,ty12b,&
                    &tx11c,tx12c,ty11c,ty12c)
    call readlubgh(ntgh,trgh,tx11g,tx12g,ty11g,ty12g,ty11h,ty12h)
    call readlubm(ntm,trm,txm,tym,tzm)


    readflag = 1
    end if

    x11a=0
    x12a=0
    y11a=0
    y12a=0
    y11b=0
    y12b=0
    x11c=0
    x12c=0
    y11c=0
    y12c=0
    x11g=0
    x12g=0
    y11g=0
    y12g=0
    y11h=0
    y12h=0
    xm=0
    ym=0
    zm=0

    do i=1,ntabc-1
            if((trabc(i).le.r).and.(trabc(i+1).gt.r)) then
            x11a=tx11a(i)+(r-trabc(i))*(tx11a(i+1)-tx11a(i))/(trabc(i+1)-trabc(i))
            x12a=tx12a(i)+(r-trabc(i))*(tx12a(i+1)-tx12a(i))/(trabc(i+1)-trabc(i))
            y11a=ty11a(i)+(r-trabc(i))*(ty11a(i+1)-ty11a(i))/(trabc(i+1)-trabc(i))
            y12a=ty12a(i)+(r-trabc(i))*(ty12a(i+1)-ty12a(i))/(trabc(i+1)-trabc(i))
            y11b=ty11b(i)+(r-trabc(i))*(ty11b(i+1)-ty11b(i))/(trabc(i+1)-trabc(i))
            y12b=ty12b(i)+(r-trabc(i))*(ty12b(i+1)-ty12b(i))/(trabc(i+1)-trabc(i))
            x11c=tx11c(i)+(r-trabc(i))*(tx11c(i+1)-tx11c(i))/(trabc(i+1)-trabc(i))
            x12c=tx12c(i)+(r-trabc(i))*(tx12c(i+1)-tx12c(i))/(trabc(i+1)-trabc(i))
            y11c=ty11c(i)+(r-trabc(i))*(ty11c(i+1)-ty11c(i))/(trabc(i+1)-trabc(i))
            y12c=ty12c(i)+(r-trabc(i))*(ty12c(i+1)-ty12c(i))/(trabc(i+1)-trabc(i))
            exit
            end if
    end do

    do i=1,ntgh-1
            if((trgh(i).le.r).and.(trgh(i+1).gt.r)) then
            x11g=tx11g(i)+(r-trgh(i))*(tx11g(i+1)-tx11g(i))/(trgh(i+1)-trgh(i))
            x12g=tx12g(i)+(r-trgh(i))*(tx12g(i+1)-tx12g(i))/(trgh(i+1)-trgh(i))
            y11g=ty11g(i)+(r-trgh(i))*(ty11g(i+1)-ty11g(i))/(trgh(i+1)-trgh(i))
            y12g=ty12g(i)+(r-trgh(i))*(ty12g(i+1)-ty12g(i))/(trgh(i+1)-trgh(i))
            y11h=ty11h(i)+(r-trgh(i))*(ty11h(i+1)-ty11h(i))/(trgh(i+1)-trgh(i))
            y12h=ty12h(i)+(r-trgh(i))*(ty12h(i+1)-ty12h(i))/(trgh(i+1)-trgh(i))
            exit
            end if
    end do

    do i=1,ntm-1
            if((trm(i).le.r).and.(trm(i+1).gt.r)) then
            xm=(txm(i)+(r-trm(i))*(txm(i+1)-txm(i))/(trm(i+1)-trm(i)))/2.0_cp
            ym=(tym(i)+(r-trm(i))*(tym(i+1)-tym(i))/(trm(i+1)-trm(i)))/2.0_cp
            zm=(tzm(i)+(r-trm(i))*(tzm(i+1)-tzm(i))/(trm(i+1)-trm(i)))/2.0_cp
            exit
            end if
    end do

    end subroutine lubtabcalc

    !*****************************************************

    subroutine readlubabc(ntabc,trabc,tx11a,tx12a,ty11a,ty12a,ty11b,ty12b,&
    &tx11c,tx12c,ty11c,ty12c)
    implicit none
    integer, intent(in) :: ntabc
    real(kind=cp),dimension(ntabc),intent(out) :: trabc
    real(kind=cp),dimension(ntabc),intent(out) :: tx11a
    real(kind=cp),dimension(ntabc),intent(out) :: tx12a
    real(kind=cp),dimension(ntabc),intent(out) :: ty11a
    real(kind=cp),dimension(ntabc),intent(out) :: ty12a
    real(kind=cp),dimension(ntabc),intent(out) :: ty11b
    real(kind=cp),dimension(ntabc),intent(out) :: ty12b
    real(kind=cp),dimension(ntabc),intent(out) :: tx11c
    real(kind=cp),dimension(ntabc),intent(out) :: tx12c
    real(kind=cp),dimension(ntabc),intent(out) :: ty11c
    real(kind=cp),dimension(ntabc),intent(out) :: ty12c

    integer :: i
    character(len=80) :: line
    open(unit = 1001,file = 'lubdat/r2babc.dat',status = 'old')
    read(1001,'(a80)') line

    do i = 1, ntabc
    read( 1001, * ) trabc(i),tx11a(i),tx12a(i),ty11a(i),ty12a(i),&
            & ty11b(i),ty12b(i),tx11c(i),tx12c(i),ty11c(i),ty12c(i)
    end do

     close(unit=1001)
    end subroutine readlubabc

    subroutine readlubgh(ntgh,trgh,tx11g,tx12g,ty11g,ty12g,ty11h,ty12h)
    implicit none
    integer, intent(in) :: ntgh
    real(kind=cp),dimension(ntgh),intent(out) ::trgh
    real(kind=cp),dimension(ntgh),intent(out) ::tx11g
    real(kind=cp),dimension(ntgh),intent(out) ::tx12g
    real(kind=cp),dimension(ntgh),intent(out) ::ty11g
    real(kind=cp),dimension(ntgh),intent(out) ::ty12g
    real(kind=cp),dimension(ntgh),intent(out) ::ty11h
    real(kind=cp),dimension(ntgh),intent(out) ::ty12h
    integer :: i 
    character(len=80) :: line
    open( unit = 1002, file = 'lubdat/r2bgh.dat', status = 'old' )
    read(1002, '(a80)') line
    do i = 1, ntgh
    read( 1002, * ) trgh(i),tx11g(i),tx12g(i),ty11g(i),ty12g(i),ty11h(i),ty12h(i)
    end do

     close(unit=1002)
    end subroutine readlubgh

    subroutine readlubm(ntm,trm,txm,tym,tzm)
    implicit none
    integer, intent(in) :: ntm
    real(kind=cp),dimension(ntm),intent(out) ::trm
    real(kind=cp),dimension(ntm),intent(out) ::txm
    real(kind=cp),dimension(ntm),intent(out) ::tym
    real(kind=cp),dimension(ntm),intent(out) ::tzm
    integer :: i
    character(len=80) :: line
    open( unit = 1003, file = 'lubdat/r2bm.dat', status = 'old' )
    read(1003,'(a80)') line
    do i=1,ntm
    read(1003,*) trm(i),txm(i),tym(i),tzm(i)
    end do

     close(unit=1003)

    end subroutine readlubm

  end module hydro_lub_mod
