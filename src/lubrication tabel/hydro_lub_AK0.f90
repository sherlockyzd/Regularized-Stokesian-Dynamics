!******************************************************
!******************************************************
!******************************************************

subroutine lubmxcalc(ntotal,P_pos,rfu,rfe,rse)
use prutil
use TENSORS,only:Y21,Y22,Y2,front
USE method
use CONFIG, only:RADII
implicit none
integer, intent(in) :: ntotal
real(kind=cp), dimension(3,ntotal), intent(in) :: P_pos
real(kind=cp), dimension(6*ntotal,6*ntotal), intent(out) :: rfu
real(kind=cp), dimension(6*ntotal,5*ntotal), intent(out) :: rfe
real(kind=cp), dimension(5*ntotal,5*ntotal), intent(out) :: rse

real(kind=cp) :: r,rx(3)
real(kind=cp) :: x11a,x12a,y11a,y12a,y11b,y12b,x11c,x12c,y11c,y12c
real(kind=cp) :: x11g,x12g,y11g,y12g,y11h,y12h,x11m,x12m,y11m,y12m,z11m,z12m
real(kind=cp), dimension(3) :: e
real(kind=cp), dimension(2,2,3,3) :: alub
real(kind=cp), dimension(2,2,3,3) :: blub
real(kind=cp), dimension(2,2,3,3) :: club
real(kind=cp), dimension(2,2,3,3,3) :: glubtensor
real(kind=cp), dimension(2,2,3,3,3) :: hlubtensor

real(kind=cp), dimension(2,2,5,3) :: glubtensor_Y2,hlubtensor_Y2

real(kind=cp), dimension(2,2,3,3,3,3) :: mlubtensor
real(kind=cp), dimension(2,2,5,5) :: mmx
!real(kind=cp), dimension(ntotal,3) :: p_pos

real(kind=cp), dimension(6*ntotal,6*ntotal) :: APPl
real(kind=cp), dimension(5*ntotal,6*ntotal) :: AQPl
!real(kind=cp), dimension(5*ntotal,6*ntotal) :: APQl
real(kind=cp), dimension(5*ntotal,5*ntotal) :: AQQl
!real(kind=cp), dimension(3*ntotal,3*ntotal) :: APPl_B

integer :: alpha,beta,NN
integer :: i,j,k,l,m,n

real*8 aTemp,V0(3,3,3),V1(5,3),V2(5,5),r_ave,s


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
        rx=p_pos(:,alpha)-p_pos(:,beta)
        r=sqrt(sum(rx*rx))
      !if(IsPeriod) then
      !  CALL PER_SKEW_CORRECTION(rx,r)
      !endif
        alub=0.0_8
        blub=0.0_8
        club=0.0_8
        glubtensor=0.0_8
        hlubtensor=0.0_8
        mlubtensor=0.0_8
        r_ave=0.5_8*(radii(alpha)+radii(beta))
        s=r/r_ave      
        if(s.lt.2.0_8) then
                write(*,*) 'Contact Error, r='
                write(*,*) s,alpha,beta
                !stop
                cycle
        else if(s.gt.4.5_8) then
                cycle                
        else
        ! initialize the lub2bmx matrix
        ! build the 2body lubrication resistance matrix lub2bmx(22*22) 
        call lubtabcalc(s,x11a,x12a,y11a,y12a,y11b,y12b,x11c,x12c,y11c,y12c,&
              & x11g,x12g,y11g,y12g,y11h,y12h,x11m,x12m,y11m,y12m,z11m,z12m)


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
        mlubtensor(1,1,i,j,k,l)=x11m*(3.0_8*(e(i)*e(j)-kd(i,j)/3.0_8)*(e(k)*e(l)&
        & -kd(k,l)/3.0_8)/2.0_8)+y11m*(e(i)*kd(j,l)*e(k)+e(j)*kd(i,l)*e(k)&
        & +e(i)*kd(j,k)*e(l)+e(j)*kd(i,k)*e(l)-4.0_8*e(i)*e(j)*e(k)*e(l))&
        & /2.0_8+z11m*(kd(i,k)*kd(j,l)+kd(j,k)*kd(i,l)-kd(i,j)*kd(k,l)& 
        & +e(i)*e(j)*kd(k,l)+kd(i,j)*e(k)*e(l)-e(i)*kd(j,l)*e(k)&
        & -e(j)*kd(i,l)*e(k)-e(i)*kd(j,k)*e(l)-e(j)*kd(i,k)*e(l)&
        & +e(i)*e(j)*e(k)*e(l))/2.0_8

        mlubtensor(1,2,i,j,k,l)=x12m*(3.0_8*(e(i)*e(j)-kd(i,j)/3.0_8)*(e(k)*e(l)&
        & -kd(k,l)/3.0_8)/2.0_8)+y12m*(e(i)*kd(j,l)*e(k)+e(j)*kd(i,l)*e(k)&
        & +e(i)*kd(j,k)*e(l)+e(j)*kd(i,k)*e(l)-4.0_8*e(i)*e(j)*e(k)*e(l))&
        & /2.0_8+z12m*(kd(i,k)*kd(j,l)+kd(j,k)*kd(i,l)-kd(i,j)*kd(k,l)& 
        & +e(i)*e(j)*kd(k,l)+kd(i,j)*e(k)*e(l)-e(i)*kd(j,l)*e(k)&
        & -e(j)*kd(i,l)*e(k)-e(i)*kd(j,k)*e(l)-e(j)*kd(i,k)*e(l)&
        & +e(i)*e(j)*e(k)*e(l))/2.0_8
        end forall
      
        forall(i=1:3,j=1:3)
        alub(2,2,i,j)=alub(1,1,i,j) !x22a=x11a,y22a=y11a
        alub(2,1,i,j)=alub(1,2,i,j) !x12a=x21a,y12a=y21a
        blub(2,2,i,j)= blub(1,1,i,j) !y11b=y22b
        blub(2,1,i,j)=-blub(1,2,i,j) !y12b=-y21b
        club(2,2,i,j)=club(1,1,i,j) !x11c=x22c
        club(2,1,i,j)=club(1,2,i,j) !x21c=x12c
        end forall

        forall(i=1:3,j=1:3,k=1:3)
        glubtensor(2,2,i,j,k)= glubtensor(1,1,i,j,k) ! x22g=-x11g, y22g=-y11g
        glubtensor(2,1,i,j,k)=-glubtensor(1,2,i,j,k) ! x21g=-x12g, y21g=-y12g
        hlubtensor(2,2,i,j,k)=hlubtensor(1,1,i,j,k) !  y22h=y11h 
        hlubtensor(2,1,i,j,k)=hlubtensor(1,2,i,j,k) ! y21h=y12h
        end forall



        forall(i=1:3,j=1:3,k=1:3,l=1:3)
        mlubtensor(2,2,i,j,k,l)=mlubtensor(1,1,i,j,k,l) ! x22g=-x11g, y22g=-y11g
        mlubtensor(2,1,i,j,k,l)=mlubtensor(1,2,i,j,k,l) ! x21g=-x12g, y21g=-y12g
        end forall

       ! call DD_transY2(mmx,mlubtensor)



       do m=1,2
        do n=1,2
           
          V0=glubtensor(m,n,1:3,1:3,1:3)
          call RTD_front_transY2(V1,V0)
          glubtensor_Y2(m,n,1:5,1:3)=V1
           
          V0=hlubtensor(m,n,1:3,1:3,1:3)
          call RTD_front_transY2(V1,V0)
          hlubtensor_Y2(m,n,1:5,1:3)=V1
        enddo
      enddo
           
       do m=1,2
        do n=1,2

          call DD_transY2(V2,mlubtensor(m,n,1:3,1:3,1:3,1:3))
          mmx(m,n,1:5,1:5)=V2

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

       
!#ifdef fronttest
       AQPl(5*(alpha-1)+1:5*(alpha-1)+5,3*(alpha-1)+1:3*(alpha-1)+3) = &
       AQPl(5*(alpha-1)+1:5*(alpha-1)+5,3*(alpha-1)+1:3*(alpha-1)+3) + glubtensor_Y2(1,1,1:5,1:3)
       AQPl(5*(beta-1)+1:5*(beta-1)+5,3*(beta-1)+1:3*(beta-1)+3) = &
       AQPl(5*(beta-1)+1:5*(beta-1)+5,3*(beta-1)+1:3*(beta-1)+3) + glubtensor_Y2(2,2,1:5,1:3)

       AQPl(5*(alpha-1)+1:5*(alpha-1)+5,3*NN+3*(alpha-1)+1:3*NN+3*(alpha-1)+3) = &
       AQPl(5*(alpha-1)+1:5*(alpha-1)+5,3*NN+3*(alpha-1)+1:3*NN+3*(alpha-1)+3) + hlubtensor_Y2(1,1,1:5,1:3)
       AQPl(5*(beta-1)+1:5*(beta-1)+5,3*NN+3*(beta-1)+1:3*NN+3*(beta-1)+3) = &
       AQPl(5*(beta-1)+1:5*(beta-1)+5,3*NN+3*(beta-1)+1:3*NN+3*(beta-1)+3) + hlubtensor_Y2(2,2,1:5,1:3)
!#endif

#ifdef backttest
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
       AQQl(5*(alpha-1)+1:5*(alpha-1)+5,5*(alpha-1)+1:5*(alpha-1)+5)+mmx(1,1,1:5,1:5)
       AQQl(5*(beta-1)+1:5*(beta-1)+5,5*(beta-1)+1:5*(beta-1)+5) = &
       AQQl(5*(beta-1)+1:5*(beta-1)+5,5*(beta-1)+1:5*(beta-1)+5)+mmx(2,2,1:5,1:5)
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


!#ifdef fronttest  
       AQPl(5*(alpha-1)+1:5*(alpha-1)+5,3*(beta-1)+1:3*(beta-1)+3) = &
       AQPl(5*(alpha-1)+1:5*(alpha-1)+5,3*(beta-1)+1:3*(beta-1)+3) + glubtensor_Y2(1,2,1:5,1:3)
       AQPl(5*(beta-1)+1:5*(beta-1)+5,3*(alpha-1)+1:3*(alpha-1)+3) = &
       AQPl(5*(beta-1)+1:5*(beta-1)+5,3*(alpha-1)+1:3*(alpha-1)+3) + glubtensor_Y2(2,1,1:5,1:3)

       AQPl(5*(alpha-1)+1:5*(alpha-1)+5,3*NN+3*(beta-1)+1:3*NN+3*(beta-1)+3) = &
       AQPl(5*(alpha-1)+1:5*(alpha-1)+5,3*NN+3*(beta-1)+1:3*NN+3*(beta-1)+3) + hlubtensor_Y2(1,2,1:5,1:3)
       AQPl(5*(beta-1)+1:5*(beta-1)+5,3*NN+3*(alpha-1)+1:3*NN+3*(alpha-1)+3) = &
       AQPl(5*(beta-1)+1:5*(beta-1)+5,3*NN+3*(alpha-1)+1:3*NN+3*(alpha-1)+3) + hlubtensor_Y2(2,1,1:5,1:3)
!#endif
     
#ifdef backtest  
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
        AQQl(5*(alpha-1)+1:5*(alpha-1)+5,5*(beta-1)+1:5*(beta-1)+5) + mmx(1,2,1:5,1:5)

        AQQl(5*(beta-1)+1:5*(beta-1)+5,5*(alpha-1)+1:5*(alpha-1)+5) = &
        AQQl(5*(beta-1)+1:5*(beta-1)+5,5*(alpha-1)+1:5*(alpha-1)+5)+mmx(2,1,1:5,1:5)
       
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



subroutine lubtabcalc(r,x11a,x12a,y11a,y12a,y11b,y12b,x11c,x12c,y11c,y12c,&
 & x11g,x12g,y11g,y12g,y11h,y12h,x11m,x12m,y11m,y12m,z11m,z12m)
use prutil
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
real(kind=cp), intent(out) :: x11m
real(kind=cp), intent(out) :: y11m
real(kind=cp), intent(out) :: z11m
real(kind=cp), intent(out) :: x12m
real(kind=cp), intent(out) :: y12m
real(kind=cp), intent(out) :: z12m

integer, parameter :: nt=38
real(kind=cp),dimension(nt) , save :: tr
real(kind=cp),dimension(nt) , save:: tx11a
real(kind=cp),dimension(nt) , save:: tx12a
real(kind=cp),dimension(nt) , save:: ty11a
real(kind=cp),dimension(nt) , save:: ty12a
real(kind=cp),dimension(nt) , save:: ty11b
real(kind=cp),dimension(nt) , save:: ty12b
real(kind=cp),dimension(nt) , save:: tx11c
real(kind=cp),dimension(nt) , save:: tx12c
real(kind=cp),dimension(nt) , save:: ty11c
real(kind=cp),dimension(nt) , save:: ty12c
real(kind=cp),dimension(nt) , save:: tx11g
real(kind=cp),dimension(nt) , save:: tx12g
real(kind=cp),dimension(nt) , save:: ty11g
real(kind=cp),dimension(nt) , save:: ty12g
real(kind=cp),dimension(nt) , save:: ty11h
real(kind=cp),dimension(nt) , save:: ty12h
real(kind=cp),dimension(nt) , save:: tx11m
real(kind=cp),dimension(nt) , save:: tx12m
real(kind=cp),dimension(nt) , save:: ty11m
real(kind=cp),dimension(nt) , save:: ty12m
real(kind=cp),dimension(nt) , save:: tz11m
real(kind=cp),dimension(nt) , save:: tz12m

integer :: i
real(kind=cp):: r_tr


integer, save :: readflag=0
if (readflag.eq.0) then
call readlub(nt,tr,tx11a,tx12a,ty11a,ty12a,ty11b,ty12b,&
          & tx11c,tx12c,ty11c,ty12c,tx11g,tx12g,ty11g,ty12g, &
          & ty11h,ty12h,tx11m,tx12m,ty11m,ty12m,tz11m,tz12m)
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
x11m=0
x12m=0
y11m=0
y12m=0
z11m=0
z12m=0

do i=1,nt-1
        if((tr(i).le.r).and.(tr(i+1).gt.r)) then
        r_tr=(r-tr(i))/(tr(i+1)-tr(i))    
        x11a=tx11a(i)+r_tr*(tx11a(i+1)-tx11a(i))
        x12a=tx12a(i)+r_tr*(tx12a(i+1)-tx12a(i))

        y11a=ty11a(i)+r_tr*(ty11a(i+1)-ty11a(i))
        y12a=ty12a(i)+r_tr*(ty12a(i+1)-ty12a(i))
 
        y11b=ty11b(i)+r_tr*(ty11b(i+1)-ty11b(i))
        y12b=ty12b(i)+r_tr*(ty12b(i+1)-ty12b(i))

        x11c=tx11c(i)+r_tr*(tx11c(i+1)-tx11c(i))
        x12c=tx12c(i)+r_tr*(tx12c(i+1)-tx12c(i))

        y11c=ty11c(i)+r_tr*(ty11c(i+1)-ty11c(i))
        y12c=ty12c(i)+r_tr*(ty12c(i+1)-ty12c(i))

        x11g=tx11g(i)+r_tr*(tx11g(i+1)-tx11g(i))
        x12g=tx12g(i)+r_tr*(tx12g(i+1)-tx12g(i))

        y11g=ty11g(i)+r_tr*(ty11g(i+1)-ty11g(i))
        y12g=ty12g(i)+r_tr*(ty12g(i+1)-ty12g(i))

        y11h=ty11h(i)+r_tr*(ty11h(i+1)-ty11h(i))
        y12h=ty12h(i)+r_tr*(ty12h(i+1)-ty12h(i))

        x11m=tx11m(i)+r_tr*(tx11m(i+1)-tx11m(i))
        x12m=tx12m(i)+r_tr*(tx12m(i+1)-tx12m(i))

        y11m=ty11m(i)+r_tr*(ty11m(i+1)-ty11m(i))
        y12m=ty12m(i)+r_tr*(ty12m(i+1)-ty12m(i))

        z11m=tz11m(i)+r_tr*(tz11m(i+1)-tz11m(i))
        z12m=tz12m(i)+r_tr*(tz12m(i+1)-tz12m(i))
        exit
        end if
end do
end subroutine lubtabcalc

!*****************************************************

subroutine readlub(nt,tr,tx11a,tx12a,ty11a,ty12a,ty11b,ty12b,&
& tx11c,tx12c,ty11c,ty12c,tx11g,tx12g,ty11g,ty12g, &
          & ty11h,ty12h,tx11m,tx12m,ty11m,ty12m,tz11m,tz12m)
!readlub(nt,tr,tx11a,tx12a,ty11a,ty12a,ty11b,ty12b,&
!          & tx11c,tx12c,ty11c,ty12c,tx11g,tx12g,ty11g,ty12g, &
!         & ty11h,ty12h,tx11m,tx12m,ty11m,ty12m,tz11m,tz12m)
use prutil
implicit none
integer, intent(in) :: nt
real(kind=cp),dimension(nt),intent(out) :: tr
real(kind=cp),dimension(nt),intent(out) :: tx11a
real(kind=cp),dimension(nt),intent(out) :: tx12a
real(kind=cp),dimension(nt),intent(out) :: ty11a
real(kind=cp),dimension(nt),intent(out) :: ty12a
real(kind=cp),dimension(nt),intent(out) :: ty11b
real(kind=cp),dimension(nt),intent(out) :: ty12b
real(kind=cp),dimension(nt),intent(out) :: tx11c
real(kind=cp),dimension(nt),intent(out) :: tx12c
real(kind=cp),dimension(nt),intent(out) :: ty11c
real(kind=cp),dimension(nt),intent(out) :: ty12c
real(kind=cp),dimension(nt),intent(out) :: tx11g
real(kind=cp),dimension(nt),intent(out) :: tx12g
real(kind=cp),dimension(nt),intent(out) :: ty11g
real(kind=cp),dimension(nt),intent(out) :: ty12g
real(kind=cp),dimension(nt),intent(out) :: ty11h
real(kind=cp),dimension(nt),intent(out) :: ty12h
real(kind=cp),dimension(nt),intent(out) :: tx11m
real(kind=cp),dimension(nt),intent(out) :: tx12m
real(kind=cp),dimension(nt),intent(out) :: ty11m
real(kind=cp),dimension(nt),intent(out) :: ty12m
real(kind=cp),dimension(nt),intent(out) :: tz11m
real(kind=cp),dimension(nt),intent(out) :: tz12m

integer :: i
real(kind=cp):: trtest,lamdatest,gamatest
character(len=80) :: line
open(unit = 1001,file = 'lubdat/scalars_general_resistance_text_d.txt',status = 'old')
read(1001,'(a80)') line
read(1001,'(a80)') line

do i = 1, nt
read( 1001, * ) tr(i),lamdatest,gamatest,tx11a(i),ty11a(i),&
            & ty11b(i),tx11c(i),ty11c(i), &
            & tx11g(i),ty11g(i),ty11h(i), &
            & tx11m(i),ty11m(i),tz11m(i)

read( 1001, * ) trtest,lamdatest,gamatest,tx12a(i),ty12a(i),&
            & ty12b(i),tx12c(i),ty12c(i), &
            & tx12g(i),ty12g(i),ty12h(i), &
            & tx12m(i),ty12m(i),tz12m(i)
end do

 close(unit=1001)
end subroutine readlub


