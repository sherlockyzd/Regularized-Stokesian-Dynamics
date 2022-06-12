! WRITES THE CURRENT CONFIGURATION ON FILE

  PROGRAM MAIN
  use tensors

  implicit none
  integer::i,iLb,NLB,Nlambda,j
  REAL*8 :: grmobmx(11*2,11*2)
  REAL*8 :: Resist(11*2,11*2)
  REAL*8 ::p_pos(3,2),RADII(2),ei(3)
  LOGICAL::THERE
  REAL*8, ALLOCATABLE :: lambda(:),Lb0_list(:)

  write(*,*) 'start-------------'
  CALL INIT_KDirac
  CALL INIT_EPS
  call INIT_Y2()
  !write(*,*) 'Y2=',Y2
  call INIT_II()



  INQUIRE( FILE='lambda.dat', EXIST=THERE ) 
  IF ( THERE ) THEN
    OPEN(10,FILE='lambda.dat')
    READ(10,*) Nlambda
    ALLOCATE( lambda(Nlambda) )
    do i=1,Nlambda
    READ(10,*) lambda(i)
    enddo
    CLOSE(10)
  ELSE
    WRITE(0,*) "error: no control_file.dat file "  &
    // "present in this folder"
    WRITE(0,*) " Please provide one"
    CALL EXIT()
  endif


  NLB=60
  ALLOCATE( Lb0_list(NLB) )
  Lb0_list(1:10)=(/(i,i=1,10)/)*1.0d0*1e-5
  Lb0_list(11:15)=(/(i,i=1,5)/)*2.0d0*1e-4
  Lb0_list(16:20)=(/(i,i=1,5)/)*2.0d0*1e-3
  Lb0_list(21:25)=(/(i,i=1,5)/)*2.0d0*1e-2
  Lb0_list(26:50)=(/(i,i=1,25)/)*0.1d0+0.1d0
  Lb0_list(51:60)=(/(i,i=1,10)/)*0.5d0+2.6d0
  !Lb0_list(36:45)=(/(i,i=1,10)/)*0.05d0+1.d0
  Lb0_list=Lb0_list+2.0_8
  !write(*,*) Lb0_list

  open(20,file='tablefar.txt')
  write(20,'(24A24)') 'lambda','s^','X11A','X12A','Y11A','Y12A','Y11B','Y12B',  &
              & 'X11C','X12C','Y11C','Y12C','X11G','X12G', & 
              & 'Y11G','Y12G','Y11H','Y12H','X11M','X12M', &
              & 'Y11M','Y12M','Z11M','Z12M'
  ei=0.0_8
  ei(1)=1.0_8
  call Init_L(ei)
  !Nlambda=1
  lam_Cal: do i=1,Nlambda
  !NLB=1
    LB_size: do iLb =1,NLB !100

    RADII=1.0_8
    RADII(2)=lambda(i)

    p_pos=0.0_8
    p_pos(1,2)=Lb0_list(iLb)*(RADII(1)+RADII(2))/2.0_8
    
    do j=1,2
    write(*,*) 'conf=',p_pos(:,j)
    enddo
    write(*,*) 'RADII=',RADII(:)
    CALL GRPERY_MOB(grmobmx,p_pos,RADII)

    Resist=grmobmx
    call MATREV(Resist,22)
    call FarfieldScalar(Resist)
    write(20,"(24ES24.15)") lambda(i),Lb0_list(iLb)-2.0_8,X11A,X12A,Y11A,Y12A,Y11B,Y12B,&
    X11C,X12C,Y11C,Y12C,X11G,X12G,Y11G,Y12G,Y11H,Y12H,X11M,X12M,Y11M,Y12M,Z11M,Z12M

    
    


    enddo LB_size
  enddo lam_Cal

  DEALLOCATE(LB0_list)  
  DEALLOCATE(lambda)

  CLOSE(20) ! xyz file for vmd
  write(*,*)'end mob!!'

    
  END

