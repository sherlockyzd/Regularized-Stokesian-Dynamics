module PrtclNBS_Munjiza
  use PrtclFill
  use PrtclTypeDef
  implicit none
  private
  integer::mbox,nPair,mPair
  real(RK)::xst_cs,yst_cs,zst_cs,cell_len
  type PrtclPair_info
    integer::pid,pjd
    logical::IsHaveGhost
    type(real4)::posI,posJ
  end type PrtclPair_info
  type(PrtclPair_info),dimension(:),allocatable::PrtclPair
  
  type(integer3),dimension(:),allocatable:: box_index ! integer coordinate of box
  integer,dimension(:),allocatable :: NextX
  integer,dimension(:),allocatable :: NextY
  integer,dimension(:),allocatable :: NextZ
  integer,dimension(:),allocatable :: HeadY
  integer,dimension(:),allocatable :: HeadX
  integer,dimension(:),allocatable :: HeadX0
  integer,dimension(:,:),allocatable :: HeadZ  ! head list for iy, current row 
  integer,dimension(:,:),allocatable :: HeadZ0 ! head list for (iy-1), lower row
    
  type::NBS_Munjiza
    integer:: nx         ! number of divisions in x direction
    integer:: ny         ! number of divisions in y direction
    integer:: nz         ! number of divisions in z direction
    integer:: num_Cnsv_cntct = 0 !number of conservative contacts in the broad search phase
  contains
    procedure:: Init_NBSM
    procedure:: BuildZList0
        
    ! performing contact search (includes all steps)
    procedure:: ContactSearch=> NBSM_ContactSearch
    procedure:: LoopNBSMask  => NBSM_LoopNBSMask
  end type NBS_Munjiza
  type(NBS_Munjiza),public:: m_NBS_Munjiza
  public::  NBS_Munjiza_Finalize,nPair,PrtclPair
contains
  !******************************************************************
  ! Initializing NBS_Munjiza object
  !******************************************************************
  subroutine Init_NBSM(this)
    implicit none
    class(NBS_Munjiza)::this

    ! locals
    type(integer3)::numCell
    real(RK)::xed_cs,yed_cs,zed_cs
    integer::iErr1,iErr2,iErr3,iErr4,iErr5,iErr6,iErr7,iErr8,iErr9,iErrSum

    cell_len= LenForCS
    xst_cs = simDomainMin%x -cell_len*1.05_RK
    yst_cs = simDomainMin%y -cell_len*1.05_RK
    zst_cs = simDomainMin%z -cell_len*1.05_RK
    xed_cs = simDomainMax%x +cell_len*1.05_RK
    yed_cs = simDomainMax%y +cell_len*1.05_RK 
    zed_cs = simDomainMax%z +cell_len*1.05_RK
    this%nx = int((xed_cs-xst_cs)/cell_len)+1
    this%ny = int((yed_cs-yst_cs)/cell_len)+1
    this%nz = int((zed_cs-zst_cs)/cell_len)+1
    numcell = integer3(this%nx, this%ny, this%nz)
    print*,"Contact search method is NBS Munjiza"
    print*,"Cell size is [m]: "// trim(num2str(cell_len))
    print*,"Number of cells considered is (x,y,z) :"//trim(num2str(numCell))
    
    mBox = nPrtcl+nGhostPrtcl
    allocate(box_index(mBox),   STAT=iErr1) 
    allocate(HeadY(this%ny),         STAT=iErr2)
    allocate(HeadX(this%nx),         STAT=iErr3)
    allocate(HeadX0(this%nx),        STAT=iErr4)
    allocate(HeadZ(0:1,0:this%nz+1), STAT=iErr5)
    allocate(HeadZ0(0:2,0:this%nz+1),STAT=iErr6)
    allocate(NextY(mBox),       STAT=iErr7)
    allocate(NextX(mBox),       STAT=iErr8)
    allocate(NextZ(mBox),       STAT=iErr9)
    iErrSum=abs(iErr1)+abs(iErr2)+abs(iErr3)+abs(iErr4)+abs(iErr5)+abs(iErr6)+abs(iErr7)+abs(iErr8)+abs(iErr9)
    if(iErrSum/=0) print*,"Init_NBSM: Allocation failed "
    HeadY = -1; NextY = -1
    HeadX = -1; NextX = -1; HeadX0 = -1
    HeadZ = -1; NextZ = -1; HeadZ0 = -1
    
    nPair=0; mPair=0
  end subroutine Init_NBSM

  !******************************************************************
  ! performing a contact search on all particles (includes all steps)
  !******************************************************************
  subroutine NBSM_ContactSearch(this)    
    implicit none
    class(NBS_Munjiza):: this      

    ! locals
    integer::ix,iy,iz

    this%num_Cnsv_cntct = 0
    if(nPrtcl < 1) return

    call clcBoxIndex()
    call BuildYList()

    HeadX0 = -1
    DO iy = 1, this%ny
      call BuildXList(iy)
      if(HeadY(iy).ne. -1 ) then
        HeadZ(0,:)  = -1
        HeadZ0(0,:) = -1
        call this%BuildZlist0(1,1)
        do ix = 1, this%nx
          call BuildZlist(ix,1)
          call this%BuildZlist0(ix,2)

          if(HeadX(ix).ne.-1) then
            do iz=1,this%nz
              call this%LoopNBSMask(iz)
            enddo
          endif
          HeadZ(0,:) = HeadZ(1,:)  ! same row, subs
          HeadZ0(0,:)= HeadZ0(1,:) ! lower row, subs
          HeadZ0(1,:)= HeadZ0(2,:)
        enddo
      endif
      HeadX0 = HeadX
    ENDDO
  end subroutine NBSM_ContactSearch

  !******************************************************************
  ! calculating integer coordinates of all boxes
  !******************************************************************
  subroutine clcBoxIndex()
    implicit none

    ! locals
    integer::i,m
    real(RK)::rpdx,rpdy,rpdz

    m=1
    do i=1, nPrtcl
      rpdx = GPrtcl_PosR(m)%x - xst_cs
      rpdy = GPrtcl_PosR(m)%y - yst_cs
      rpdz = GPrtcl_PosR(m)%z - zst_cs
      box_index(i)%x = floor(rpdx/cell_len)+1
      box_index(i)%y = floor(rpdy/cell_len)+1
      box_index(i)%z = floor(rpdz/cell_len)+1
      m=m+1
    enddo

    m=1
    do i=nPrtcl+1, mBox
      rpdx = GhostP_PosR(m)%x - xst_cs
      rpdy = GhostP_PosR(m)%y - yst_cs
      rpdz = GhostP_PosR(m)%z - zst_cs
      box_index(i)%x = floor(rpdx/cell_len)+1
      box_index(i)%y = floor(rpdy/cell_len)+1
      box_index(i)%z = floor(rpdz/cell_len)+1
      m=m+1    
    enddo
  end subroutine clcBoxIndex
    
  !******************************************************************
  ! constructing Ylist of particles
  !******************************************************************
  subroutine BuildYList()
    implicit none
    integer:: i,iy

    HeadY = -1  ! nullifying list Y
    
    do i=1,mBox
      iy=box_index(i)%y
      NextY(i) = HeadY(iy) 
      HeadY(iy)= i
    enddo
  end subroutine BuildYList

  !******************************************************************
  ! constructing Xlist of row iy 
  !******************************************************************
  subroutine BuildXList(iy)
    implicit none
    integer,intent(in)::iy ! row index
    integer::n,ix

    ! nullifying the xlist of current row but keeps the previous raw
    HeadX= -1
        
    n = HeadY(iy)
    do while (n .ne. -1)
      ix = box_index(n)%x
      NextX(n) = HeadX(ix)
      HeadX(ix)= n  
      n = NextY(n)
    enddo
  end subroutine BuildXList

  !*********************************************************************
  !   Constructing the ZList of column ix, ix-1, or ix+1 depending on the value of m
  !*********************************************************************
  subroutine BuildZList(ix,m)
    implicit none
    integer,intent(in) :: ix ! col index
    integer,intent(in) :: m  ! the column location with resect to ix
    integer:: n,iz
 
    ! nullifying the zlist of current col ix, but keeps the previous and next cols
    ! 0: previous (left) column, 1: current column, 2: next (right) column
    HeadZ(m,:) = -1
    n = HeadX(ix+m-1) ! reading from current row iy
    do while ( n .ne. -1 )
      iz = box_index(n)%z
      NextZ(n) = HeadZ(m,iz)
      HeadZ(m,iz) = n
      n = NextX(n)
    enddo
  end subroutine BuildZList
    
  !*********************************************************************
  !   Constructing the ZList0 of column ix, ix-1, or ix+1 depending on the value of m
  !*********************************************************************
  subroutine BuildZList0(this,ix,m)
    implicit none
    class(NBS_Munjiza):: this
    integer,intent(in):: ix ! column index
    integer,intent(in):: m  ! the column location with respect to ix
    integer:: n,iz,ixm
        
    ! 0: previous (left) column, 1: current column, 2: next (right) column
    HeadZ0(m,:) = -1
    ixm=ix+m-1
    if(ixm>this%nx) return
   
    n = HeadX0(ixm) ! reading from the row below iy (or iy-1)
    do while ( n .ne. -1 )
      iz = box_index(n)%z
      NextZ(n) = HeadZ0(m,iz)
      HeadZ0(m,iz) = n
      n = NextX(n)
    enddo
  end subroutine BuildZList0

  !**********************************************************************
  ! finding contacts between particles in the target cell and particles in
  ! cells determined by NBS mask.  
  !**********************************************************************
  subroutine NBSM_LoopNBSMask(this, iz)
    implicit none
    class(NBS_Munjiza) this
    integer,intent(in)  :: iz
    integer m, n, i, lx
  
    m = HeadZ(1,iz)
    DO WHILE(m.ne.-1)
      IF(m>nPrtcl) THEN    !==================================

        !over particles in the same cell but not the same particle (to prevent self-contact)
        n = NextZ(m)
        DO WHILE(n.ne.-1)
          if(n<nPrtcl+1) then
            call FineSearch2(n,m)
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          endif
          n = NextZ(n)
        ENDDO

        ! over particles in (ix, iy , iz-1)
        n = HeadZ(1,iz-1)
        DO WHILE (n.ne.-1)
          if(n<nPrtcl+1) then
            call FineSearch2(n,m)
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          endif
          n = NextZ(n)
        ENDDO

        ! over particles in all cells located at (ix-1) and (iy)
        do i = -1,1
          n = HeadZ(0,iz+i)
          DO WHILE(n.ne.-1)
            if(n<nPrtcl+1) then
              call FineSearch2(n,m)
              this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
            endif
            n = NextZ(n)
          ENDDO
        enddo
                        
        ! over particles in all 9 cells located at row (iy-1)
        do lx = 0,2
          do i=-1,1
            n = HeadZ0(lx,iz+i)
            DO WHILE (n.ne.-1)
              if(n<nPrtcl+1) then
                call FineSearch2(n,m)
                this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
              endif
              n = NextZ(n)
            ENDDO
          enddo
        enddo

      ELSE                             !==================================
       !over particles in the same cell but not the same particle (to prevent self-contact)
        n = NextZ(m)
        DO WHILE(n.ne.-1)
          if(n>nPrtcl) then
            call FineSearch2(m,n)
          else
            call FineSearch1(m,n)
          endif
          this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          n = NextZ(n)
        ENDDO

        ! over particles in (ix, iy , iz-1)
        n = HeadZ(1,iz-1)
        DO WHILE (n.ne.-1)
          if(n>nPrtcl) then
            call FineSearch2(m,n)
          else
            call FineSearch1(m,n)
          endif
          this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          n = NextZ(n)
        ENDDO

        ! over particles in all cells located at (ix-1) and (iy)
        do i = -1,1
          n = HeadZ(0,iz+i)
          DO WHILE(n.ne.-1)
            if(n>nPrtcl) then
              call FineSearch2(m,n)
            else
              call FineSearch1(m,n)
            endif
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
            n = NextZ(n)
          ENDDO
        enddo
                        
        ! over particles in all 9 cells located at row (iy-1)
        do lx = 0,2
          do i=-1,1
            n = HeadZ0(lx,iz+i)
            DO WHILE (n.ne.-1)
              if(n>nPrtcl) then
                call FineSearch2(m,n)
              else
                call FineSearch1(m,n)
              endif
              this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
              n = NextZ(n)
            ENDDO
          enddo
        enddo
      ENDIF                            !==================================
      m = NextZ(m)
    ENDDO
  end subroutine NBSM_LoopNBSMask

  !********************************************************************** 
  ! particle fine search(Moving-particles with Moving-particles)
  !**********************************************************************    
  subroutine FineSearch1(pid1,pid2)
    implicit none
    integer,intent(in):: pid1,pid2

    ! locals
    integer::pid1Temp,pid2Temp
    real(RK)::dx,dy,dz,dr,d2sum,dr2,drlub,ovrlp

    pid1Temp= pid1
    pid2Temp= pid2
    dr= GPrtcl_PosR(pid1Temp)%w + GPrtcl_PosR(pid2Temp)%w
    drlub= dr+ ExtraLenForCS
    dr2= drlub*drlub

    dx= GPrtcl_PosR(pid1Temp)%x - GPrtcl_PosR(pid2Temp)%x
    d2sum=dx*dx;             if(d2sum>dr2) return
    dy= GPrtcl_PosR(pid1Temp)%y - GPrtcl_PosR(pid2Temp)%y
    d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
    dz= GPrtcl_PosR(pid1Temp)%z - GPrtcl_PosR(pid2Temp)%z
    d2sum=dz*dz+d2sum;       if(d2sum>dr2) return
    ovrlp = dr-sqrt(d2sum)        

    nPair=nPair+1
    if(nPair>mPair) call reallocate_PrtclPair(nPair)
    PrtclPair(nPair)%IsHaveGhost=.false.
    if(GPrtcl_id(pid1Temp) < GPrtcl_id(pid2Temp)) then
      PrtclPair(nPair)%pid =pid1Temp
      PrtclPair(nPair)%pjd =pid2Temp
      PrtclPair(nPair)%posI=GPrtcl_PosR(pid1Temp)
      PrtclPair(nPair)%posJ=GPrtcl_PosR(pid2Temp)
    else
      PrtclPair(nPair)%pid =pid2Temp
      PrtclPair(nPair)%pjd =pid1Temp
      PrtclPair(nPair)%posI=GPrtcl_PosR(pid2Temp)
      PrtclPair(nPair)%posJ=GPrtcl_PosR(pid1Temp)  
    endif
    
    !if(ovrlp>=zero) then
      ! this is a convention, the lower id should be the first item in the contact pair (particle & particle)
      !if(GPrtcl_id(pid1Temp) < GPrtcl_id(pid2Temp) ) then
        !call GPPW_CntctList%AddContactPP(pid1Temp,pid2Temp,ovrlp)
      !else
        !call GPPW_CntctList%AddContactPP(pid2Temp,pid1Temp,ovrlp)
      !endif
    !else
      !call GPPW_CntctList%AddLubForcePP(pid1Temp,pid2Temp,-ovrlp)
    !endif
  end subroutine FineSearch1

  !********************************************************************** 
  ! particle fine search (Moving-particles with Ghost-particles)
  !**********************************************************************    
  subroutine FineSearch2(pid1,pid2)
    implicit none
    integer,intent(in):: pid1,pid2

    ! locals
    integer::gid,pid1Temp
    real(RK)::dx,dy,dz,dr,d2sum,dr2,drlub,ovrlp

    pid1Temp= pid1
    gid     = pid2-nPrtcl
    dr= GPrtcl_PosR(pid1Temp)%w + GhostP_PosR(gid)%w
    drlub= dr+ ExtraLenForCS
    dr2= drlub*drlub

    dx= GPrtcl_PosR(pid1Temp)%x - GhostP_PosR(gid)%x
    d2sum=dx*dx;             if(d2sum>dr2) return
    dy= GPrtcl_PosR(pid1Temp)%y - GhostP_PosR(gid)%y
    d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
    dz= GPrtcl_PosR(pid1Temp)%z - GhostP_PosR(gid)%z
    d2sum=dz*dz+d2sum;       if(d2sum>dr2) return
    ovrlp = dr-sqrt(d2sum)
    
    nPair=nPair+1
    if(nPair>mPair) call reallocate_PrtclPair(nPair)
    PrtclPair(nPair)%IsHaveGhost=.true.
    PrtclPair(nPair)%pid =pid1Temp
    PrtclPair(nPair)%pjd =GhostP_id(gid)
    PrtclPair(nPair)%posI=GPrtcl_PosR(pid1Temp)
    PrtclPair(nPair)%posJ=GhostP_PosR(gid)
    
    !if(ovrlp>=zero) then
      !call GPPW_CntctList%AddContactPPG(pid1Temp,gid,ovrlp)
    !else
      !call GPPW_CntctList%AddLubForcePPG(pid1Temp,gid,-ovrlp) 
    !endif
  end subroutine FineSearch2

  !********************************************************************** 
  ! NBS_Munjiza_Finalize
  !**********************************************************************   
  subroutine NBS_Munjiza_Finalize()
    implicit none
    if(allocated(box_index))deallocate(box_index)
    if(allocated(NextX))deallocate(NextX)
    if(allocated(NextY))deallocate(NextY)
    if(allocated(NextZ))deallocate(NextZ)
    if(allocated(HeadY))deallocate(HeadY)
    if(allocated(HeadX))deallocate(HeadX)
    if(allocated(HeadX0))deallocate(HeadX0)
    if(allocated(HeadZ))deallocate(HeadZ)
    if(allocated(HeadZ0))deallocate(HeadZ0)
  end subroutine NBS_Munjiza_Finalize
  
  !********************************************************************** 
  ! reallocate_PrtclPair
  !**********************************************************************
  subroutine reallocate_PrtclPair(ng)
    implicit none
    integer,intent(in)::ng
    
    ! locals
    integer::sizep,sizen
    type(PrtclPair_info),dimension(:),allocatable::PrtclPairVec
    
    sizep= mPair
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= max(sizen,ng+1)  
    mPair= sizen

    if(sizep>0) then
      call move_alloc(PrtclPair,PrtclPairVec)
      allocate(PrtclPair(sizen))
      PrtclPair(1:sizep)=PrtclPairVec
      deallocate(PrtclPairVec) 
    else
      allocate(PrtclPair(sizen))
    endif
  end subroutine reallocate_PrtclPair
end module PrtclNBS_Munjiza
