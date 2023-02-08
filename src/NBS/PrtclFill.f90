module PrtclFill
  use PrtclTypeDef
  implicit none
  public
  integer::nPrtcl=0,nGhostPrtcl=0
  real(RK)::LenForCS,ExtraLenForCS
  type(real3)::simDomainMin,simDomainMax,simLen
  integer,dimension(:),allocatable::GPrtcl_id
  integer,dimension(:),allocatable::GhostP_id
  type(real4),dimension(:),allocatable::GPrtcl_PosR
  type(real4),dimension(:),allocatable::GhostP_PosR

  integer,private::mGhostPrtcl=0
  real(RK),private:: xst0_cs,xst1_cs,xst2_cs
  real(RK),private:: xed0_cs,xed1_cs,xed2_cs
  real(RK),private:: yst0_cs,yst1_cs,yst2_cs
  real(RK),private:: yed0_cs,yed1_cs,yed2_cs
  real(RK),private:: zst0_cs,zst1_cs,zst2_cs
  real(RK),private:: zed0_cs,zed1_cs,zed2_cs
contains

  !**********************************************************************
  ! Prtcl_Fill_Init
  !**********************************************************************
  subroutine Prtcl_Fill_Init(nP,PrtclPosR,simDomMin,simDomMax,ExtraLenCS,zpShift,zmShift)
    implicit none
    integer,intent(in)::nP
    real(RK),dimension(4,nP),intent(in)::PrtclPosR
    real(RK),intent(in)::ExtraLenCS,zpShift,zmShift
    real(RK),dimension(3),intent(in)::simDomMin,simDomMax
    
    ! locals
    integer::pid,ngp
    real(RK)::px,py,pz,maxDiam=0.0_RK
    integer,dimension(:),allocatable:: IntVec
    type(real4),dimension(:),allocatable::Real4Vec
    
    nPrtcl=nP
    if(nPrtcl<1) return
    simDomainMin=simDomMin
    simDomainMax=simDomMax
    simLen=simDomainMax-simDomainMin
    allocate(GPrtcl_PosR(nPrtcl),GPrtcl_id(nPrtcl))
    do pid=1,nPrtcl
      GPrtcl_id(pid)=pid
      GPrtcl_PosR(pid)%x=PrtclPosR(1,pid)
      GPrtcl_PosR(pid)%y=PrtclPosR(2,pid)      
      GPrtcl_PosR(pid)%z=PrtclPosR(3,pid)
      GPrtcl_PosR(pid)%w=PrtclPosR(4,pid)
      if(GPrtcl_PosR(pid)%w>maxDiam)maxDiam=GPrtcl_PosR(pid)%w
    enddo
    maxDiam=maxDiam*2.0_RK
    ExtraLenForCS=ExtraLenCS
    LenForCS=maxDiam+ExtraLenCS
    mGhostPrtcl=nPrtcl
    allocate(GhostP_PosR(mGhostPrtcl),GhostP_id(mGhostPrtcl))
    
    if(1.02_RK*LenForCS>min(min(simLen%x,simLen%y),simLen%z) ) then
      print*,'Prtcl_Fill_Init: so big Diameter+ExtraLenCS'
      stop
    endif
    xst0_cs = simDomainMin%x-LenForCS
    xst1_cs = simDomainMin%x
    xst2_cs = simDomainMin%x+LenForCS
    xed0_cs = simDomainMax%x-LenForCS
    xed1_cs = simDomainMax%x
    xed2_cs = simDomainMax%x+LenForCS

    yst0_cs = simDomainMin%y-LenForCS
    yst1_cs = simDomainMin%y
    yst2_cs = simDomainMin%y+LenForCS    
    yed0_cs = simDomainMax%y-LenForCS
    yed1_cs = simDomainMax%y
    yed2_cs = simDomainMax%y+LenForCS 

    zst0_cs = simDomainMin%z-LenForCS
    zst1_cs = simDomainMin%z
    zst2_cs = simDomainMin%z+LenForCS    
    zed0_cs = simDomainMax%z-LenForCS
    zed1_cs = simDomainMax%z
    zed2_cs = simDomainMax%z+LenForCS
    
    nGhostPrtcl=0
    ! step1: Handle y-dir
    do pid=1,nPrtcl
      py=GPrtcl_PosR(nGhostPrtcl)%y
      if(py<=yst2_cs) then
        nGhostPrtcl=nGhostPrtcl+1
        if(nGhostPrtcl>mGhostPrtcl) call reallocate_ghost_for_Cntct(nGhostPrtcl)
        GhostP_id(nGhostPrtcl)  = GPrtcl_id(pid)
        GhostP_PosR(nGhostPrtcl)= GPrtcl_PosR(pid)
        GhostP_PosR(nGhostPrtcl)%y = py+simLen%y
      endif
      if(py>=yed0_cs) then
        nGhostPrtcl=nGhostPrtcl+1
        if(nGhostPrtcl>mGhostPrtcl) call reallocate_ghost_for_Cntct(nGhostPrtcl)
        GhostP_id(nGhostPrtcl)  = GPrtcl_id(pid)
        GhostP_PosR(nGhostPrtcl)= GPrtcl_PosR(pid)
        GhostP_PosR(nGhostPrtcl)%y = py-simLen%y
      endif
    enddo
    
    ! step2: send to xp_axis, and receive from xm_dir
    ngp=nGhostPrtcl
    do pid=1,ngp    ! consider the previous ghost particle firstly
      px=GhostP_PosR(pid)%x
      if(px>=xed0_cs ) then
        nGhostPrtcl=nGhostPrtcl+1
        if(nGhostPrtcl>mGhostPrtcl) call reallocate_ghost_for_Cntct(nGhostPrtcl)
        GhostP_id(nGhostPrtcl)     = GhostP_id(pid)
        GhostP_PosR(nGhostPrtcl)   = GhostP_PosR(pid)
        GhostP_PosR(nGhostPrtcl)%x = px-simLen%x
      endif
    enddo
    do pid=1,nPrtcl
      px=GPrtcl_PosR(nGhostPrtcl)%x
      if(px>=xed0_cs) then
        nGhostPrtcl=nGhostPrtcl+1
        if(nGhostPrtcl>mGhostPrtcl) call reallocate_ghost_for_Cntct(nGhostPrtcl)
        GhostP_id(nGhostPrtcl)     = GPrtcl_id(pid)
        GhostP_PosR(nGhostPrtcl)   = GPrtcl_PosR(nGhostPrtcl)
        GhostP_PosR(nGhostPrtcl)%x = px-simLen%x
      endif
    enddo

    ! step3: send to xm_axis, and receive from xp_dir
    ngp=nGhostPrtcl
    do pid=1,ngp    ! consider the previous ghost particle firstly
      px=GhostP_PosR(pid)%x
      if(px<=xst2_cs .and. px>=xst1_cs) then
        nGhostPrtcl=nGhostPrtcl+1
        if(nGhostPrtcl>mGhostPrtcl) call reallocate_ghost_for_Cntct(nGhostPrtcl)
        GhostP_id(nGhostPrtcl)     = GhostP_id(pid)
        GhostP_PosR(nGhostPrtcl)   = GhostP_PosR(pid)
        GhostP_PosR(nGhostPrtcl)%x = px+simLen%x
      endif
    enddo
    do pid=1,nPrtcl
      px=GPrtcl_PosR(nGhostPrtcl)%x
      if(px<=xst2_cs .and. px>=xst1_cs) then
        nGhostPrtcl=nGhostPrtcl+1
        if(nGhostPrtcl>mGhostPrtcl) call reallocate_ghost_for_Cntct(nGhostPrtcl)
        GhostP_id(nGhostPrtcl)     = GPrtcl_id(pid)
        GhostP_PosR(nGhostPrtcl)   = GPrtcl_PosR(nGhostPrtcl)
        GhostP_PosR(nGhostPrtcl)%x = px+simLen%x
      endif
    enddo

    ! step4: send to zp_axis, and receive from zm_dir
    ngp=nGhostPrtcl
    do pid=1,ngp    ! consider the previous ghost particle firstly
      pz=GhostP_PosR(pid)%z
      if(pz>=zed0_cs ) then
        nGhostPrtcl=nGhostPrtcl+1
        if(nGhostPrtcl>mGhostPrtcl) call reallocate_ghost_for_Cntct(nGhostPrtcl)
        GhostP_id(nGhostPrtcl)     = GhostP_id(pid)
        GhostP_PosR(nGhostPrtcl)   = GhostP_PosR(pid)
        GhostP_PosR(nGhostPrtcl)%z = pz-simLen%z
      endif
    enddo
    do pid=1,nPrtcl
      pz=GPrtcl_PosR(nGhostPrtcl)%z
      if(pz>=zed0_cs) then
        nGhostPrtcl=nGhostPrtcl+1
        if(nGhostPrtcl>mGhostPrtcl) call reallocate_ghost_for_Cntct(nGhostPrtcl)
        GhostP_id(nGhostPrtcl)     = GPrtcl_id(pid)
        GhostP_PosR(nGhostPrtcl)   = GPrtcl_PosR(nGhostPrtcl)
        GhostP_PosR(nGhostPrtcl)%z = pz-simLen%z
      endif
    enddo  

    ! step5: send to zm_axis, and receive from zp_dir
    ngp=nGhostPrtcl
    do pid=1,ngp    ! consider the previous ghost particle firstly
      pz=GhostP_PosR(pid)%z
      if(pz<=zst2_cs .and. pz>=zst1_cs) then
        nGhostPrtcl=nGhostPrtcl+1
        if(nGhostPrtcl>mGhostPrtcl) call reallocate_ghost_for_Cntct(nGhostPrtcl)
        GhostP_id(nGhostPrtcl)     = GhostP_id(pid)
        GhostP_PosR(nGhostPrtcl)   = GhostP_PosR(pid)
        GhostP_PosR(nGhostPrtcl)%z = pz+simLen%z
      endif
    enddo
    do pid=1,nPrtcl
      pz=GPrtcl_PosR(nGhostPrtcl)%z
      if(pz<=zst2_cs .and. pz>=zst1_cs) then
        nGhostPrtcl=nGhostPrtcl+1
        if(nGhostPrtcl>mGhostPrtcl) call reallocate_ghost_for_Cntct(nGhostPrtcl)
        GhostP_id(nGhostPrtcl)     = GPrtcl_id(pid)
        GhostP_PosR(nGhostPrtcl)   = GPrtcl_PosR(nGhostPrtcl)
        GhostP_PosR(nGhostPrtcl)%z = pz+simLen%z
      endif
    enddo
    
    if(nGhostPrtcl>0) then
      call move_alloc(GhostP_id,IntVec)
      allocate(GhostP_id(nGhostPrtcl))
      GhostP_id=IntVec(1:nGhostPrtcl)
      deallocate(IntVec)
      
      call move_alloc(GhostP_PosR,Real4Vec)
      allocate(GhostP_PosR(nGhostPrtcl))
      GhostP_PosR=Real4Vec(1:nGhostPrtcl)
      deallocate(Real4Vec)
    else
      deallocate(GhostP_id,GhostP_PosR)
    endif
  end subroutine Prtcl_Fill_Init

  !**********************************************************************
  ! reallocate_ghost_for_Cntct
  !**********************************************************************  
  subroutine reallocate_ghost_for_Cntct(ng)
    implicit none
    integer,intent(in)::ng

    ! locals
    integer:: sizep,sizen
    integer,dimension(:),allocatable:: IntVec
    type(real4),dimension(:),allocatable::Real4Vec

    sizep= mGhostPrtcl
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= min(sizen,nPrtcl)
    sizen= max(sizen,ng+1)
    mGhostPrtcl=sizen  ! NOTE HERE, sometimes mGhostPrtcl CAN bigger than DEM_Opt%numPrtcl

    ! ======= integer vector part =======
    call move_alloc(GhostP_id,IntVec)
    allocate(GhostP_id(sizen))
    GhostP_id(1:sizep)=IntVec
    deallocate(IntVec)

    ! ======= real4 vercor part =======
    call move_alloc(GhostP_PosR,Real4Vec)
    allocate(GhostP_PosR(sizen))
    GhostP_PosR(1:sizep)=Real4Vec
    deallocate(Real4Vec)
  end subroutine reallocate_ghost_for_Cntct
  
  !**********************************************************************
  ! Prtcl_Fill_Finalize
  !**********************************************************************
  subroutine Prtcl_Fill_Finalize()
    implicit none
    if(allocated(GPrtcl_id))  deallocate(GPrtcl_id)
    if(allocated(GhostP_id))  deallocate(GhostP_id)
    if(allocated(GPrtcl_PosR))deallocate(GPrtcl_PosR)
    if(allocated(GhostP_PosR))deallocate(GhostP_PosR)
  end subroutine Prtcl_Fill_Finalize
end module PrtclFill
