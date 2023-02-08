module DEM_gz
  use MPI
  use Prtcl_TypeDef
  use Prtcl_decomp_2d
  use Prtcl_DEMSystem
  use Prtcl_Parameters
  use Prtcl_IOAndVisu
  use Prtcl_LogInfo
  use Prtcl_Variables
#ifdef DEM_DEBUG
  use Prtcl_CL_and_CF
#endif
  implicit none
  private 
  
  logical::IsInitial=.false.
  public::DEM_calc_initial_position
contains

subroutine DEM_calc_initial_position(CONF,LB,L_mag,RADII,NN,itimeDEM)
  implicit none
  integer,intent(in)::NN,itimeDEM
  real(RK),intent(in)::LB(3),RADII(NN),L_mag
  real(RK),intent(inout)::CONF(3,NN)
  
  ! locals
  character(len=64)::chDEMPrm
  integer :: ierror,i
  
  if(.not.IsInitial) then
   call MPI_INIT(ierror)
   call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
   IsInitial=.true.
  endif

  ! read DEM options
  if(command_argument_count()/=1 .and. nrank==0) write(*,*)'command argument wrong!'
  call get_command_argument(1,chDEMPrm)
  call DEM_opt%ReadDEMOption(chDEMPrm)
  
  !==================== 2022-01-13 Zheng Gong
  DEM_opt%SimDomain_min = zero_r3
  DEM_opt%SimDomain_max = real3(LB(1),LB(2),LB(3))*1.0D-3
  DEM_opt%numPrtcl =NN
  !==================== 2022-01-13 Zheng Gong
  
  call DEM_decomp%Init_DECOMP(chDEMPrm)
  call DEM%Initialize(chDEMPrm,CONF,RADII,L_mag) ! Topest level initialing for DEM body
#ifdef DEM_DEBUG
  call GPPW_CntctList%printCL(DEM_opt%ifirst-1)
#endif

  !call DEM_IO%dump_visu(DEM_opt%ifirst-1)
  print*, nrank,GPrtcl_list%nlocal,GPrtcl_list%mlocalFix
  do i= 1,itimeDEM !DEM_opt%ifirst, DEM_opt%ilast
    call DEM%iterate(i)
  enddo
  !call DEM_IO%Final_visu()


  do i=1,NN
    CONF(1,i)=GPrtcl_PosR(i)%x*1.0D+3
    CONF(2,i)=GPrtcl_PosR(i)%y*1.0D+3
    CONF(3,i)=GPrtcl_PosR(i)%z*1.0D+3
  enddo

  call DEM%Finalize()
  if(nrank==0)call DEMLogInfo%OutInfo("Good job! DEM finished successfully!",1)
  if(nrank==0)call DEMLogInfo%CloseFile()
  !call MPI_FINALIZE(ierror)
end subroutine DEM_calc_initial_position
end module DEM_gz
