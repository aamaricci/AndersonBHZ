! MODEL Hamiltonian is: H + \eta*\Gamma_\eta
! H in k-space reads:
! |     h^{2x2}(k)              &            0d0              |
! |         0d0                 &        [h^{2x2}]*(-k)       |
!
!
! h^{2x2}(k):=
!
! | m-(Cos{kx}+Cos{ky})         & \lambda*(Sin{kx}-i*Sin{ky}) |
! | \lambda*(Sin{kx}+i*Sin{ky}) & -m+(Cos{kx}+Cos{ky})        |
!
!\eta is a random variable uniformly distributed within -W:W
! acting in the channel defined by the Gamma matrix \Gamma_\eta
! \eta=0 => N disorder
! \eta=1 => Tz disorder
! \eta=2 => Sz disorder
program anderson_bhz_2d
  USE COMMON
  USE LCM_SQUARE
  implicit none

  integer,parameter                         :: Norb=2,Nspin=2
  real(8)                                   :: z2,sp_chern(2),bT
  !Solution
  real(8),dimension(:),allocatable          :: Ev
  complex(8),dimension(:,:),allocatable     :: H
  real(8),dimension(:,:,:),allocatable      :: Nii
  real(8),dimension(:),allocatable          :: Tzii,Szii
  complex(8),dimension(:,:,:,:),allocatable :: Gf
  integer                                   :: ilat,iorb,ispin,io,i,ii,id
  logical                                   :: bool

  call init_parallel()  
  !  
  !Read input:
  call parse_cmd_variable(inputFILE,"inputFILE",default="inputABHZ.conf")
  call parse_cmd_variable(idumFILE,"idumFILE",default="list_idum")
  call parse_input_variable(Nx,"Nx",inputFILE,default=10)
  call parse_input_variable(Wdis,"WDIS",inputFILE,default=0d0)
  call parse_input_variable(idum,"IDUM",inputFILE,default=1234567)
  call parse_input_variable(disorder_type,"DISORDER_TYPE",inputFILE,default=0)
  call parse_input_variable(bhz_pbc,"BHZ_PBC",inputFILE,default=.true.)
  call parse_input_variable(mh,"MH",inputFILE,default=3.d0)
  call parse_input_variable(lambda,"LAMBDA",inputFILE,default=0.3d0)  
  call parse_input_variable(mu,"MU",inputFILE,default=0.d0)
  call parse_input_variable(temp,"temp",inputFILE,default=0.0001d0)
  call parse_input_variable(wmix,"WMIX",inputFILE,default=0.5d0)
  call parse_input_variable(Lfreq,"Lfreq",inputFILE,default=1024)
  call parse_input_variable(wmin,"WMIN",inputFILE,default=-5d0)
  call parse_input_variable(wmax,"wmax",inputFILE,default= 5d0)
  call parse_input_variable(ieta,"iETA",inputFILE,default=4.d-2)
  call parse_input_variable(p_field,"P_FIELD",inputFILE,default=0.01d0)
  call parse_input_variable(it_error,"IT_ERROR",inputFILE,default=1d-5)
  call parse_input_variable(maxiter,"MAXITER",inputFILE,default=100)
  call parse_input_variable(with_lcm,"WITH_lcm",inputFILE,default=.false.)
  call parse_input_variable(with_mats_gf,"WITH_MATS_GF",inputFILE,default=.false.)
  call parse_input_variable(with_real_gf,"WITH_REAL_GF",inputFILE,default=.false.)
  call parse_input_variable(Nblock,"NBLOCK",inputFILE,default=4)
  if(MPImaster)call save_input_file(inputFILE)
  if(MPImaster)call print_input()
  !

  !Save variables into DMFT_TOOLS memory pool
  bT = 1d0/temp
  call add_ctrl_var(bT,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(mu,"xmu")
  call add_ctrl_var(wmin,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(ieta,"eps")



  !SETUP COMMON DIMENSION
  Ny   = Nx
  Nk   = Nx*Ny
  Nlat = Nx*Ny
  !
  Nso  = Nspin*Norb
  Nlso = Nlat*Nso
  Nocc = Nlso/2

  !SETUP THE GAMMA MATRICES:
  call setup_GammaMatrices()


  !Set the basis vectors square lattice
  call TB_set_ei([1d0,0d0],[0d0,1d0])
  call TB_set_bk([pi2,0d0],[0d0,pi2])



  !SOLVE THE HOMOGENOUS PROBLEM:  
  if(MPImaster)write(*,*) "Solve homogeneous model with Nk="//str(Nk)
  allocate(Hk(Nso,Nso,Nk))
  call TB_build_model(Hk,hk_model,Nso,[Nx,Ny])
  z2 = hk_to_spin_Chern(Hk,[Nx,Ny],spin=1)
  if(MPImaster)write(*,*)"get spin Chern UP:",z2
  z2 = hk_to_spin_Chern(Hk,[Nx,Ny],spin=2)
  if(MPImaster)write(*,*)"get spin Chern DW:",z2
  z2 = hk_to_spin_Chern(Hk,[Nx,Ny])
  if(MPImaster)write(*,*)"get Z2:",z2
  if(MPImaster)write(*,*)""
  if(MPImaster)write(*,*)""
  deallocate(Hk)



  inquire(file=str(idumFILE),exist=bool)
  if(bool)then
     Nidum=file_length(str(idumFILE))
     allocate(list_idum(Nidum))
     open(unit=100,file=str(idumFILE))
     do i=1,Nidum
        read(100,*)id,list_idum(i)
     enddo
     close(100)
  else
     Nidum=1
     allocate(list_idum(1))
     list_idum(1)=idum
  endif


  call getcwd(here)


  do ii=1,Nidum
     idum=list_idum(ii)
     !
     dir="IDUM_"//str(idum)
     call create_dir(str(dir))
     call chdir(str(dir))
     !
     !####################################
     !< Build up disorder:
     call start_timer()
     call setup_Abhz()
     !Solve the A_BHZ (non-interacting):
     call solve_Anderson_bhz()
     if(MPImaster)then
        call save_array("sz_"//str(idum)//".dat",Szii)
        call save_array("tz_"//str(idum)//".dat",Tzii)
        call save_array("n_l1s1_"//str(idum)//".dat",Nii(:,1,1))
        call save_array("n_l1s2_"//str(idum)//".dat",Nii(:,1,2))
        call save_array("n_l2s1_"//str(idum)//".dat",Nii(:,2,1))
        call save_array("n_l2s2_"//str(idum)//".dat",Nii(:,2,2))
     endif
     call push_Bloch(H,Ev)
     !Get topological info:
     sp_chern(1) = single_point_spin_chern(spin=1)
     sp_chern(2) = single_point_spin_chern(spin=2)
     if(MPImaster)call save_array("z2.dat",(sp_chern(1)-sp_chern(2))/2d0 )
     if(MPImaster)call save_array("spin_chern.dat",sp_chern)
     if(MPImaster)print*,"spin_Chern UP,DW:",sp_chern
     !
     if(with_lcm)then
        call pbc_local_spin_chern_marker(spin=1,lcm=LsCM)
        if(MPImaster)call splot3d("PBC_Local_SpinChern_Marker.dat",dble(arange(1,Nx)),dble(arange(1,Ny)),LsCM)
     endif
     !< Get GF if required
     if(.not.allocated(Gf))allocate(Gf(Nlat,Nso,Nso,Lfreq))
     if(with_mats_gf)then
        call get_gf(Gf,'mats')
        if(MPImaster)call write_gf(gf_reshape(Gf,Nspin,Norb,Nlat),"Gloc_"//str(idum),'mats',iprint=4,itar=.true.)
     endif
     !
     if(with_real_gf)then
        call get_gf(Gf,'real')
        if(MPImaster)call write_gf(gf_reshape(Gf,Nspin,Norb,Nlat),"Gloc_"//str(idum),'real',iprint=4,itar=.true.)
     endif
     deallocate(Gf)
     !< Free memory:
     call free_Abhz()
     call stop_timer("IDUM: "//str(idum))
     if(MPImaster)write(*,*)""
     !####################################
     !
     call chdir(str(here))
  end do


  call end_parallel()


contains


  subroutine solve_Anderson_BHZ()
    complex(8),dimension(Nlso,NLso) :: fE
    complex(8),dimension(Nlso,NLso) :: Rho
    integer                         :: io,jo,ilat,iorb,ispin
    !
    if(.not.allocated(Nii))allocate(Nii(Nlat,Nspin,Norb))
    if(.not.allocated(Szii))allocate(Szii(Nlat))
    if(.not.allocated(Tzii))allocate(Tzii(Nlat))
    if(.not.allocated(H))allocate(H(Nlso,Nlso))
    if(.not.allocated(Ev))allocate(Ev(Nlso))
    !
    if(MPImaster)call start_timer("Solve Anderson BHZ")
    !
    H = Hij(:,:,1)
#ifdef _SCALAPACK
    if(MpiSize==1)then
       call eigh(H,Ev)
    else
       call p_eigh(H,Ev,Nblock)
    endif
#else
    call eigh(H,Ev)
#endif
    !
    !Actual solution:
    fE  = one*diag(fermi(Ev,1d0/temp))
    Rho = ( H.mx.fE ).mx.conjg(transpose(H))
    ! Rho = matmul(H , matmul(fE, conjg(transpose(H))) )
    !
    !Get observables:
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,iorb=1:Norb)
       io = iorb+(ispin-1)*Norb+(ilat-1)*Nspin*Norb
       Nii(ilat,ispin,iorb) = Rho(io,io)
    enddo
    do ilat=1,Nlat
       Szii(ilat) = 0.5d0*sum(Nii(ilat,1,:)) - 0.5d0*sum(Nii(ilat,2,:)) !N_up - N_dw
       Tzii(ilat) = 0.5d0*sum(Nii(ilat,:,1)) - 0.5d0*sum(Nii(ilat,:,2)) !N_1  - N_2
    enddo
    !
    if(MPImaster)call stop_timer()
  end subroutine solve_Anderson_BHZ


end program anderson_bhz_2d








