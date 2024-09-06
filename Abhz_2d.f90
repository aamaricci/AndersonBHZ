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
  real(8)                                   :: z2,sp_chern(2)
  !Solution
  real(8),dimension(:),allocatable          :: Ev
  complex(8),dimension(:,:),allocatable     :: H
  real(8),dimension(:,:,:),allocatable      :: Nii
  real(8),dimension(:),allocatable          :: Tzii,Szii
  complex(8),dimension(:,:,:,:),allocatable :: Gf
  integer                                   :: ilat,iorb,ispin,io


  call init_MPI()
  call init_BLACS()
  master = get_master_BLACS()


  
  !Read input:
  call parse_cmd_variable(inputFILE,"inputFILE",default="inputABHZ.conf")
  call parse_input_variable(Nx,"Nx",inputFILE,default=10)
  call parse_input_variable(Wdis,"WDIS",inputFILE,default=0d0)
  call parse_input_variable(idum,"IDUM",inputFILE,default=1234567)
  call parse_input_variable(disorder_type,"DISORDER_TYPE",inputFILE,default=0)
  call parse_input_variable(bhz_pbc,"BHZ_PBC",inputFILE,default=.true.)
  call parse_input_variable(mh,"MH",inputFILE,default=3.d0)
  call parse_input_variable(lambda,"LAMBDA",inputFILE,default=0.3d0)  
  call parse_input_variable(xmu,"XMU",inputFILE,default=0.d0)
  call parse_input_variable(beta,"BETA",inputFILE,default=1000.d0)
  call parse_input_variable(wmix,"WMIX",inputFILE,default=0.5d0)
  call parse_input_variable(Lfreq,"Lfreq",inputFILE,default=1024)
  call parse_input_variable(wmin,"WMIN",inputFILE,default=-5d0)
  call parse_input_variable(wmax,"wmax",inputFILE,default= 5d0)
  call parse_input_variable(eps,"EPS",inputFILE,default=4.d-2)
  call parse_input_variable(sb_field,"SB_FIELD",inputFILE,default=0.01d0)
  call parse_input_variable(it_error,"IT_ERROR",inputFILE,default=1d-5)
  call parse_input_variable(maxiter,"MAXITER",inputFILE,default=100)
  call parse_input_variable(with_lcm,"WITH_lcm",inputFILE,default=.false.)
  call parse_input_variable(with_mats_gf,"WITH_MATS_GF",inputFILE,default=.false.)
  call parse_input_variable(with_real_gf,"WITH_REAL_GF",inputFILE,default=.false.)
  call parse_input_variable(Nblock,"NBLOCK",inputFILE,default=4)
  call save_input_file(inputFILE)
  call print_input()
  !

  !Save variables into DMFT_TOOLS memory pool
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wmin,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(eps,"eps")


  
  !SETUP COMMON DIMENSION
  Ny   = Nx
  Nk   = Nx*Ny
  Nlat = Nx*Ny
  !
  Nso  = Nspin*Norb
  Nlso = Nlat*Nso
  Nocc = Nlso/2

  !SETUP THE GAMMA MATRICES:
  gamma0 = kron( pauli_sigma_0, pauli_tau_0)
  gammaX = kron( pauli_sigma_z, pauli_tau_x)
  gammaY = kron( pauli_sigma_0,-pauli_tau_y)
  gamma5 = kron( pauli_sigma_0, pauli_tau_z)
  gammaS = kron( pauli_sigma_z, pauli_tau_0)


  !Set the basis vectors square lattice
  call TB_set_ei([1d0,0d0],[0d0,1d0])
  call TB_set_bk([pi2,0d0],[0d0,pi2])



  !SOLVE THE HOMOGENOUS PROBLEM:  
  if(master)write(*,*) "Solve homogeneous model with Nk="//str(Nk)
  allocate(Hk(Nso,Nso,Nk))
  call TB_build_model(Hk,hk_model,Nso,[Nx,Ny])
  z2 = hk_to_spin_Chern(Hk,[Nx,Ny],spin=1)
  if(master)write(*,*)"get spin Chern UP:",z2
  z2 = hk_to_spin_Chern(Hk,[Nx,Ny],spin=2)
  if(master)write(*,*)"get spin Chern DW:",z2
  z2 = hk_to_spin_Chern(Hk,[Nx,Ny])
  if(master)write(*,*)"get Z2:",z2
  if(master)write(*,*)""
  if(master)write(*,*)""
  deallocate(Hk)


  !< Build up disorder:
  call setup_Abhz()



  !Solve the A_BHZ (non-interacting):
  call solve_Anderson_bhz()

  if(master)then
     call save_array("Ebhz.restart",Ev)
     call save_array("sz_"//str(idum)//".dat",Szii)
     call save_array("tz_"//str(idum)//".dat",Tzii)
     call save_array("n_l1s1_"//str(idum)//".dat",Nii(:,1,1))
     call save_array("n_l2s1_"//str(idum)//".dat",Nii(:,1,2))
     call save_array("n_l2s1_"//str(idum)//".dat",Nii(:,2,1))
     call save_array("n_l2s2_"//str(idum)//".dat",Nii(:,2,2))
  endif


  call push_Bloch(H,Ev)

  !Get topological info:
  sp_chern(1) = single_point_spin_chern(spin=1)
  sp_chern(2) = single_point_spin_chern(spin=2)
  if(master)call save_array("spin_chern.dat",sp_chern)
  if(master)print*,"spin_Chern UP,DW:",sp_chern
  !
  if(with_lcm)then
     call pbc_local_spin_chern_marker(spin=1,lcm=LsCM)
     if(master)call splot3d("PBC_Local_SpinChern_Marker.dat",dble(arange(1,Nx)),dble(arange(1,Ny)),LsCM)
  endif


  !< Get GF:
  if(with_mats_gf)then
     allocate(Gf(Nlat,Nso,Nso,Lfreq))
     call get_gf(Gf,'mats')
     if(master)call write_gf(gf_reshape(Gf,Nspin,Norb,Nlat),"Gloc_"//str(idum),'mats',iprint=5,itar=.true.)
     deallocate(Gf)
  endif

  if(with_real_gf)then
     allocate(Gf(Nlat,Nso,Nso,Lfreq))
     call get_gf(Gf,'real')
     if(master)call write_gf(gf_reshape(Gf,Nspin,Norb,Nlat),"Gloc_"//str(idum),'real',iprint=5,itar=.true.)
     deallocate(Gf)
  endif


  call finalize_BLACS()


contains


  subroutine solve_Anderson_BHZ()
#ifdef _SCALAPACK
    complex(8),dimension(Nlso,NLso) :: fE
#else
    real(8),dimension(Nlso)         :: fE
#endif
    complex(8),dimension(Nlso,NLso) :: Rho
    integer                         :: io,jo,ilat,iorb,ispin
    !
    if(.not.allocated(Nii))allocate(Nii(Nlat,Nspin,Norb))
    if(.not.allocated(Szii))allocate(Szii(Nlat))
    if(.not.allocated(Tzii))allocate(Tzii(Nlat))
    if(.not.allocated(H))allocate(H(Nlso,Nlso))
    if(.not.allocated(Ev))allocate(Ev(Nlso))
    !
    if(master)call start_timer("Solve Anderson BHZ")
    !
    H = Hij(:,:,1)
#ifdef _SCALAPACK
    call p_eigh(H,Ev,Nblock)
#else
    call eigh(H,Ev)
#endif
    !
    !Actual solution:
#ifdef _SCALAPACK
    fE  = one*diag(fermi(Ev,beta))
    Rho = ( H.px.fE ).px.conjg(transpose(H))
#else
    fE  = fermi(Ev,beta)
    Rho = matmul(H , matmul(diag(fE), conjg(transpose(H))) )
#endif
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
    if(master)call stop_timer()
  end subroutine solve_Anderson_BHZ


end program anderson_bhz_2d








