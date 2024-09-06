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
program mf_anderson_bhz_2d
  USE COMMON
  USE LCM_SQUARE
  implicit none
  integer,parameter                         :: Norb=2,Nspin=2
  real(8)                                   :: Uloc,Jh
  real(8)                                   :: z2,sp_chern(2)
  real(8),dimension(:),allocatable          :: Ev
  complex(8),dimension(:,:),allocatable     :: H
  real(8),dimension(:,:,:),allocatable      :: Nii
  real(8),dimension(:),allocatable          :: Tzii,Szii
  complex(8),dimension(:,:,:,:),allocatable :: Gf
  integer                                   :: ilat,iorb,ispin,io
  logical                                   :: converged,iexist
  integer                                   :: Iter,Nsuccess=2
  real(8),dimension(:,:),allocatable        :: params,params_prev

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
  call parse_input_variable(uloc,"ULOC",inputFILE,default=0.5d0)
  call parse_input_variable(Jh,"JH",inputFILE,default=0.1d0)
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


  !Start MF HERE:
  !>Read OR INIT MF params
  allocate(params(Nlat,2),params_prev(Nlat,2))
  do ilat=1,Nlat
     params(ilat,:)= [sb_field,sb_field]   ![Tz,Sz]
  enddo
  inquire(file="params.restart",exist=iexist)
  if(iexist)then
     call read_array("params.restart",params)     
     params(:,2)=params(:,2)+sb_field
  endif
  !
  !> Mean-Field cycle with linear mixing
  if(master)call save_array("params.init",params)
  converged=.false. ; iter=0
  do while(.not.converged.AND.iter<maxiter)
     iter=iter+1
     call start_loop(iter,maxiter,"MF-loop")
     !
     !call symmetrize_params(params)
     call solve_MF_bhz(iter,params)
     if(master)then
        call save_array("Ebhz.dat",Ev)
        call save_array("sz_"//str(idum)//".dat",Szii)
        call save_array("tz_"//str(idum)//".dat",Tzii)
        call save_array("n_l1s1_"//str(idum)//".dat",Nii(:,1,1))
        call save_array("n_l2s1_"//str(idum)//".dat",Nii(:,1,2))
        call save_array("n_l2s1_"//str(idum)//".dat",Nii(:,2,1))
        call save_array("n_l2s2_"//str(idum)//".dat",Nii(:,2,2))
     endif
     !
     if(iter>1)params = wmix*params + (1d0-wmix)*params_prev
     params_prev = params
     !
     converged = check_convergence_local(params,it_error,1,maxiter) 
     !
     call end_loop
  end do
  if(master)call save_array("params.restart",params)


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



  !call finalize_MPI()
  call finalize_BLACS()




contains




  subroutine solve_MF_bhz(iter,a)
    integer                                 :: iter
    real(8),dimension(Nlat,2),intent(inout) :: a
#ifdef _SCALAPACK
    complex(8),dimension(Nlso,NLso)         :: fE
#else
    real(8),dimension(Nlso)                 :: fE
#endif
    complex(8),dimension(Nlso,NLso)         :: Rho
    integer                                 :: io,jo,ilat,iorb,ispin
    !
    !
    if(.not.allocated(Nii))allocate(Nii(Nlat,Nspin,Norb))
    if(.not.allocated(Szii))allocate(Szii(Nlat))
    if(.not.allocated(Tzii))allocate(Tzii(Nlat))
    if(.not.allocated(H))allocate(H(Nlso,Nlso))
    if(.not.allocated(Ev))allocate(Ev(Nlso))
    !
    H = Hij(:,:,1) + mf_Hij_correction(a)
    !
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
       Tzii(ilat) = 0.5d0*sum(Nii(ilat,:,1)) - 0.5d0*sum(Nii(ilat,:,2)) !N_1  - N_2
       Szii(ilat) = 0.5d0*sum(Nii(ilat,1,:)) - 0.5d0*sum(Nii(ilat,2,:)) !N_up - N_dw
    enddo
    a(:,1) = Tzii
    a(:,2) = Szii
    if(master)call stop_timer()
  end subroutine solve_MF_bhz


  function mf_Hij_correction(a) result(HijMF)
    real(8),dimension(Nlat,2)               :: a
    complex(8),dimension(Nlat,Nso,Nso)      :: Htmp
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: HijMF
    integer                                 :: ilat,io,jo,i,j
    HijMF=zero
    do ilat=1,Nlat
       Htmp(ilat,:,:) = -a(ilat,1)*(Uloc-5d0*Jh)/2d0*Gamma5 -a(ilat,2)*(Uloc+Jh)/2d0*GammaS
       do io=1,Nso
          do jo=1,Nso
             i = io + (ilat-1)*Nso
             j = jo + (ilat-1)*Nso
             HijMF(i,j) = Htmp(ilat,io,jo)
          enddo
       enddo
    enddo
  end function mf_Hij_correction





  function get_Hmf(a,Hij) result(Hmat)
    real(8),dimension(Nlat,2)                :: a
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Hij
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Hmat
    Hmat = Hij + mf_Hij_correction(a)
  end function get_Hmf



  subroutine symmetrize_params(a)
    real(8),dimension(Nlat,Nso) :: a
    integer                     :: ilat
    do ilat=2,Nlat
       a(ilat,:) = a(1,:)
    enddo
  end subroutine symmetrize_params



end program mf_anderson_bhz_2d







! !Build up the real-space Hamiltonian thru FT:"
! allocate(Hij(Nlat*Nso,Nlat*Nso,1)) ; Hij = zero
! call start_timer
! do ilat=1,Nlat
!    vecRi = Rgrid(ilat,:)
!    do jlat=1,Nlat
!       vecRj = Rgrid(jlat,:)
!       !
!       Htmp = zero
!       do ik=1,Nktot
!          vecK = Kgrid(ik,:)
!          arg=dot_product(vecK,vecRj-vecRi)
!          Htmp(:,:)= Htmp(:,:) + exp(xi*arg)*hk_model(vecK,Nso)/Nktot
!       enddo
!       !
!       do io=1,Nso
!          i = io + (ilat-1)*Nso
!          do jo=1,Nso
!             j = jo + (jlat-1)*Nso
!             !
!             Hij(i,j,1) = Htmp(io,jo)
!             !
!          enddo
!       enddo
!       !
!    enddo
!    call eta(ilat,Nlat)
! enddo
! where(abs(Hij)<1.d-6)Hij=zero
! call stop_timer


! print*,Nlat*Nso
! allocate(Evals(Nlat*Nso))
! allocate(rhoH(Nlat*Nso,Nlat*Nso))
! call eigh(Hij(:,:,1),Evals)
! rhoDiag = fermi(Evals,beta)
! rhoH    = matmul(Hij(:,:,1) , matmul(diag(rhoDiag), conjg(transpose(Hij(:,:,1)))) )

! ilat=1
! do io=1,Nso
!    dens(io) = dreal(rhoH(io+(ilat-1)*Nso,io+(ilat-1)*Nso))
! enddo
! write(*,"(A,10F14.9)")"Occupations =",(dens(io),io=1,Nso),sum(dens)
! open(10,file="Robservables.nint")
! write(10,"(10F20.12)")(dens(io),io=1,Nso),sum(dens)
! close(10)


