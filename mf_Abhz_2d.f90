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
  real(8),dimension(:),allocatable          :: Tzii,Szii,N1ii,N2ii
  complex(8),dimension(:,:,:,:),allocatable :: Gf
  integer                                   :: ilat,iorb,ispin,io
  logical                                   :: converged,iexist
  integer                                   :: Iter,Nsuccess=2
  real(8),dimension(:,:),allocatable        :: params,params_prev


  call init_parallel()
  !  
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
  call parse_input_variable(mu,"MU",inputFILE,default=0.d0)
  call parse_input_variable(temp,"temp",inputFILE,default=0.0001d0)
  call parse_input_variable(wmix,"WMIX",inputFILE,default=0.5d0)
  call parse_input_variable(Lfreq,"Lfreq",inputFILE,default=1024)
  call parse_input_variable(wmin,"WMIN",inputFILE,default=-5d0)
  call parse_input_variable(wmax,"wmax",inputFILE,default= 5d0)
  call parse_input_variable(ieta,"IETA",inputFILE,default=4.d-2)
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
  !
  !Save variables into DMFT_TOOLS memory pool
  call add_ctrl_var(1d0/temp,"BETA")
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


  !< Build up disorder:
  call setup_Abhz()


  !Start MF HERE:
  !>Read OR INIT MF params
  allocate(params(Nlat,2),params_prev(Nlat,2))
  do ilat=1,Nlat
     params(ilat,:)= [p_field,p_field]   ![Tz,Sz]
  enddo
  inquire(file="params.restart",exist=iexist)
  if(iexist)then
     call read_array("params.restart",params)     
     params(:,2)=params(:,2)+p_field
  endif
  !
  !> Mean-Field cycle with linear mixing
  if(MPImaster)call save_array("params.init",params)
  converged=.false. ; iter=0
  do while(.not.converged.AND.iter<maxiter)
     iter=iter+1
     call start_loop(iter,maxiter,"MF-loop")
     !
     !call symmetrize_params(params)
     call solve_MF_bhz(iter,params)
     if(MPImaster)then
        write(*,*)"E(Tz), sd(Tz):",get_mean(Tzii),get_sd(Tzii)
        write(*,*)"E(Sz), sd(Sz):",get_mean(Szii),get_sd(Szii)
        write(*,*)"E(N1), sd(N1):",get_mean(N1ii),get_sd(N1ii)
        write(*,*)"E(N2), sd(N2):",get_mean(N2ii),get_sd(N2ii)
        call save_array("Ebhz.dat",Ev)
        call save_array("tz_"//str(idum)//".dat",Tzii)
        call save_array("sz_"//str(idum)//".dat",Szii)
        call save_array("n1_"//str(idum)//".dat",N1ii)
        call save_array("n2_"//str(idum)//".dat",N2ii)
     endif
     !
     if(iter>1)params = wmix*params + (1d0-wmix)*params_prev
     params_prev = params
     !
     converged = check_convergence_local(params,it_error,1,maxiter) 
     !
     call end_loop
  end do
  if(MPImaster)call save_array("params.restart",params)


  call push_Bloch(H,Ev)

  !Get topological info:
  sp_chern(1) = single_point_spin_chern(spin=1)
  sp_chern(2) = single_point_spin_chern(spin=2)
  if(MPImaster)call save_array("spin_chern.dat",sp_chern)
  if(MPImaster)print*,"spin_Chern UP,DW:",sp_chern
  !
  if(with_lcm)then
     call pbc_local_spin_chern_marker(spin=1,lcm=LsCM)
     if(MPImaster)call splot3d("PBC_Local_SpinChern_Marker.dat",dble(arange(1,Nx)),dble(arange(1,Ny)),LsCM)
  endif


  !< Get GF:
  if(with_mats_gf)then
     allocate(Gf(Nlat,Nso,Nso,Lfreq))
     call get_gf(Gf,'mats')
     if(MPImaster)call write_gf(gf_reshape(Gf,Nspin,Norb,Nlat),"Gloc_"//str(idum),'mats',iprint=5,itar=.true.)
     deallocate(Gf)
  endif

  if(with_real_gf)then
     allocate(Gf(Nlat,Nso,Nso,Lfreq))
     call get_gf(Gf,'real')
     if(MPImaster)call write_gf(gf_reshape(Gf,Nspin,Norb,Nlat),"Gloc_"//str(idum),'real',iprint=5,itar=.true.)
     deallocate(Gf)
  endif



  call end_parallel()




contains




  subroutine solve_MF_bhz(iter,a)
    integer                                 :: iter
    real(8),dimension(Nlat,2),intent(inout) :: a
    complex(8),dimension(Nlso,NLso)         :: fE
    complex(8),dimension(Nlso,NLso)         :: Rho
    integer                                 :: io,jo,ilat,iorb,ispin
    !
    !
    if(.not.allocated(Nii))allocate(Nii(Nlat,Nspin,Norb))
    if(.not.allocated(Szii))allocate(Szii(Nlat))
    if(.not.allocated(Tzii))allocate(Tzii(Nlat))
    if(.not.allocated(N1ii))allocate(N1ii(Nlat))
    if(.not.allocated(N2ii))allocate(N2ii(Nlat))
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
       N1ii(ilat) = sum(Nii(ilat,:,1))
       N2ii(ilat) = sum(Nii(ilat,:,2))
       Tzii(ilat) = 0.5d0*(N1ii(ilat) - N2ii(ilat))
       Szii(ilat) = 0.5d0*sum(Nii(ilat,1,:)) - 0.5d0*sum(Nii(ilat,2,:)) !N_up - N_dw
    enddo
    a(:,1) = Tzii
    a(:,2) = Szii
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






