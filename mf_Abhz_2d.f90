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
  real(8)                                   :: z2,sp_chern(2),bT,error,error_prev
  !Solution
  real(8),dimension(:),allocatable          :: Ev
  complex(8),dimension(:,:),allocatable     :: H
  real(8),dimension(:,:,:),allocatable      :: Nii
  real(8),dimension(:),allocatable          :: Tzii,Szii,N1ii,N2ii
  complex(8),dimension(:,:,:,:),allocatable :: Gf
  integer                                   :: ilat,iorb,ispin,io,ii,id
  logical                                   :: converged,iexist,ioidum,post_processing
  integer                                   :: Iter,Nsuccess=2,ix,iy
  real(8),dimension(:,:),allocatable        :: params,params_prev
  character(len=100) :: pfile

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
  call parse_input_variable(with_lcm,"WITH_LCM",inputFILE,default=.false.)
  call parse_input_variable(with_real_gf,"WITH_REAL_GF",inputFILE,default=.false.)
  call parse_input_variable(Nblock,"NBLOCK",inputFILE,default=4)
  if(MPImaster)call save_input_file(inputFILE)
  if(MPImaster)call print_input()
  !
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
  z2 = hk_to_spin_Chern(Hk,[Nx,Ny])
  if(MPImaster)write(*,*)"get Z2:",z2
  if(MPImaster)write(*,*)""
  if(MPImaster)write(*,*)""
  deallocate(Hk)

  post_processing=with_real_gf.OR.with_lcm


  inquire(file=str(idumFILE),exist=iexist)
  if(iexist)then
     Nidum = file_length(str(idumFILE))
     allocate(list_idum(Nidum))
     open(unit=100,file=str(idumFILE))
     do ii=1,Nidum
        read(100,*)id,list_idum(ii)
     enddo
     close(100)
  else
     Nidum = 1
     allocate(list_idum(Nidum))
     list_idum = idum
  endif

  call getcwd(here)



  !####################################
  do ii=1,Nidum
     idum  = list_idum(ii)
     dir   = "IDUM_"//str(idum)
     pfile = "params_"//str(idum)
     !
     !Create IDUM directory and ENTER it     
     if(MpiMaster)call system("mkdir -p "//str(dir))
     if(MpiMaster)call system("cp -fv "//str(pfile)//".restart "//str(dir)//"/")
     call chdir(str(dir))
     !
     !< Build up disorder:
     if(MpiMaster)call start_timer()
     call setup_Abhz()
     !
     !
     if(post_processing)then
        !< POST-PROCESSING: Get GF or LCM 
        call post_process_MFbhz()
     else
        !>MF Solution:
        !read params
        allocate(params(Nlat,2),params_prev(Nlat,2))
        if(MPImaster)then
           do ilat=1,Nlat
              params(ilat,:)= [p_field,p_field*afm_sign(ilat)]
           enddo
           inquire(file=str(pfile)//".restart",exist=iexist)
           if(iexist)then
              call read_array(str(pfile)//".restart",params)     
              params(:,2)=params(:,2)+p_field*afm_sign(ilat)
           endif
        endif
        call Bcast_MPI(MPI_COMM_WORLD,params)
        !
        !> Mean-Field cycle with linear mixing
        if(MPImaster)call save_array(str(pfile)//".init",params)
        converged=.false. ; iter=0
        error=1d0
        error_prev=error
        do while(.not.converged.AND.iter<maxiter)
           iter=iter+1
           call start_loop(iter,maxiter,"MF-loop")
           !
           call solve_MF_Abhz(iter,params)           
           !
           if(iter>1)params = wmix*params + (1d0-wmix)*params_prev
           params_prev = params
           !
           converged = check_convergence_local(params,it_error,1,maxiter,oerr=error)
           if(error > error_prev)then
              wmix=wmix*(0.95d0) !decrease wmix by 5%
              if(MPImaster)write(*,*)"Reduced wmix:",wmix
           endif
           error_prev=error
           call end_loop
        end do
        if(MPImaster)call save_array(str(pfile)//".restart",params)
        if(MpiMaster)call system("cp -fv "//str(pfile)//".restart "//str(here)//"/")
        !
        call push_Bloch(H,Ev)
        !
        !Get topological info:
        sp_chern(1) = single_point_spin_chern(spin=1)
        sp_chern(2) = single_point_spin_chern(spin=2)
        if(MPImaster)call save_array("z2.dat",(sp_chern(1)-sp_chern(2))/2d0 )
        if(MPImaster)call save_array("spin_chern.dat",sp_chern)
        if(MpiMaster)call save_array("Ebhz.dat",Ev)
        if(MpiMaster)call save_array("tz.dat",Tzii)
        if(MpiMaster)call save_array("sz.dat",Szii)
        if(MpiMaster)call save_array("afm.dat",afm_vec())
        if(MpiMaster)call save_array("Etz.dat",get_mean(Tzii))
        if(MpiMaster)call save_array("Esz.dat",get_AFMmean(Szii))
        if(MPImaster)print*,"spin_Chern UP,DW:",sp_chern
     endif
     !
     !< Free memory:
     call chdir(str(here))
     call free_Abhz()
     deallocate(params,params_prev)
     if(MpiMaster)call stop_timer("IDUM: "//str(idum))
     if(MPImaster)write(*,*)""
  enddo
  !####################################


  call end_parallel()


contains


  subroutine solve_MF_Abhz(iter,a)
    integer                                 :: iter
    real(8),dimension(Nlat,2),intent(inout) :: a
    complex(8),dimension(Nlso,NLso)         :: fE
    complex(8),dimension(Nlso,NLso)         :: Rho
    integer                                 :: io,jo,ilat,iorb,ispin
    real(8)                                 :: eSz
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
    if(MPImaster)then
       write(*,*)"E(Tz):",get_mean(Tzii)
       write(*,*)"E(Sz):",get_AFMmean(Szii)
    endif
  end subroutine solve_MF_Abhz




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



  subroutine post_process_MFbhz()
    !>Read MF params
    allocate(params(Nlat,2),params_prev(Nlat,2))
    if(MPImaster)then
       inquire(file=str(pfile)//".restart",exist=iexist)
       if(.not.iexist)stop "Can not read params.restart file: can not do post-processing"
       call read_array(str(pfile)//".restart",params)     
    endif
    call Bcast_MPI(MPI_COMM_WORLD,params)
    !
    !Solve MF problem once with restart parameters
    call solve_MF_Abhz(1,params)
    call push_Bloch(H,Ev)
    !
    if(with_real_gf)then
       if(.not.allocated(Gf))allocate(Gf(Nlat,Nso,Nso,Lfreq))
       call get_gf(Gf,'real')
       if(MPImaster)call write_gf(gf_reshape(Gf,Nspin,Norb,Nlat),"Gloc",'real',iprint=4,itar=.true.)
       if(allocated(Gf))deallocate(Gf)
    endif
    !
    if(with_lcm)then
       if(bhz_pbc)then
          call pbc_local_spin_chern_marker(spin=1,lcm=LsCM)
          if(MPImaster)call splot3d("PBC_Local_SpinChern_Marker.dat",dble(arange(1,Nx)),dble(arange(1,Ny)),LsCM)
       else
          call obc_local_spin_chern_marker(spin=1,lcm=LsCM)
          if(MpiMaster)call splot3d("OBC_Local_SpinChern_Marker.dat",dble(arange(1,Nx)),dble(arange(1,Ny)),LsCM)
       endif
    endif
  end subroutine post_process_MFbhz



end program mf_anderson_bhz_2d






