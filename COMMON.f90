MODULE COMMON
#ifdef _SCALAPACK
  USE SCIFOR, operator(.mx.) => operator(.px.)
  USE DMFT_TOOLS
  USE MPI
#else
  USE SCIFOR, operator(.mx.) => operator(.x.)
  USE DMFT_TOOLS
#endif
  implicit none

  integer                                   :: Nso=0
  integer                                   :: Nlat=0
  integer                                   :: Nx=0
  integer                                   :: Ny=0
  integer                                   :: Nlso=0
  integer                                   :: Nocc=0
  integer                                   :: Nblock=0
  !
  integer                                   :: Nk
  integer                                   :: Nkpath,Npts,Nidum
  complex(8),dimension(4,4)                 :: Gamma0
  complex(8),dimension(4,4)                 :: Gamma5
  complex(8),dimension(4,4)                 :: GammaX
  complex(8),dimension(4,4)                 :: GammaY
  complex(8),dimension(4,4)                 :: GammaS
  !Model parameters:
  character(len=20)                         :: inputFILE
  character(len=20)                         :: idumFILE
  real(8)                                   :: Wdis
  integer                                   :: Idum
  integer                                   :: disorder_type
  logical                                   :: bhz_pbc
  real(8)                                   :: mh
  real(8)                                   :: lambda
  real(8)                                   :: wmix
  integer                                   :: Lfreq
  real(8)                                   :: wmin,wmax
  real(8)                                   :: temp
  real(8)                                   :: mu,ieta
  real(8)                                   :: p_field
  real(8)                                   :: it_error
  integer                                   :: MaxIter
  logical                                   :: with_lcm
  logical                                   :: with_mats_gf
  logical                                   :: with_real_gf
  !Random 
  real(8),allocatable,dimension(:)          :: erandom
  integer,allocatable,dimension(:)          :: list_idum
  character(len=255)                        :: here
  character(len=20)                         :: dir
  !Hamiltonian
  complex(8),dimension(:,:,:),allocatable   :: Hk
  integer,dimension(:,:),allocatable        :: Links
  complex(8),dimension(:,:,:),allocatable   :: Hij
  complex(8),dimension(:,:,:,:),allocatable :: Hlat
  complex(8),dimension(:,:),allocatable     :: Sz
  real(8),dimension(:,:),allocatable        :: LsCM
  !
  complex(8),dimension(:,:),allocatable     :: U
  complex(8),dimension(:,:),allocatable     :: PSzP
  real(8),dimension(:),allocatable          :: E
  real(8),dimension(:),allocatable          :: Epsp
  !MPI:
  logical                                   :: MpiStatus=.false.
  logical                                   :: MpiMaster=.true.
  integer                                   :: MpiRank=0
  integer                                   :: MpiSize=1
  integer                                   :: Mpi_Ierr



contains


  !SETUP THE GAMMA MATRICES:
  subroutine setup_GammaMatrices()
    gamma0 = kron( pauli_sigma_0, pauli_tau_0)
    gammaX = kron( pauli_sigma_z, pauli_tau_x)
    gammaY = kron( pauli_sigma_0,-pauli_tau_y)
    gamma5 = kron( pauli_sigma_0, pauli_tau_z)
    gammaS = kron( pauli_sigma_z, pauli_tau_0)
  end subroutine setup_GammaMatrices



  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: ek
    real(8)                   :: kx,ky
    complex(8),dimension(N,N) :: hk
    if(N/=4)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    ek = -1d0*(cos(kx)+cos(ky))
    Hk = (Mh+ek)*Gamma5 + lambda*sin(kx)*GammaX + lambda*sin(ky)*GammaY
  end function hk_model



  function ts_model(link,Nso) result(Hts)
    integer                       :: link
    integer                       :: Nso
    complex(8),dimension(Nso,Nso) :: Hts
    select case(link)
    case (0) !LOCAL PART
       Hts =  Mh*Gamma5
    case (1) !RIGHT HOPPING
       Hts = -0.5d0*Gamma5 + xi*0.5d0*lambda*GammaX
    case (2) !UP HOPPING
       Hts = -0.5d0*Gamma5 + xi*0.5d0*lambda*GammaY
    case (3) !LEFT HOPPING
       Hts = -0.5d0*Gamma5 - xi*0.5d0*lambda*GammaX
    case (4) !DOWN HOPPING
       Hts = -0.5d0*Gamma5 - xi*0.5d0*lambda*GammaY
    case default 
       stop "ts_model ERROR: link != {0,...,4}"
    end select
  end function ts_model





  subroutine setup_Abhz(verbose)
    logical,optional :: verbose
    logical :: verbose_
    logical :: bool
    integer :: ilat
    verbose_=.true.;if(present(verbose))verbose_=verbose    
    if(allocated(Hlat))deallocate(Hlat)
    allocate(Hlat(Nso,Nso,Nlat,Nlat))
    if(allocated(Hij))deallocate(Hij)
    allocate(Hij(Nlso,Nlso,1))
    if(allocated(erandom))deallocate(erandom)
    allocate(erandom(Nlat))
    call check_dimension("setup_disorder")
    !
    if(verbose_.AND.MPImaster)write(*,*)
    if(verbose_.AND.MPImaster)write(*,*)"Solve IDUM:",idum
    if(verbose_.AND.MPImaster)write(*,*)"Nlso      :",Nlso
    !
    if(allocated(Links))deallocate(Links)
    allocate(Links(4,2))
    Links(1,:) = [1 ,0]
    Links(2,:) = [0 ,1]
    Links(3,:) = [-1,0]
    Links(4,:) = [0,-1]
    call TB_build_model(Hlat,ts_model,Nso,[Nx,Nx],Links,pbc=bhz_pbc)
    !
    call mersenne_init(idum)
    call mt_random(erandom)
    erandom=(2d0*erandom-1d0)*Wdis/2d0
    !
    inquire(file='erandom_'//str(idum)//'.restart',exist=bool)
    if(bool)then
       if(file_length('erandom_'//str(idum)//'.restart')/=Nlat)&
            stop "setup_disorder error: size(erandom_"//str(idum)//".restart) != Nlat"
       call read_array('erandom_'//str(idum)//'.restart',erandom)
    endif
    if(verbose_.AND.MPImaster)call save_array('erandom_'//str(idum)//'.used',erandom)
    !
    do ilat=1,Nlat
       select case(disorder_type)
       case default;stop "setup_disorder error: disorder_type [0:2]"
       case(0)                  !N
          Hlat(:,:,ilat,ilat) = Hlat(:,:,ilat,ilat) + erandom(ilat)*Gamma0
       case(1)                  !Tz
          Hlat(:,:,ilat,ilat) = Hlat(:,:,ilat,ilat) + erandom(ilat)*Gamma5
       case(2)                  !Sz
          Hlat(:,:,ilat,ilat) = Hlat(:,:,ilat,ilat) + erandom(ilat)*GammaS
       end select
    enddo
    !
    !Reshape:
    Hij(:,:,1) = reshape_rank4_to_rank2(Hlat,Nso,Nlat)
    !
    if(allocated(Sz))deallocate(Sz)
    allocate(Sz(Nlso,Nlso))
    Sz = kron(eye(Nlat),GammaS)
    !
    !
  end subroutine setup_Abhz


  subroutine free_Abhz
    if(allocated(Hlat))deallocate(Hlat)
    if(allocated(Hij))deallocate(Hij)
    if(allocated(erandom))deallocate(erandom)
    if(allocated(Links))deallocate(Links)
    if(allocated(Sz))deallocate(Sz)
    if(allocated(U))deallocate(U)
    if(allocated(E))deallocate(E)
    if(allocated(PSzP))deallocate(PSzP)
    if(allocated(Epsp))deallocate(Epsp)
  end subroutine free_Abhz

    
  

  !------------------------------------------------------------------
  ! Checks whether all the dimensions of the systems have been set
  ! properly in the calling procedure.
  !------------------------------------------------------------------
  subroutine check_dimension(caller)
    character(len=*) :: caller
#ifdef _SCALAPACK
    if(Nblock==0)stop str(caller)//" error: Nblock not set"
#endif
    if(Nlso==0)stop str(caller)//" error: Nlso not set"
    if(Nso==0)stop str(caller)//" error: Nso not set"
    if(Nx==0)stop str(caller)//" error: Nx not set"
    if(Ny==0)stop str(caller)//" error: Ny not set"
    if(Nlat==0)stop str(caller)//" error: Nlat not set"
    if(Nocc==0)stop str(caller)//" error: Nocc not set"
    if(Nlat/=Nx*Ny)stop str(caller)//" error: Nlat != Nx*Ny"
    if(.not.allocated(Hlat)) stop str(caller)//" error: Hlat not allocated"
    if(allocated(Hlat))call assert_shape(Hlat,[Nso,Nso,Nlat,Nlat],"check_dimension","Hlat")
    if(allocated(Hij))call assert_shape(Hij,[Nlso,Nlso,1],"check_dimension","Hij")
    if(allocated(U))call assert_shape(U,[Nlso,Nlso],"check_dimension","U")
    if(allocated(E))call assert_shape(E,[Nlso],"check_dimension","E")
    if(allocated(Sz))call assert_shape(Sz,[Nlso,Nlso],"check_dimension","Sz")
    if(allocated(PSzP))call assert_shape(PSzP,[Nocc,Nocc],"check_dimension","PSzP")
    if(allocated(Epsp))call assert_shape(Epsp,[Nocc],"check_dimension","Epsp")
  end subroutine check_dimension




  !------------------------------------------------------------------
  ! Push Bloch states and levels to internal arrays U, E
  !------------------------------------------------------------------
  subroutine push_Bloch(Evec,Evals)
    complex(8),dimension(:,:) :: Evec
    real(8),dimension(:)      :: Evals
    call check_dimension("push_Bloch")
    if(allocated(U).AND.allocated(E))stop "push_Bloch warning: U, E allocated"
    call assert_shape(Evec,[Nlso,Nlso],"push_Bloch","Evec")
    call assert_shape(Evals,[Nlso]    ,"push_Bloch","Evals")
    if(allocated(U))deallocate(U)
    if(allocated(E))deallocate(E)
    allocate(U, source=Evec)
    allocate(E, source=Evals)
    if(MPImaster)call save_array("E.dat",E)
    if(MPImaster)call check_Egap("push_Bloch")
  end subroutine push_Bloch




  subroutine check_Pgap(caller)
    character(len=*) :: caller
    integer          :: N
    real(8)          :: Ep,Em,Pgap
    N    = int(Nocc/2)
    Ep   = Epsp(N+1)
    Em   = Epsp(N)
    Pgap = Ep - Em
    if(MPImaster)write(*,"(A,G21.13)")"PSzP gap=",Pgap
    !
    if(Pgap<1d-12)then
       print*, str(caller)//" warning: closing of the PSzP spectrum"
    elseif(Ep*Em>0d0)then
       print*, str(caller)//" warning: PSzP spectrum not symmetric"
    else
       return
    endif
  end subroutine check_Pgap

  


  subroutine check_Egap(caller)
    character(len=*) :: caller
    integer          :: N
    real(8)          :: Ep,Em,Egap
    N    = Nocc
    Ep   = E(N+1)
    Em   = E(N)
    Egap = Ep - Em
    if(MPImaster)write(*,"(A,G21.13)")"E gap=",Egap
    !
    if(Egap<1d-12)then
       print*,str(caller)//" warning: closing of the H spectrum"
    elseif(Ep*Em>0d0)then
       print*,str(caller)//" warning: H spectrum not symmetric"
    else
       return
    endif
  end subroutine check_Egap






  subroutine get_gf(Gloc,axis)
    complex(8),dimension(Nlat,Nso,Nso,Lfreq) :: Gloc
    character(len=*)                         :: axis
    complex(8),dimension(Lfreq)              :: wfreq
    complex(8),dimension(Nlso,Lfreq)         :: csi
    complex(8),dimension(Nlat,Nso,Nso,Lfreq) :: Gtmp
    integer                                  :: i,ilat,io,jo,is,js
    !
    call check_dimension("get_gf")
    !
    if(MPImaster)write(*,"(A)")"Get local Green's function, axis:"//str(axis)
    wfreq = build_frequency_array(axis)
    !
    !Allocate and setup the Matsubara freq.
    forall(i=1:Lfreq)csi(:,i)=one/(wfreq(i)+mu-E(:))
    !
    Gloc = zero
    if(MPImaster)call start_timer
    do i=1+MPIrank,Lfreq,MPIsize
       do concurrent(ilat=1:Nlat,io=1:Nso,jo=1:Nso)
          is = io + (ilat-1)*Nso
          js = jo + (ilat-1)*Nso
          Gloc(ilat,io,jo,i) = sum(U(is,:)*conjg(U(js,:))*csi(:,i))!can use matmul
       enddo
    enddo
#ifdef _SCALAPACK
    call MPI_AllReduce(MPI_IN_PLACE,Gloc,size(Gloc),MPI_Double_Complex,MPI_Sum,MPI_COMM_WORLD,mpi_ierr)
#endif
    if(MPImaster)call stop_timer
  end subroutine get_gf



  function build_frequency_array(axis) result(wfreq)
    character(len=*)                    :: axis
    complex(8),dimension(:),allocatable :: wfreq
    if(allocated(wfreq))deallocate(wfreq)
    allocate(wfreq(Lfreq))
    !
    select case(to_lower(axis))
    case default;
       stop "build_frequency_array ERROR: axis undefined. axis=[matsubara,realaxis]"
    case("matsubara","mats","m")
       wfreq = dcmplx(0d0,pi*temp*(2*arange(1,Lfreq)-1))
    case("realaxis","real","r")
       wfreq = dcmplx(linspace(wmin,wmax,Lfreq),ieta)
    end select
    return
  end function build_frequency_array




  function reshape_rank4_to_rank2(MatNN,Nso,Nlat) result(Kmat)
    integer                                 :: Nso,Nlat
    complex(8),dimension(Nso,Nso,Nlat,Nlat) :: MatNN
    complex(8),dimension(Nso*Nlat,Nso*Nlat) :: Kmat
    integer                                 :: iso,jso,i,j,ii,jj
    do concurrent(i=1:Nlat,j=1:Nlat,iso=1:Nso,jso=1:Nso)
       ii = iso + (i-1)*Nso
       jj = jso + (j-1)*Nso
       Kmat(ii,jj) = MatNN(iso,jso,i,j)
    enddo
  end function reshape_rank4_to_rank2

  function reshape_rank2_to_rank4(Kmat,Nso,Nlat) result(MatNN)
    integer                                 :: Nso,Nlat
    complex(8),dimension(Nso*Nlat,Nso*Nlat) :: Kmat
    complex(8),dimension(Nso,Nso,Nlat,Nlat) :: MatNN
    integer                                 :: iso,jso,i,j,ii,jj
    do concurrent(i=1:Nlat,j=1:Nlat,iso=1:Nso,jso=1:Nso)
       ii = iso + (i-1)*Nso
       jj = jso + (j-1)*Nso
       MatNN(iso,jso,i,j)  =  Kmat(ii,jj)
    enddo
  end function reshape_rank2_to_rank4



  function reshape_rank3_to_rank5(MatNN,Nspin,Norb,Nlat) result(Kmat)
    integer                                          :: Nspin,Norb,Nlat
    complex(8),dimension(Nlat,Nspin*Norb,Nspin*Norb) :: MatNN
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: Kmat
    integer                                          :: ilat
    integer                                          :: ispin,jspin
    integer                                          :: iorb,jorb
    integer                                          :: io,jo
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       io = iorb + (ispin-1)*Norb
       jo = jorb + (jspin-1)*Norb
       Kmat(ilat,ispin,jspin,iorb,jorb) = MatNN(ilat,io,jo)
    enddo
  end function reshape_rank3_to_rank5

  function reshape_rank5_to_rank3(Kmat,Nspin,Norb,Nlat) result(MatNN)
    integer                                          :: Nspin,Norb,Nlat
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: Kmat
    complex(8),dimension(Nlat,Nspin*Norb,Nspin*Norb) :: MatNN
    integer                                          :: ilat
    integer                                          :: ispin,jspin
    integer                                          :: iorb,jorb
    integer                                          :: io,jo
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       io = iorb + (ispin-1)*Norb
       jo = jorb + (jspin-1)*Norb
       MatNN(ilat,io,jo) = Kmat(ilat,ispin,jspin,iorb,jorb)
    enddo
  end function reshape_rank5_to_rank3



  function gf_reshape(MatNN,Nspin,Norb,Nlat) result(Kmat)
    integer                                                :: Nspin,Norb,Nlat
    complex(8),dimension(Nlat,Nspin*Norb,Nspin*Norb,Lfreq) :: MatNN
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,Lfreq) :: Kmat
    integer                                                :: i
    do i=1,Lfreq
       Kmat(:,:,:,:,:,i) = reshape_rank3_to_rank5(MatNN(:,:,:,i),Nspin,Norb,Nlat)
    enddo
  end function gf_reshape



  subroutine init_parallel
#ifdef _SCALAPACK
    call init_BLACS()           !includes init_MPI
    MPImaster = get_master_BLACS()
    MPIrank   = get_rank_BLACS()
    MPISize   = get_size_BLACS()
#endif
  end subroutine init_parallel


  subroutine end_parallel
#ifdef _SCALAPACK
    call finalize_BLACS()
    MPImaster = .false.
    MPIrank   = 0
    MPISize   = 1
#endif
  end subroutine end_parallel

END MODULE COMMON
