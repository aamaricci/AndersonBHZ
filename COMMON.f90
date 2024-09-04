MODULE COMMON
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none


  integer                                   :: Nso=0
  integer                                   :: Nlat=0
  integer                                   :: Nx=0
  integer                                   :: Ny=0
  integer                                   :: Nlso=0
  integer                                   :: Nocc=0
  !
  integer                                   :: Nk
  integer                                   :: Nkpath,Npts
  complex(8),dimension(4,4)                 :: Gamma0
  complex(8),dimension(4,4)                 :: Gamma5
  complex(8),dimension(4,4)                 :: GammaX
  complex(8),dimension(4,4)                 :: GammaY
  complex(8),dimension(4,4)                 :: GammaS
  !Model parameters:
  character(len=20)                         :: inputFILE
  real(8)                                   :: Wdis
  integer                                   :: Idum
  integer                                   :: disorder_type
  logical                                   :: bhz_pbc
  real(8)                                   :: mh
  real(8)                                   :: lambda
  real(8)                                   :: xmu
  real(8)                                   :: beta
  real(8)                                   :: wmix
  integer                                   :: Lfreq
  real(8)                                   :: wmin,wmax
  real(8)                                   :: eps
  real(8)                                   :: sb_field
  real(8)                                   :: it_error
  integer                                   :: MaxIter
  integer                                   :: Iter,Nsuccess=2
  logical                                   :: converged
  logical                                   :: with_lcm
  logical                                   :: with_mats_gf
  logical                                   :: with_real_gf
  real(8)                                   :: Uloc,Jh,JU
  !Mean-Fields:
  real(8)                                   :: mfSz,mfTz,mfRz,mfNt
  !Random 
  real(8),allocatable,dimension(:)          :: erandom
  !Hamiltonian
  complex(8),dimension(:,:,:),allocatable   :: Hk
  integer,dimension(:,:),allocatable        :: Links
  complex(8),dimension(:,:,:),allocatable   :: Hij
  complex(8),dimension(:,:,:,:),allocatable :: Hlat
  complex(8),dimension(:,:),allocatable     :: Sz
  real(8),dimension(:,:),allocatable        :: LsCM


  
contains



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



  subroutine setup_Abhz()
    logical :: bool
    integer :: ilat
    if(allocated(Hlat))deallocate(Hlat)
    allocate(Hlat(Nso,Nso,Nlat,Nlat))
    if(allocated(Hij))deallocate(Hij)
    allocate(Hij(Nlso,Nlso,1))
    if(allocated(erandom))deallocate(erandom)
    allocate(erandom(Nlat))
    call check_dimension("setup_disorder")
    !
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
    call save_array('erandom_'//str(idum)//'.used',erandom)
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
    !Reshape:
    Hij(:,:,1) = reshape_rank4_to_rank2(Hlat,Nso,Nlat)
    !
    open(99,file="list_idum.dat",access='append')
    write(99,*)idum
    close(99)
    !
    allocate(Sz(Nlso,Nlso))
    Sz = kron(eye(Nlat),GammaS)
    !
  end subroutine setup_Abhz





  !------------------------------------------------------------------
  ! Checks whether all the dimensions of the systems have been set
  ! properly in the calling procedure.
  !------------------------------------------------------------------
  subroutine check_dimension(caller)
    character(len=*) :: caller
    if(Nlso==0)stop str(caller)//" error: Nlso not set"
    if(Nso==0)stop str(caller)//" error: Nso not set"
    if(Nx==0)stop str(caller)//" error: Nx not set"
    if(Ny==0)stop str(caller)//" error: Ny not set"
    if(Nlat==0)stop str(caller)//" error: Nlat not set"
    if(Nocc==0)stop str(caller)//" error: Nocc not set"
    if(Nlat/=Nx*Ny)stop str(caller)//" error: Nlat != Nx*Ny"
    if(.not.allocated(Hlat)) stop str(caller)//" error: Hlat not allocated"
    call assert_shape(Hlat,[Nso,Nso,Nlat,Nlat],"check_dimension","Hlat")
    call assert_shape(Hij,[Nlso,Nlso,1],"check_dimension","Hij")

    ! stop str(caller)//" error: "
    ! stop str(caller)//" error: "
  end subroutine check_dimension




  subroutine get_gf(U,E,Gloc,axis)
    complex(8),dimension(Nlso,Nlso),intent(in)             :: U
    real(8),dimension(Nlso),intent(in)                     :: E
    complex(8),dimension(Nlat,Nso,Nso,Lfreq),intent(inout) :: Gloc
    character(len=*)                                       :: axis
    complex(8),dimension(Lfreq)                            :: wfreq
    complex(8),dimension(Nlso,Lfreq)                       :: csi
    integer                                                :: i,ilat,io,jo,is,js
    !
    !
    write(*,"(A)")"Get local Green's function, axis:"//str(axis)
    wfreq = build_frequency_array(axis)
    !
    !Allocate and setup the Matsubara freq.
    forall(i=1:Lfreq)csi(:,i)=one/(wfreq(i)+xmu-E(:))
    !
    call start_timer
    do concurrent(ilat=1:Nlat,io=1:Nso,jo=1:Nso,i=1:Lfreq)
       is = io + (ilat-1)*Nso
       js = jo + (ilat-1)*Nso
       Gloc(ilat,io,jo,i) = sum(U(is,:)*conjg(U(js,:))*csi(:,i))
    enddo
    call stop_timer
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
       wfreq = dcmplx(0d0,pi/beta*(2*arange(1,Lfreq)-1))
    case("realaxis","real","r")
       wfreq = dcmplx(linspace(wmin,wmax,Lfreq),eps)
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



END MODULE COMMON
