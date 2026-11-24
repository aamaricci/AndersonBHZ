MODULE COMMON
  USE SCIFOR, operator(.mx.) => operator(.x.)
  USE DMFT_TOOLS
  implicit none

  integer                              :: Nso=0
  integer                              :: Nlat=0
  integer                              :: Nx=0
  integer                              :: Ny=0
  integer                              :: Nlso=0
  !
  complex(8),dimension(4,4),public     :: Gamma0
  complex(8),dimension(4,4),public     :: Gamma5
  complex(8),dimension(4,4),public     :: GammaX
  complex(8),dimension(4,4),public     :: GammaY
  complex(8),dimension(4,4),public     :: GammaS
  !
  integer                              :: Nk
  integer                              :: Nkpath,Npts,Nidum
  !Model parameters:
  character(len=512)                   :: inputFILE
  character(len=512)                   :: idumFILE
  character(len=512)                   :: prefix="."
  real(8)                              :: Wdis
  integer                              :: Idum
  integer                              :: disorder_type=0
  integer                              :: spin_blocks=0
  character(len=16)                    :: method="lapack"
  character(len=16)                    :: cheb_kind="jackson"
  integer                              :: cheb_order=256
  logical                              :: bhz_pbc
  real(8)                              :: mh
  real(8)                              :: lambda
  real(8)                              :: wmix
  integer                              :: Lfreq
  real(8)                              :: wmin,wmax
  real(8)                              :: temp
  real(8)                              :: mu,ieta
  real(8)                              :: p_field
  real(8)                              :: it_error
  integer                              :: MaxIter
  logical                              :: with_lcm

  real(8)                              :: Uloc,Jh,sb_field
  !Random
  integer,allocatable,dimension(:)     :: list_idum
  character(len=255)                   :: here
  character(len=20)                    :: dir
  real(8),dimension(:,:),allocatable   :: LsCM
  !
  logical,private                      :: wd_active=.false.
  integer,private                      :: wd_seed=0
  character(len=:),allocatable,private :: previous_wd
  character(len=:),allocatable,private :: active_wd
  character(len=:),allocatable,private :: wd_prefix


contains

  subroutine parse_input_file(filename,default)
    ! Read INPUTFILE from the command line without list-directed
    ! character parsing. This preserves absolute paths, whose slash
    ! would otherwise terminate a list-directed Fortran READ.
    character(len=*),intent(out) :: filename
    character(len=*),intent(in)  :: default
    character(len=1024)          :: argument,key,value
    integer                      :: i,separator,n
    !
    filename = default
    do i=1,command_argument_count()
      call get_command_argument(i,argument)
      separator = index(argument,"=")
      if(separator<2)cycle
      key = to_upper(adjustl(argument(:separator-1)))
      if(trim(key)/="INPUTFILE")cycle
      value = adjustl(argument(separator+1:))
      n = len_trim(value)
      if(n>=2)then
        if((value(1:1)=="'".and.value(n:n)=="'").or. &
          (value(1:1)=='"'.and.value(n:n)=='"')) &
          value = value(2:n-1)
      endif
      if(len_trim(value)<1)error stop "parse_input_file: empty INPUTFILE"
      filename = trim(value)
      write(0,*)"Variable INPUTFILE updated to "//trim(filename)
    enddo
  end subroutine parse_input_file



  ! Read and validate the complete input contract for one AndersonBHZ
  ! driver. INPUTFILE is supplied by the caller and may come from the
  ! command line. METHOD is
  ! selected from its first letter: L means LAPACK and C Chebychev.
  subroutine read_input(InputFile)
    character(len=*),intent(in)  :: InputFile
    !
    ! Common model, disorder and representation controls.
    call parse_input_variable(prefix,"PREFIX",InputFile,default=".")
    call parse_input_variable(Nx,"NX",InputFile,default=10)
    call parse_input_variable(Wdis,"WDIS",InputFile,default=0d0)
    call parse_input_variable(idum,"IDUM",InputFile,default=1234567)
    call parse_input_variable(disorder_type,"DISORDER_TYPE",InputFile,default=0)
    call parse_input_variable(bhz_pbc,"BHZ_PBC",InputFile,default=.true.)
    call parse_input_variable(mh,"MH",InputFile,default=3d0)
    call parse_input_variable(lambda,"LAMBDA",InputFile,default=0.3d0)
    call parse_input_variable(spin_blocks,"SPIN_BLOCKS",InputFile,default=2,&
      comment="Spin representation: 0, 1 or 2")
    !
    ! Driver-specific physical and output controls.
    call parse_input_variable(Uloc,"ULOC",InputFile,default=1d0)
    call parse_input_variable(Jh,"JH",InputFile,default=0.125d0)
    call parse_input_variable(sb_field,"SB_FIELD",InputFile,default=0.01d0)
    call parse_input_variable(with_lcm,"WITH_LCM",InputFile,default=.false.)
    !
    ! Solver and Chebychev controls.
    call parse_input_variable(method,"METHOD",InputFile,default="lapack",&
      comment="Solver method:Lapack, Chebychev")
    call parse_input_variable(cheb_order,"CHEB_ORDER",InputFile,default=256)
    call parse_input_variable(cheb_kind,"CHEB_KIND",InputFile,default="jackson")
    !
    ! Chemical potential, frequency grid and convergence controls.
    call parse_input_variable(mu,"MU",InputFile,default=0d0)
    call parse_input_variable(ieta,"IETA",InputFile,default=4d-2)
    call parse_input_variable(Lfreq,"LFREQ",InputFile,default=1024)
    call parse_input_variable(wmin,"WMIN",InputFile,default=-10d0)
    call parse_input_variable(wmax,"WMAX",InputFile,default=10d0)
    call parse_input_variable(wmix,"WMIX",InputFile,default=0.5d0)
    call parse_input_variable(it_error,"IT_ERROR",InputFile,default=1d-5)
    call parse_input_variable(maxiter,"MAXITER",InputFile,default=100)
    !
    !
    ! Normalize the method from its first letter only.
    select case(to_lower(str(method(1:1))))
     case("l");method = "lapack"
     case("c");method = "chebychev"
     case default;error stop "read_input: METHOD must start with L or C"
    end select
    cheb_kind = to_lower(trim(cheb_kind))
    !
    ! Validate individual values and coupled selections.
    call validate_spin_selection("read_input")
    if(Nx<2)error stop "read_input: lattice size must be at least 2"
    if(Wdis<0d0)error stop "read_input: WDIS must be nonnegative"
    if(lambda<0d0)error stop "read_input: LAMBDA must be nonnegative"
    if(cheb_order<1)error stop "read_input: CHEB_ORDER must be positive"
    if(str(cheb_kind)/="none".and.str(cheb_kind)/="jackson") &
      error stop "read_input: CHEB_KIND must be none or jackson"
    if(wmax<=wmin)error stop "read_input: require WMAX > WMIN"
    if(ieta<=0d0)error stop "read_input: IETA or EPS must be positive"
    if(wmix<=0d0.or.wmix>1d0)error stop "read_input: require 0 < WMIX <= 1"
    if(it_error<=0d0.or.maxiter<1)error stop "read_input: invalid convergence controls"
    if(Lfreq<1)error stop "read_input: LFREQ must be positive"
    if(Uloc<0d0.or.Jh<0d0)error stop "read_input: ULOC and JH must be nonnegative"
    !
    ! Derived dimensions and common numerical controls.
    Ny   = Nx
    Nk   = Nx*Ny
    Nlat = Nk
    Nso  = 4
    Nlso = Nlat*Nso
    temp = 0d0
    call add_ctrl_var(2,"NORB")
    call add_ctrl_var(2,"NSPIN")
    call add_ctrl_var(mu,"xmu")
    call add_ctrl_var(wmin,"wini")
    call add_ctrl_var(wmax,"wfin")
    call add_ctrl_var(ieta,"eps")
    call print_input()
    call setup_GammaMatrices()
  end subroutine read_input


  ! Initialize the local spin-orbital operators for the BHZ basis.
  subroutine setup_GammaMatrices()
    gamma0 = kron( pauli_sigma_0, pauli_tau_0)
    gammaX = kron( pauli_sigma_z, pauli_tau_x)
    gammaY = kron( pauli_sigma_0,-pauli_tau_y)
    gamma5 = kron( pauli_sigma_0, pauli_tau_z)
    gammaS = kron( pauli_sigma_z, pauli_tau_0)
  end subroutine setup_GammaMatrices



  ! Validate the global spin representation selected by the input.
  ! SPIN_BLOCKS=1 solves the up sector and obtains down by time
  ! reversal, SPIN_BLOCKS=2 solves both sectors independently, and
  ! SPIN_BLOCKS=0 uses the complete 4x4 spin-orbital Hamiltonian.
  ! The one-sector representation is invalid for longitudinal spin
  ! disorder because DISORDER_TYPE=2 breaks time-reversal symmetry.
  subroutine validate_spin_selection(caller)
    character(len=*),intent(in) :: caller
    if(.not.any(spin_blocks==[0,1,2]))error stop trim(caller)//": SPIN_BLOCKS must be 0, 1 or 2"
    if(disorder_type<0.or.disorder_type>2)error stop trim(caller)//": DISORDER_TYPE must be in [0,2]"
    if(spin_blocks==1.and.disorder_type==2)error stop trim(caller)//": Sz disorder requires SPIN_BLOCKS=0 or 2"
  end subroutine validate_spin_selection




  ! Enter PREFIX/seed_<IDUM> and retain the caller's absolute working
  ! directory. Repeated calls for the active seed are harmless. Call
  ! reset_seed_wd before changing IDUM or PREFIX.
  subroutine set_seed_wd()
    character(len=4096)          :: location
    character(len=:),allocatable :: destination
    integer                      :: status
    !
    if(len_trim(prefix)==0) &
      error stop "set_seed_wd: PREFIX must be a nonempty directory"
    if(wd_active)then
      if(idum/=wd_seed.or.trim(prefix)/=wd_prefix) &
        error stop "set_seed_wd: reset the working directory first"
      call getcwd(location,status)
      if(status/=0.or.trim(location)/=active_wd) &
        error stop "set_seed_wd: active working directory changed"
      return
    endif
    !
    call getcwd(location,status)
    if(status/=0)error stop "set_seed_wd: cannot read working directory"
    previous_wd = trim(location)
    destination = trim(prefix)//"/seed_"//str(idum)
    call make_output_directory(destination)
    call chdir(destination,status)
    if(status/=0)error stop "set_seed_wd: cannot enter seed directory"
    call getcwd(location,status)
    if(status/=0)error stop "set_seed_wd: cannot resolve seed directory"
    active_wd = trim(location)
    wd_prefix = trim(prefix)
    wd_seed   = idum
    wd_active = .true.
  end subroutine set_seed_wd

  ! Restore the absolute working directory saved by set_seed_wd.
  subroutine reset_seed_wd()
    integer :: status
    !
    if(.not.wd_active)return
    call chdir(previous_wd,status)
    if(status/=0)error stop "reset_seed_wd: cannot restore working directory"
    wd_active = .false.
    wd_seed   = 0
    deallocate(previous_wd,active_wd,wd_prefix)
  end subroutine reset_seed_wd

  ! Return the active seed directory, or its path relative to the caller.
  function seed_directory() result(folder)
    character(len=:),allocatable :: folder
    !
    if(wd_active)then
      if(idum/=wd_seed.or.trim(prefix)/=wd_prefix) &
        error stop "seed_directory: active seed and PREFIX disagree"
      folder = active_wd
    else
      folder = trim(prefix)//"/seed_"//str(idum)
    endif
  end function seed_directory








  function file_basename(path) result(name)
    character(len=*),intent(in)      :: path
    character(len=:),allocatable     :: name
    integer                          :: separator
    !
    separator = scan(trim(path),"/",back=.true.)
    name = path(separator+1:len_trim(path))
  end function file_basename


  subroutine make_output_directory(path)
    ! Create PATH recursively. Quote it as a single shell argument.
    character(len=*),intent(in)  :: path
    character(len=:),allocatable :: quoted
    integer                      :: i,status,command_status
    !
    quoted = "'"
    do i=1,len_trim(path)
      if(path(i:i)=="'")then
        quoted = quoted//"'"//achar(34)//"'"//achar(34)//"'"
      else
        quoted = quoted//path(i:i)
      endif
    enddo
    quoted = quoted//"'"
    call execute_command_line("mkdir -p "//quoted,exitstat=status, &
      cmdstat=command_status)
    if(command_status/=0.or.status/=0) &
      error stop "make_output_directory: cannot create output directory"
  end subroutine make_output_directory




  ! Return the zero-temperature Fermi weight Theta(-X). As an
  ! ELEMENTAL function it accepts both scalar and array arguments.
  elemental function fermi_zero(x) result(fermi_weight)
    real(8),intent(in) :: x
    real(8)            :: fermi_weight
    !
    if(x<0d0)then
      fermi_weight = 1d0
    elseif(x>0d0)then
      fermi_weight = 0d0
    else
      fermi_weight = 0.5d0
    endif
  end function fermi_zero





  function xcoord(ilat) result(ix)
    integer :: ilat
    integer :: ix
    ix = mod(ilat,Nx);if(ix==0)ix=Nx
  end function xcoord

  function ycoord(ilat) result(iy)
    integer :: ilat
    integer :: iy
    iy = (ilat-1)/Nx+1
  end function ycoord

  function afm_sign(ilat) result(afm)
    integer :: ilat
    integer :: Afm
    integer :: ix,iy
    ix  = xcoord(ilat)
    iy  = ycoord(ilat)
    Afm = (-1)**(ix+(iy-1))
  end function afm_sign

  function afm_vec() result(afm)
    real(8),dimension(Nlat) :: Afm
    integer                 :: ilat
    do ilat=1,Nlat
      Afm(ilat) = afm_sign(ilat)
    enddo
  end function afm_vec

  function get_AFMmean(Szii) result(eSz)
    real(8),dimension(Nlat) :: Szii
    real(8)                 :: eSz
    integer,dimension(Nlat) :: Afm
    integer :: ilat
    do ilat=1,Nlat
      Afm(ilat) = afm_sign(ilat)
    enddo
    eSz = get_mean(Szii*Afm)
  end function get_AFMmean





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
    if(Nlat/=Nx*Ny)stop str(caller)//" error: Nlat != Nx*Ny"
  end subroutine check_dimension




  subroutine read_list_idum(idumFILE,Nidum)
    character(len=*) :: idumFILE
    integer          :: Nidum,unit,id,ii
    logical          :: iexist
    inquire(file=str(idumFILE),exist=iexist)
    if(iexist)then
      Nidum = file_length(str(idumFILE))
      allocate(list_idum(Nidum))
      open(unit=free_unit(unit),file=str(idumFILE))
      do ii=1,Nidum
        read(unit,*)id,list_idum(ii)
      enddo
      close(unit)
    else
      Nidum = 1
      allocate(list_idum(Nidum))
      list_idum = idum
    endif
  end subroutine read_list_idum





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


  subroutine add_to(vec,vals)
    real(8),dimension(:),allocatable,intent(inout) :: vec
    real(8),dimension(:),intent(in)                :: vals
    real(8)                                        :: val
    real(8),dimension(:),allocatable               :: tmp
    integer                                        :: i,n,m
    !
    m = size(vals)
    if (allocated(vec)) then
      n = size(vec)
      allocate(tmp(n+m))
      tmp(:n) = vec
      call move_alloc(tmp,vec)
      n = n + m
    else
      n = m
      allocate(vec(n))
    end if
    !
    !Put val as last m entries:
    vec(n-m+1:n) = vals
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine add_to




end module COMMON
