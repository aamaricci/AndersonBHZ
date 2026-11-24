MODULE COMMON
  USE SCIFOR, operator(.mx.) => operator(.x.)
  USE DMFT_TOOLS
  implicit none

  integer                                   :: Nso=0
  integer                                   :: Nlat=0
  integer                                   :: Nx=0
  integer                                   :: Ny=0
  integer                                   :: Nlso=0
  !
  integer                                   :: Nk
  integer                                   :: Nkpath,Npts,Nidum
  complex(8),dimension(4,4)                 :: Gamma0
  complex(8),dimension(4,4)                 :: Gamma5
  complex(8),dimension(4,4)                 :: GammaX
  complex(8),dimension(4,4)                 :: GammaY
  complex(8),dimension(4,4)                 :: GammaS
  !Model parameters:
  character(len=512)                        :: inputFILE
  character(len=512)                        :: idumFILE
  character(len=512)                        :: prefix="."
  real(8)                                   :: Wdis
  integer                                   :: Idum
  integer                                   :: disorder_type=0
  integer                                   :: spin_blocks=0
  character(len=16)                         :: method="lapack"
  character(len=16)                         :: cheb_kind="jackson"
  integer                                   :: cheb_order=256
  real(8)                                   :: cheb_tol=1d-12
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
  logical                                   :: with_spectrum
  integer                                   :: spectrum_order=1024
  integer                                   :: spectrum_points=2048
  real(8)                                   :: Uloc,Jh,sb_field
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
  logical,private                           :: wd_active=.false.
  integer,private                           :: wd_seed=0
  integer,private                           :: disorder_seed=0
  character(len=:),allocatable,private      :: previous_wd
  character(len=:),allocatable,private      :: active_wd
  character(len=:),allocatable,private      :: wd_prefix


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




  subroutine read_input(InputFile,driver)
    ! Read and validate the complete input contract for one AndersonBHZ
    ! driver. INPUTFILE is supplied by the caller and may come from the
    ! command line. DRIVER accepts "Abhz" or "mf_Abhz". METHOD is
    ! selected from its first letter: L means LAPACK and C Chebychev.
    character(len=*),intent(in)  :: InputFile
    character(len=*),intent(in)  :: driver
    character(len=:),allocatable :: driver_
    character(len=1)             :: method_initial
    !
    driver_ = to_lower(trim(driver))
    !
    ! Common model, disorder and representation controls.
    call parse_input_variable(prefix,"PREFIX",InputFile,default=".")
    if(driver_=="mf_abhz")then
      call parse_input_variable(Nx,"NKX",InputFile,default=25)
    else
      call parse_input_variable(Nx,"NX",InputFile,default=10)
    endif
    call parse_input_variable(Wdis,"WDIS",InputFile,default=0d0)
    call parse_input_variable(idum,"IDUM",InputFile,default=1234567)
    call parse_input_variable(disorder_type,"DISORDER_TYPE", &
      InputFile,default=0)
    call parse_input_variable(bhz_pbc,"BHZ_PBC",InputFile, &
      default=.true.)
    call parse_input_variable(mh,"MH",InputFile,default=3d0)
    call parse_input_variable(lambda,"LAMBDA",InputFile,default=0.3d0)
    call parse_input_variable(spin_blocks,"SPIN_BLOCKS",InputFile, &
      default=merge(2,1,driver_=="mf_abhz"), &
      comment="Spin representation: 0, 1 or 2")
    !
    ! Solver and Chebychev controls.
    call parse_input_variable(method,"METHOD",InputFile,default="lapack")
    call parse_input_variable(cheb_order,"CHEB_ORDER",InputFile, &
      default=256)
    call parse_input_variable(cheb_kind,"CHEB_KIND",InputFile, &
      default="jackson")
    if(driver_=="mf_abhz") &
      call parse_input_variable(cheb_tol,"CHEB_TOL",InputFile, &
      default=1d-12)
    !
    ! Chemical potential, frequency grid and convergence controls.
    if(driver_=="mf_abhz")then
      call parse_input_variable(mu,"XMU",InputFile,default=0d0)
      call parse_input_variable(ieta,"EPS",InputFile,default=4d-2)
      call parse_input_variable(wmin,"WMIN",InputFile,default=-10d0)
      call parse_input_variable(wmax,"WMAX",InputFile,default=10d0)
    else
      call parse_input_variable(mu,"MU",InputFile,default=0d0)
      call parse_input_variable(ieta,"IETA",InputFile,default=4d-2)
      call parse_input_variable(Lfreq,"LFREQ",InputFile,default=1024)
      call parse_input_variable(wmin,"WMIN",InputFile,default=-5d0)
      call parse_input_variable(wmax,"WMAX",InputFile,default=5d0)
    endif
    call parse_input_variable(wmix,"WMIX",InputFile,default=0.5d0)
    call parse_input_variable(it_error,"IT_ERROR",InputFile,default=1d-5)
    call parse_input_variable(maxiter,"MAXITER",InputFile,default=100)
    !
    ! Driver-specific physical and output controls.
    select case(driver_)
    case("abhz")
      call parse_input_variable(p_field,"P_FIELD",InputFile,default=0.01d0)
      call parse_input_variable(with_lcm,"WITH_LCM",InputFile, &
        default=.false.)
      call parse_input_variable(with_real_gf,"WITH_REAL_GF",InputFile, &
        default=.false.)
      with_spectrum = .false.
    case("mf_abhz")
      call parse_input_variable(Uloc,"ULOC",InputFile,default=1d0)
      call parse_input_variable(Jh,"JH",InputFile,default=0.125d0)
      call parse_input_variable(sb_field,"SB_FIELD",InputFile, &
        default=0.01d0)
      call parse_input_variable(with_spectrum,"WITH_SPECTRUM", &
        InputFile,default=.false.)
      call parse_input_variable(spectrum_points,"L",InputFile, &
        default=2048)
      call parse_input_variable(spectrum_order,"SPECTRUM_ORDER", &
        InputFile,default=1024)
      with_lcm     = .false.
      with_real_gf = .false.
    case default
      error stop "read_input: DRIVER must be Abhz or mf_Abhz"
    end select
    !
    ! Normalize the method from its first letter only.
    method = to_lower(adjustl(method))
    if(len_trim(method)<1)error stop "read_input: empty METHOD"
    method_initial = method(1:1)
    select case(method_initial)
    case("l")
      method = "lapack"
    case("c")
      method = "chebychev"
    case default
      error stop "read_input: METHOD must start with L or C"
    end select
    cheb_kind = to_lower(trim(cheb_kind))
    !
    ! Validate individual values and coupled selections.
    call validate_spin_selection("read_input")
    if(Nx<2)error stop "read_input: lattice size must be at least 2"
    if(Wdis<0d0)error stop "read_input: WDIS must be nonnegative"
    if(lambda<0d0)error stop "read_input: LAMBDA must be nonnegative"
    if(cheb_order<1)error stop "read_input: CHEB_ORDER must be positive"
    if(cheb_tol<=0d0)error stop "read_input: CHEB_TOL must be positive"
    if(trim(cheb_kind)/="none".and.trim(cheb_kind)/="jackson") &
      error stop "read_input: CHEB_KIND must be none or jackson"
    if(wmax<=wmin)error stop "read_input: require WMAX > WMIN"
    if(ieta<=0d0)error stop "read_input: IETA or EPS must be positive"
    if(wmix<=0d0.or.wmix>1d0) &
      error stop "read_input: require 0 < WMIX <= 1"
    if(it_error<=0d0.or.maxiter<1) &
      error stop "read_input: invalid convergence controls"
    if(driver_=="abhz")then
      if(Lfreq<1)error stop "read_input: LFREQ must be positive"
      if(method=="chebychev".and.(with_lcm.or.with_real_gf)) &
        error stop "read_input: LCM and real GF require LAPACK"
    else
      if(Uloc<0d0.or.Jh<0d0) &
        error stop "read_input: ULOC and JH must be nonnegative"
      if(spectrum_points<2.or.spectrum_order<1) &
        error stop "read_input: invalid spectrum controls"
    endif
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
  end subroutine read_input




  subroutine validate_spin_selection(caller)
    ! Validate the global spin representation selected by the input.
    ! SPIN_BLOCKS=1 solves the up sector and obtains down by time
    ! reversal, SPIN_BLOCKS=2 solves both sectors independently, and
    ! SPIN_BLOCKS=0 uses the complete 4x4 spin-orbital Hamiltonian.
    ! The one-sector representation is invalid for longitudinal spin
    ! disorder because DISORDER_TYPE=2 breaks time-reversal symmetry.
    character(len=*),intent(in) :: caller
    !
    if(.not.any(spin_blocks==[0,1,2])) &
      error stop trim(caller)//": SPIN_BLOCKS must be 0, 1 or 2"
    if(disorder_type<0.or.disorder_type>2) &
      error stop trim(caller)//": DISORDER_TYPE must be in [0,2]"
    if(spin_blocks==1.and.disorder_type==2) &
      error stop trim(caller)// &
        ": Sz disorder requires SPIN_BLOCKS=0 or 2"
  end subroutine validate_spin_selection


  subroutine set_seed_wd()
    ! Enter PREFIX/seed_<IDUM> and retain the caller's absolute working
    ! directory. Repeated calls for the active seed are harmless. Call
    ! reset_seed_wd before changing IDUM or PREFIX.
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


  subroutine reset_seed_wd()
    ! Restore the absolute working directory saved by set_seed_wd.
    integer :: status
    !
    if(.not.wd_active)return
    call chdir(previous_wd,status)
    if(status/=0) &
      error stop "reset_seed_wd: cannot restore working directory"
    wd_active = .false.
    wd_seed   = 0
    deallocate(previous_wd,active_wd,wd_prefix)
  end subroutine reset_seed_wd


  function seed_directory() result(folder)
    ! Return the active seed directory, or its path relative to the caller.
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
  end function Afm_Vec

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





  subroutine setup_Abhz(verbose,spH,spin)
    ! Build the square BHZ Hamiltonian for IDUM (serial). With SPH,
    ! SPIN=+1/-1 returns one 2*Nlat spin sector and SPIN=0 returns the
    ! complete 4*Nlat Hamiltonian. SPIN defaults to zero. Without SPH,
    ! populate the dense COMMON arrays for legacy callers.
    ! VERBOSE defaults to true and controls logging and disorder output.
    ! A matching erandom_<IDUM>.restart overrides the generated disorder.
    logical,optional,intent(in)                  :: verbose
    type(sparse_zmatrix_csr),optional,intent(out) :: spH
    integer,optional,intent(in)                  :: spin
    logical                                      :: verbose_,bool,new_disorder
    integer                                      :: ilat
    integer                                      :: spin_
    complex(8)                                   :: channel(4,4)
    !
    call check_dimension("setup_Abhz")
    spin_ = 0
    if(present(spin))spin_=spin
    if(abs(spin_)>1) &
      error stop "setup_Abhz: SPIN must be -1, 0 or +1"
    if(.not.present(spH).and.spin_/=0) &
      error stop "setup_Abhz: dense legacy setup requires SPIN=0"
    select case(disorder_type)
    case(0)
       channel = Gamma0
    case(1)
       channel = Gamma5
    case(2)
       channel = GammaS
    case default
       error stop "setup_Abhz: DISORDER_TYPE must be in [0,2]"
    end select
    verbose_ = .true.
    if(present(verbose))verbose_=verbose
    !
    new_disorder = .not.allocated(erandom)
    if(new_disorder)then
       allocate(erandom(Nlat))
       call mersenne_init(idum)
       call mt_random(erandom)
       erandom = (2d0*erandom-1d0)*Wdis/2d0
       inquire(file='erandom_'//str(idum)//'.restart',exist=bool)
       if(bool)then
          if(file_length('erandom_'//str(idum)//'.restart')/=Nlat) &
            error stop "setup_Abhz: invalid disorder restart length"
          call read_array('erandom_'//str(idum)//'.restart',erandom)
       endif
       disorder_seed = idum
    elseif(disorder_seed/=idum)then
       error stop "setup_Abhz: call free_Abhz before changing IDUM"
    endif
    if(verbose_.and.new_disorder)then
       write(*,*) "Solve IDUM:",idum," Nlso:",Nlso
       call save_array('erandom_'//str(idum)//'.used',erandom)
    endif
    !
    if(allocated(Links))deallocate(Links)
    allocate(Links(4,2))
    Links(1,:) = [1,0]
    Links(2,:) = [0,1]
    Links(3,:) = [-1,0]
    Links(4,:) = [0,-1]
    if(present(spH))then
       if(spin_==0)then
          call build_Abhz_csr(spH,channel)
       else
          call build_Abhz_block_csr(spH,spin_)
       endif
    else
       if(allocated(Hlat))deallocate(Hlat)
       if(allocated(Hij))deallocate(Hij)
       if(allocated(Sz))deallocate(Sz)
       allocate(Hlat(Nso,Nso,Nlat,Nlat),Hij(Nlso,Nlso,1))
       call TB_build_model(Hlat,ts_model,Nso,[Nx,Ny],Links,pbc=bhz_pbc)
       do ilat=1,Nlat
          Hlat(:,:,ilat,ilat) = Hlat(:,:,ilat,ilat) + erandom(ilat)*channel
       enddo
       Hij(:,:,1) = reshape_rank4_to_rank2(Hlat,Nso,Nlat)
       Sz = kron(eye(Nlat),GammaS)
    endif
  end subroutine setup_Abhz


  function ts_model_block(link,spin) result(Hts)
    ! Return one conserved-spin 2x2 hopping block in orbital space.
    integer,intent(in)          :: link,spin
    complex(8),dimension(2,2)   :: Hts
    !
    select case(link)
    case(0)
       Hts = Mh*pauli_tau_z
    case(1)
       Hts = -0.5d0*pauli_tau_z &
         +xi*0.5d0*lambda*spin*pauli_tau_x
    case(2)
       Hts = -0.5d0*pauli_tau_z &
         -xi*0.5d0*lambda*pauli_tau_y
    case(3)
       Hts = -0.5d0*pauli_tau_z &
         -xi*0.5d0*lambda*spin*pauli_tau_x
    case(4)
       Hts = -0.5d0*pauli_tau_z &
         +xi*0.5d0*lambda*pauli_tau_y
    case default
       error stop "ts_model_block: LINK must be in [0,4]"
    end select
  end function ts_model_block


  subroutine build_Abhz_block_csr(spH,spin)
    ! Build one 2*Nlat spin sector. New CSR columns are inserted in
    ! increasing order and coincident periodic links are accumulated.
    type(sparse_zmatrix_csr),intent(inout) :: spH
    integer,intent(in)                     :: spin
    complex(8)                             :: blocks(2,2,0:4),block(2,2)
    complex(8)                             :: disorder_block(2,2)
    integer                                :: coords(Nlat,2)
    integer                                :: target(0:4),r(2)
    integer                                :: ilat,jlat,link,io,jo,i,j
    !
    select case(disorder_type)
    case(0)
       disorder_block = pauli_tau_0
    case(1)
       disorder_block = pauli_tau_z
    case(2)
       disorder_block = spin*pauli_tau_0
    case default
       error stop "build_Abhz_block_csr: invalid DISORDER_TYPE"
    end select
    call TB_build_CoordGrid([Nx,Ny],coords)
    call spH%init(2*Nlat,2*Nlat)
    do ilat=1,Nlat
       target(0) = ilat
       blocks(:,:,0) = ts_model_block(0,spin) &
         +erandom(ilat)*disorder_block
       do link=1,4
          r = coords(ilat,:) + Links(link,:)
          if(bhz_pbc)r = modulo(r-1,[Nx,Ny])+1
          target(link) = TB_find_IndxCoord(r,[Nx,Ny])
          blocks(:,:,link) = ts_model_block(link,spin)
       enddo
       do while(any(target>0))
          jlat = minval(target,mask=target>0)
          block = zero
          do link=0,4
             if(target(link)/=jlat)cycle
             block = block + blocks(:,:,link)
             target(link) = 0
          enddo
          do io=1,2
             i = io + (ilat-1)*2
             do jo=1,2
                j = jo + (jlat-1)*2
                if(block(io,jo)==zero)cycle
                call spH%fast_insert(block(io,jo),i,j)
             enddo
          enddo
       enddo
    enddo
  end subroutine build_Abhz_block_csr


  subroutine assemble_Abhz_full(spHup,spHdown,spHfull)
    ! Assemble the site-major 4*Nlat Hamiltonian from two 2*Nlat
    ! sectors. The spin-off-diagonal blocks are exactly zero.
    type(sparse_zmatrix_csr),intent(in)  :: spHup,spHdown
    type(sparse_zmatrix_csr),intent(out) :: spHfull
    integer                             :: spin,ib,jb,ifull,jfull
    integer                             :: ilat,jlat,iorb,jorb,pos
    !
    if(spHup%Nrow/=2*Nlat.or.spHdown%Nrow/=2*Nlat) &
      error stop "assemble_Abhz_full: invalid spin-block dimensions"
    call spHfull%init(4*Nlat,4*Nlat)
    do spin=1,2
       do ib=1,2*Nlat
          ilat  = (ib-1)/2+1
          iorb  = modulo(ib-1,2)+1
          ifull = iorb + (spin-1)*2 + (ilat-1)*4
          if(spin==1)then
             do pos=1,spHup%row(ib)%Size
                jb = spHup%row(ib)%cols(pos)
                jlat = (jb-1)/2+1
                jorb = modulo(jb-1,2)+1
                jfull = jorb + (spin-1)*2 + (jlat-1)*4
                call spHfull%fast_insert( &
                  spHup%row(ib)%vals(pos),ifull,jfull)
             enddo
          else
             do pos=1,spHdown%row(ib)%Size
                jb = spHdown%row(ib)%cols(pos)
                jlat = (jb-1)/2+1
                jorb = modulo(jb-1,2)+1
                jfull = jorb + (spin-1)*2 + (jlat-1)*4
                call spHfull%fast_insert( &
                  spHdown%row(ib)%vals(pos),ifull,jfull)
             enddo
          endif
       enddo
    enddo
  end subroutine assemble_Abhz_full


  subroutine build_Abhz_csr(spH,channel)
    ! Insert each physical column once and in increasing order.
    ! For Nx=2 with PBC, opposite links can reach the same site.
    type(sparse_zmatrix_csr),intent(inout) :: spH
    complex(8),intent(in)                 :: channel(4,4)
    complex(8)                            :: blocks(Nso,Nso,0:4),block(Nso,Nso)
    integer                               :: coords(Nlat,2),target(0:4),r(2)
    integer                               :: ilat,jlat,link,io,jo,i,j
    !
    call TB_build_CoordGrid([Nx,Ny],coords)
    call spH%init(Nlso,Nlso)
    do ilat=1,Nlat
       target(0) = ilat
       blocks(:,:,0) = ts_model(0,Nso) + erandom(ilat)*channel
       do link=1,4
          r = coords(ilat,:) + Links(link,:)
          if(bhz_pbc)r = modulo(r-1,[Nx,Ny])+1
          target(link) = TB_find_IndxCoord(r,[Nx,Ny])
          blocks(:,:,link) = ts_model(link,Nso)
       enddo
       do while(any(target>0))
          jlat = minval(target,mask=target>0)
          block = zero
          do link=0,4
             if(target(link)/=jlat)cycle
             block = block + blocks(:,:,link)
             target(link) = 0
          enddo
          do io=1,Nso
             i = io + (ilat-1)*Nso
             do jo=1,Nso
                j = jo + (jlat-1)*Nso
                if(block(io,jo)==zero)cycle
                call spH%fast_insert(block(io,jo),i,j)
             enddo
          enddo
       enddo
    enddo
  end subroutine build_Abhz_csr


  subroutine free_Abhz
    if(allocated(Hlat))deallocate(Hlat)
    if(allocated(Hij))deallocate(Hij)
    if(allocated(erandom))deallocate(erandom)
    if(allocated(Links))deallocate(Links)
    if(allocated(Sz))deallocate(Sz)
    disorder_seed = 0
  end subroutine free_Abhz




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
    if(allocated(Hlat))call assert_shape(Hlat,[Nso,Nso,Nlat,Nlat],"check_dimension","Hlat")
    if(allocated(Hij))call assert_shape(Hij,[Nlso,Nlso,1],"check_dimension","Hij")
    if(allocated(Sz))call assert_shape(Sz,[Nlso,Nlso],"check_dimension","Sz")
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




END MODULE COMMON
