MODULE COMMON
  USE SCIFOR, operator(.mx.) => operator(.x.)
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
  logical,private                           :: wd_active=.false.
  integer,private                           :: wd_seed=0
  integer,private                           :: disorder_seed=0
  character(len=:),allocatable,private      :: previous_wd
  character(len=:),allocatable,private      :: active_wd
  character(len=:),allocatable,private      :: wd_prefix


contains


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


  subroutine cheb_fermi_coefficients(scale,shift,fermi_level,coeff,kind)
    ! Return the Chebyshev coefficients of Theta(fermi_level-E) on
    ! [shift-scale,shift+scale] and activate them in SciFortran.
    ! The series takes the value 1/2 at E=fermi_level. KIND accepts
    ! "none" or "jackson" and defaults to "jackson".
    real(8),intent(in)                    :: scale,shift,fermi_level
    real(8),dimension(0:),intent(out)     :: coeff
    character(len=*),optional,intent(in)  :: kind
    character(len=16)                     :: kernel
    real(8)                               :: xfermi,theta,angle,g
    integer                               :: m,order
    !
    if(scale<=0d0) &
      error stop "cheb_fermi_coefficients: SCALE must be positive"
    order  = ubound(coeff,1)
    kernel = "jackson"
    if(present(kind))kernel=to_lower(trim(kind))
    select case(kernel)
    case("none","n")
      kernel = "none"
    case("jackson","j")
      kernel = "jackson"
    case default
      error stop "cheb_fermi_coefficients: KIND must be none or jackson"
    end select
    !
    coeff  = 0d0
    xfermi = (fermi_level-shift)/scale
    if(xfermi<=-1d0)then
      call cheb_set(scale,shift,coeff)
      return
    elseif(xfermi>=1d0)then
      coeff(0) = 1d0
      call cheb_set(scale,shift,coeff)
      return
    endif
    !
    angle   = acos(xfermi)
    coeff(0) = (pi-angle)/pi
    do m=1,order
      coeff(m) = -2d0*sin(m*angle)/(pi*m)
    enddo
    if(kernel=="jackson")then
      theta = pi/(order+2)
      do m=0,order
        g = ((order+2-m)*cos(m*theta) + &
          sin(m*theta)/tan(theta))/(order+2)
        coeff(m) = coeff(m)*g
      enddo
    endif
    call cheb_set(scale,shift,coeff)
  end subroutine cheb_fermi_coefficients


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
    if(allocated(U))deallocate(U)
    if(allocated(E))deallocate(E)
    if(allocated(PSzP))deallocate(PSzP)
    if(allocated(Epsp))deallocate(Epsp)
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
    if(Nocc==0)stop str(caller)//" error: Nocc not set"
    if(Nlat/=Nx*Ny)stop str(caller)//" error: Nlat != Nx*Ny"
    if(allocated(Hlat))call assert_shape(Hlat,[Nso,Nso,Nlat,Nlat],"check_dimension","Hlat")
    if(allocated(Hij))call assert_shape(Hij,[Nlso,Nlso,1],"check_dimension","Hij")
    if(allocated(U))call assert_shape(U,[Nlso,Nlso],"check_dimension","U")
    if(allocated(E))call assert_shape(E,[Nlso],"check_dimension","E")
    if(allocated(Sz))call assert_shape(Sz,[Nlso,Nlso],"check_dimension","Sz")
    if(allocated(PSzP))call assert_shape(PSzP,[Nocc,Nocc],"check_dimension","PSzP")
    if(allocated(Epsp))call assert_shape(Epsp,[Nocc],"check_dimension","Epsp")
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
    call save_array("E.dat",E)
    call check_Egap("push_Bloch")
  end subroutine push_Bloch




  subroutine check_Pgap(caller)
    character(len=*) :: caller
    integer          :: N
    real(8)          :: Ep,Em,Pgap
    N    = int(Nocc/2)
    Ep   = Epsp(N+1)
    Em   = Epsp(N)
    Pgap = Ep - Em
    write(*,"(A,G21.13)")"PSzP gap=",Pgap
    call save_array("PSzPgap.dat",Pgap)
    !
    if(Pgap<1d-12)then
       write(*,*)str(caller)//" warning: closing of the PSzP spectrum"
    elseif(Ep*Em>0d0)then
       write(*,*)str(caller)//" warning: PSzP spectrum not symmetric"
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
    write(*,"(A,G21.13)")"E gap=",Egap
    call save_array("Egap.dat",Egap)
    !
    ! if(Egap<1d-12)then
    !    write(*,*)str(caller)//" warning: closing of the H spectrum"
    ! elseif(Ep*Em>0d0)then
    !    write(*,*)str(caller)//" warning: H spectrum not symmetric"
    ! endif
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
    write(*,"(A)")"Get local Green's function, axis:"//str(axis)
    wfreq = build_frequency_array(axis)
    !
    !Allocate and setup the Matsubara freq.
    forall(i=1:Lfreq)csi(:,i)=one/(wfreq(i)+mu-E(:))
    !
    Gloc = zero
    call start_timer
    do i=1,Lfreq
       do concurrent(ilat=1:Nlat,io=1:Nso,jo=1:Nso)
          is = io + (ilat-1)*Nso
          js = jo + (ilat-1)*Nso
          Gloc(ilat,io,jo,i) = sum(U(is,:)*conjg(U(js,:))*csi(:,i))!can use matmul
       enddo
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
