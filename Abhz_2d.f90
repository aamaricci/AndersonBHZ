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
!\eta is a random variable uniformly distributed within -W/2:W/2
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
  real(8),dimension(:),allocatable          :: dens_up,dens_down
  real(8),dimension(:),allocatable          :: dens_full
  real(8),dimension(:),allocatable          :: Ev_up,Ev_down
  complex(8),dimension(:,:),allocatable     :: U_up,U_down
  integer                                   :: ilat,iorb,io
  integer                                   :: spin_blocks
  type(sparse_zmatrix_csr),target           :: spHup,spHdown
  type(sparse_zmatrix_csr),target           :: spHfull
  type(sparse_zmatrix_csr),pointer          :: active_spH=>null()
  integer                                   :: cheb_order
  character(len=32)                         :: method,cheb_kind


  !
  !Read input:
  call parse_cmd_variable(inputFILE,"inputFILE",default="inputABHZ.conf")
  call parse_input_variable(prefix,"PREFIX",inputFILE,default=".")
  call parse_input_variable(Nx,"Nx",inputFILE,default=10)
  call parse_input_variable(Wdis,"WDIS",inputFILE,default=0d0)
  call parse_input_variable(idum,"IDUM",inputFILE,default=1234567)
  call parse_input_variable(disorder_type,"DISORDER_TYPE",inputFILE,default=0)
  call parse_input_variable(bhz_pbc,"BHZ_PBC",inputFILE,default=.true.)
  call parse_input_variable(mh,"MH",inputFILE,default=3.d0)
  call parse_input_variable(lambda,"LAMBDA",inputFILE,default=0.3d0)
  call parse_input_variable(mu,"MU",inputFILE,default=0.d0)
  call parse_input_variable(wmix,"WMIX",inputFILE,default=0.5d0)
  call parse_input_variable(Lfreq,"Lfreq",inputFILE,default=1024)
  call parse_input_variable(wmin,"WMIN",inputFILE,default=-5d0)
  call parse_input_variable(wmax,"wmax",inputFILE,default= 5d0)
  call parse_input_variable(ieta,"iETA",inputFILE,default=4.d-2)
  call parse_input_variable(p_field,"P_FIELD",inputFILE,default=0.01d0)
  call parse_input_variable(it_error,"IT_ERROR",inputFILE,default=1d-5)
  call parse_input_variable(maxiter,"MAXITER",inputFILE,default=100)
  call parse_input_variable(with_lcm,"WITH_lcm",inputFILE,default=.false.)
  call parse_input_variable(with_real_gf,"WITH_REAL_GF",inputFILE,default=.false.)
  call parse_input_variable(spin_blocks,"SPIN_BLOCKS",inputFILE, &
    default=1,comment="Spin representation: 1, 2 or 4")
  call parse_input_variable(method,"METHOD",inputFILE,default="lapack")
  call parse_input_variable(cheb_order,"CHEB_ORDER",inputFILE,default=256)
  call parse_input_variable(cheb_kind,"CHEB_KIND",inputFILE, &
    default="jackson")
  method = to_lower(trim(method))
  if(method/="lapack".and.method/="chebyshev") &
    error stop "Abhz: METHOD must be lapack or chebyshev"
  if(.not.any(spin_blocks==[1,2,4])) &
    error stop "Abhz: SPIN_BLOCKS must be 1, 2 or 4"
  if(spin_blocks==1.and.disorder_type==2) &
    error stop "Abhz: Sz disorder requires SPIN_BLOCKS=2 or 4"
  if(Nx<2)error stop "Abhz: require NX >= 2"
  if(method=="chebyshev")then
    if(cheb_order<1)error stop "Abhz: CHEB_ORDER must be positive"
    if(with_lcm.or.with_real_gf) &
      error stop "Abhz: LCM and real GF require METHOD=lapack"
  endif
  call print_input()
  !

  !Save variables into DMFT_TOOLS memory pool
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
  write(*,*) "Solve homogeneous model with Nk="//str(Nk)
  allocate(Hk(Nso,Nso,Nk))
  call TB_build_model(Hk,hk_model,Nso,[Nx,Ny])
  z2 = hk_to_spin_Chern(Hk,[Nx,Ny])
  write(*,*)"get Z2:",z2
  write(*,*)""
  write(*,*)""
  deallocate(Hk)

  !SETUP THE OUTPUT DIRECTORY
  call set_seed_wd()
  call save_input_file(file_basename(inputFILE))

  ! Build and solve the selected spin representation.
  call start_timer()
  select case(spin_blocks)
  case(1)
    call setup_Abhz(spH=spHup,spin=1)
    spHdown = spHup%conjg()
    call solve_sector(spHup,dens_up,U_up,Ev_up)
    allocate(dens_down,source=dens_up)
    if(method=="lapack")then
      allocate(Ev_down,source=Ev_up)
      allocate(U_down,source=conjg(U_up))
      call assemble_Abhz_full(spHup,spHdown,spHfull)
      call assemble_full_eigensystem()
    endif
  case(2)
    call setup_Abhz(spH=spHup,spin=1)
    call setup_Abhz(verbose=.false.,spH=spHdown,spin=-1)
    call solve_sector(spHup,dens_up,U_up,Ev_up)
    call solve_sector(spHdown,dens_down,U_down,Ev_down)
    if(method=="lapack")then
      call assemble_Abhz_full(spHup,spHdown,spHfull)
      call assemble_full_eigensystem()
    endif
  case(4)
    call setup_Abhz(spH=spHfull,spin=0)
    call solve_sector(spHfull,dens_full,H,Ev)
    allocate(dens_up(2*Nlat),dens_down(2*Nlat))
    do ilat=1,Nlat
      do iorb=1,Norb
        dens_up(iorb+(ilat-1)*Norb) = &
          dens_full(iorb+(ilat-1)*Nso)
        dens_down(iorb+(ilat-1)*Norb) = &
          dens_full(iorb+Norb+(ilat-1)*Nso)
      enddo
    enddo
  end select
  call write_local_observables()
  !
  if(method=="lapack")then
    call push_Bloch(H,Ev)
    deallocate(H,Ev)
    Sz = kron(eye(Nlat),GammaS)
    !
    !Get topological info:
    sp_chern(1) = single_point_spin_chern(spin=1)
    sp_chern(2) = single_point_spin_chern(spin=2)
    call save_array("z2.dat",(sp_chern(1)-sp_chern(2))/2d0 )
    call save_array("spin_chern.dat",sp_chern)
    write(*,*)"spin_Chern UP,DW:",sp_chern

    if(with_real_gf)then
      if(.not.allocated(Gf))allocate(Gf(Nlat,Nso,Nso,Lfreq))
      call get_gf(Gf,'real')
      call write_gf(gf_reshape(Gf,Nspin,Norb,Nlat), &
        "Gloc_"//str(idum),'real',iprint=4,itar=.true.)
      if(allocated(Gf))deallocate(Gf)
    endif
    !
    if(with_lcm)then
      if(bhz_pbc)then
        call pbc_local_spin_chern_marker(spin=1,lcm=LsCM)
        call splot3d("PBC_Local_SpinChern_Marker.dat", &
          dble(arange(1,Nx)),dble(arange(1,Ny)),LsCM)
      else
        call obc_local_spin_chern_marker(spin=1,lcm=LsCM)
        call splot3d("OBC_Local_SpinChern_Marker.dat", &
          dble(arange(1,Nx)),dble(arange(1,Ny)),LsCM)
      endif
    endif
  endif

  !< Free memory:
  call free_seed_solution()
  if(spHup%status)call spHup%free()
  if(spHdown%status)call spHdown%free()
  if(spHfull%status)call spHfull%free()
  call free_Abhz()
  call reset_seed_wd()
  call stop_timer("IDUM: "//str(idum))
  write(*,*)""



contains




  subroutine solve_sector(spH,density,evec,eval)
    ! Compute the diagonal of Theta(mu-H). LAPACK optionally returns
    ! the complete ordered eigensystem; Chebyshev keeps H sparse.
    type(sparse_zmatrix_csr),target,intent(inout) :: spH
    real(8),dimension(:),allocatable,intent(out)  :: density
    complex(8),dimension(:,:),allocatable, &
      optional,intent(out)                       :: evec
    real(8),dimension(:),allocatable, &
      optional,intent(out)                       :: eval
    complex(8),dimension(:,:),allocatable        :: dense
    real(8),dimension(:),allocatable             :: levels,fermi_weight
    integer                                      :: n,io
    !
    n = spH%Nrow
    allocate(density(n))
    call start_timer("Solve Anderson BHZ sector N="//str(n))
    select case(method)
    case("lapack")
      if(.not.present(evec).or..not.present(eval)) &
        error stop "solve_sector: LAPACK requires EVECS and EVALS"
      allocate(dense(n,n),levels(n),fermi_weight(n))
      dense = as_matrix(spH)
      call eigh(dense,levels)
      where(levels<mu)
        fermi_weight = 1d0
      elsewhere(levels>mu)
        fermi_weight = 0d0
      elsewhere
        fermi_weight = 0.5d0
      endwhere
      do io=1,n
        density(io) = sum(fermi_weight*abs(dense(io,:))**2)
      enddo
      allocate(evec,source=dense)
      allocate(eval,source=levels)
      deallocate(dense,levels,fermi_weight)
    case("chebyshev")
      active_spH => spH
      call cheb_density(spH,density)
      nullify(active_spH)
    end select
    call stop_timer()
  end subroutine solve_sector


  subroutine cheb_density(spH,density)
    type(sparse_zmatrix_csr),intent(in) :: spH
    real(8),intent(out)                 :: density(:)
    complex(8),allocatable             :: seed(:)
    complex(8)                         :: value(1),diag
    real(8)                             :: coeff(0:cheb_order),scale,shift
    real(8)                             :: lower,upper,center,radius
    integer                             :: n,io,info
    !
    n = spH%Nrow
    if(size(density)/=n) &
      error stop "cheb_density: invalid DENSITY size"
    lower = huge(1d0)
    upper = -huge(1d0)
    do io=1,n
      diag = spH%get(io,io)
      center = real(diag,8)
      radius = sum(abs(spH%row(io)%vals))-abs(diag)
      lower = min(lower,center-radius)
      upper = max(upper,center+radius)
    enddo
    shift = 0.5d0*(upper+lower)
    scale = max(0.5d0*(upper-lower),1d-12)/0.95d0
    call cheb_fermi_coefficients(scale,shift,mu,coeff,cheb_kind)
    !
    allocate(seed(n),source=zero)
    do io=1,n
      seed(io) = one
      call cheb_projection(MatVec_H,seed,[io],value,info)
      if(info/=0) &
        error stop "Abhz: Chebyshev occupation projection failed"
      density(io) = real(value(1),8)
      seed(io) = zero
    enddo
    deallocate(seed)
    write(*,*) "Chebyshev: degree, scale, shift =", &
      cheb_order,scale,shift
  end subroutine cheb_density


  subroutine assemble_full_eigensystem()
    ! Embed both spin-sector eigensystems in the site-major 4x4 basis
    ! and sort all columns by increasing energy for LCM_SQUARE.
    complex(8),allocatable :: raw_vectors(:,:)
    real(8),allocatable    :: raw_levels(:)
    integer,allocatable    :: order(:)
    integer                :: nblock,ib,ifull,ilat,iorb,istate
    !
    nblock = 2*Nlat
    if(size(Ev_up)/=nblock.or.size(Ev_down)/=nblock) &
      error stop "assemble_full_eigensystem: invalid eigenvalue size"
    allocate(raw_vectors(Nlso,Nlso),source=zero)
    allocate(raw_levels(Nlso),order(Nlso))
    raw_levels(1:nblock) = Ev_up
    raw_levels(nblock+1:Nlso) = Ev_down
    do ib=1,nblock
      ilat = (ib-1)/Norb+1
      iorb = modulo(ib-1,Norb)+1
      ifull = iorb + (ilat-1)*Nso
      raw_vectors(ifull,1:nblock) = U_up(ib,:)
      ifull = iorb + Norb + (ilat-1)*Nso
      raw_vectors(ifull,nblock+1:Nlso) = U_down(ib,:)
    enddo
    call sort(raw_levels,order)
    allocate(Ev,source=raw_levels)
    allocate(H(Nlso,Nlso))
    do istate=1,Nlso
      H(:,istate) = raw_vectors(:,order(istate))
    enddo
    deallocate(raw_vectors,raw_levels,order)
  end subroutine assemble_full_eigensystem


  subroutine write_local_observables()
    ! Combine the two orbital spin sectors into site-resolved output.
    ! For SPIN_BLOCKS=1, DENS_DOWN is the TR partner of DENS_UP.
    ! The reported Tz and Sz include a factor of one half.
    !
    allocate(Nii(Nlat,Nspin,Norb))
    do ilat=1,Nlat
      do iorb=1,Norb
        io = iorb + (ilat-1)*Norb
        Nii(ilat,1,iorb) = dens_up(io)
        Nii(ilat,2,iorb) = dens_down(io)
      enddo
    enddo
    Szii = 0.5d0*(sum(Nii(:,1,:),dim=2)-sum(Nii(:,2,:),dim=2))
    Tzii = 0.5d0*(sum(Nii(:,:,1),dim=2)-sum(Nii(:,:,2),dim=2))
    write(*,*)"E(Tz), sd(Tz):",get_mean(Tzii),get_sd(Tzii)
    write(*,*)"E(Sz), sd(Sz):",get_mean(Szii),get_sd(Szii)
    call save_array("sz_"//str(idum)//".dat",Szii)
    call save_array("tz_"//str(idum)//".dat",Tzii)
  end subroutine write_local_observables


  subroutine free_seed_solution()
    ! Release all arrays owned by one disorder realization. COMMON
    ! eigensystem arrays are released separately by free_Abhz.
    !
    if(associated(active_spH))nullify(active_spH)
    if(allocated(dens_up))deallocate(dens_up)
    if(allocated(dens_down))deallocate(dens_down)
    if(allocated(dens_full))deallocate(dens_full)
    if(allocated(Ev_up))deallocate(Ev_up)
    if(allocated(Ev_down))deallocate(Ev_down)
    if(allocated(U_up))deallocate(U_up)
    if(allocated(U_down))deallocate(U_down)
    if(allocated(H))deallocate(H)
    if(allocated(Ev))deallocate(Ev)
    if(allocated(Nii))deallocate(Nii)
    if(allocated(Tzii))deallocate(Tzii)
    if(allocated(Szii))deallocate(Szii)
    if(allocated(Gf))deallocate(Gf)
  end subroutine free_seed_solution


  subroutine MatVec_H(Nloc,vin,vout)
    integer    :: Nloc
    complex(8) :: vin(Nloc),vout(Nloc)
    if(.not.associated(active_spH)) &
      error stop "MatVec_H: no active sparse Hamiltonian"
    vout = active_spH%dot(vin)
  end subroutine MatVec_H
end program anderson_bhz_2d
