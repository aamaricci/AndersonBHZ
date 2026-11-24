!> Rebuild one free or converged mean-field Anderson-BHZ Hamiltonian
!> and write compact site/disorder spectral statistics.
program spectrum_bhz_2d_Anderson
  USE COMMON
  USE BHZ_BUILD
  USE BHZ_SPECTRA
  implicit none

  character(len=16)                :: spectrum_source
  character(len=16)                :: spectrum_backend
  integer                          :: spectrum_order
  real(8)                          :: spectrum_floor
  real(8)                          :: spectrum_negative_tolerance
  real(8),allocatable              :: params(:,:)
  type(sparse_zmatrix_csr),target  :: spHup,spHdown,spHfull
  type(sparse_zmatrix_csr),target  :: spHmf_up,spHmf_down,spHmf_full
  logical                          :: mean_field

  ! Read the common Hamiltonian contract and spectral controls.
  call parse_input_file(inputFILE,"inputABHZ.conf")
  call read_input(inputFILE)
  call parse_input_variable(spectrum_source,"SPECTRUM_SOURCE", &
    inputFILE,default="free",comment="Free or mean_field")
  call parse_input_variable(spectrum_backend,"SPECTRUM_METHOD", &
    inputFILE,default=method,comment="Lapack or Chebychev")
  call parse_input_variable(spectrum_order,"SPECTRUM_ORDER", &
    inputFILE,default=max(cheb_order,1024))
  call parse_input_variable(spectrum_floor,"SPECTRUM_RHO_FLOOR", &
    inputFILE,default=1d-300)
  call parse_input_variable(spectrum_negative_tolerance, &
    "SPECTRUM_NEGATIVE_TOLERANCE",inputFILE,default=1d-12)
  call normalize_source(spectrum_source,mean_field)
  call configure_Abhz_spectra(spectrum_backend,spectrum_order, &
    spectrum_floor,spectrum_negative_tolerance)

  ! Recreate exactly the seed-local Hamiltonian.
  call set_seed_wd()
  call save_input_file("spectrum_"//file_basename(inputFILE))
  call TB_set_ei([1d0,0d0],[0d0,1d0])
  call TB_set_bk([pi2,0d0],[0d0,pi2])
  select case(spin_blocks)
   case(1)
    call setup_Abhz(spH=spHup,spin=1)
   case(2)
    call setup_Abhz(spH=spHup,spin=1)
    call setup_Abhz(spH=spHdown,spin=-1,verbose=.false.)
   case(0)
    call setup_Abhz(spH=spHfull,spin=0)
  end select

  ! The MF spectrum is one additional solve at the final mixed fields.
  if(mean_field)then
    call read_converged_parameters(params)
    select case(spin_blocks)
     case(1)
      spHmf_up = build_mf_Abhz(spHup,params,spin=1)
     case(2)
      spHmf_up   = build_mf_Abhz(spHup,params,spin=1)
      spHmf_down = build_mf_Abhz(spHdown,params,spin=-1)
     case(0)
      spHmf_full = build_mf_Abhz(spHfull,params,spin=0)
    end select
  endif

  ! Reduce all local spectra without persisting site-frequency arrays.
  if(mean_field)then
    call write_Abhz_spectrum_statistics("mean_field", &
      spHmf_up,spHmf_down,spHmf_full)
  else
    call write_Abhz_spectrum_statistics("free",spHup,spHdown,spHfull)
  endif

  ! Release seed-local state.
  if(spHup%status)call spHup%free()
  if(spHdown%status)call spHdown%free()
  if(spHfull%status)call spHfull%free()
  if(spHmf_up%status)call spHmf_up%free()
  if(spHmf_down%status)call spHmf_down%free()
  if(spHmf_full%status)call spHmf_full%free()
  if(allocated(params))deallocate(params)
  call free_Abhz()
  call reset_seed_wd()

contains


  ! Normalize the source name and expose whether an MF restart is used.
  subroutine normalize_source(source,is_mean_field)
    character(len=*),intent(inout) :: source
    logical,intent(out)            :: is_mean_field
    !
    source = to_lower(trim(source))
    select case(trim(source))
     case("free","noninteracting","non_interacting")
      source        = "free"
      is_mean_field = .false.
     case("mf","meanfield","mean_field")
      source        = "mean_field"
      is_mean_field = .true.
     case default
      error stop "spectrum_Abhz: invalid SPECTRUM_SOURCE"
    end select
  end subroutine normalize_source


  ! Read final mixed fields only after proving that the corresponding
  ! MF history ended below IT_ERROR. No symmetry-breaking field is
  ! added during post-processing.
  subroutine read_converged_parameters(a)
    real(8),allocatable,intent(out) :: a(:,:)
    character(len=1024)             :: line
    integer                         :: unit,status,iteration,last_iteration
    real(8)                         :: error,last_error
    logical                         :: exists
    !
    inquire(file="params.restart",exist=exists)
    if(.not.exists) &
      error stop "spectrum_Abhz: params.restart not found"
    inquire(file="mf_history.dat",exist=exists)
    if(.not.exists) &
      error stop "spectrum_Abhz: mf_history.dat not found"
    !
    last_iteration = 0
    last_error     = huge(1d0)
    open(newunit=unit,file="mf_history.dat",status="old",action="read")
    do
      read(unit,'(A)',iostat=status)line
      if(status/=0)exit
      line = adjustl(line)
      if(len_trim(line)==0.or.line(1:1)=="#")cycle
      read(line,*,iostat=status)iteration,error
      if(status/=0)cycle
      last_iteration = iteration
      last_error     = error
    enddo
    close(unit)
    if(last_iteration<1) &
      error stop "spectrum_Abhz: empty mean-field history"
    if(last_error>=it_error) &
      error stop "spectrum_Abhz: mean-field solution did not converge"
    !
    allocate(a(Nlat,2))
    call read_array("params.restart",a)
    write(*,*)"Spectrum MF restart: iteration, error =", &
      last_iteration,last_error
  end subroutine read_converged_parameters


end program spectrum_bhz_2d_Anderson
