MODULE BHZ_SPECTRA
  USE COMMON
  implicit none
  private

  integer,parameter                     :: Nchannel=5
  character(len=16),parameter           :: channel_name(Nchannel) = [ &
    character(len=16) :: "up_orb1","up_orb2","down_orb1", &
    "down_orb2","total"]

  type :: spectrum_statistics
    integer                              :: count=0
    real(8),allocatable                  :: mean(:),m2(:)
    real(8),allocatable                  :: mean_log(:),m2_log(:)
    real(8),allocatable                  :: raw_mean(:),raw_min(:)
    integer,allocatable                  :: negative_count(:)
    integer,allocatable                  :: floor_count(:)
  end type spectrum_statistics

  character(len=16)                     :: spectrum_method="lapack"
  integer                               :: spectrum_order=1024
  real(8)                               :: rho_floor=1d-300
  real(8)                               :: negative_tolerance=1d-12
  type(sparse_zmatrix_csr),pointer      :: active_spH=>null()

  public :: configure_Abhz_spectra
  public :: write_Abhz_spectrum_statistics

contains


  ! Configure the post-processing backend. METHOD accepts LAPACK or
  ! CHEBYCHEV. ORDER is used only by CHEBYCHEV. RHO_FLOOR enters only
  ! log(rho), while NEGATIVE_TOLERANCE controls the diagnostic count;
  ! negative spectral values are clipped to zero for linear moments.
  subroutine configure_Abhz_spectra(method,order,floor,tolerance)
    character(len=*),intent(in) :: method
    integer,intent(in)          :: order
    real(8),intent(in)          :: floor,tolerance
    !
    select case(to_lower(str(method(1:1))))
     case("l");spectrum_method="lapack"
     case("c");spectrum_method="chebychev"
     case default
      error stop "configure_Abhz_spectra: invalid METHOD"
    end select
    if(order<1) &
      error stop "configure_Abhz_spectra: ORDER must be positive"
    if(floor<=0d0) &
      error stop "configure_Abhz_spectra: FLOOR must be positive"
    if(tolerance<0d0) &
      error stop "configure_Abhz_spectra: invalid TOLERANCE"
    spectrum_order      = order
    rho_floor           = floor
    negative_tolerance  = tolerance
  end subroutine configure_Abhz_spectra


  ! Compute reduced local spectral statistics for the selected spin
  ! representation. For SPIN_BLOCKS=1/2, SPHUP is required and
  ! SPHDOWN is additionally required in mode 2. For SPIN_BLOCKS=0,
  ! SPHFULL is required. Only O(Nchannel*Lfreq) statistics and one
  ! local spectrum are retained; no site-resolved spectrum is written.
  subroutine write_Abhz_spectrum_statistics(source,spHup,spHdown,spHfull)
    character(len=*),intent(in)                       :: source
    type(sparse_zmatrix_csr),target,optional,intent(in) :: spHup
    type(sparse_zmatrix_csr),target,optional,intent(in) :: spHdown
    type(sparse_zmatrix_csr),target,optional,intent(in) :: spHfull
    type(spectrum_statistics)                         :: stats(Nchannel)
    real(8)                                           :: omega(Lfreq)
    integer                                           :: ichannel
    !
    ! Check the matrix contract before numerical work.
    select case(spin_blocks)
     case(1)
      if(.not.present(spHup)) &
        error stop "write_Abhz_spectrum_statistics: missing SPHUP"
     case(2)
      if(.not.present(spHup).or..not.present(spHdown)) &
        error stop "write_Abhz_spectrum_statistics: missing spin block"
     case(0)
      if(.not.present(spHfull)) &
        error stop "write_Abhz_spectrum_statistics: missing SPHFULL"
     case default
      error stop "write_Abhz_spectrum_statistics: invalid SPIN_BLOCKS"
    end select
    !
    omega = linspace(wmin,wmax,Lfreq)
    do ichannel=1,Nchannel
      call initialize_statistics(stats(ichannel))
    enddo
    !
    ! Evaluate one local four-channel spectrum at a time.
    select case(spectrum_method)
     case("lapack")
      call spectra_lapack(stats,omega,spHup,spHdown,spHfull)
     case("chebychev")
      call spectra_chebychev(stats,omega,spHup,spHdown,spHfull)
     case default
      error stop "write_Abhz_spectrum_statistics: invalid backend"
    end select
    !
    ! Persist only the sufficient per-frequency statistics.
    do ichannel=1,Nchannel
      if(any(stats(ichannel)%negative_count>0)) &
        write(*,*)"Spectrum warning: significant negative LDOS in ", &
          trim(channel_name(ichannel))
    enddo
    call write_statistics_files(stats,omega)
    call write_spectrum_metadata(source)
    do ichannel=1,Nchannel
      call free_statistics(stats(ichannel))
    enddo
  end subroutine write_Abhz_spectrum_statistics


  ! Diagonalize the selected matrix representation and reduce each
  ! local diagonal resolvent immediately after it is formed.
  subroutine spectra_lapack(stats,omega,spHup,spHdown,spHfull)
    type(spectrum_statistics),intent(inout)             :: stats(:)
    real(8),intent(in)                                  :: omega(:)
    type(sparse_zmatrix_csr),optional,intent(in)         :: spHup
    type(sparse_zmatrix_csr),optional,intent(in)         :: spHdown
    type(sparse_zmatrix_csr),optional,intent(in)         :: spHfull
    complex(8),allocatable                              :: Uup(:,:),Udown(:,:)
    complex(8),allocatable                              :: Ufull(:,:)
    real(8),allocatable                                 :: Eup(:),Edown(:),Efull(:)
    real(8)                                             :: rho(Nchannel,Lfreq)
    integer                                             :: ilat,iorb,row
    !
    select case(spin_blocks)
     case(1,2)
      call dense_eigensystem(spHup,Uup,Eup)
      if(spin_blocks==2)call dense_eigensystem(spHdown,Udown,Edown)
      do ilat=1,Nlat
        do iorb=1,2
          row = iorb+(ilat-1)*2
          call lapack_local_spectrum(Uup(row,:),Eup,omega,rho(iorb,:))
          if(spin_blocks==1)then
            rho(iorb+2,:) = rho(iorb,:)
          else
            call lapack_local_spectrum(Udown(row,:),Edown,omega, &
              rho(iorb+2,:))
          endif
        enddo
        rho(5,:) = sum(rho(1:4,:),dim=1)
        call update_all_statistics(stats,rho)
      enddo
     case(0)
      call dense_eigensystem(spHfull,Ufull,Efull)
      do ilat=1,Nlat
        do iorb=1,4
          row = iorb+(ilat-1)*4
          call lapack_local_spectrum(Ufull(row,:),Efull,omega, &
            rho(iorb,:))
        enddo
        rho(5,:) = sum(rho(1:4,:),dim=1)
        call update_all_statistics(stats,rho)
      enddo
    end select
    !
    if(allocated(Uup))deallocate(Uup,Eup)
    if(allocated(Udown))deallocate(Udown,Edown)
    if(allocated(Ufull))deallocate(Ufull,Efull)
  end subroutine spectra_lapack


  ! Materialize one dense copy and return its ordered Hermitian
  ! eigensystem. The caller owns the returned arrays.
  subroutine dense_eigensystem(spH,evec,eval)
    type(sparse_zmatrix_csr),intent(in)       :: spH
    complex(8),allocatable,intent(out)        :: evec(:,:)
    real(8),allocatable,intent(out)           :: eval(:)
    integer                                   :: n
    !
    n = spH%Nrow
    if(n<1.or.spH%Ncol/=n) &
      error stop "dense_eigensystem: invalid SPH"
    allocate(evec(n,n),eval(n))
    evec = as_matrix(spH)
    call eigh(evec,eval)
  end subroutine dense_eigensystem


  ! Evaluate a positive Lorentz-broadened diagonal spectral function
  ! from one row of a complete eigensystem.
  subroutine lapack_local_spectrum(weights,eval,omega,rho)
    complex(8),intent(in) :: weights(:)
    real(8),intent(in)    :: eval(:),omega(:)
    real(8),intent(out)   :: rho(:)
    integer               :: iw
    !
    if(size(weights)/=size(eval).or.size(rho)/=size(omega)) &
      error stop "lapack_local_spectrum: invalid shape"
    do iw=1,size(omega)
      rho(iw) = sum(abs(weights)**2*ieta/pi/ &
        ((omega(iw)+mu-eval)**2+ieta**2))
    enddo
  end subroutine lapack_local_spectrum


  ! Compute one diagonal Chebyshev moment vector at a time. Coefficients
  ! are retained only for the common frequency grid and one spin sector.
  subroutine spectra_chebychev(stats,omega,spHup,spHdown,spHfull)
    type(spectrum_statistics),intent(inout)              :: stats(:)
    real(8),intent(in)                                   :: omega(:)
    type(sparse_zmatrix_csr),target,optional,intent(in)   :: spHup
    type(sparse_zmatrix_csr),target,optional,intent(in)   :: spHdown
    type(sparse_zmatrix_csr),target,optional,intent(in)   :: spHfull
    complex(8),allocatable                               :: coeff_up(:,:)
    complex(8),allocatable                               :: coeff_down(:,:)
    complex(8),allocatable                               :: coeff_full(:,:)
    real(8)                                              :: scale_up,shift_up
    real(8)                                              :: scale_down,shift_down
    real(8)                                              :: scale_full,shift_full
    real(8)                                              :: rho(Nchannel,Lfreq)
    integer                                              :: ilat,iorb,row
    !
    select case(spin_blocks)
     case(1,2)
      call prepare_chebychev(spHup,omega,scale_up,shift_up,coeff_up)
      if(spin_blocks==2) &
        call prepare_chebychev(spHdown,omega,scale_down,shift_down, &
          coeff_down)
      do ilat=1,Nlat
        do iorb=1,2
          row = iorb+(ilat-1)*2
          call chebychev_local_spectrum(spHup,row,scale_up,shift_up, &
            coeff_up,rho(iorb,:))
          if(spin_blocks==1)then
            rho(iorb+2,:) = rho(iorb,:)
          else
            call chebychev_local_spectrum(spHdown,row,scale_down, &
              shift_down,coeff_down,rho(iorb+2,:))
          endif
        enddo
        rho(5,:) = sum(rho(1:4,:),dim=1)
        call update_all_statistics(stats,rho)
        call eta(ilat,Nlat)
      enddo
     case(0)
      call prepare_chebychev(spHfull,omega,scale_full,shift_full, &
        coeff_full)
      do ilat=1,Nlat
        do iorb=1,4
          row = iorb+(ilat-1)*4
          call chebychev_local_spectrum(spHfull,row,scale_full, &
            shift_full,coeff_full,rho(iorb,:))
        enddo
        rho(5,:) = sum(rho(1:4,:),dim=1)
        call update_all_statistics(stats,rho)
        call eta(ilat,Nlat)
      enddo
    end select
    !
    if(allocated(coeff_up))deallocate(coeff_up)
    if(allocated(coeff_down))deallocate(coeff_down)
    if(allocated(coeff_full))deallocate(coeff_full)
  end subroutine spectra_chebychev


  ! Estimate safe spectral bounds and tabulate resolvent coefficients
  ! for the complete output grid. The coefficient table is independent
  ! of the system size.
  subroutine prepare_chebychev(spH,omega,scale,shift,coeff)
    type(sparse_zmatrix_csr),intent(in) :: spH
    real(8),intent(in)                  :: omega(:)
    real(8),intent(out)                 :: scale,shift
    complex(8),allocatable,intent(out)  :: coeff(:,:)
    integer                             :: iw
    !
    call spectrum_rescaling(spH,scale,shift)
    allocate(coeff(0:spectrum_order,size(omega)))
    do iw=1,size(omega)
      call resolvent_coefficients(dcmplx(omega(iw)+mu,ieta), &
        scale,shift,coeff(:,iw))
    enddo
    write(*,*)"Spectrum Chebychev: degree, scale, shift =", &
      spectrum_order,scale,shift
  end subroutine prepare_chebychev


  ! Return the exact truncated Chebyshev coefficients of
  ! 1/(Z-E), E=SHIFT+SCALE*x. With zeta=(Z-SHIFT)/SCALE,
  ! c_0=1/(SCALE*sqrt(zeta^2-1)) and
  ! c_m=2*q^m/(SCALE*sqrt(zeta^2-1)), q=zeta-sqrt(zeta^2-1).
  ! The square-root sign is selected so that abs(q)<1.
  subroutine resolvent_coefficients(z,scale,shift,coeff)
    complex(8),intent(in)  :: z
    real(8),intent(in)     :: scale,shift
    complex(8),intent(out) :: coeff(0:)
    complex(8)             :: zeta,root,q,power
    integer                :: order
    !
    if(scale<=0d0.or.size(coeff)<1) &
      error stop "resolvent_coefficients: invalid input"
    zeta = (z-shift)/scale
    root = sqrt(zeta*zeta-one)
    q    = zeta-root
    if(abs(q)>1d0)then
      root = -root
      q    = zeta-root
    endif
    if(abs(q)>=1d0) &
      error stop "resolvent_coefficients: invalid resolvent branch"
    coeff(0) = one/(scale*root)
    power = one
    do order=1,ubound(coeff,1)
      power        = power*q
      coeff(order) = 2d0*power/(scale*root)
    enddo
  end subroutine resolvent_coefficients


  ! Return one diagonal resolvent on the configured frequency grid.
  subroutine chebychev_local_spectrum(spH,row,scale,shift,coeff,rho)
    type(sparse_zmatrix_csr),target,intent(in) :: spH
    integer,intent(in)                         :: row
    real(8),intent(in)                         :: scale,shift
    complex(8),intent(in)                      :: coeff(0:,:)
    real(8),intent(out)                        :: rho(:)
    complex(8),allocatable                     :: seed(:)
    complex(8)                                 :: moment(1,0:spectrum_order)
    complex(8)                                 :: gloc(size(rho))
    integer                                    :: info
    !
    if(row<1.or.row>spH%Nrow.or.size(coeff,2)/=size(rho)) &
      error stop "chebychev_local_spectrum: invalid input"
    allocate(seed(spH%Nrow),source=zero)
    seed(row) = one
    active_spH => spH
    call cheb_set(scale,shift)
    call cheb_moments(MatVec_spectrum,seed,[row],moment,info)
    nullify(active_spH)
    if(info/=0) &
      error stop "chebychev_local_spectrum: moment recurrence failed"
    gloc = matmul(moment(1,:),coeff)
    rho  = -dimag(gloc)/pi
    deallocate(seed)
  end subroutine chebychev_local_spectrum


  ! Gershgorin bounds with a five-percent rescaling margin.
  subroutine spectrum_rescaling(spH,scale,shift)
    type(sparse_zmatrix_csr),intent(in) :: spH
    real(8),intent(out)                 :: scale,shift
    complex(8)                          :: diagonal
    real(8)                             :: lower,upper,center,radius
    integer                             :: row
    !
    if(spH%Nrow<1.or.spH%Ncol/=spH%Nrow) &
      error stop "spectrum_rescaling: invalid SPH"
    lower = huge(1d0)
    upper = -huge(1d0)
    do row=1,spH%Nrow
      diagonal = spH%get(row,row)
      center   = dreal(diagonal)
      radius   = sum(abs(spH%row(row)%vals))-abs(diagonal)
      lower    = min(lower,center-radius)
      upper    = max(upper,center+radius)
    enddo
    shift = 0.5d0*(upper+lower)
    scale = max(0.5d0*(upper-lower),1d-12)/0.95d0
  end subroutine spectrum_rescaling


  ! Sparse matrix-vector callback for the active spectral sector.
  subroutine MatVec_spectrum(nloc,vin,vout)
    integer    :: nloc
    complex(8) :: vin(nloc),vout(nloc)
    !
    if(.not.associated(active_spH)) &
      error stop "MatVec_spectrum: no active Hamiltonian"
    vout = active_spH%dot(vin)
  end subroutine MatVec_spectrum


  ! Allocate and zero one channel accumulator.
  subroutine initialize_statistics(stats)
    type(spectrum_statistics),intent(out) :: stats
    !
    allocate(stats%mean(Lfreq),source=0d0)
    allocate(stats%m2(Lfreq),source=0d0)
    allocate(stats%mean_log(Lfreq),source=0d0)
    allocate(stats%m2_log(Lfreq),source=0d0)
    allocate(stats%raw_mean(Lfreq),source=0d0)
    allocate(stats%raw_min(Lfreq),source=huge(1d0))
    allocate(stats%negative_count(Lfreq),source=0)
    allocate(stats%floor_count(Lfreq),source=0)
    stats%count = 0
  end subroutine initialize_statistics


  ! Update the five channel accumulators with one site's spectra.
  subroutine update_all_statistics(stats,rho)
    type(spectrum_statistics),intent(inout) :: stats(:)
    real(8),intent(in)                      :: rho(:,:)
    integer                                 :: ichannel
    !
    if(size(stats)/=Nchannel.or.size(rho,1)/=Nchannel.or. &
      size(rho,2)/=Lfreq) &
      error stop "update_all_statistics: invalid shape"
    do ichannel=1,Nchannel
      call update_statistics(stats(ichannel),rho(ichannel,:))
    enddo
  end subroutine update_all_statistics


  ! Apply a vector Welford update. Negative values are diagnosed and
  ! clipped to zero for rho moments. RHO_FLOOR is used only for logs.
  subroutine update_statistics(stats,raw_rho)
    type(spectrum_statistics),intent(inout) :: stats
    real(8),intent(in)                      :: raw_rho(:)
    real(8)                                 :: rho(size(raw_rho))
    real(8)                                 :: logrho(size(raw_rho))
    real(8)                                 :: delta(size(raw_rho))
    integer                                 :: new_count
    !
    if(size(raw_rho)/=Lfreq) &
      error stop "update_statistics: invalid spectrum size"
    new_count = stats%count+1
    stats%negative_count = stats%negative_count &
      +merge(1,0,raw_rho < -negative_tolerance)
    stats%raw_min = min(stats%raw_min,raw_rho)
    delta = raw_rho-stats%raw_mean
    stats%raw_mean = stats%raw_mean+delta/new_count
    !
    rho = max(raw_rho,0d0)
    stats%floor_count = stats%floor_count+merge(1,0,rho<rho_floor)
    logrho = log(max(rho,rho_floor))
    delta = rho-stats%mean
    stats%mean = stats%mean+delta/new_count
    stats%m2 = stats%m2+delta*(rho-stats%mean)
    delta = logrho-stats%mean_log
    stats%mean_log = stats%mean_log+delta/new_count
    stats%m2_log = stats%m2_log+delta*(logrho-stats%mean_log)
    stats%count = new_count
  end subroutine update_statistics


  ! Write compact statistics and diagnostics with numbered columns.
  subroutine write_statistics_files(stats,omega)
    type(spectrum_statistics),intent(in) :: stats(:)
    real(8),intent(in)                   :: omega(:)
    character(len=128)                   :: filename
    integer                              :: ichannel,iw,unit
    !
    do ichannel=1,Nchannel
      filename = "spectrum_"//trim(channel_name(ichannel))//"_stats.dat"
      open(newunit=unit,file=trim(filename),status="replace")
      write(unit,'(A)') &
        "# Reduced local spectral statistics; population M2"
      write(unit,'(A)') &
        "# 1.omega 2.count 3.mean_rho 4.M2_rho 5.mean_logrho 6.M2_logrho"
      do iw=1,Lfreq
        write(unit,'(ES26.17E3,1X,I0,1X,4(ES26.17E3,1X))') &
          omega(iw),stats(ichannel)%count, &
          stats(ichannel)%mean(iw),stats(ichannel)%m2(iw), &
          stats(ichannel)%mean_log(iw),stats(ichannel)%m2_log(iw)
      enddo
      close(unit)
      !
      filename = "spectrum_"//trim(channel_name(ichannel))// &
        "_diagnostics.dat"
      open(newunit=unit,file=trim(filename),status="replace")
      write(unit,'(A)') &
        "# Positivity and logarithmic-floor diagnostics"
      write(unit,'(A)') &
        "# 1.omega 2.negative_count 3.raw_min 4.floor_count 5.raw_mean"
      do iw=1,Lfreq
        write(unit,'(ES26.17E3,1X,I0,1X,ES26.17E3,1X,I0,1X,ES26.17E3)') &
          omega(iw),stats(ichannel)%negative_count(iw), &
          stats(ichannel)%raw_min(iw), &
          stats(ichannel)%floor_count(iw), &
          stats(ichannel)%raw_mean(iw)
      enddo
      close(unit)
    enddo
  end subroutine write_statistics_files


  ! Record all controls needed to reject incompatible seed statistics.
  subroutine write_spectrum_metadata(source)
    character(len=*),intent(in) :: source
    integer                     :: unit
    !
    open(newunit=unit,file="spectrum_metadata.dat",status="replace")
    write(unit,'(A)')"format spectrum_stats_v1"
    write(unit,'(A)')"source "//trim(to_lower(source))
    write(unit,'(A)')"method "//trim(spectrum_method)
    write(unit,'(A,I0)')"seed ",idum
    write(unit,'(A,I0)')"nlat ",Nlat
    write(unit,'(A,I0)')"lfreq ",Lfreq
    write(unit,'(A,ES26.17E3)')"wmin ",wmin
    write(unit,'(A,ES26.17E3)')"wmax ",wmax
    write(unit,'(A,ES26.17E3)')"mu ",mu
    write(unit,'(A,ES26.17E3)')"ieta ",ieta
    write(unit,'(A,I0)')"spectrum_order ",spectrum_order
    write(unit,'(A,ES26.17E3)')"rho_floor ",rho_floor
    write(unit,'(A,ES26.17E3)') &
      "negative_tolerance ",negative_tolerance
    write(unit,'(A,I0)')"spin_blocks ",spin_blocks
    write(unit,'(A)') &
      "channels up_orb1 up_orb2 down_orb1 down_orb2 total"
    close(unit)
  end subroutine write_spectrum_metadata


  ! Release one channel accumulator.
  subroutine free_statistics(stats)
    type(spectrum_statistics),intent(inout) :: stats
    !
    if(allocated(stats%mean))deallocate(stats%mean,stats%m2)
    if(allocated(stats%mean_log)) &
      deallocate(stats%mean_log,stats%m2_log)
    if(allocated(stats%raw_mean)) &
      deallocate(stats%raw_mean,stats%raw_min)
    if(allocated(stats%negative_count)) &
      deallocate(stats%negative_count,stats%floor_count)
    stats%count = 0
  end subroutine free_statistics


end module BHZ_SPECTRA
