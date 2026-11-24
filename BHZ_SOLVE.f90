MODULE BHZ_SOLVE
  USE SCIFOR
  USE DMFT_TOOLS
  USE COMMON
  implicit none
  private

  complex(8),dimension(:,:),allocatable,public :: U,U_up,U_down
  complex(8),dimension(:,:),allocatable,public :: Uocc_up,Uocc_down
  complex(8),dimension(:,:),allocatable,public :: PSzP
  real(8),dimension(:),allocatable,public      :: E,Ev_up,Ev_down
  real(8),dimension(:),allocatable,public      :: Epsp
  integer,public                               :: Nocc=0
  integer,public                               :: Nocc_up=0,Nocc_down=0

  public :: solve_Abhz
  public :: free_Abhz_solution
  public :: push_Bloch
  public :: count_occupied_levels
  public :: check_Pgap
  public :: check_Egap
  public :: get_gf

  type(sparse_zmatrix_csr),pointer :: active_spH=>null()

contains


  subroutine cheb_fermi_coefficients(scale,shift,fermi_level,coeff,kind)
    ! Return coefficients of Theta(FERMI_LEVEL-E) and activate them in
    ! SciFortran. KIND accepts NONE or JACKSON and defaults to JACKSON.
    real(8),intent(in)                   :: scale,shift,fermi_level
    real(8),dimension(0:),intent(out)    :: coeff
    character(len=*),optional,intent(in) :: kind
    character(len=16)                    :: kernel
    real(8)                              :: xfermi,theta,angle,g
    integer                              :: m,order
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
      error stop "cheb_fermi_coefficients: invalid KIND"
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
    angle    = acos(xfermi)
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


  subroutine solve_Abhz(spH,density,evec,eval)
    ! Compute diag[Theta(MU-SPH)] at zero temperature. SPH remains
    ! owned by the caller. LAPACK uses a temporary dense copy and may
    ! optionally return the complete ordered eigensystem in EVEC and
    ! EVAL. Chebychev keeps SPH sparse and returns only DENSITY.
    type(sparse_zmatrix_csr),target,intent(in)   :: spH
    real(8),dimension(:),allocatable,intent(out) :: density
    complex(8),dimension(:,:),allocatable, &
      optional,intent(out)                       :: evec
    real(8),dimension(:),allocatable, &
      optional,intent(out)                       :: eval
    complex(8),dimension(:,:),allocatable        :: dense
    real(8),dimension(:),allocatable             :: levels,fermi_weight
    integer                                      :: n,io
    !
    ! Check the Hamiltonian and optional eigensystem contract.
    n = spH%Nrow
    if(n<1.or.spH%Ncol/=n) &
      error stop "solve_Abhz: SPH must be a nonempty square matrix"
    if(present(evec).neqv.present(eval)) &
      error stop "solve_Abhz: EVEC and EVAL must be both present"
    !
    ! Solve with the globally selected representation.
    allocate(density(n))
    call start_timer("Solve Anderson BHZ sector N="//str(n))
    select case(method)
    case("lapack")
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
      if(present(evec))then
        allocate(evec,source=dense)
        allocate(eval,source=levels)
      endif
      deallocate(dense,levels,fermi_weight)
    case("chebychev")
      active_spH => spH
      call cheb_density(spH,density)
      nullify(active_spH)
    case default
      error stop "solve_Abhz: invalid METHOD"
    end select
    call stop_timer()
  end subroutine solve_Abhz


  subroutine cheb_density(spH,density)
    ! Compute the diagonal of the Chebychev Fermi projector for SPH.
    type(sparse_zmatrix_csr),intent(in) :: spH
    real(8),intent(out)                 :: density(:)
    complex(8),allocatable             :: seed(:)
    complex(8)                         :: value(1),diag
    real(8)                            :: coeff(0:cheb_order)
    real(8)                            :: scale,shift,lower,upper
    real(8)                            :: center,radius
    integer                            :: n,io,info
    !
    ! Check dimensions and estimate spectral bounds by Gershgorin.
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
    ! Evaluate each diagonal projector element with a basis seed.
    allocate(seed(n),source=zero)
    do io=1,n
      seed(io) = one
      call cheb_projection(MatVec_Abhz,seed,[io],value,info)
      if(info/=0) &
        error stop "cheb_density: Chebychev projection failed"
      density(io) = real(value(1),8)
      seed(io) = zero
    enddo
    deallocate(seed)
    write(*,*) "Chebychev: degree, scale, shift =", &
      cheb_order,scale,shift
  end subroutine cheb_density


  subroutine MatVec_Abhz(Nloc,vin,vout)
    integer    :: Nloc
    complex(8) :: vin(Nloc),vout(Nloc)
    !
    if(.not.associated(active_spH)) &
      error stop "MatVec_Abhz: no active sparse Hamiltonian"
    vout = active_spH%dot(vin)
  end subroutine MatVec_Abhz


  subroutine push_Bloch(Evec,Evals)
    ! Store a complete ordered eigensystem in the module state.
    complex(8),dimension(:,:),intent(in) :: Evec
    real(8),dimension(:),intent(in)      :: Evals
    !
    call check_dimension("push_Bloch")
    if(size(Evec,1)/=Nlso.or.size(Evec,2)/=Nlso) &
      error stop "push_Bloch: invalid EVEC shape"
    if(size(Evals)/=Nlso) &
      error stop "push_Bloch: invalid EVALS size"
    if(allocated(U).or.allocated(E)) &
      error stop "push_Bloch: eigensystem already allocated"
    allocate(U,source=Evec)
    allocate(E,source=Evals)
    Nocc = count_occupied_levels(E,"full")
    call save_array("E.dat",E)
    call check_Egap("push_Bloch")
  end subroutine push_Bloch


  function count_occupied_levels(eval,sector) result(noccupied)
    ! Count zero-temperature levels below MU. A state at MU is
    ! reported and excluded because it does not define a projector.
    real(8),dimension(:),intent(in) :: eval
    character(len=*),intent(in)     :: sector
    integer                         :: noccupied
    real(8),parameter               :: mu_tol=1d-12
    !
    if(any(abs(eval-mu)<=mu_tol)) &
      write(*,*)"occupied levels warning: level at MU in ",trim(sector)
    noccupied = count(eval<mu-mu_tol)
    if(noccupied<1.or.noccupied>=size(eval)) &
      error stop "count_occupied_levels: empty or complete sector"
  end function count_occupied_levels


  subroutine check_Pgap(caller)
    character(len=*),intent(in) :: caller
    integer                     :: nminus,nplus
    real(8)                     :: Ep,Em,Pgap
    real(8),parameter           :: spin_tol=1d-10
    !
    if(.not.allocated(Epsp)) &
      error stop "check_Pgap: projected-spin spectrum is not set"
    nminus = count(Epsp < -spin_tol)
    nplus  = count(Epsp >  spin_tol)
    if(nminus+nplus/=Nocc)then
      write(*,*)str(caller)//" warning: closing of the PSzP spectrum"
      return
    endif
    if(nminus<1.or.nplus<1)then
      write(*,*)str(caller)//" warning: empty PSzP sector"
      return
    endif
    Em   = maxval(Epsp,mask=Epsp<0d0)
    Ep   = minval(Epsp,mask=Epsp>0d0)
    Pgap = Ep-Em
    write(*,"(A,G21.13)")"PSzP gap=",Pgap
    call save_array("PSzPgap.dat",Pgap)
    if(Pgap<2d0*spin_tol) &
      write(*,*)str(caller)//" warning: closing of the PSzP spectrum"
  end subroutine check_Pgap


  subroutine check_Egap(caller)
    character(len=*),intent(in) :: caller
    real(8)                     :: Ep,Em,Egap
    !
    if(.not.allocated(E)) &
      error stop "check_Egap: spectrum is not set"
    if(Nocc<1.or.Nocc>=size(E)) &
      error stop "check_Egap: invalid occupation"
    Ep   = E(Nocc+1)
    Em   = E(Nocc)
    Egap = Ep-Em
    write(*,"(A,G21.13)")"E gap=",Egap
    call save_array("Egap.dat",Egap)
  end subroutine check_Egap


  subroutine get_gf(Gloc,axis)
    ! Evaluate the local Green function from the stored eigensystem.
    complex(8),dimension(Nlat,Nso,Nso,Lfreq),intent(out) :: Gloc
    character(len=*),intent(in)                          :: axis
    complex(8),dimension(Lfreq)                         :: wfreq
    complex(8),dimension(Nlso,Lfreq)                    :: csi
    integer                                             :: i,ilat
    integer                                             :: io,jo,is,js
    !
    call check_dimension("get_gf")
    if(.not.allocated(U).or..not.allocated(E)) &
      error stop "get_gf: complete eigensystem is not set"
    write(*,"(A)")"Get local Green's function, axis:"//str(axis)
    wfreq = build_frequency_array(axis)
    forall(i=1:Lfreq)csi(:,i)=one/(wfreq(i)+mu-E(:))
    Gloc = zero
    call start_timer
    do i=1,Lfreq
      do concurrent(ilat=1:Nlat,io=1:Nso,jo=1:Nso)
        is = io+(ilat-1)*Nso
        js = jo+(ilat-1)*Nso
        Gloc(ilat,io,jo,i) = &
          sum(U(is,:)*conjg(U(js,:))*csi(:,i))
      enddo
    enddo
    call stop_timer
  end subroutine get_gf


  subroutine free_Abhz_solution()
    ! Release all eigensystems and projector data owned by this module.
    if(allocated(U))deallocate(U)
    if(allocated(U_up))deallocate(U_up)
    if(allocated(U_down))deallocate(U_down)
    if(allocated(Uocc_up))deallocate(Uocc_up)
    if(allocated(Uocc_down))deallocate(Uocc_down)
    if(allocated(E))deallocate(E)
    if(allocated(Ev_up))deallocate(Ev_up)
    if(allocated(Ev_down))deallocate(Ev_down)
    if(allocated(PSzP))deallocate(PSzP)
    if(allocated(Epsp))deallocate(Epsp)
    Nocc      = 0
    Nocc_up   = 0
    Nocc_down = 0
    nullify(active_spH)
  end subroutine free_Abhz_solution


END MODULE BHZ_SOLVE
