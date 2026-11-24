MODULE BHZ_SOLVE
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

  public :: get_Abhz_EigenSolution
  public :: compose_Abhz_EigenSolution
  public :: free_Abhz_EigenSolution
  !
  public :: push_Bloch
  !
  public :: count_occupied_levels
  !
  public :: check_Pgap
  public :: check_Egap
  !
  public :: setup_fermi_chebychev
  public :: apply_fermi_chebychev


  public :: get_gf

  type(sparse_zmatrix_csr),pointer :: active_spH=>null()

contains


  ! Compute diag[Theta(MU-SPH)] at zero temperature.
  ! SPH remains owned by the caller.
  ! LAPACK uses a temporary dense copy and may optionally
  ! return the complete ordered eigensystem in EVEC and EVAL.
  ! CHEBYCHEV keeps SPH sparse and returns only DENSITIES.
  subroutine get_Abhz_EigenSolution(spH,density,evec,eval)
    type(sparse_zmatrix_csr),target,intent(in)                 :: spH
    real(8),dimension(:),allocatable,intent(out)               :: density
    complex(8),dimension(:,:),allocatable,optional,intent(out) :: evec
    real(8),dimension(:),allocatable,optional,intent(out)      :: eval
    complex(8),dimension(:,:),allocatable                      :: dense
    real(8),dimension(:),allocatable                           :: levels
    integer                                                    :: n,io
    !
    ! Check the Hamiltonian and optional eigensystem contract.
    n = spH%Nrow
    !
    if(n<1.or.spH%Ncol/=n)error stop "solve_Abhz: SPH must be a nonempty square matrix"
    if(present(evec).neqv.present(eval)) error stop "solve_Abhz: EVEC and EVAL must be both present"
    !
    ! Solve with the globally selected representation.
    allocate(density(n))
    call start_timer("Solve Anderson BHZ sector N="//str(n))
    !
    select case(method)
     case("lapack")
      allocate(dense(n,n),levels(n))
      dense = as_matrix(spH)
      call eigh(dense,levels)
      do io=1,n
        density(io) = sum(fermi_zero(levels-mu)*abs(dense(io,:))**2)
      enddo
      if(present(evec))then
        allocate(evec,source=dense)
      endif
      if(present(eval))then
        allocate(eval,source=levels)
      endif
      deallocate(dense,levels)
      !
     case("chebychev")
      active_spH => spH
      call cheb_solve(spH,density)
      nullify(active_spH)
      !
     case default
      error stop "solve_Abhz: invalid METHOD"
    end select
    call stop_timer()
  end subroutine get_Abhz_EigenSolution





  ! Assemble the two globally stored spin-sector eigensystems in the
  ! complete site-major basis and sort the combined spectrum into U,E.
  ! For SPIN_BLOCKS=1 the caller first reconstructs the down sector by
  ! time reversal; no diagonalization is performed here.
  subroutine compose_Abhz_EigenSolution()
    complex(8),allocatable :: raw_vectors(:,:)
    real(8),allocatable    :: raw_levels(:)
    integer,allocatable    :: order(:)
    integer                :: nblock,ib,ifull,ilat,iorb,istate
    !
    nblock = 2*Nlat
    if(.not.allocated(U_up).or..not.allocated(U_down)) &
      error stop "rebuild_Abhz_solution: missing eigenvectors"
    if(.not.allocated(Ev_up).or..not.allocated(Ev_down)) &
      error stop "rebuild_Abhz_solution: missing eigenvalues"
    if(size(Ev_up)/=nblock.or.size(Ev_down)/=nblock) &
      error stop "rebuild_Abhz_solution: invalid eigenvalue size"
    !
    allocate(raw_vectors(Nlso,Nlso),source=zero)
    allocate(raw_levels(Nlso),order(Nlso))
    raw_levels(1:nblock)       = Ev_up
    raw_levels(nblock+1:Nlso)  = Ev_down
    do ib=1,nblock
      ilat = (ib-1)/2+1
      iorb = modulo(ib-1,2)+1
      ifull = iorb+(ilat-1)*4
      raw_vectors(ifull,1:nblock) = U_up(ib,:)
      ifull = iorb+2+(ilat-1)*4
      raw_vectors(ifull,nblock+1:Nlso) = U_down(ib,:)
    enddo
    call sort(raw_levels,order)
    if(allocated(U))deallocate(U)
    if(allocated(E))deallocate(E)
    allocate(U(Nlso,Nlso),E(Nlso))
    E = raw_levels
    do istate=1,Nlso
      U(:,istate) = raw_vectors(:,order(istate))
    enddo
    call save_array("E.dat",E)
    deallocate(raw_vectors,raw_levels,order)
  end subroutine compose_Abhz_EigenSolution



  ! Release all dense eigensystems, occupied bases and projected-spin
  ! data owned by BHZ_SOLVE. The sparse Hamiltonians remain owned by
  ! the caller.
  subroutine free_Abhz_EigenSolution()
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
  end subroutine free_Abhz_EigenSolution











  ! Compute the diagonal density of the Chebychev approximation to the
  ! zero-temperature Fermi projector. Each diagonal entry is obtained
  ! by applying the same sparse polynomial filter to one basis vector.
  subroutine cheb_solve(spH,density)
    type(sparse_zmatrix_csr),intent(in) :: spH
    real(8),intent(out)                 :: density(:)
    complex(8),allocatable              :: seed(:)
    complex(8)                          :: value(1)
    integer                             :: n,io,info
    !
    ! Prepare the polynomial used for all diagonal elements.
    n = spH%Nrow
    if(size(density)/=n) error stop "cheb_density: invalid DENSITY size"
    !
    call setup_fermi_chebychev(spH)
    !
    allocate(seed(n),source=zero)
    do io=1,n
      seed(io) = one
      call cheb_projection(MatVec_Abhz,seed,[io],value,info)
      if(info/=0) &
        error stop "cheb_density: Chebychev projection failed"
      density(io) = real(value(1),8)
      seed(io) = zero
    enddo
  end subroutine cheb_solve


  ! Activate the T=0 filter for the square CSR matrix SPH. Uses MU,
  ! CHEB_ORDER and CHEB_KIND from input. SciFortran stores one active
  ! expansion per process. Prepare it again when switching matrices.
  subroutine setup_fermi_chebychev(spH)
    type(sparse_zmatrix_csr),intent(in) :: spH
    complex(8)                         :: diag
    real(8)                            :: coeff(0:cheb_order)
    real(8)                            :: scale,shift,lower,upper
    real(8)                            :: center,radius
    integer                            :: io
    !
    if(spH%Nrow<1.or.spH%Nrow/=spH%Ncol) &
      error stop "setup_fermi_chebychev: invalid SPH shape"
    lower = huge(1d0)
    upper = -huge(1d0)
    do io=1,spH%Nrow
      diag   = spH%get(io,io)
      center = dreal(diag)
      radius = sum(abs(spH%row(io)%vals))-abs(diag)
      lower  = min(lower,center-radius)
      upper  = max(upper,center+radius)
    enddo
    shift = 0.5d0*(upper+lower)
    scale = max(0.5d0*(upper-lower),1d-12)/0.95d0
    call cheb_fermi_coefficients(scale,shift,mu,coeff,cheb_kind)
    !
    write(*,*) "Chebychev: degree, scale, shift =", &
      cheb_order,scale,shift
  end subroutine setup_fermi_chebychev


  ! Return W=P_M V, with V, W and ROWS of length SPH%Nrow and ROWS
  ! equal to [1,...,SPH%Nrow]. The caller owns SPH and must have
  ! activated its expansion via setup_fermi_chebychev. Serial within
  ! each process, since the callback uses shared SciFortran state.
  subroutine apply_fermi_chebychev(spH,rows,v,w)
    type(sparse_zmatrix_csr),target,intent(in) :: spH
    integer,intent(in)                         :: rows(:)
    complex(8),intent(in)                      :: v(:)
    complex(8),intent(out)                     :: w(:)
    integer                                    :: n,info
    !
    n = spH%Nrow
    if(size(v)/=n.or.size(w)/=n.or.size(rows)/=n) &
      error stop "apply_fermi_chebychev: invalid vector size"
    active_spH => spH
    call cheb_projection(MatVec_Abhz,v,rows,w,info)
    nullify(active_spH)
    if(info/=0) &
      error stop "apply_fermi_chebychev: Chebychev projection failed"
  end subroutine apply_fermi_chebychev



  ! Build the Chebychev coefficients of Theta(FERMI_LEVEL-E) and
  ! activate them in SciFortran. KIND selects no damping or the Jackson
  ! kernel and defaults to JACKSON.
  subroutine cheb_fermi_coefficients(scale,shift,fermi_level,coeff,kind)
    real(8),intent(in)                   :: scale,shift,fermi_level
    real(8),dimension(0:),intent(out)    :: coeff
    character(len=*),optional,intent(in) :: kind
    character(len=16)                    :: kernel
    real(8)                              :: xfermi,theta,angle,g
    integer                              :: m,order
    !
    if(scale<=0d0)error stop "cheb_fermi_coefficients: SCALE must be positive"
    order  = ubound(coeff,1)
    kernel = "jackson"
    if(present(kind))kernel=to_lower(trim(kind))
    select case(kernel)
     case("none","n")   ;kernel = "none"
     case("jackson","j");kernel = "jackson"
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
    !
    angle    = acos(xfermi)
    coeff(0) = (pi-angle)/pi
    do m=1,order
      coeff(m) = -2d0*sin(m*angle)/(pi*m)
    enddo
    if(kernel=="jackson")then
      theta = pi/(order+2)
      do m=0,order
        g = ((order+2-m)*cos(m*theta) + sin(m*theta)/tan(theta))/(order+2)
        coeff(m) = coeff(m)*g
      enddo
    endif
    call cheb_set(scale,shift,coeff)
  end subroutine cheb_fermi_coefficients



  ! Apply the currently active sparse BHZ Hamiltonian to one vector.
  ! This private callback connects SciFortran's Chebychev recurrence to
  ! the caller-owned CSR matrix selected by the surrounding operation.
  subroutine MatVec_Abhz(Nloc,vin,vout)
    integer    :: Nloc
    complex(8) :: vin(Nloc),vout(Nloc)
    !
    if(.not.associated(active_spH)) &
      error stop "MatVec_Abhz: no active sparse Hamiltonian"
    vout = active_spH%dot(vin)
  end subroutine MatVec_Abhz











  ! Store a complete, ordered dense eigensystem in the module state,
  ! determine its T=0 occupation and check the single-particle gap.
  ! This is the Lapack/full-matrix entry point used by legacy callers.
  subroutine push_Bloch(Evec,Evals)
    complex(8),dimension(:,:),intent(in) :: Evec
    real(8),dimension(:),intent(in)      :: Evals
    !
    call check_dimension("push_Bloch")
    if(size(Evec,1)/=Nlso.or.size(Evec,2)/=Nlso)error stop "push_Bloch: invalid EVEC shape"
    if(size(Evals)/=Nlso) error stop "push_Bloch: invalid EVALS size"
    if(allocated(U).or.allocated(E)) error stop "push_Bloch: eigensystem already allocated"
    !
    allocate(U,source=Evec)
    allocate(E,source=Evals)
    Nocc = count_occupied_levels(E,"full")
    !
    call save_array("E.dat",E)
    !
    call check_Egap("push_Bloch")
    !
  end subroutine push_Bloch





  ! Count zero-temperature levels strictly below MU in an ordered or
  ! unordered spectrum. A state at MU is reported and excluded because
  ! it does not define an unambiguous occupied-state projector.
  function count_occupied_levels(eval,sector) result(noccupied)
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





  ! Inspect the spectrum of the occupied-space operator P Sz P and
  ! write its gap. A closing is reported as a warning because it makes
  ! the two projected-spin topological sectors ill-defined.
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


  ! Compute and write the dense single-particle gap across MU from the
  ! globally stored ordered spectrum E and occupation Nocc.
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






  ! Evaluate the local Green function by the spectral representation of
  ! the globally stored complete eigensystem. This is a dense Lapack
  ! post-processing path and is not part of the Chebychev LCM workflow.
  subroutine get_gf(Gloc,axis)
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




end module BHZ_SOLVE
