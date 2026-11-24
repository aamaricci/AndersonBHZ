MODULE BHZ_LCM
  USE SCIFOR
  USE DMFT_TOOLS
  USE COMMON
  USE BHZ_SOLVE
  implicit none
  private
  public :: hk_to_Chern
  public :: hk_to_spin_Chern
  public :: spin_chern_numbers
  public :: local_spin_chern_marker


contains


  function spin_chern_numbers() result(chern)
    ! Return the physical up/down Chern numbers for the globally
    ! selected spin representation. For SPIN_BLOCKS=1 or 2 the
    ! occupied sector bases are prepared internally from U_UP/DOWN and
    ! EV_UP/DOWN. SPIN_BLOCKS=0 uses the projected-spin sectors of the
    ! complete globally stored eigensystem U,E.
    real(8) :: chern(2)
    !
    select case(spin_blocks)
    case(1)
      call setup_occupied_bases()
      chern(1) = single_point_chern_basis(Uocc_up)
      chern(2) = -chern(1)
    case(2)
      call setup_occupied_bases()
      chern(1) = single_point_chern_basis(Uocc_up)
      chern(2) = single_point_chern_basis(Uocc_down)
    case(0)
      if(.not.allocated(U).or..not.allocated(E)) &
        error stop "spin_chern_numbers: full eigensystem is not set"
      chern(1) = single_point_spin_chern(spin=2)
      chern(2) = single_point_spin_chern(spin=1)
    case default
      error stop "spin_chern_numbers: invalid SPIN_BLOCKS"
    end select
  end function spin_chern_numbers




  subroutine local_spin_chern_marker(marker_up,marker_down)
    ! Return the physical up/down local Chern markers using the global
    ! boundary and spin representation. Occupied bases are prepared
    ! internally. The caller never handles UOCC explicitly.
    real(8),dimension(:,:),allocatable,intent(out) :: marker_up
    real(8),dimension(:,:),allocatable,intent(out) :: marker_down
    !
    select case(spin_blocks)
    case(1)
      call setup_occupied_bases()
      if(bhz_pbc)then
        call pbc_local_chern_marker_basis(Uocc_up,2,marker_up)
      else
        call obc_local_chern_marker_basis(Uocc_up,2,marker_up)
      endif
      allocate(marker_down,source=-marker_up)
    case(2)
      call setup_occupied_bases()
      if(bhz_pbc)then
        call pbc_local_chern_marker_basis(Uocc_up,2,marker_up)
        call pbc_local_chern_marker_basis(Uocc_down,2,marker_down)
      else
        call obc_local_chern_marker_basis(Uocc_up,2,marker_up)
        call obc_local_chern_marker_basis(Uocc_down,2,marker_down)
      endif
    case(0)
      if(bhz_pbc)then
        call pbc_local_spin_chern_marker(spin=2,lcm=marker_up)
        call pbc_local_spin_chern_marker(spin=1,lcm=marker_down)
      else
        call obc_local_spin_chern_marker(spin=2,lcm=marker_up)
        call obc_local_spin_chern_marker(spin=1,lcm=marker_down)
      endif
    case default
      error stop "local_spin_chern_marker: invalid SPIN_BLOCKS"
    end select
  end subroutine local_spin_chern_marker




  subroutine setup_occupied_bases()
    ! Build the global zero-temperature occupied bases from the
    ! globally stored sector eigensystems. A state at MU is reported
    ! and excluded because it does not define an idempotent projector.
    if(spin_blocks==0)return
    if(spin_blocks==1.and.allocated(Uocc_up))return
    if(spin_blocks==2.and.allocated(Uocc_up).and. &
      allocated(Uocc_down))return
    if(.not.allocated(U_up).or..not.allocated(Ev_up)) &
      error stop "setup_occupied_bases: up eigensystem is not set"
    call set_occupied_basis(U_up,Ev_up,"up",Uocc_up,Nocc_up)
    if(spin_blocks==1)then
      Nocc_down = Nocc_up
    else
      if(.not.allocated(U_down).or..not.allocated(Ev_down)) &
        error stop "setup_occupied_bases: down eigensystem is not set"
      call set_occupied_basis(U_down,Ev_down,"down", &
        Uocc_down,Nocc_down)
    endif
    Nocc = Nocc_up+Nocc_down
  end subroutine setup_occupied_bases




  subroutine set_occupied_basis(evec,eval,sector,Uocc,noccupied)
    complex(8),dimension(:,:),intent(in)              :: evec
    real(8),dimension(:),intent(in)                   :: eval
    character(len=*),intent(in)                       :: sector
    complex(8),dimension(:,:),allocatable,intent(out) :: Uocc
    integer,intent(out)                               :: noccupied
    real(8),parameter                                 :: mu_tol=1d-12
    integer                                           :: n
    !
    n = size(eval)
    if(any(shape(evec)/=[n,n])) &
      error stop "set_occupied_basis: invalid eigenvector shape"
    if(n>1.and.any(eval(2:)<eval(:n-1))) &
      error stop "set_occupied_basis: eigenvalues are not ordered"
    if(any(abs(eval-mu)<=mu_tol)) &
      write(*,*)"occupied levels warning: level at MU in ",trim(sector)
    noccupied = count(eval<mu-mu_tol)
    if(noccupied<1.or.noccupied>=n) &
      error stop "set_occupied_basis: empty or complete sector"
    allocate(Uocc(n,noccupied),source=evec(:,1:noccupied))
    write(*,*)"Occupied states [",trim(sector),"]:",noccupied
  end subroutine set_occupied_basis




  !####################################################################
  !                 H(k) ==>   C   and  C_\sigma
  !####################################################################
  !+------------------------------------------------------------------+
  !get Chern number of H(k) using discretized Berry flux in the BZ
  !+------------------------------------------------------------------+
  function hk_to_Chern(Hk,Nkvec,plot_berry) result(chern)
    complex(8),intent(in),dimension(:,:,:)    :: Hk    ![Nlso][Nlso][Nktot]
    integer,intent(in),dimension(2)           :: Nkvec ![Nk1][Nk2]: prod(Nkvec)=Nktot
    logical,optional                          :: plot_berry
    real(8)                                   :: Chern
    !
    integer                                   :: Nocc
    integer                                   :: Nktot
    integer                                   :: Nkx,Nky
    integer                                   :: i,j
    integer                                   :: ikx,iky
    integer                                   :: ikxP,ikyP
    integer                                   :: ik,iocc
    complex(8),dimension(:,:),allocatable     :: Eigvec ![Nlso][Nlso]
    real(8),dimension(:),allocatable          :: Eigval ![Nlso]
    complex(8),dimension(:,:,:),allocatable   :: Gmat
    complex(8),dimension(:,:,:,:),allocatable :: BlochStates ![Nkx][Nky][Noccupied][Nlso]
    complex(8),dimension(4)                   :: Ulink
    real(8),dimension(:,:),allocatable        :: BerryCurvature
    real(8)                                   :: berry_phase
    integer                                   :: unit
    logical                                   :: bool
    !
    bool = .false.;if(present(plot_berry))bool=plot_berry
    !
    Nocc  = Nso/2
    Nkx   = Nkvec(1)
    Nky   = Nkvec(2)
    Nktot = product(Nkvec)
    call assert_shape(Hk,[Nso,Nso,Nktot],"Get_Chern_NUmber","Hk")
    if(Nkx*Nky/=Nktot)stop "ERROR Get_Chern_Number: Nktot = prod(Nkvec)"
    !
    !
    !1. Get the Bloch states from H(:,:,k)
    allocate(Eigvec(Nso,Nso))
    allocate(Eigval(Nso))
    allocate(BlochStates(Nkx,Nky,Nocc,Nso))
    allocate(BerryCurvature(Nkx,Nky))
    allocate(Gmat(4,Nocc,Nocc))
    ik=0
    do ikx=1,Nkx
      do iky=1,Nky
        ik=ik+1
        Eigvec = Hk(:,:,ik)
        call eigh(Eigvec,Eigval)
        do iocc=1,Nocc
          BlochStates(ikx,iky,iocc,:) = Eigvec(:,iocc)
        enddo
      enddo
    enddo
    deallocate(Eigvec,Eigval)
    !
    !
    !2. Evaluate the Berry Curvature
    chern=0d0
    do ikx= 1, Nkx
      ikxP = modulo(ikx,Nkx) + 1
      do iky= 1, Nky
        ikyP = modulo(iky,Nky) + 1
        !
        if(Nocc==1)then
          Ulink(1) = dot_product(BlochStates(ikx,iky,1,:)  , BlochStates(ikx,ikyP,1,:))
          Ulink(2) = dot_product(BlochStates(ikx,ikyP,1,:) , BlochStates(ikxP,ikyP,1,:))
          Ulink(3) = dot_product(BlochStates(ikxP,ikyP,1,:), BlochStates(ikxP,iky,1,:))
          Ulink(4) = dot_product(BlochStates(ikxP,iky,1,:) , BlochStates(ikx,iky,1,:))
        else
          do concurrent(i=1:Nocc,j=1:Nocc)
            gmat(1,i,j) = dot_product(BlochStates(ikx,iky,i,:)  , BlochStates(ikx,ikyP,j,:))
            gmat(2,i,j) = dot_product(BlochStates(ikx,ikyP,i,:) , BlochStates(ikxP,ikyP,j,:))
            gmat(3,i,j) = dot_product(BlochStates(ikxP,ikyP,i,:), BlochStates(ikxP,iky,j,:))
            gmat(4,i,j) = dot_product(BlochStates(ikxP,iky,i,:) , BlochStates(ikx,iky,j,:))
          enddo
          Ulink(1) = det(gmat(1,:,:))
          Ulink(2) = det(gmat(2,:,:))
          Ulink(3) = det(gmat(3,:,:))
          Ulink(4) = det(gmat(4,:,:))
        endif
        !
        berry_phase             = dimag(zlog( product(Ulink(:))  ))
        chern                   = chern + berry_phase/pi2
        BerryCurvature(ikx,iky) = berry_phase!*one_over_area=Nkx/pi2*Nkx/pi2
        !
      enddo
    enddo
    !
    open(unit=free_unit(unit),file="Hk_to_Chern.dat")
    write(unit,*)chern
    close(unit)
    !
    if(bool)&
      call splot3d("Berry_Curvature.dat",&
      linspace(0d0,pi2,Nkx,iend=.false.),&
      linspace(0d0,pi2,Nky,iend=.false.),&
      BerryCurvature(:Nkx,:Nky))
  end function hk_to_Chern

  !+------------------------------------------------------------------+
  !get spin Chern number of H(k) using discretized Berry flux in the BZ
  !+------------------------------------------------------------------+
  function hk_to_spin_Chern(Hk,Nkvec,spin,berry) result(sp_chern)
    complex(8),intent(in),dimension(:,:,:)      :: Hk    ![Nlso][Nlso][Nktot]
    integer,intent(in),dimension(2)             :: Nkvec ![Nk1][Nk2]: prod(Nkvec)=Nktot
    integer,optional                            :: spin
    real(8),dimension(:,:),allocatable,optional :: berry
    real(8)                                     :: sp_chern
    !
    integer                                     :: spin_
    integer                                     :: Nocc
    integer                                     :: Nktot
    integer                                     :: Nkx,Nky
    integer                                     :: i,j,ispin
    integer                                     :: ikx,iky
    integer                                     :: ikxP,ikyP
    integer                                     :: ik,iocc
    complex(8),dimension(:,:),allocatable       :: Eigvec ![Nlso][Nlso]
    real(8),dimension(:),allocatable            :: Eigval ![Nlso]
    complex(8),dimension(:,:,:,:),allocatable   :: BlochStates ![Nkx][Nky][Noccupied][Nlso]
    complex(8),dimension(4)                     :: Ulink
    real(8)                                     :: berry_phase,chern(2)
    integer                                     :: unit
    logical                                     :: bool
    !
    spin_ = 0;if(present(spin))spin_=spin
    Nocc  = Nso/2
    Nkx   = Nkvec(1)
    Nky   = Nkvec(2)
    Nktot = product(Nkvec)
    call assert_shape(Hk,[Nso,Nso,Nktot],"Get_Chern_NUmber","Hk")
    if(Nkx*Nky/=Nktot)stop "ERROR Get_Chern_Number: Nktot = prod(Nkvec)"
    !
    !
    !1. Get the Bloch states from H(:,:,k)
    allocate(Eigvec(Nso,Nso))
    allocate(Eigval(Nso))
    allocate(BlochStates(Nkx,Nky,Nocc,Nso))
    if(present(berry))then
      if(allocated(berry))deallocate(berry)
      allocate(berry(Nkx,Nky))
    endif
    !
    ik=0
    do ikx=1,Nkx
      do iky=1,Nky
        ik=ik+1
        Eigvec = Hk(:,:,ik)
        call eigh(Eigvec,Eigval)
        do iocc=1,Nocc
          BlochStates(ikx,iky,iocc,:) = Eigvec(:,iocc)
        enddo
      enddo
    enddo
    deallocate(Eigvec,Eigval)
    !
    !
    !2. Evaluate the Berry Curvature
    chern=0d0
    do ikx= 1, Nkx
      ikxP = modulo(ikx,Nkx) + 1
      do iky= 1, Nky
        ikyP = modulo(iky,Nky) + 1
        !
        select case(spin_)
         case default;stop "hk_to_spin_Chern error: spin != {0,1,2}: 0=Z2/Both, 1=UP, 2=DW"
         case(0)
          do ispin=1,2
            Ulink(1)    = dot_product(BlochStates(ikx,iky,ispin,:)  , BlochStates(ikx,ikyP,ispin,:))
            Ulink(2)    = dot_product(BlochStates(ikx,ikyP,ispin,:) , BlochStates(ikxP,ikyP,ispin,:))
            Ulink(3)    = dot_product(BlochStates(ikxP,ikyP,ispin,:), BlochStates(ikxP,iky,ispin,:))
            Ulink(4)    = dot_product(BlochStates(ikxP,iky,ispin,:) , BlochStates(ikx,iky,ispin,:))
            berry_phase = -dimag(zlog( product(Ulink(:))  ))/pi2
            chern(ispin)= chern(ispin) + berry_phase
          enddo
          sp_chern = 0.5d0*(chern(1)-chern(2))
          open(unit=free_unit(unit),file="Hk_to_z2.dat")
          write(unit,*)sp_chern
          close(unit)
         case(1,2)
          ispin=spin_
          Ulink(1)    = dot_product(BlochStates(ikx,iky,ispin,:)  , BlochStates(ikx,ikyP,ispin,:))
          Ulink(2)    = dot_product(BlochStates(ikx,ikyP,ispin,:) , BlochStates(ikxP,ikyP,ispin,:))
          Ulink(3)    = dot_product(BlochStates(ikxP,ikyP,ispin,:), BlochStates(ikxP,iky,ispin,:))
          Ulink(4)    = dot_product(BlochStates(ikxP,iky,ispin,:) , BlochStates(ikx,iky,ispin,:))
          berry_phase = dimag(zlog( product(Ulink(:))  ))/pi2
          chern(ispin)= chern(ispin) + berry_phase
          if(present(berry))&
            berry(ikx,iky) = berry_phase!*one_over_area=Nkx/pi2*Nkx/pi2
          sp_chern = chern(ispin)
          open(unit=free_unit(unit),file="Hk_to_spin_Chern_s"//str(spin_)//".dat")
          write(unit,*)sp_chern
          close(unit)
        end select
      enddo
    enddo
    !
    !
  end function hk_to_spin_Chern







  !####################################################################
  !              Single Point C and C_\sigma
  !####################################################################
  !+------------------------------------------------------------------+
  ! Evaluate the Chern number of an explicitly supplied occupied
  ! basis using the single-point approximation. UOCC has shape
  ! [Nlat*Ninternal,Noccupied] and may represent either one conserved
  ! spin sector or the complete occupied space.
  !+------------------------------------------------------------------+
  function single_point_chern_basis(Uocc,method) result(chern)
    complex(8),dimension(:,:),intent(in)     :: Uocc
    character(len=*),optional,intent(in)     :: method
    real(8)                                  :: chern
    !
    character(len=1)                         :: method_
    real(8),dimension(2)                     :: b1,b2
    complex(8),dimension(:,:),allocatable    :: Ub1,Ub2,Umb1,Umb2
    complex(8),dimension(:,:),allocatable    :: Vb1,Vb2,Vmb1,Vmb2
    complex(8)                               :: sum_chern
    integer                                  :: i,nrow,nocc
    !
    method_='s'
    if(present(method))method_=method
    nrow = size(Uocc,1)
    nocc = size(Uocc,2)
    if(Nlat<1)error stop "single_point_chern_basis: invalid NLAT"
    if(modulo(nrow,Nlat)/=0.or.nocc<1.or.nocc>nrow) &
      error stop "single_point_chern_basis: invalid UOCC shape"
    !
    call TB_get_bk(b1,b2)
    call start_timer("single_point_chern_basis")
    Ub1 = periodic_gauge(Uocc,b1)
    Ub2 = periodic_gauge(Uocc,b2)
    Vb1 = dual_state_basis(Uocc,Ub1)
    Vb2 = dual_state_basis(Uocc,Ub2)
    !
    sum_chern = zero
    if(to_lower(method_)=='s')then
      Umb1  = periodic_gauge(Uocc,-b1)
      Umb2  = periodic_gauge(Uocc,-b2)
      Vmb1  = dual_state_basis(Uocc,Umb1)
      Vmb2  = dual_state_basis(Uocc,Umb2)
      do i=1,nocc
        sum_chern = sum_chern + &
          dot_product(Vmb1(:,i)-Vb1(:,i), &
          Vmb2(:,i)-Vb2(:,i))/4d0
      enddo
    else
      do i=1,nocc
        sum_chern = sum_chern + dot_product(Vb1(:,i),Vb2(:,i))
      enddo
    endif
    call stop_timer()
    chern = dimag(sum_chern)/pi
  end function single_point_chern_basis




  !+------------------------------------------------------------------+
  ! Evaluate the OBC local Chern marker from an explicitly supplied
  ! occupied basis. NINTERNAL is the number of local components in
  ! UOCC, e.g. two orbitals for one conserved spin sector.
  !+------------------------------------------------------------------+
  subroutine obc_local_chern_marker_basis(Uocc,ninternal,lcm)
    complex(8),dimension(:,:),intent(in)       :: Uocc
    integer,intent(in)                         :: ninternal
    real(8),dimension(:,:),allocatable, &
      intent(out)                              :: lcm
    !
    complex(8),dimension(:,:),allocatable      :: P,Xcomm_P,Ycomm_P
    real(8),dimension(:,:),allocatable         :: R
    complex(8),dimension(:,:),allocatable      :: Q2
    real(8),dimension(:,:,:,:),allocatable     :: Q4
    integer                                    :: nrow,ix,iy,ilat
    !
    nrow = size(Uocc,1)
    if(ninternal<1.or.nrow/=Nlat*ninternal) &
      error stop "obc_local_chern_marker_basis: invalid NINTERNAL"
    if(size(Uocc,2)<1.or.size(Uocc,2)>nrow) &
      error stop "obc_local_chern_marker_basis: invalid UOCC shape"
    !
    call start_timer("obc_local_chern_marker_basis")
    allocate(P(nrow,nrow),Xcomm_P(nrow,nrow),Ycomm_P(nrow,nrow))
    P = Uocc.mx.transpose(conjg(Uocc))
    R = TB_build_Rcoord([Nx,Ny],to_home=.false., &
      ninternal=ninternal)
    Xcomm_P = (diag(one*R(:,1)).mx.P) - &
      (P.mx.diag(one*R(:,1)))
    Ycomm_P = (diag(one*R(:,2)).mx.P) - &
      (P.mx.diag(one*R(:,2)))
    Q2 = -4d0*pi*dimag((P.mx.Xcomm_P).mx.Ycomm_P)
    Q4 = reshape_rank2_to_rank4(Q2,ninternal,Nlat)
    !
    allocate(lcm(Nx,Ny),source=0d0)
    do ix=1,Nx
      do iy=1,Ny
        ilat = ix+(iy-1)*Nx
        lcm(ix,iy) = trace(Q4(:,:,ilat,ilat))
      enddo
    enddo
    call stop_timer()
  end subroutine obc_local_chern_marker_basis




  !+------------------------------------------------------------------+
  ! Evaluate the PBC local Chern marker from an explicitly supplied
  ! occupied basis. The overlap inversion is restricted to the actual
  ! number of occupied states in UOCC.
  !+------------------------------------------------------------------+
  subroutine pbc_local_chern_marker_basis(Uocc,ninternal,lcm,method)
    complex(8),dimension(:,:),intent(in)       :: Uocc
    integer,intent(in)                         :: ninternal
    real(8),dimension(:,:),allocatable, &
      intent(out)                              :: lcm
    character(len=*),optional,intent(in)       :: method
    !
    character(len=1)                           :: method_
    real(8),dimension(2)                       :: b1,b2
    complex(8),dimension(:,:),allocatable      :: Ub1,Ub2,Umb1,Umb2
    complex(8),dimension(:,:),allocatable      :: Vb1,Vb2,Vmb1,Vmb2
    complex(8),dimension(:,:),allocatable      :: P,Pgs
    complex(8),dimension(:,:),allocatable      :: Pb1,Pb2,Pmb1,Pmb2
    complex(8),dimension(:,:),allocatable      :: Chern_Q2
    real(8),dimension(:,:,:,:),allocatable     :: Chern_Q4
    integer                                    :: nrow,ix,iy,ilat
    !
    method_='s'
    if(present(method))method_=method
    nrow = size(Uocc,1)
    if(ninternal<1.or.nrow/=Nlat*ninternal) &
      error stop "pbc_local_chern_marker_basis: invalid NINTERNAL"
    if(size(Uocc,2)<1.or.size(Uocc,2)>nrow) &
      error stop "pbc_local_chern_marker_basis: invalid UOCC shape"
    !
    call TB_get_bk(b1,b2)
    call start_timer("pbc_local_chern_marker_basis")
    Ub1 = periodic_gauge(Uocc,b1)
    Ub2 = periodic_gauge(Uocc,b2)
    Vb1 = dual_state_basis(Uocc,Ub1)
    Vb2 = dual_state_basis(Uocc,Ub2)
    Pgs = Uocc.mx.transpose(conjg(Uocc))
    Pb1 = Vb1.mx.transpose(conjg(Vb1))
    Pb2 = Vb2.mx.transpose(conjg(Vb2))
    P   = (Pb1.mx.Pb2) - (Pb2.mx.Pb1)
    !
    if(to_lower(method_)=='s')then
      Umb1  = periodic_gauge(Uocc,-b1)
      Umb2  = periodic_gauge(Uocc,-b2)
      Vmb1  = dual_state_basis(Uocc,Umb1)
      Vmb2  = dual_state_basis(Uocc,Umb2)
      Pmb1  = Vmb1.mx.transpose(conjg(Vmb1))
      Pmb2  = Vmb2.mx.transpose(conjg(Vmb2))
      ! Symmetric differences require all four commutators.
      P = ((Pb1.mx.Pb2) - (Pb2.mx.Pb1)) + &
        ((Pmb1.mx.Pmb2) - (Pmb2.mx.Pmb1)) - &
        ((Pb1.mx.Pmb2) - (Pmb2.mx.Pb1)) - &
        ((Pmb1.mx.Pb2) - (Pb2.mx.Pmb1))
      P = P/4d0
    endif
    !
    Chern_Q2 = dimag(P.mx.Pgs)/pi2*Nlat
    Chern_Q4 = reshape_rank2_to_rank4(Chern_Q2,ninternal,Nlat)
    allocate(lcm(Nx,Ny),source=0d0)
    do ix=1,Nx
      do iy=1,Ny
        ilat = ix+(iy-1)*Nx
        lcm(ix,iy) = trace(Chern_Q4(:,:,ilat,ilat))
      enddo
    enddo
    call stop_timer()
  end subroutine pbc_local_chern_marker_basis




  !+------------------------------------------------------------------+
  !Evaluate Chern number using single point approximation
  !Ref. Ceresoli-Resta (2007)
  !<https://journals.aps.org/prb/abstract/10.1103/PhysRevB.76.012405>
  !+------------------------------------------------------------------+
  function single_point_chern(method) result(sp_chern)
    character(len=*),optional                 :: method
    real(8)                                   :: sp_chern
    !
    character(len=1)                          :: method_
    real(8),dimension(2)                      :: b1,b2
    complex(8),dimension(:,:),allocatable     :: Ub1,Ub2,Umb1,Umb2
    complex(8),dimension(:,:),allocatable     :: Vb1,Vb2,Vmb1,Vmb2
    complex(8)                                :: sum_chern
    integer                                   :: i,j
    !
    method_='s' ;if(present(method))method_=method
    !
    call check_dimension("single_point_chern")
    call assert_shape(U,[size(U,1),Nlso],"single_point_chern","U")
    !
    call TB_get_bk(b1,b2)
    !
    call start_timer("single_point_chern")
    !
    Ub1 = periodic_gauge(U,b1)
    Ub2 = periodic_gauge(U,b2)
    !
    Vb1 = dual_state(U,Ub1)
    Vb2 = dual_state(U,Ub2)
    !
    sum_chern = zero
    if(to_lower(method_)=='s')then
      Umb1 = periodic_gauge(U,-b1)
      Umb2 = periodic_gauge(U,-b2)
      Vmb1 = dual_state(U,Umb1)
      Vmb2 = dual_state(U,Umb2)
      do i=1,Nocc
        sum_chern = sum_chern + dot_product((Vmb1(:,i)-Vb1(:,i)),(Vmb2(:,i)-Vb2(:,i)))/4d0
      enddo
    else
      do i=1,Nocc
        sum_chern = sum_chern + dot_product(Vb1(:,i),Vb2(:,i))
      enddo
    end if
    !
    call stop_timer()
    !
    sp_chern = dimag(sum_chern)/pi
    !
  end function single_point_chern



  !+------------------------------------------------------------------+
  !Evaluate spin Chern number using single point approximation
  ! Ref. Favata-Marrazzo (2023)
  !<https://iopscience.iop.org/article/10.1088/2516-1075/acba6f/meta>
  !+------------------------------------------------------------------+
  ! function single_point_spin_chern(U,E,Sz,spin,method,gap) result(sp_chern)
  function single_point_spin_chern(spin,method) result(sp_chern)
    integer                                    :: spin
    character(len=*),optional                  :: method
    real(8)                                    :: sp_chern
    !
    character(len=1)                           :: method_
    real(8),dimension(2)                       :: b1,b2
    complex(8),dimension(:,:),allocatable      :: Q
    complex(8),dimension(:,:),allocatable      :: Qb1,Qb2,Qmb1,Qmb2
    complex(8),dimension(:,:),allocatable      :: Vb1,Vb2,Vmb1,Vmb2
    complex(8)                                 :: sum_chern
    integer                                    :: i,m,N,a
    !
    method_='s' ;if(present(method))method_=method
    !
    if(spin<1.OR.spin>2)stop "single_point_spin_chern error: spin < 1 OR spin > 2"
    call check_dimension("single_point_chern")
    call assert_shape(U,[size(U,1),Nlso],"single_point_chern","U")
    !
    call TB_get_bk(b1,b2)
    !
    if(.not.allocated(PSzP))call eigh_PSzP()
    call projected_spin_sector(spin,a,N)
    !
    call start_timer("single_point_spinChern_s"//str(spin))
    !
    !|q_i0> = sum_m=1,Nocc q_i(m)|u_m0>
    ! q     = U x P_{a'b'} [Nlso,Nocc]
    allocate(Q(Nlso,Nocc))
    Q = (U(:,1:Nocc).mx.PSzP)
    !
    Qb1 = periodic_gauge(Q,b1)
    Qb2 = periodic_gauge(Q,b2)
    !
    Vb1 = dual_state(Q,Qb1,spin)
    Vb2 = dual_state(Q,Qb2,spin)
    !
    sum_chern = zero
    if(to_lower(method_)=='s')then
      Qmb1 = periodic_gauge(Q,-b1)
      Qmb2 = periodic_gauge(Q,-b2)
      Vmb1 = dual_state(Q,Qmb1,spin)
      Vmb2 = dual_state(Q,Qmb2,spin)
      do i=1,N
        sum_chern = sum_chern + &
          dot_product(Vmb1(:,i)-Vb1(:,i), &
          Vmb2(:,i)-Vb2(:,i))/4d0
      enddo
    else
      do i=1,N
        sum_chern = sum_chern + dot_product(Vb1(:,i),Vb2(:,i))
      enddo
    end if
    !
    call stop_timer()
    !
    sp_chern = dimag(sum_chern)/pi
    !
  end function single_point_spin_chern





  !####################################################################
  !    Finite System/OBC local/space resolved C and C_\sigma markers
  !####################################################################
  !+------------------------------------------------------------------+
  !Evaluate the local Chern marker
  !Ref. Bianco-Resta (2011)
  !<https://doi.org/10.1103/PhysRevB.84.241106>
  !+------------------------------------------------------------------+
  subroutine obc_local_chern_marker(lcm)
    real(8),dimension(:,:),allocatable     :: lcm
    !
    complex(8),dimension(:,:),allocatable  :: P,Xcomm_P,Ycomm_P
    real(8),dimension(:,:),allocatable     :: R
    real(8),dimension(:,:,:,:),allocatable :: Q4
    integer                                :: ix,iy,i,j,ilat,jlat,io,jo
    !
    call check_dimension("OBC_local_chern_marker")
    call assert_shape(U,[size(U,1),Nlso],"OBC_local_chern_marker","H")
    !
    allocate(P(Nlso,Nlso))
    allocate(Xcomm_P(Nlso,Nlso))
    allocate(Ycomm_P(Nlso,Nlso))
    !
    call start_timer("obc_local_ChernMarker")
    !
    P = ( U(:,1:Nocc).mx.transpose(conjg(U(:,1:Nocc))) )
    !
    allocate(R(Nlso,2))
    R       = TB_build_Rcoord([Nx,Ny],to_home=.false.)
    Xcomm_P = (diag(one*R(:,1)).mx.P) - (P.mx.diag(one*R(:,1)))
    Ycomm_P = (diag(one*R(:,2)).mx.P) - (P.mx.diag(one*R(:,2)))
    P       = -4*pi*dimag(((P.mx.Xcomm_P).mx.Ycomm_P))
    !
    allocate(Q4(Nso,Nso,Nlat,Nlat))
    Q4 = reshape_rank2_to_rank4(P,Nso,Nlat)
    !
    if(allocated(lcm))deallocate(lcm)
    allocate(lcm(Nx,Ny))
    lcm = 0d0
    do ix = 1,Nx
      do iy = 1,Ny
        ilat = ix + (iy-1)*Nx
        lcm(ix,iy) = trace(Q4(:,:,ilat,ilat))
      enddo
    enddo
    call stop_timer()
    !
  end subroutine obc_local_chern_marker


  !+------------------------------------------------------------------+
  !Evaluate the local spin Chern marker
  !Ref. Bau'-Marrazzo (2024b)
  ! <https://arxiv.org/abs/2404.04598>`
  ! as used also in
  !Ref. Amaricci et al (2017)
  !<https://journals.aps.org/prb/abstract/10.1103/PhysRevB.95.205120>
  !+------------------------------------------------------------------+
  subroutine obc_local_spin_chern_marker(spin,lcm)
    integer                                    :: spin
    real(8),dimension(:,:),allocatable         :: lcm
    !
    complex(8),dimension(:,:),allocatable      :: Q
    complex(8),dimension(Nlso,Nlso)            :: P,Xcomm_P,Ycomm_P
    real(8),dimension(Nlso,2)                  :: R
    complex(8),dimension(Nlso,Nlso)            :: Q2
    real(8),dimension(Nso,Nso,Nlat,Nlat)       :: Q4
    integer                                    :: i,ix,iy,ilat,m,N,a
    if(spin<1.OR.spin>2)stop "OBC_local_spin_chern_marker error: spin < 1 OR spin > 2"
    !
    call check_dimension("OBC_local_spin_chern_marker")
    call assert_shape(U,[size(U,1),Nlso],"OBC_local_spin_chern_marker","U")
    !
    if(.not.allocated(PSzP))call eigh_PSzP()
    call projected_spin_sector(spin,a,N)
    !
    call start_timer("obc_local_spinChernMarker_s"//str(spin))
    !
    allocate(Q(Nlso,Nocc))
    Q = ( U(:,1:Nocc).mx. PSzP )
    !
    !GS projectors P_-, P_+
    select case(spin)
     case(1);P = Q(:,1:N).mx.transpose(conjg(Q(:,1:N)))
     case(2);P = Q(:,a+1:a+N).mx.transpose(conjg(Q(:,a+1:a+N)))
    end select
    !
    R       = TB_build_Rcoord([Nx,Ny],to_home=.false.)
    Xcomm_P = (diag(one*R(:,1)).mx.P) - (P.mx.diag(one*R(:,1)))
    Ycomm_P = (diag(one*R(:,2)).mx.P) - (P.mx.diag(one*R(:,2)))
    Q2      = -4*pi*dimag(((P.mx.Xcomm_P).mx.Ycomm_P))
    Q4      = reshape_rank2_to_rank4(Q2,Nso,Nlat)
    !
    if(allocated(lcm))deallocate(lcm)
    allocate(lcm(Nx,Ny))
    lcm = 0d0
    do ix = 1,Nx
      do iy = 1,Ny
        ilat = ix + (iy-1)*Nx
        lcm(ix,iy) = trace(Q4(:,:,ilat,ilat))
      enddo
    enddo
    !
    call stop_timer()
    !
  end subroutine obc_local_spin_chern_marker





  !####################################################################
  !      Supercell/PBC local/space resolved C and C_\sigma markers
  !####################################################################
  !+------------------------------------------------------------------+
  !Evaluate the local spin Chern marker
  !Ref. Bau'-Marrazzo (2024b)
  ! <https://arxiv.org/abs/2404.04598>`
  !+------------------------------------------------------------------+
  subroutine pbc_local_chern_marker(lcm,method)
    real(8),dimension(:,:),allocatable         :: lcm
    character(len=*),optional                  :: method
    !
    character(len=1)                           :: method_
    real(8),dimension(2)                       :: b1,b2
    complex(8),dimension(Nlso,Nlso)            :: Ub1,Ub2,Umb1,Umb2
    complex(8),dimension(Nlso,Nlso)            :: Vb1,Vb2,Vmb1,Vmb2
    complex(8),dimension(Nlso,Nlso)            :: P,Pgs,Pb1,Pb2,Pmb1,Pmb2
    complex(8),dimension(Nlso,Nlso)            :: Chern_Q2
    real(8),dimension(Nso,Nso,Nlat,Nlat)       :: Chern_Q4
    integer                                    :: i,ix,iy,ilat
    !
    method_='s' ;if(present(method))method_=method
    !
    call check_dimension("PBC_local_chern_marker")
    call assert_shape(U,[size(U,1),Nlso],"PBC_local_chern_marker","H")
    !
    call TB_get_bk(b1,b2)
    !
    call start_timer("pbc_local_ChernMarker")
    !
    Ub1 = periodic_gauge(U,b1)
    Ub2 = periodic_gauge(U,b2)
    Vb1 = dual_state(U,Ub1)
    Vb2 = dual_state(U,Ub2)
    !
    Pgs = ( U(:,1:Nocc).mx.transpose(conjg(U(:,1:Nocc))) )
    Pb1 = ( Vb1(:,1:Nocc).mx.transpose(conjg(Vb1(:,1:Nocc))) )
    Pb2 = ( Vb2(:,1:Nocc).mx.transpose(conjg(Vb2(:,1:Nocc))) )
    P   = (Pb1.mx.Pb2) - (Pb2.mx.Pb1)
    !
    if(to_lower(method_)=='s')then
      Umb1 = periodic_gauge(U,-b1)
      Umb2 = periodic_gauge(U,-b2)
      Vmb1 = dual_state(U,Umb1)
      Vmb2 = dual_state(U,Umb2)
      Pmb1 = ( Vmb1(:,1:Nocc).mx.transpose(conjg(Vmb1(:,1:Nocc))) )
      Pmb2 = ( Vmb2(:,1:Nocc).mx.transpose(conjg(Vmb2(:,1:Nocc))) )
      ! Symmetric differences require all four commutators.
      P    =  ((Pb1.mx.Pb2) - (Pb2.mx.Pb1)) + &
        ((Pmb1.mx.Pmb2) - (Pmb2.mx.Pmb1)) - &
        ((Pb1 .mx.Pmb2) - (Pmb2.mx.Pb1) ) - &
        ((Pmb1.mx.Pb2)  - (Pb2 .mx.Pmb1))
      P    = P/4d0
    endif
    !
    Chern_Q2 = dimag((P.mx.Pgs))/pi2*Nlat
    !
    Chern_Q4 = reshape_rank2_to_rank4(Chern_Q2,Nso,Nlat)
    !
    if(allocated(lcm))deallocate(lcm)
    allocate(lcm(Nx,Ny))
    lcm = 0d0
    do ix = 1,Nx
      do iy = 1,Ny
        ilat = ix + (iy-1)*Nx
        lcm(ix,iy) = trace(Chern_Q4(:,:,ilat,ilat))
      enddo
    enddo
    !
    call stop_timer()
    !
  end subroutine pbc_local_chern_marker


  !+------------------------------------------------------------------+
  !Evaluate the local spin Chern marker
  !Ref. Bau'-Marrazzo (2024b)
  ! <https://arxiv.org/abs/2404.04598>`
  !+------------------------------------------------------------------+
  subroutine pbc_local_spin_chern_marker(spin,lcm,method)
    integer                                    :: spin
    real(8),dimension(:,:),allocatable         :: lcm
    character(len=*),optional                  :: method
    !
    character(len=1)                           :: method_
    real(8),dimension(2)                       :: b1,b2
    complex(8),dimension(:,:),allocatable      :: Q
    complex(8),dimension(:,:),allocatable      :: Ub1,Ub2,Umb1,Umb2
    complex(8),dimension(:,:),allocatable      :: Vb1,Vb2,Vmb1,Vmb2
    complex(8),dimension(Nlso,Nlso)            :: Pb1,Pb2,Pmb1,Pmb2
    complex(8),dimension(Nlso,Nlso)            :: P,Pgs
    complex(8),dimension(Nlso,Nlso)            :: Chern_Q2
    real(8),dimension(Nso,Nso,Nlat,Nlat)       :: Chern_Q4
    integer                                    :: i,ix,iy,ilat,m,N,a
    !
    method_='s' ;if(present(method))method_=method
    !
    if(spin<1.OR.spin>2)stop "PBC_local_spin_chern_marker error: spin < 1 OR spin > 2"
    !
    call check_dimension("PBC_local_spin_chern_marker")
    call assert_shape(U,[size(U,1),Nlso],"PBC_local_spin_chern_marker","U")
    !
    !Start Building LCM_+,-
    call TB_get_bk(b1,b2)
    !
    if(.not.allocated(PSzP))call eigh_PSzP()
    call projected_spin_sector(spin,a,N)
    !
    !
    call start_timer("pbc_local_spinChernMarker_s"//str(spin))
    !
    allocate(Q(Nlso,Nocc))!;Q=zero
    Q = ( U(:,1:Nocc).mx. PSzP )
    !
    !GS projectors Pgs_-, Pgs_+
    select case(spin)
     case(1);Pgs = Q(:,1:N) .mx. transpose(conjg(Q(:,1:N)))
     case(2);Pgs = Q(:,a+1:a+N).mx. &
       transpose(conjg(Q(:,a+1:a+N)))
    end select
    !
    Ub1 = periodic_gauge(Q,b1)
    Ub2 = periodic_gauge(Q,b2)
    Vb1 = dual_state(Q,Ub1,spin)
    Vb2 = dual_state(Q,Ub2,spin)
    Pb1 = Vb1(:,1:N).mx.transpose(conjg(Vb1(:,1:N)))
    Pb2 = Vb2(:,1:N).mx.transpose(conjg(Vb2(:,1:N)))
    P   = (Pb1.mx.Pb2) - (Pb2.mx.Pb1)
    !
    if(to_lower(method_)=='s')then
      Umb1 = periodic_gauge(Q,-b1)
      Umb2 = periodic_gauge(Q,-b2)
      Vmb1 = dual_state(Q,Umb1,spin)
      Vmb2 = dual_state(Q,Umb2,spin)
      Pmb1 = (Vmb1(:,1:N).mx.transpose(conjg(Vmb1(:,1:N))) )
      Pmb2 = (Vmb2(:,1:N).mx.transpose(conjg(Vmb2(:,1:N))) )
      ! Symmetric differences require all four commutators.
      P    =  ((Pb1.mx.Pb2) - (Pb2.mx.Pb1)) + &
        ((Pmb1.mx.Pmb2) - (Pmb2.mx.Pmb1)) - &
        ((Pb1.mx.Pmb2) - (Pmb2.mx.Pb1) )    - &
        ((Pmb1.mx.Pb2)  - (Pb2.mx.Pmb1))
      P    = P/4d0
    endif
    !
    !
    Chern_Q2 = dimag((P.mx.Pgs))/pi2*Nlat
    !
    Chern_Q4 = reshape_rank2_to_rank4(Chern_Q2,Nso,Nlat)
    !
    if(allocated(lcm))deallocate(lcm)
    allocate(lcm(Nx,Ny))
    lcm = 0d0
    do ix = 1,Nx
      do iy = 1,Ny
        ilat = ix + (iy-1)*Nx
        lcm(ix,iy) = trace(Chern_Q4(:,:,ilat,ilat))
      enddo
    enddo
    !
    call stop_timer()
    !
  end subroutine pbc_local_spin_chern_marker











  !##################################################################
  !                 COMPUTATIONAL ROUTINES
  !##################################################################
  !------------------------------------------------------------------
  !Build the periodic gauge:  |u_j(b)> = e^{+ib.R}|u_j(0)>
  !using Bloch eigenstates (U) and reciprocal basis states b
  !(b can be any k-point indeed)
  !------------------------------------------------------------------
  function Periodic_Gauge(U,b) result(Ub)
    complex(8),dimension(:,:),intent(in)  :: U
    real(8),dimension(2),intent(in)       :: b
    complex(8),dimension(:,:),allocatable :: Ub
    real(8),dimension(:,:),allocatable    :: R
    real(8),dimension(:),allocatable      :: RdotB
    complex(8),dimension(:),allocatable   :: Varg
    integer                               :: i,j,N,nrow,ninternal
    !
    nrow = size(U,1)
    N    = size(U,2)
    if(Nlat<1.or.modulo(nrow,Nlat)/=0) &
      error stop "Periodic_Gauge: invalid U row count"
    ninternal = nrow/Nlat
    !
    !Build the position array:
    R = TB_build_Rcoord([Nx,Ny],ninternal=ninternal)
    allocate(RdotB(nrow),Varg(nrow))
    forall(i=1:nrow)RdotB(i) = dot_product(R(i,:),b)
    Varg = dcmplx(cos(RdotB),sin(RdotB))
    !
    if(allocated(Ub))deallocate(Ub)
    allocate(Ub(nrow,N));Ub=zero
    forall(j=1:N)Ub(:,j) = Varg(:)*U(:,j)
  end function Periodic_Gauge




  !------------------------------------------------------------------
  ! Build dual states for all columns of an explicitly supplied
  ! occupied basis. U and Ub have shape [Nrow,Noccupied].
  !------------------------------------------------------------------
  function Dual_State_Basis(U,Ub) result(Vb)
    complex(8),dimension(:,:),intent(in)  :: U
    complex(8),dimension(:,:),intent(in)  :: Ub
    complex(8),dimension(:,:),allocatable :: Vb
    complex(8),dimension(:,:),allocatable :: S
    integer                               :: nrow,nocc
    !
    nrow = size(U,1)
    nocc = size(U,2)
    if(any(shape(Ub)/=[nrow,nocc])) &
      error stop "Dual_State_Basis: U and UB shapes differ"
    if(nocc<1.or.nocc>nrow) &
      error stop "Dual_State_Basis: invalid occupied basis"
    allocate(S(nocc,nocc))
    S = transpose(conjg(U)).mx.Ub
    call inv(S)
    allocate(Vb(nrow,nocc))
    Vb = Ub.mx.S
  end function Dual_State_Basis



  !------------------------------------------------------------------
  ! Build the dual state implementing the parallel transport:
  !  !|\~{u}_j(b)⟩ =\sum_{a=1,...,Nocc} S^{−1}_{aj}(b)|u_a(b)⟩
  !
  ! Use Bloch eigenstates (U) and Periodic gage (Ub).
  ! Spin=0,1,2. 0=spinless, 1,2=up,dw
  !
  ! S_mn = <u_m0|u_nb> = <u_m0|e^{-ib.R}|u_n0>
  !------------------------------------------------------------------
  function Dual_State(U,Ub,spin) result(Vb)
    complex(8),dimension(:,:)             :: U
    complex(8),dimension(:,:)             :: Ub
    integer,optional                      :: spin
    complex(8),dimension(:,:),allocatable :: Vb
    complex(8),dimension(:,:),allocatable :: S
    integer                               :: spin_
    integer                               :: N,a
    !
    spin_ = 0;if(present(spin))spin_=spin
    if(spin_<0.OR.spin_>2)stop "dual_state error: spin < 0 OR spin > 2"
    !
    N = Nocc
    a = 0
    if(spin_>0)call projected_spin_sector(spin_,a,N)
    !
    allocate(S(N,N));S=zero
    S = ( transpose(conjg(U(:,a+1:a+N))).mx. Ub(:,a+1:a+N) )
    !
    call inv(S)
    !
    if(allocated(Vb))deallocate(Vb)
    allocate(Vb(Nlso,N))
    Vb = ( Ub(:,a+1:a+N).mx.S)
    !
  end function Dual_State




  !------------------------------------------------------------------
  ! Return the column offset and rank of one spectral sector of
  ! Uocc^+ Sz Uocc. EIGH orders EPSP increasingly, so SPIN=1 is the
  ! negative sector and SPIN=2 is the positive sector. The ranks need
  ! not be equal when time-reversal symmetry is broken.
  !------------------------------------------------------------------
  subroutine projected_spin_sector(spin,offset,nsector)
    integer,intent(in)  :: spin
    integer,intent(out) :: offset,nsector
    real(8),parameter   :: spin_tol=1d-10
    integer             :: nminus,nplus,nzero
    !
    if(spin<1.or.spin>2) &
      error stop "projected_spin_sector: SPIN must be 1 or 2"
    if(.not.allocated(Epsp).or.size(Epsp)/=Nocc) &
      error stop "projected_spin_sector: EPSP is unavailable"
    nminus = count(Epsp < -spin_tol)
    nplus  = count(Epsp >  spin_tol)
    nzero  = Nocc-nminus-nplus
    if(nzero>0) &
      error stop "projected_spin_sector: projected spin gap is closed"
    select case(spin)
    case(1)
      offset  = 0
      nsector = nminus
    case(2)
      offset  = nminus
      nsector = nplus
    end select
    if(nsector<1) &
      error stop "projected_spin_sector: empty spin sector"
  end subroutine projected_spin_sector




  !------------------------------------------------------------------
  ! Build the Bloch ground state projected representation of the
  ! spin operators Sz: PSzP = P^+ . Sz. P
  ! where P = \sum_occ|gs><gs| is the projection onto occupied states
  ! PSzP_{a'b'}  = [U^+]_{a'a} Sz_{ab} U_{b,b'}
  !------------------------------------------------------------------
  subroutine eigh_PSzP()
    if(allocated(PSzP))deallocate(PSzP)
    if(allocated(Epsp))deallocate(Epsp)
    allocate(PSzP(Nocc,Nocc))
    allocate(Epsp(Nocc))
    PSzP = conjg(transpose(U(:,1:Nocc))).mx. (Sz.mx.U(:,1:Nocc))
    call eigh(PSzP,Epsp)
    call save_array("Epszp.dat",Epsp)
    call check_Pgap("eigh_PSzP")
    !
  end subroutine eigh_PSzP
  !
  function PSzP_Matrix() result(PSzP)
    complex(8),dimension(:,:),allocatable :: PSzP
    PSzP = ( conjg(transpose(U(:,1:Nocc))).mx. (Sz.mx.U(:,1:Nocc)) )
  end function PSzP_Matrix








  !------------------------------------------------------------------
  ! Get the chemical potential: unused. Can be implemented using fzero
  !------------------------------------------------------------------
  function chemical_potential(e,beta,nocc) result(mu)
    real(8),dimension(:) :: e
    real(8)              :: beta
    integer              :: nocc
    real(8)              :: mu
    real(8)              :: a, b, n,err
    integer              :: info,iter
    integer,parameter    :: Niter=200
    real(8),parameter    :: eps=1d-7
    a = minval(e); b = maxval(e)
    do iter=1,200
      mu = 0.5d0*(a + b)
      n  = sum(fermi(e-mu, beta))
      err= n - nocc
      if(err<0d0)then
        a = mu
      else
        b = mu
      endif
      print*,iter,err
      if(abs(err)<eps)return
    enddo
    stop "ERROR chemical_potential: failed after 200 iterations"
  end function chemical_potential




  !------------------------------------------------------------------
  ! Smearing the Fermi step function using temperature cutoff
  !------------------------------------------------------------------
  function smearing(state,u,evals,beta,mu) result(weight)
    complex(8),dimension(:)   :: state
    complex(8),dimension(:,:) :: u
    real(8),dimension(:)      :: evals
    real(8)                   :: beta
    real(8)                   :: mu
    real(8)                   :: weight,ww,fw
    integer                   :: i,N
    N = size(u,2)
    weight = 1d0
    if(beta<1d3)then
      weight  = 0d0
      do i=1,N
        ww = dot_product(state, u(:,i))
        fw = fermi(evals(i)-mu,beta)
        weight = weight+fw*abs(ww)**2
      enddo
    endif
  end function smearing








  !------------------------------------------------------------------
  ! Returns the coordinates of the orbitals in the supercell.
  ! >> THIS WORKS ONLY FOR THE SQUARE LATTICE SO FAR <<
  !------------------------------------------------------------------
  function TB_build_Rcoord(Nvec,to_home,ninternal) result(Rcoord)
    integer,dimension(:)               :: Nvec     ![Nx,Ny]
    logical,optional                   :: to_home
    integer,optional                   :: ninternal
    real(8),dimension(:,:),allocatable :: Rcoord   ![Nx*Ny*Nso,D=2]
    integer                            :: Nlat     !
    real(8),dimension(:,:),allocatable :: Rgrid    ![Nx*Ny,2]
    integer                            :: Ndim,Ngrd
    real(8),dimension(2)               :: e1,e2
    integer                            :: igrd,io,i,nint
    integer                            :: ir,ix,iy,iz,Nr(3),D
    logical :: rescale
    !
    rescale = .true. ;if(present(to_home))rescale=to_home
    nint = Nso
    if(present(ninternal))nint=ninternal
    if(nint<1)error stop "TB_build_Rcoord: invalid NINTERNAL"
    !
    Ndim = size(Nvec)
    Ngrd = product(Nvec)
    !
    if(allocated(Rcoord))deallocate(Rcoord)
    allocate(Rcoord(Ngrd*nint,Ndim))
    !
    D = size(Nvec)
    Nr=1
    do ir=1,D
      Nr(ir)=Nvec(ir)
    enddo
    if(product(Nr)/=product(Nvec))stop "TB_build_Rcoord ERROR: product(Nrvec) != product(Nr)"
    !
    call TB_get_ei(e1,e2)
    !
    ir=0
    do iy=1,Nr(2)
      do ix=1,Nr(1)
        do io=1,nint
          ir=ir+1
          Rcoord(ir,:) = dble(ix)*e1(1:D) + dble(iy)*e2(1:D)
          if(rescale)Rcoord(ir,:) = dble(ix)/Nr(1)*e1(1:D) + dble(iy)/Nr(2)*e2(1:D)
        enddo
      enddo
    enddo
  end function TB_build_Rcoord








END MODULE BHZ_LCM










 ! !------------------------------------------------------------------
 ! ! calcola il numero di chern localmente nello spazio reale,
 ! ! alla bianco & resta, per strip obc su x e pbc su y.
 ! ! prende in input i proiettori p e q sulle bande occupate complessive,
 ! ! che sono matrici nx*nky x nx*nky
 ! ! get_local_chern is the unsymmetric version requiring less memory
 ! ! c(r) = -4 pi (2/v_uc) im tr (x_p y_q)
 ! !------------------------------------------------------------------
 ! subroutine local_chern_marker(U,lcm)
 !   complex(8),dimension(:,:),intent(in)      :: U    ![Nlso,Nlso]
 !   real(8),dimension(:,:),allocatable        :: lcm  ![Nx,Ny]
 !   !
 !   real(8),dimension(:,:),allocatable        :: R    ![Nlso,D=2]
 !   complex(8),dimension(:,:),allocatable     :: P,Q  ![Nlso,Nlso]
 !   complex(8),dimension(:,:,:,:),allocatable :: Q4   ![Nso,Nso,Nlat,Nlat]
 !   integer                                   :: i,j,ilat
 !   real(8)                                   :: Rx1,Ry2
 !   !
 !   call check_dimension("local_chern_marker")
 !   call assert_shape(U,[size(U,1),Nlso],"local_chern_marker","H")
 !   !
 !   allocate(P(Nlso,Nlso),Q(Nlso,Nlso))
 !   P = zero
 !   Q = zero
 !   !
 !   do j=1,Nocc
 !      P = P + outerprod(U(:,j),conjg(U(:,j)))
 !   enddo
 !   Q = zeye(Nlso) - P
 !   !
 !   allocate(R(Nlso,2))
 !   R = TB_build_Rcoord([Nx,Ny],to_home=.false.)
 !   do concurrent(i=1:Nlso,j=1:Nlso)
 !      rx1 = R(i,1)
 !      ry2 = R(j,2)
 !      Q(i,j) = rx1*ry2*Q(i,j)
 !   enddo
 !   !
 !   Q = matmul(P, matmul(Q,P))
 !   !
 !   allocate(Q4(Nso,Nso,Nlat,Nlat))
 !   Q4 = reshape_rank2_to_rank4(Q,Nso,Nlat)
 !   !
 !   if(allocated(lcm))deallocate(lcm)
 !   allocate(lcm(Nx,Ny))
 !   lcm=0d0
 !   do i = 1,Nx
 !      do j = 1,Ny
 !         ilat = i + (j-1)*Nx
 !         lcm(i,j) = lcm(i,j) + dimag(trace(Q4(:,:,ilat,ilat)) )*4*pi
 !      enddo
 !   enddo
 !   !
 !   return
 ! end subroutine local_chern_marker
