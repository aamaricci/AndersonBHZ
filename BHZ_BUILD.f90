MODULE BHZ_BUILD
  USE COMMON
  implicit none
  private

  complex(8),dimension(:,:,:),allocatable,public   :: Hk
  complex(8),dimension(:,:,:),allocatable,public   :: Hij
  complex(8),dimension(:,:,:,:),allocatable,public :: Hlat
  complex(8),dimension(:,:),allocatable,public     :: Sz
  integer,dimension(:,:),allocatable               :: Links
  real(8),dimension(:),allocatable                 :: erandom
  integer                                          :: disorder_seed=0

  public :: hk_model
  !
  public :: setup_Abhz
  public :: free_Abhz
  !
  public :: compose_Abhz_spH
  public :: build_mf_Abhz



contains


  ! Return the mean-field Hamiltonian associated with the converged
  ! local fields A(:,1:2)=[Tz,Sz]. SPH0 is one 2*Nlat spin block for
  ! SPIN=+1/-1, or the complete 4*Nlat matrix for SPIN=0. The returned
  ! CSR matrix owns its storage; SPH0 and A are left unchanged.
  function build_mf_Abhz(spH0,a,spin) result(spHmf)
    type(sparse_zmatrix_csr),intent(in) :: spH0
    real(8),intent(in)                  :: a(:,:)
    integer,intent(in)                  :: spin
    type(sparse_zmatrix_csr)            :: spHmf
    type(sparse_zmatrix_csr)            :: correction
    complex(8),allocatable              :: block(:,:)
    real(8)                             :: orbital_coupling
    real(8)                             :: spin_coupling
    integer                             :: ilat,nblock
    !
    ! Check dimensions before allocating the correction.
    if(abs(spin)>1) &
      error stop "build_mf_Abhz: SPIN must be -1, 0 or +1"
    if(size(a,1)/=Nlat.or.size(a,2)/=2) &
      error stop "build_mf_Abhz: invalid A shape"
    nblock = merge(Nso,2,spin==0)
    if(spH0%Nrow/=Nlat*nblock.or.spH0%Ncol/=spH0%Nrow) &
      error stop "build_mf_Abhz: invalid SPH0 shape"
    !
    ! Build the local Hartree-Fock contribution.
    orbital_coupling = (Uloc-5d0*Jh)/2d0
    spin_coupling    = (Uloc+Jh)/2d0
    allocate(block(nblock,nblock))
    call correction%init(Nlat*nblock,Nlat*nblock)
    do ilat=1,Nlat
      if(spin==0)then
        block = -a(ilat,1)*orbital_coupling*Gamma5 &
          -a(ilat,2)*spin_coupling*GammaS
      else
        block = -a(ilat,1)*orbital_coupling*pauli_tau_z &
          -spin*a(ilat,2)*spin_coupling*pauli_tau_0
      endif
      call insert_local_block(correction,block,ilat)
    enddo
    !
    ! Keep both sparse operands alive while materializing their sum.
    spHmf = spH0 + correction
    call correction%free()
    deallocate(block)
  end function build_mf_Abhz


  ! Insert one dense block on site ILAT in a caller-owned CSR matrix.
  subroutine insert_local_block(spH,block,ilat)
    type(sparse_zmatrix_csr),intent(inout) :: spH
    complex(8),intent(in)                  :: block(:,:)
    integer,intent(in)                     :: ilat
    integer                                :: i,j,io,jo,nblock
    !
    if(size(block,1)/=size(block,2)) &
      error stop "insert_local_block: BLOCK must be square"
    if(ilat<1.or.ilat>Nlat) &
      error stop "insert_local_block: invalid ILAT"
    nblock = size(block,1)
    do io=1,nblock
      i = io+(ilat-1)*nblock
      do jo=1,nblock
        if(block(io,jo)==zero)cycle
        j = jo+(ilat-1)*nblock
        call spH%fast_insert(block(io,jo),i,j)
      enddo
    enddo
  end subroutine insert_local_block



  ! Return the clean 4x4 Bloch Hamiltonian H(k) of the square BHZ
  ! model. This dense momentum-space form is used by translationally
  ! invariant calculations; disorder is introduced only in real space.
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







  ! Build the square BHZ Hamiltonian for IDUM (serial).
  ! With SPH:
  ! SPIN=+1/-1 returns one 2*Nlat spin sector and SPIN=0 returns the
  ! complete 4*Nlat Hamiltonian. SPIN defaults to zero. Without SPH,
  ! populate the dense COMMON arrays for legacy callers.
  ! VERBOSE defaults to true and controls logging and disorder output.
  ! A matching erandom_<IDUM>.restart overrides the generated disorder.
  subroutine setup_Abhz(verbose,spH,spin)
    logical,optional,intent(in)                   :: verbose
    type(sparse_zmatrix_csr),optional,intent(out) :: spH
    integer,optional,intent(in)                   :: spin
    logical                                       :: verbose_,bool,new_disorder
    integer                                       :: ilat
    integer                                       :: spin_
    complex(8)                                    :: channel(4,4)
    !
    verbose_ = .true. ; if(present(verbose))verbose_=verbose
    spin_ = 0         ; if(present(spin))spin_=spin
    !
    call check_dimension("setup_Abhz")
    !
    if(abs(spin_)>1) error stop "setup_Abhz: SPIN must be -1, 0 or +1"
    if(.not.present(spH).and.spin_/=0)error stop "setup_Abhz: dense legacy setup requires SPIN=0"
    !
    select case(disorder_type)
     case(0);channel = Gamma0
     case(1);channel = Gamma5
     case(2);channel = GammaS
     case default
      error stop "setup_Abhz: DISORDER_TYPE must be in [0,2]"
    end select
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
      else
        ! The realized field is the compact authoritative input for
        ! later spectral post-processing. It also protects old runs
        ! from future changes in the random-number call sequence.
        inquire(file='erandom_'//str(idum)//'.used',exist=bool)
        if(bool)then
          if(file_length('erandom_'//str(idum)//'.used')/=Nlat) &
            error stop "setup_Abhz: invalid disorder used length"
          call read_array('erandom_'//str(idum)//'.used',erandom)
        endif
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
    !
    if(present(spH))then
      if(spin_==0)then
        call build_Abhz_sparse_full(spH,channel)
      else
        call build_Abhz_sparse_block(spH,spin_)
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


  ! Release the dense model, lattice links, spin operator and disorder
  ! realization allocated by setup_Abhz. Call this before changing IDUM.
  subroutine free_Abhz
    if(allocated(Hlat))deallocate(Hlat)
    if(allocated(Hij))deallocate(Hij)
    if(allocated(erandom))deallocate(erandom)
    if(allocated(Links))deallocate(Links)
    if(allocated(Sz))deallocate(Sz)
    disorder_seed = 0
  end subroutine free_Abhz








  ! Insert each physical column once and in increasing order.
  ! For Nx=2 with PBC, opposite links can reach the same site.
  subroutine build_Abhz_sparse_full(spH,channel)
    type(sparse_zmatrix_csr),intent(inout) :: spH
    complex(8),intent(in)                  :: channel(4,4)
    complex(8)                             :: blocks(Nso,Nso,0:4),block(Nso,Nso)
    integer                                :: coords(Nlat,2),target(0:4),r(2)
    integer                                :: ilat,jlat,link,io,jo,i,j
    !
    call TB_build_CoordGrid([Nx,Ny],coords)
    !
    call spH%init(Nlso,Nlso)
    !
    do ilat=1,Nlat
      target(0) = ilat
      blocks(:,:,0) = ts_model(0,Nso) + erandom(ilat)*channel
      do link=1,4
        r = coords(ilat,:) + Links(link,:)
        if(bhz_pbc)r = modulo(r-1,[Nx,Ny])+1
        target(link) = TB_find_IndxCoord(r,[Nx,Ny])
        blocks(:,:,link) = ts_model(link,Nso)
      enddo
      !
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
        !
      enddo
      !
    enddo
  end subroutine build_Abhz_sparse_full

  ! Build one 2*Nlat spin sector. New CSR columns are inserted in
  ! increasing order and coincident periodic links are accumulated.
  subroutine build_Abhz_sparse_block(spH,spin)
    type(sparse_zmatrix_csr),intent(inout) :: spH
    integer,intent(in)                     :: spin
    complex(8)                             :: blocks(2,2,0:4),block(2,2)
    complex(8)                             :: disorder_block(2,2)
    integer                                :: coords(Nlat,2)
    integer                                :: target(0:4),r(2)
    integer                                :: ilat,jlat,link,io,jo,i,j
    !
    select case(disorder_type)
     case(0);disorder_block = pauli_tau_0
     case(1);disorder_block = pauli_tau_z
     case(2);disorder_block = spin*pauli_tau_0
     case default
      error stop "build_Abhz_sparse_block: invalid DISORDER_TYPE"
    end select
    !
    call TB_build_CoordGrid([Nx,Ny],coords)
    !
    call spH%init(2*Nlat,2*Nlat)
    !
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
  end subroutine build_Abhz_sparse_block


  ! Return the clean 4x4 on-site term or directed nearest-neighbour
  ! hopping selected by LINK. This dense block is used to assemble the
  ! full real-space BHZ Hamiltonian.
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

  ! Return one conserved-spin 2x2 hopping block in orbital space.
  function ts_model_block(link,spin) result(Hts)
    integer,intent(in)       :: link,spin
    complex(8),dimension(2,2) :: Hts
    !
    select case(link)
     case(0) !LOCAL PART
      Hts = Mh*pauli_tau_z
     case(1) !RIGHT HOPPING
      Hts = -0.5d0*pauli_tau_z + xi*0.5d0*lambda*spin*pauli_tau_x
     case(2) !UP HOPPING
      Hts = -0.5d0*pauli_tau_z - xi*0.5d0*lambda*pauli_tau_y
     case(3) !LEFT HOPPING
      Hts = -0.5d0*pauli_tau_z - xi*0.5d0*lambda*spin*pauli_tau_x
     case(4) !DOWN HOPPING
      Hts = -0.5d0*pauli_tau_z + xi*0.5d0*lambda*pauli_tau_y
     case default
      error stop "ts_model_block: LINK must be in [0,4]"
    end select
  end function ts_model_block













  ! Assemble a sparse site-major 4*Nlat Hamiltonian from the two
  ! conserved-spin CSR blocks. This is an exact block embedding: the
  ! spin-off-diagonal matrix elements are set to zero.
  function compose_Abhz_spH(spHup,spHdown) result(spHfull)
    type(sparse_zmatrix_csr),intent(in)  :: spHup,spHdown
    type(sparse_zmatrix_csr)             :: spHfull
    integer                              :: spin,ib,jb,ifull,jfull
    integer                              :: ilat,jlat,iorb,jorb,pos
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
  end function compose_Abhz_spH








end module BHZ_BUILD
