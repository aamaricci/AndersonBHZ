!> Solve the 2D BHZ + Kanamori n.n + Disorder model using Mean-Field
program MF_bhz_2d_Anderson
  USE COMMON
  USE BHZ_SOLVE
  implicit none
  integer,parameter                             :: Norb=2,Nspin=2
  integer                                       :: io,ispin,iorb,ilat
  type(sparse_zmatrix_csr)                      :: spHij_up,spHij_down
  type(sparse_zmatrix_csr)                      :: spHij_full,spHmf
  integer                                       :: Iter,Nsuccess=2
  logical                                       :: iexist,converged
  real(8),dimension(:,:),allocatable            :: params,params_prev

  real(8)                                       :: mf_error

  complex(8)                                    :: spectrum_z

  call parse_input_file(inputFILE,"inputBHZ.conf")
  call read_input(inputFILE,"mf_Abhz")

  !SETUP THE GAMMA MATRICES:
  call setup_GammaMatrices()

  !SETUP THE OUTPUT DIRECTORY AND THE SPARSE BHZ HAMILTONIAN
  call set_seed_wd()
  call save_input_file(file_basename(inputFILE))
  call TB_set_ei([1d0,0d0],[0d0,1d0])
  call TB_set_bk([pi2,0d0],[0d0,pi2])
  select case(spin_blocks)
  case(1)
    call setup_Abhz(spH=spHij_up,spin=1)
  case(2)
    call setup_Abhz(spH=spHij_up,spin=1)
    call setup_Abhz(verbose=.false.,spH=spHij_down,spin=-1)
  case(0)
    call setup_Abhz(spH=spHij_full,spin=0)
  end select
  !
  !
  !Start MF HERE:
  allocate(params(Nlat,2))
  allocate(params_prev(Nlat,2))
  do ilat=1,Nlat
    params(ilat,:) = [sb_field,sb_field]   ![Tz,Sz]
  enddo
  if(spin_blocks==1)params(:,2)=0d0
  inquire(file="params.restart",exist=iexist)
  if(iexist)then
    call read_array("params.restart",params)
    if(spin_blocks==1)then
      params(:,2) = 0d0
    else
      params(:,2) = params(:,2)+sb_field
    endif
  endif
  call save_array("params.init",params)
  !
  !
  !DO MEAN-FIELD CYCLE
  open(100,file="tz_sz.dat",status="replace")
  open(101,file="dens_1.dat",status="replace")
  open(102,file="dens_2.dat",status="replace")
  open(103,file="mf_history.dat",status="replace")
  write(103,'(A)') "# 1.iter 2.error 3.mean_Tz 4.mean_Sz"
  converged = .false. ; iter=0
  do while(.not.converged.AND.iter<maxiter)
    iter = iter+1
    call start_loop(iter,maxiter,"MF-loop")
    !
!   call symmetrize_params(params)
    call solve_MF_bhz(params)
    !
    if(iter>1)params = wmix*params + (1d0-wmix)*params_prev
    params_prev = params
    !
    converged = check_convergence_local(params,it_error,nsuccess, &
      maxiter,oerr=mf_error)
    write(103,*)iter,mf_error,sum(params(:,1))/Nlat, &
      sum(params(:,2))/Nlat
    !
    call end_loop
  end do
  call save_array("params.restart",params)
  close(100)
  close(101)
  close(102)
  close(103)
  if(mf_error>=it_error) stop 1
  !
  !
  if(with_spectrum) call write_spectrum()
  if(spHij_up%status)call spHij_up%free()
  if(spHij_down%status)call spHij_down%free()
  if(spHij_full%status)call spHij_full%free()
  if(spHmf%status)call spHmf%free()
  call free_Abhz_solution()
  call free_Abhz()
  call reset_seed_wd()




contains




  subroutine solve_MF_bhz(a)
    real(8),dimension(Nlat,2),intent(inout) :: a
    real(8),dimension(:),allocatable        :: dens_up,dens_down
    real(8),dimension(:),allocatable        :: dens_full
    real(8)                                 :: dens(Nlat,Nspin,Norb)
    integer                                 :: iorb,ispin
    !
    dens = 0d0
    select case(spin_blocks)
    case(1)
      call set_mf_hamiltonian(spHij_up,a,spin=1)
      call solve_Abhz(spHmf,dens_up)
      allocate(dens_down,source=dens_up)
    case(2)
      call set_mf_hamiltonian(spHij_up,a,spin=1)
      call solve_Abhz(spHmf,dens_up)
      call set_mf_hamiltonian(spHij_down,a,spin=-1)
      call solve_Abhz(spHmf,dens_down)
    case(0)
      call set_mf_hamiltonian(spHij_full,a,spin=0)
      call solve_Abhz(spHmf,dens_full)
    end select

    if(spin_blocks/=0)then
      do concurrent(ilat=1:Nlat,iorb=1:Norb)
        io = iorb+(ilat-1)*Norb
        dens(ilat,1,iorb) = dens_up(io)
        dens(ilat,2,iorb) = dens_down(io)
      enddo
    else
      do concurrent(ilat=1:Nlat,ispin=1:Nspin,iorb=1:Norb)
        io = iorb+(ispin-1)*Norb+(ilat-1)*Nspin*Norb
        dens(ilat,ispin,iorb) = dens_full(io)
      enddo
    endif
    !
    do ilat=1,Nlat
      a(ilat,1) = 0.5d0*sum(dens(ilat,:,1)) &
        -0.5d0*sum(dens(ilat,:,2))
      a(ilat,2) = 0.5d0*sum(dens(ilat,1,:)) &
        -0.5d0*sum(dens(ilat,2,:))
    enddo
    if(spin_blocks==1)a(:,2)=0d0
    !
    print*," Tz=",sum(a(:,1))/Nlat," Sz=",sum(a(:,2))/Nlat
    rewind(100);rewind(101);rewind(102)
    write(100,'(A)') "# 1.Tz 2.Sz"
    write(101,'(A)') "# 1.n_up_orb1 2.n_down_orb1"
    write(102,'(A)') "# 1.n_up_orb2 2.n_down_orb2"
    do ilat=1,Nlat
      write(100,*)a(ilat,1),a(ilat,2)
      write(101,*)(dens(ilat,ispin,1),ispin=1,Nspin)
      write(102,*)(dens(ilat,ispin,2),ispin=1,Nspin)
    enddo
    if(allocated(dens_up))deallocate(dens_up)
    if(allocated(dens_down))deallocate(dens_down)
    if(allocated(dens_full))deallocate(dens_full)
  end subroutine solve_MF_bhz





  subroutine set_mf_hamiltonian(spH0,a,spin)
    ! Build the current mean-field Hamiltonian for one spin sector
    ! (SPIN=+1/-1) or for the complete site-major basis (SPIN=0).
    type(sparse_zmatrix_csr),intent(in)     :: spH0
    real(8),dimension(Nlat,2),intent(in)    :: a
    integer,intent(in)                     :: spin
    !
    if(abs(spin)>1) &
      error stop "set_mf_hamiltonian: SPIN must be -1, 0 or +1"
    if(spHmf%status)call spHmf%free()
    spHmf = spH0 + mf_Hij_correction(a,spin)
  end subroutine set_mf_hamiltonian





  function mf_Hij_correction(a,spin) result(spdHmf)
    ! Return the local mean-field correction in one 2x2 spin sector
    ! (SPIN=+1/-1) or in the complete 4x4 basis (SPIN=0).
    real(8),dimension(Nlat,2),intent(in) :: a
    integer,intent(in)                  :: spin
    type(sparse_zmatrix_csr)            :: spdHmf
    complex(8),allocatable              :: dHmf(:,:)
    real(8)                             :: orbital_coupling,spin_coupling
    integer                             :: ilat,nblock
    !
    if(abs(spin)>1) &
      error stop "mf_Hij_correction: SPIN must be -1, 0 or +1"
    orbital_coupling = (Uloc-5d0*Jh)/2d0
    spin_coupling    = (Uloc+Jh)/2d0
    nblock = merge(Nso,Norb,spin==0)
    allocate(dHmf(nblock,nblock))
    call spdHmf%init(Nlat*nblock,Nlat*nblock)
    do ilat=1,Nlat
      if(spin==0)then
        dHmf = -a(ilat,1)*orbital_coupling*Gamma5 &
          -a(ilat,2)*spin_coupling*GammaS
      else
        dHmf = -a(ilat,1)*orbital_coupling*pauli_tau_z &
          -spin*a(ilat,2)*spin_coupling*pauli_tau_0
      endif
      call sp_insert_block(spdHmf,dHmf,ilat,ilat)
    enddo
    deallocate(dHmf)
  end function mf_Hij_correction




  ! Post-processing uses the final mixed fields in params.restart.
  !We could use DMFT_TOOLS for this...
  subroutine write_spectrum()
    real(8)                             :: omega(spectrum_points)
    complex(8),allocatable             :: Gdiag(:,:,:,:),Gmean(:,:,:)
    complex(8),allocatable             :: Gup(:,:),Gdown(:,:),Gfull(:,:)
    integer                             :: iw,iso
    !
    omega = linspace(wmin,wmax,spectrum_points)
    allocate(Gdiag(Nlat,Nso,Nso,spectrum_points))
    allocate(Gmean(Nso,Nso,spectrum_points))
    Gdiag = zero
    !
    select case(spin_blocks)
    case(1)
      call set_mf_hamiltonian(spHij_up,params,spin=1)
      call current_sector_spectrum(omega,Gup)
      allocate(Gdown,source=Gup)
    case(2)
      call set_mf_hamiltonian(spHij_up,params,spin=1)
      call current_sector_spectrum(omega,Gup)
      call set_mf_hamiltonian(spHij_down,params,spin=-1)
      call current_sector_spectrum(omega,Gdown)
    case(0)
      call set_mf_hamiltonian(spHij_full,params,spin=0)
      call current_sector_spectrum(omega,Gfull)
    end select

    if(spin_blocks/=0)then
      do ilat=1,Nlat
        do iorb=1,Norb
          io = iorb+(ilat-1)*Norb
          Gdiag(ilat,iorb,iorb,:) = Gup(io,:)
          iso = iorb+Norb
          Gdiag(ilat,iso,iso,:) = Gdown(io,:)
        enddo
      enddo
    else
      do ilat=1,Nlat
        do iso=1,Nso
          io = iso+(ilat-1)*Nso
          Gdiag(ilat,iso,iso,:) = Gfull(io,:)
        enddo
      enddo
    endif
    !
    Gmean = sum(Gdiag,1)/Nlat
    call write_gf(dcmplx(omega,ieta),"real")
    call write_gf(Gdiag,"G_ii",axis='real',iprint=1)
    call write_gf(Gmean,"G",axis='real',iprint=1)
    if(allocated(Gup))deallocate(Gup)
    if(allocated(Gdown))deallocate(Gdown)
    if(allocated(Gfull))deallocate(Gfull)
    deallocate(Gdiag,Gmean)
  end subroutine write_spectrum





  subroutine current_sector_spectrum(omega,Gsector)
    ! Return the diagonal resolvent of the active sparse mean-field
    ! sector on OMEGA. The result has shape [spHmf%Nrow,size(OMEGA)].
    real(8),dimension(:),intent(in)                 :: omega
    complex(8),dimension(:,:),allocatable,intent(out) :: Gsector
    complex(8),dimension(:,:),allocatable           :: Hmf,moments
    complex(8),dimension(:),allocatable             :: Vseed
    complex(8),allocatable                          :: local_mom(:,:)
    complex(8)                                      :: coeff(0:spectrum_order)
    real(8),dimension(:),allocatable                :: Evals
    real(8)                                         :: scale,shift
    integer                                         :: n,iw,io,info,nop,numop
    !
    n = spHmf%Nrow
    if(n<1)error stop "current_sector_spectrum: empty Hamiltonian"
    allocate(Gsector(n,size(omega)))
    select case(method)
    case("lapack")
      allocate(Hmf(n,n),Evals(n))
      Hmf = as_matrix(spHmf)
      call eigh(Hmf,Evals)
      do iw=1,size(omega)
        do io=1,n
          Gsector(io,iw) = sum(abs(Hmf(io,:))**2/ &
            (dcmplx(omega(iw)+mu,ieta)-Evals))
        enddo
      enddo
      deallocate(Hmf,Evals)
    case("chebychev")
      call cheb_rescaling(spHmf,scale,shift)
      call cheb_set(scale,shift)
      allocate(Vseed(n),local_mom(1,0:spectrum_order))
      allocate(moments(n,0:spectrum_order))
      numop = 0
      do io=1,n
        Vseed     = zero
        Vseed(io) = one
        call cheb_moments(MatVec_Hmf,Vseed,[io],local_mom,info,nop)
        if(info/=0) &
          error stop "current_sector_spectrum: moment recurrence failed"
        moments(io,:) = local_mom(1,:)
        numop         = numop+nop
      enddo
      do iw=1,size(omega)
        spectrum_z = dcmplx(omega(iw)+mu,ieta)
        call cheb_coefficients(spectral_resolvent,scale,shift,coeff, &
          kind="none",tol=cheb_tol,activate=.false.)
        Gsector(:,iw) = matmul(moments,coeff)
      enddo
      deallocate(Vseed,local_mom,moments)
      write(*,*) "Spectrum: degree, eta, scale, shift, MatVec calls =", &
        spectrum_order,ieta,scale,shift,numop
    end select
  end subroutine current_sector_spectrum





  subroutine sp_insert_block(self,a,ilat,jlat)
    class(sparse_zmatrix_csr),intent(inout) :: self
    complex(8),dimension(:,:),intent(in)    :: a
    integer,intent(in)                      :: ilat,jlat
    integer                                 :: i,j,io,jo,nblock
    !
    if(size(a,1)/=size(a,2)) &
      error stop "sp_insert_block: A must be square"
    nblock = size(a,1)
    do io=1,nblock
      do jo=1,nblock
        i = io + (ilat-1)*nblock
        j = jo + (jlat-1)*nblock
        if(a(io,jo)==zero)cycle
        call self%insert(a(io,jo),i,j)
      enddo
    enddo
  end subroutine sp_insert_block







  function spectral_resolvent(E) result(G)
    real(8),intent(in)                          :: E
    complex(8)                                  :: G
    G = one/(spectrum_z-E)
  end function spectral_resolvent

  subroutine MatVec_Hmf(Nloc,vin,vout)
    integer                                     :: Nloc
    complex(8)                                  :: vin(Nloc),vout(Nloc)
    vout = spHmf%dot(vin)
  end subroutine MatVec_Hmf
  subroutine cheb_rescaling(H,scale,shift)
    type(sparse_zmatrix_csr) :: H
    real(8),intent(out)      :: scale,shift
    complex(8)               :: C
    real(8)                  :: lower,upper,center,radius
    integer                  :: io
    lower =  huge(1d0)
    upper = -huge(1d0)
    do io=1,H%Nrow
      C      = H%get(io,io)
      center = dreal( C )
      radius = sum(abs(H%row(io)%vals))-abs(C)
      lower  = min(lower,center-radius)
      upper  = max(upper,center+radius)
    enddo
    shift = 0.5d0*(upper+lower)
    scale = max(0.5d0*(upper-lower),1d-12)/0.95d0
  end subroutine cheb_rescaling
  subroutine symmetrize_params(a)
    real(8),dimension(Nlat,2)                   :: a
    integer                                     :: ilat
    do ilat=2,Nlat
      a(ilat,:) = a(1,:)
    enddo
  end subroutine symmetrize_params



end program MF_bhz_2d_Anderson
