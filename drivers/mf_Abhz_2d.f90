!> Solve the 2D BHZ + Kanamori n.n + Disorder model using Mean-Field
program MF_bhz_2d_Anderson
  USE COMMON
  USE BHZ_BUILD
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


  call parse_input_file(inputFILE,"inputBHZ.conf")
  call read_input(inputFILE)

  !SETUP THE GAMMA MATRICES:

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
    call setup_Abhz(spH=spHij_down,spin=-1)
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
  if(spHij_up%status)call spHij_up%free()
  if(spHij_down%status)call spHij_down%free()
  if(spHij_full%status)call spHij_full%free()
  if(spHmf%status)call spHmf%free()
  call free_Abhz_EigenSolution()
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
      call get_Abhz_EigenSolution(spHmf,dens_up)
      allocate(dens_down,source=dens_up)
    case(2)
      call set_mf_hamiltonian(spHij_up,a,spin=1)
      call get_Abhz_EigenSolution(spHmf,dens_up)
      call set_mf_hamiltonian(spHij_down,a,spin=-1)
      call get_Abhz_EigenSolution(spHmf,dens_down)
    case(0)
      call set_mf_hamiltonian(spHij_full,a,spin=0)
      call get_Abhz_EigenSolution(spHmf,dens_full)
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
      a(ilat,1) = 0.5d0*sum(dens(ilat,:,1)) - 0.5d0*sum(dens(ilat,:,2))
      a(ilat,2) = 0.5d0*sum(dens(ilat,1,:)) - 0.5d0*sum(dens(ilat,2,:))
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
    spHmf = build_mf_Abhz(spH0,a,spin)
  end subroutine set_mf_hamiltonian





  ! Mean-field construction is shared with spectral post-processing
  ! through BHZ_BUILD:build_mf_Abhz.




  ! Local sparse blocks are inserted inside build_mf_Abhz.







  subroutine symmetrize_params(a)
    real(8),dimension(Nlat,2)                   :: a
    integer                                     :: ilat
    do ilat=2,Nlat
      a(ilat,:) = a(1,:)
    enddo
  end subroutine symmetrize_params



end program MF_bhz_2d_Anderson
