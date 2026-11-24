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
  USE BHZ_SOLVE
  USE BHZ_LCM
  USE BHZ_OBSERVABLES
  implicit none
  integer,parameter                         :: Norb=2,Nspin=2
  real(8)                                   :: z2
  !Solution
  real(8),dimension(:),allocatable          :: dens_up,dens_down
  real(8),dimension(:),allocatable          :: dens_full
  integer                                   :: ilat,iorb
  type(sparse_zmatrix_csr),target           :: spHup,spHdown
  type(sparse_zmatrix_csr),target           :: spHfull
  !
  !Read input:
  call parse_input_file(inputFILE,"inputABHZ.conf")
  call read_input(inputFILE,"Abhz")

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
    call solve_Abhz(spHup,dens_up,U_up,Ev_up)
    allocate(dens_down,source=dens_up)
    if(method=="lapack".and.with_real_gf)then
      allocate(Ev_down,source=Ev_up)
      allocate(U_down,source=conjg(U_up))
      call assemble_Abhz_full(spHup,spHdown,spHfull)
      call assemble_Abhz_eigensystem()
    endif
  case(2)
    call setup_Abhz(spH=spHup,spin=1)
    call setup_Abhz(verbose=.false.,spH=spHdown,spin=-1)
    call solve_Abhz(spHup,dens_up,U_up,Ev_up)
    call solve_Abhz(spHdown,dens_down,U_down,Ev_down)
    if(method=="lapack".and.with_real_gf)then
      call assemble_Abhz_full(spHup,spHdown,spHfull)
      call assemble_Abhz_eigensystem()
    endif
  case(0)
    call setup_Abhz(spH=spHfull,spin=0)
    call solve_Abhz(spHfull,dens_full,U,E)
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
  call write_Abhz_density(dens_up,dens_down)
  !
  if(method=="lapack")then
    call write_Abhz_topology()
    call write_Abhz_real_gf()
    call write_Abhz_lcm()
  endif

  !< Free memory:
  call free_seed_solution()
  if(spHup%status)call spHup%free()
  if(spHdown%status)call spHdown%free()
  if(spHfull%status)call spHfull%free()
  call free_Abhz_solution()
  call free_Abhz()
  call reset_seed_wd()
  call stop_timer("IDUM: "//str(idum))
  write(*,*)""



contains




  subroutine free_seed_solution()
    ! Release arrays owned locally by one disorder realization. The
    ! module eigensystems are released by free_Abhz_solution.
    !
    if(allocated(dens_up))deallocate(dens_up)
    if(allocated(dens_down))deallocate(dens_down)
    if(allocated(dens_full))deallocate(dens_full)
  end subroutine free_seed_solution


end program anderson_bhz_2d
