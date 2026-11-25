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
  USE BHZ_BUILD
  USE BHZ_SOLVE
  USE BHZ_LCM
  USE BHZ_OBSERVABLES
  implicit none
  integer,parameter                :: Norb=2,Nspin=2
  real(8)                          :: z2
  !Solution
  real(8),dimension(:),allocatable :: dens_up,dens_down
  real(8),dimension(:),allocatable :: dens_full
  real(8),allocatable               :: marker_up(:,:),marker_down(:,:)
  integer                          :: ilat,iorb,io
  type(sparse_zmatrix_csr),target  :: spHup,spHdown
  type(sparse_zmatrix_csr),target  :: spHfull
  real(8),allocatable              :: Nii(:,:,:),Tzii(:),Szii(:)

  !
  !Read input:
  call parse_input_file(inputFILE,"ABHZ.conf")
  call read_input(inputFILE)


  !Set the basis vectors square lattice
  call TB_set_ei([1d0,0d0],[0d0,1d0])
  call TB_set_bk([pi2,0d0],[0d0,pi2])



  !SOLVE THE HOMOGENOUS PROBLEM:
  write(*,*) "Solve homogeneous model with Nk="//str(Nk)
  allocate(Hk(Nso,Nso,Nk))
  call TB_build_model(Hk,hk_model,Nso,[Nx,Ny])
  z2 = Bloch_SpinChern(Hk,[Nx,Ny])
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
    call get_Abhz_EigenSolution(spHup,dens_up,U_up,Ev_up)
    allocate(dens_down,source=dens_up)
   case(2)
    call setup_Abhz(spH=spHup,spin=1)
    call setup_Abhz(spH=spHdown,spin=-1,verbose=.false.)
    call get_Abhz_EigenSolution(spHup,dens_up,U_up,Ev_up)
    call get_Abhz_EigenSolution(spHdown,dens_down,U_down,Ev_down)
   case(0)
    call setup_Abhz(spH=spHfull,spin=0)
    call get_Abhz_EigenSolution(spHfull,dens_full,U,E)
    allocate(dens_up(2*Nlat),dens_down(2*Nlat))
    do ilat=1,Nlat
      do iorb=1,Norb
        dens_up(iorb+(ilat-1)*Norb)   = dens_full(iorb+(ilat-1)*Nso)
        dens_down(iorb+(ilat-1)*Norb) = dens_full(iorb+Norb+(ilat-1)*Nso)
      enddo
    enddo
  end select
  !
  allocate(Nii(Nlat,2,2))
  do ilat=1,Nlat
    do iorb=1,2
      io = iorb+(ilat-1)*2
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
  !
  !
  write(*,*)"Get Abhz Chern"
  call write_Abhz_topology(spHup,spHdown,spHfull,marker_up,marker_down)
  write(*,*)"Get Abhz LCM"
  call write_Abhz_lcm(spHup,spHdown,spHfull,marker_up,marker_down)

  !< Free memory:
  if(allocated(dens_up))deallocate(dens_up)
  if(allocated(dens_down))deallocate(dens_down)
  if(allocated(dens_full))deallocate(dens_full)
  if(spHup%status)call spHup%free()
  if(spHdown%status)call spHdown%free()
  if(spHfull%status)call spHfull%free()
  call free_Abhz_EigenSolution()
  call free_Abhz()
  !
  call reset_seed_wd()
  call stop_timer("IDUM: "//str(idum))
  write(*,*)""


end program anderson_bhz_2d
