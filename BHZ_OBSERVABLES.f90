MODULE BHZ_OBSERVABLES
  USE SCIFOR
  USE DMFT_TOOLS
  USE COMMON
  USE BHZ_SOLVE
  USE BHZ_LCM
  implicit none
  private

  public :: assemble_Abhz_eigensystem
  public :: write_Abhz_density
  public :: write_Abhz_topology
  public :: write_Abhz_real_gf
  public :: write_Abhz_lcm

contains


  subroutine assemble_Abhz_eigensystem()
    ! Assemble the two globally stored spin-sector eigensystems in the
    ! complete site-major basis and store the ordered result in U,E.
    ! For SPIN_BLOCKS=1 the down eigensystem must first be reconstructed
    ! by the caller through time reversal.
    complex(8),allocatable :: raw_vectors(:,:)
    real(8),allocatable    :: raw_levels(:)
    integer,allocatable    :: order(:)
    integer                :: nblock,ib,ifull,ilat,iorb,istate
    !
    nblock = 2*Nlat
    if(.not.allocated(U_up).or..not.allocated(U_down)) &
      error stop "assemble_Abhz_eigensystem: missing eigenvectors"
    if(.not.allocated(Ev_up).or..not.allocated(Ev_down)) &
      error stop "assemble_Abhz_eigensystem: missing eigenvalues"
    if(size(Ev_up)/=nblock.or.size(Ev_down)/=nblock) &
      error stop "assemble_Abhz_eigensystem: invalid eigenvalue size"
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
  end subroutine assemble_Abhz_eigensystem




  subroutine write_Abhz_density(dens_up,dens_down)
    ! Write site-resolved orbital and spin densities. DENS_UP and
    ! DENS_DOWN contain the two orbital components for every site.
    real(8),dimension(:),intent(in) :: dens_up,dens_down
    real(8),allocatable             :: Nii(:,:,:),Tzii(:),Szii(:)
    integer                         :: ilat,iorb,io
    !
    if(size(dens_up)/=2*Nlat.or.size(dens_down)/=2*Nlat) &
      error stop "write_Abhz_density: invalid sector density size"
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
  end subroutine write_Abhz_density




  subroutine write_Abhz_topology()
    ! Prepare the occupied manifolds selected by SPIN_BLOCKS, compute
    ! their physical up/down Chern numbers and write spectrum and gap
    ! information. The complete U,E eigensystem is required only for
    ! SPIN_BLOCKS=0.
    real(8) :: chern(2)
    !
    if(method/="lapack")return
    select case(spin_blocks)
    case(1,2)
      chern = spin_chern_numbers()
      call save_sector_levels()
    case(0)
      if(.not.allocated(U).or..not.allocated(E)) &
        error stop "write_Abhz_topology: full eigensystem is not set"
      Nocc = count_occupied_levels(E,"full")
      Sz = kron(eye(Nlat),GammaS)
      call save_array("E.dat",E)
      call check_Egap("write_Abhz_topology")
      chern = spin_chern_numbers()
    case default
      error stop "write_Abhz_topology: invalid SPIN_BLOCKS"
    end select
    call save_array("z2.dat",0.5d0*(chern(1)-chern(2)))
    call save_array("spin_chern.dat",chern)
    write(*,*)"spin_Chern UP,DW:",chern
  end subroutine write_Abhz_topology




  subroutine write_Abhz_real_gf()
    ! Evaluate and write the local real-axis Green function from the
    ! globally stored complete eigensystem U,E.
    complex(8),allocatable :: Gf(:,:,:,:)
    !
    if(.not.with_real_gf)return
    if(.not.allocated(U).or..not.allocated(E)) &
      error stop "write_Abhz_real_gf: complete eigensystem is not set"
    allocate(Gf(Nlat,Nso,Nso,Lfreq))
    call get_gf(Gf,"real")
    call write_gf(gf_reshape(Gf,2,2,Nlat), &
      "Gloc_"//str(idum),"real",iprint=4,itar=.true.)
  end subroutine write_Abhz_real_gf




  subroutine write_Abhz_lcm()
    ! Compute and write physical up/down local Chern markers and their
    ! spin combination. The historical marker file contains the up map.
    real(8),allocatable :: marker_up(:,:),marker_down(:,:),marker_spin(:,:)
    character(len=3)    :: boundary
    !
    if(.not.with_lcm)return
    call local_spin_chern_marker(marker_up,marker_down)
    marker_spin = 0.5d0*(marker_up-marker_down)
    boundary = merge("PBC","OBC",bhz_pbc)
    call splot3d(boundary//"_Local_Chern_Marker_up.dat", &
      dble(arange(1,Nx)),dble(arange(1,Ny)),marker_up)
    call splot3d(boundary//"_Local_Chern_Marker_down.dat", &
      dble(arange(1,Nx)),dble(arange(1,Ny)),marker_down)
    call splot3d(boundary//"_Local_SpinChern_Marker.dat", &
      dble(arange(1,Nx)),dble(arange(1,Ny)),marker_up)
    call splot3d(boundary//"_Local_SpinChern_Combination.dat", &
      dble(arange(1,Nx)),dble(arange(1,Ny)),marker_spin)
  end subroutine write_Abhz_lcm




  subroutine save_sector_levels()
    ! Save the ordered union of the spin-sector spectra and its gap.
    real(8),allocatable :: levels(:)
    integer,allocatable :: order(:)
    real(8)             :: gap
    integer             :: nblock
    !
    nblock = size(Ev_up)
    allocate(levels(2*nblock),order(2*nblock))
    levels(1:nblock) = Ev_up
    if(spin_blocks==1)then
      levels(nblock+1:) = Ev_up
    else
      if(.not.allocated(Ev_down).or.size(Ev_down)/=nblock) &
        error stop "save_sector_levels: invalid down spectrum"
      levels(nblock+1:) = Ev_down
    endif
    call sort(levels,order)
    call save_array("E.dat",levels)
    if(Nocc<1.or.Nocc>=size(levels)) &
      error stop "save_sector_levels: invalid global occupation"
    gap = levels(Nocc+1)-levels(Nocc)
    write(*,"(A,G21.13)")"E gap=",gap
    call save_array("Egap.dat",gap)
  end subroutine save_sector_levels


END MODULE BHZ_OBSERVABLES
