MODULE BHZ_OBSERVABLES
  USE COMMON
  USE BHZ_BUILD, only: Sz
  USE BHZ_SOLVE
  USE BHZ_LCM
  implicit none
  private


  public :: write_Abhz_topology
  public :: write_Abhz_lcm

contains





  ! Compute and write the global real-space Chern diagnostics for the
  ! selected solver. Lapack uses the stored eigensystem; Chebychev uses
  ! caller-owned CSR matrices and can return its PBC maps for reuse.
  ! Exact levels and the spectral gap are written only when available.
  subroutine write_Abhz_topology(spHup,spHdown,spHfull, &
    marker_up,marker_down)
    type(sparse_zmatrix_csr),optional,intent(in) :: spHup,spHdown,spHfull
    real(8),allocatable,optional,intent(out)     :: marker_up(:,:)
    real(8),allocatable,optional,intent(out)     :: marker_down(:,:)
    real(8)                                     :: chern(2),spin_chern,sem(4)
    real(8),allocatable                         :: local_chern(:)
    real(8),allocatable                         :: local_scan(:,:)
    integer                                     :: unit,i,runners(2),target_met
    !
    if(method=="lapack".and.spin_blocks==0)then
      if(.not.allocated(U))error stop "write_Abhz_topology: full eigenvectors are not set"
      if(.not.allocated(E))error stop "write_Abhz_topology: full eigenenergies are not set"
      Nocc = count_occupied_levels(E,"full")
      if(topology_method=="single_point".or. &
        (with_lcm.and.topology_method/="localizer")) &
        Sz = kron(eye(Nlat),GammaS)
      call save_array("E.dat",E)
      call check_Egap("write_Abhz_topology")
    endif
    !
    chern = LCM_SpinChern_numbers(spHup,spHdown,spHfull, &
      marker_up,marker_down,sem,local_chern,local_scan,runners)
    if(allocated(local_chern))then
      open(newunit=unit,file="localizer_eta_scan.dat",status="replace")
      write(unit,'(a)') "# 1.eta_over_eta_auto 2.eta 3.C_up 4.C_down"
      do i=1,size(local_scan,2)
        write(unit,'(4(es24.16,1x))') local_scan(:,i)
      enddo
      close(unit)
      call save_array("localizer_spin_chern.dat",local_chern)
      call save_array("localizer_chern_total.dat",sum(local_chern))
      call save_array("localizer_spin_chern_combination.dat", &
        0.5d0*(local_chern(1)-local_chern(2)))
      write(*,*) "Localizer Chern UP,DW:",local_chern
    else
      write(*,*) "Periodic localizer skipped: BHZ_PBC=F"
    endif
    !
    if(method=="lapack".and.spin_blocks/=0)call save_sector_levels()
    !
    spin_chern = 0.5d0*(chern(1)-chern(2))
    !
    !Save the physical up/down local Chern-marker maps and their half-difference
    call save_array("chern_total.dat",sum(chern))
    call save_array("spin_chern_combination.dat",spin_chern)
    call save_array("z2.dat",spin_chern)
    call save_array("spin_chern.dat",chern)
    if(topology_method=="stochastic")then
      open(newunit=unit,file="chern_sampling_sem.dat",status="replace")
      write(unit,'(a)') "# 1.runners 2.up_sem 3.down_sem 4.total_sem 5.spin_sem"
      write(unit,'(i0,4(1x,es24.16))') maxval(runners),sem
      close(unit)
      target_met = -1
      if(topology_sem_tol>0d0)target_met = merge(1,0,sem(4)<=topology_sem_tol)
      open(newunit=unit,file="chern_sampling_status.dat",status="replace")
      write(unit,'(a)') "# 1.up_runners 2.down_runners 3.spin_sem_target 4.target_met"
      write(unit,'(i0,1x,i0,1x,es24.16,1x,i0)') &
        runners,topology_sem_tol,target_met
      close(unit)
      write(*,*) "Chern sampling SEM (up, down, total, spin):",sem
      write(*,*) "Chern sampling runners (up, down), target_met:", &
        runners,target_met
    endif
    !
    write(*,*)"spin_Chern UP,DW:",chern
  end subroutine write_Abhz_topology





  ! Compute and write the physical up/down local Chern-marker maps and
  ! their half-difference. Previously returned Chebychev PBC maps are
  ! reused when present; otherwise the selected LCM kernel is invoked.
  ! The optional CSR arguments are needed only by the Chebychev path.
  subroutine write_Abhz_lcm(spHup,spHdown,spHfull,marker_up,marker_down)
    type(sparse_zmatrix_csr),optional,intent(in) :: spHup,spHdown,spHfull
    real(8),allocatable,intent(inout)            :: marker_up(:,:),marker_down(:,:)
    real(8),allocatable                          :: marker_spin(:,:)
    character(len=3)                             :: boundary
    !
    if(.not.with_lcm.or.topology_method=="localizer")return
    !
    if(.not.allocated(marker_up).or..not.allocated(marker_down)) &
      call LCM_SpinChern_marker(marker_up,marker_down,spHup,spHdown,spHfull)
    !
    marker_spin = 0.5d0*(marker_up-marker_down)
    boundary = merge("PBC","OBC",bhz_pbc)
    call splot3d(boundary//"_Local_Chern_Marker_up.dat",dble(arange(1,Nx)),dble(arange(1,Ny)),marker_up)
    call splot3d(boundary//"_Local_Chern_Marker_down.dat",dble(arange(1,Nx)),dble(arange(1,Ny)),marker_down)
    call splot3d(boundary//"_Local_SpinChern_Marker.dat",dble(arange(1,Nx)),dble(arange(1,Ny)),marker_up)
    call splot3d(boundary//"_Local_SpinChern_Combination.dat",dble(arange(1,Nx)),dble(arange(1,Ny)),marker_spin)
  end subroutine write_Abhz_lcm




  ! Merge and sort the two dense spin-sector spectra, then write the
  ! complete level list and the single-particle gap across global Nocc.
  ! In SPIN_BLOCKS=1 the time-reversed spectrum equals the up spectrum.
  subroutine save_sector_levels()
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
    if(topology_method/="single_point") &
      Nocc = count_occupied_levels(levels,"assembled")
    if(Nocc<1.or.Nocc>=size(levels)) &
      error stop "save_sector_levels: invalid global occupation"
    gap = levels(Nocc+1)-levels(Nocc)
    write(*,"(A,G21.13)")"E gap=",gap
    call save_array("Egap.dat",gap)
  end subroutine save_sector_levels


end module BHZ_OBSERVABLES
