program ed_bhz_2d_edge
  USE COMMON
  USE LCM_SQUARE
  USE EDIPACK2
  implicit none
  real(8)                                     :: z2,sp_chern(2)
  integer                                     :: iloop
  integer                                     :: ilat,iorb,jorb,ispin,jspin,io,jo,i,j
  logical                                     :: converged
  logical                                     :: spinsym
  logical                                     :: with_kinetic
  logical                                     :: with_z
  !Bath:
  integer                                     :: Nb
  real(8),allocatable,dimension(:,:)          :: Bath
  real(8),allocatable,dimension(:,:)          :: Bath_prev
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Self
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  complex(8),allocatable,dimension(:,:)         :: Gtest

  call init_parallel()

  call parse_cmd_variable(inputFILE,"inputFILE",default="inputABHZ.conf")
  call parse_input_variable(Nx,"Nx",inputFILE,default=10)
  call parse_input_variable(Wdis,"WDIS",inputFILE,default=0d0)
  call parse_input_variable(idum,"IDUM",inputFILE,default=1234567)
  call parse_input_variable(disorder_type,"DISORDER_TYPE",inputFILE,default=0)
  call parse_input_variable(bhz_pbc,"BHZ_PBC",inputFILE,default=.true.)
  call parse_input_variable(mh,"MH",inputFILE,default=1d0)
  call parse_input_variable(lambda,"LAMBDA",inputFILE,default=0.3d0)
  call parse_input_variable(spinsym,"SPINSYM",inputFILE,default=.true.)
  call parse_input_variable(wmix,"WMIX",inputFILE,default=0.5d0)
  call parse_input_variable(with_lcm,"WITH_lcm",inputFILE,default=.false.)
  call parse_input_variable(with_mats_gf,"WITH_MATS_GF",inputFILE,default=.false.)
  call parse_input_variable(with_real_gf,"WITH_REAL_GF",inputFILE,default=.false.)
  call parse_input_variable(with_kinetic,"WITH_KINETIC",inputFILE,default=.false.)
  call parse_input_variable(with_Z,"WITH_Z",inputFILE,default=.false.)
  call parse_input_variable(Nblock,"NBLOCK",inputFILE,default=4)
  !
  call ed_read_input(trim(inputFILE))


  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")


  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso=Nspin*Norb

  !SETUP COMMON DIMENSION
  Ny   = Nx
  Nk   = Nx*Ny
  Nlat = Nx*Ny
  !
  Nso  = Nspin*Norb
  Nlso = Nlat*Nso
  Nocc = Nlso/2

  !SETUP THE GAMMA MATRICES:
  call setup_GammaMatrices()


  !Allocate Functions:
  allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Weiss=zero
  allocate(Self(Nlat,Nspin,Nspin,Norb,Norb,Lmats));self=zero
  allocate(Gloc(Nlat,Nspin,Nspin,Norb,Norb,Lmats));gloc=zero
  allocate(Gtest(Nlat,Lmats));Gtest=zero


  !Set the basis vectors square lattice
  call TB_set_ei([1d0,0d0],[0d0,1d0])
  call TB_set_bk([pi2,0d0],[0d0,pi2])

  !SOLVE THE HOMOGENOUS PROBLEM:  
  if(MPImaster)write(*,*) "Solve homogeneous model with Nk="//str(Nk)
  allocate(Hk(Nso,Nso,Nk))
  call TB_build_model(Hk,hk_model,Nso,[Nx,Ny])
  z2 = hk_to_spin_Chern(Hk,[Nx,Ny],spin=1)
  if(MPImaster)write(*,*)"get spin Chern UP:",z2
  z2 = hk_to_spin_Chern(Hk,[Nx,Ny],spin=2)
  if(MPImaster)write(*,*)"get spin Chern DW:",z2
  z2 = hk_to_spin_Chern(Hk,[Nx,Ny])
  if(MPImaster)write(*,*)"get Z2:",z2
  if(MPImaster)write(*,*)""
  if(MPImaster)write(*,*)""
  deallocate(Hk)




  !< Build up disorder:
  call setup_Abhz()
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
  do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
     i = iorb + (ispin-1)*Norb + (ilat-1)*Nso
     j = jorb + (jspin-1)*Norb + (ilat-1)*Nso
     Hloc(ilat,ispin,jspin,iorb,jorb) = Hij(i,j,1)
  enddo



  !pass the total lattice H to get the local part 
  call ed_set_hloc(Hloc,Nlat)




  !=========================================  
  !POST-PROCESSING:
  !=========================================
  !Get K
  if(with_kinetic)call get_Ekin
  !
  !Get local Gf matsubara
  if(with_mats_gf)call get_gf_mats
  !
  !Get local Gf real-axis  
  if(with_real_gf)call get_gf_real
  !=========================================


  Nb=ed_get_bath_dimension()

  allocate(Bath(Nlat,Nb) )
  allocate(Bath_prev(Nlat,Nb) )
  call ed_init_solver(Bath)



  !DMFT loop:
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     !
     call ed_solve(Bath,mpi_lanc=.false.)
     call ed_get_sigma(self,Nlat,'m')
     !
     ! compute the local gf:
     call dmft_get_gloc(Hij,gloc,self,axis='m')
     !
     ! compute the Weiss field (only the Nineq ones)
     call dmft_self_consistency(gloc,self,Weiss)
     !
     ! fit baths and mix result with old baths
     call ed_chi2_fitgf(Weiss,Bath,ispin=1)
     !
     !
     !Post - fit and check error:
     call ed_spin_symmetrize_bath(Bath,save=.true.)
     if(iloop>1)Bath=wmix*Bath + (1.d0-wmix)*Bath_prev
     Bath_prev= Bath
     Gtest    = zero
     do ilat=1,Nlat
        do iorb=1,Norb
           i = iorb +  (ilat-1)*Nso 
           Gtest(ilat,:) = Gtest(ilat,:) + Weiss(ilat,1,1,iorb,iorb,:)/Norb
        enddo
     enddo
     converged = check_convergence(Gtest,dmft_error,nsuccess,nloop)
     !
     !
     call end_loop
  enddo

  !Reduce clutter
  call reduce_clutter_dmft


  !Push Sigma to topological Hamiltonian
  call push_Htop(Self)


  !Get topological info:
  sp_chern(1) = single_point_spin_chern(spin=1)
  sp_chern(2) = single_point_spin_chern(spin=2)
  if(MPImaster)call save_array("z2.dat",(sp_chern(1)-sp_chern(2))/2d0 )
  if(MPImaster)call save_array("spin_chern.dat",sp_chern)
  if(MPImaster)print*,"spin_Chern UP,DW:",sp_chern
  !
  !Get LCM: 
  if(with_lcm)then
     call pbc_local_spin_chern_marker(spin=1,lcm=LsCM)
     if(MPImaster)call splot3d("PBC_Local_SpinChern_Marker.dat",dble(arange(1,Nx)),dble(arange(1,Ny)),LsCM)
  endif


  call finalize_MPI()



contains







  subroutine push_Htop(Smats) 
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,Lmats) :: Smats
    complex(8),dimension(Nlso,Nlso)                        :: S0
    real(8),dimension(Nlso,Nlso)                           :: Z
    complex(8),dimension(Nlat,Nso,Nso)                     :: Zfoo
    real(8),dimension(Nlso)                                :: Ev
    complex(8),dimension(Nlso,Nlso)                        :: H
    integer                                                :: ilat,io,i
    !
    S0 = reshape_5to2_array(Smats(:,:,:,:,:,1))
    H  = Hij(:,:,1) + one*dreal(S0)
    !
    if(with_z)then
       do ilat=1,Nlat
          Zfoo(ilat,:,:) = select_block(ilat,S0)
          do io=1,Nso
             i = io + (ilat-1)*Nso
             Z(i,i)  = 1.d0/( 1.d0 + abs( dimag(Zfoo(ilat,io,io))/(pi/beta) ))
          enddo
       enddo
       H = matmul(sqrt(Z),matmul(H,Z))
    endif
    !
#ifdef _SCALAPACK
    if(MpiSize==1)then
       call eigh(H,Ev)
    else
       call p_eigh(H,Ev,Nblock)
    endif
#else
    call eigh(H,Ev)
#endif
    !
    call push_Bloch(H,Ev)
    !
  end subroutine push_Htop


  subroutine get_gf_mats
    !S --> Sigma(iw) real-axis:
    if(allocated(self))deallocate(self)
    if(allocated(gloc))deallocate(gloc)
    allocate(self(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(gloc(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    call ed_get_sigma(self,Nlat,'m')
    call dmft_get_gloc(Hij,gloc,self,axis='m')
    call write_gf(Gloc,"Gloc_"//str(idum),'mats',iprint=4,itar=.true.)
    call finalize_MPI()
    stop
  end subroutine get_gf_mats



  subroutine get_gf_real
    !S --> Sigma(w) real-axis:
    if(allocated(self))deallocate(self)
    if(allocated(gloc))deallocate(gloc)
    allocate(self(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(gloc(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    call ed_get_sigma(self,Nlat,'r')
    call dmft_get_gloc(Hij,gloc,self,axis='r')
    call write_gf(gloc,"Gloc_"//str(idum),'real',iprint=4,itar=.true.)
    call finalize_MPI()
    stop
  end subroutine get_gf_real



  subroutine get_Ekin
    if(allocated(self))deallocate(self)
    allocate(self(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    call ed_get_sigma(self,Nlat,'m')
    call dmft_kinetic_energy(Hij,self)
    call finalize_MPI()
    stop
  end subroutine get_Ekin


  

  function select_block(ip,Matrix) result(Vblock)
    integer                                          :: ip
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: Matrix
    complex(8),dimension(Nspin*Norb,Nspin*Norb)      :: Vblock
    integer                                          :: is,js,ispin,jspin,iorb,jorb
    Vblock=zero
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb !spin-orbit stride
       js = jorb + (jspin-1)*Norb !spin-orbit stride
       Vblock(is,js) = Matrix(ip,ispin,jspin,iorb,jorb)
    enddo
  end function select_block


  function reshape_5to2_array(fg) result(g)
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: fg
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: g
    integer                                               :: i,j,iorb,jorb,ispin,jspin
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       i = iorb + (ispin-1)*Norb + (ilat-1)*Nso
       j = jorb + (jspin-1)*Norb + (ilat-1)*Nso
       g(i,j) = fg(ilat,ispin,jspin,iorb,jorb)
    enddo
  end function reshape_5to2_array


  function reshape_2to5_array(fg) result(g)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: fg
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: g
    integer                                               :: i,j,iorb,jorb,ispin,jspin
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       i = iorb + (ispin-1)*Norb + (ilat-1)*Nso
       j = jorb + (jspin-1)*Norb + (ilat-1)*Nso
       g(ilat,ispin,jspin,iorb,jorb) = fg(i,j)
    enddo
  end function reshape_2to5_array


  function reshape_5to2_gf(fg,L) result(g)
    integer                                                 :: L
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,L)      :: fg
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L) :: g
    integer                                                 :: i,j,iorb,jorb,ispin,jspin
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       i = iorb + (ispin-1)*Norb + (ilat-1)*Nso
       j = jorb + (jspin-1)*Norb + (ilat-1)*Nso
       g(i,j,:) = fg(ilat,ispin,jspin,iorb,jorb,:)
    enddo
  end function reshape_5to2_gf


  function reshape_2to5_gf(fg,L) result(g)
    integer                                                 :: L
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L) :: fg
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,L)      :: g
    integer                                                 :: i,j,iorb,jorb,ispin,jspin
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       i = iorb + (ispin-1)*Norb + (ilat-1)*Nso
       j = jorb + (jspin-1)*Norb + (ilat-1)*Nso
       g(ilat,ispin,jspin,iorb,jorb,:) = fg(i,j,:)
    enddo
  end function reshape_2to5_gf




  subroutine reduce_clutter_dmft
    if(MpiMaster)then
       ! call file_targz(tarball="tar_gfmatrix",pattern="gfmatrix_ineq*.restart")
       call file_targz(tarball="tar_N2_correlation",pattern="N2_*.ed")
       call file_targz(tarball="tar_Sz2_correlation",pattern="Sz2_*.ed")
       call file_targz(tarball="tar_N2_correlation",pattern="N2_*.ed")
       call file_targz(tarball="tar_Z_imp",pattern="Z_*.ed")
       call file_targz(tarball="tar_Sig_imp",pattern="Sig_*.ed")
       call file_targz(tarball="tar_RDM_imp",pattern="reduced_density_matrix_*.ed")
       call file_targz(tarball="tar_energy_imp",pattern="energy_*.ed")
       call file_targz(tarball="tar_fit_weiss",pattern="fit_weiss*")
       call file_targz(tarball="tar_eigenvalues_list",pattern="eigenvalues_list*")
       call file_targz(tarball="tar_state_list",pattern="state_list*")
       call file_targz(tarball="tar_chi2fit",pattern="chi2fit_results*")       
    endif
  end subroutine reduce_clutter_dmft



end program ed_bhz_2d_edge
