subroutine mf_anderson_bhz_2d(Uloc,M,res)
  USE COMMON
  USE LCM_SQUARE
  implicit none
  real(8),intent(in)                        :: uloc,M
  integer,intent(out)                       :: res
  integer,parameter                         :: Norb=2,Nspin=2
  real(8)                                   :: JhRatio,Jh
  real(8)                                   :: z2,sp_chern(2),Esz,Etz,Ez2,threshold
  !Solution
  real(8),dimension(:),allocatable          :: Ev
  complex(8),dimension(:,:),allocatable     :: H
  real(8),dimension(:,:,:),allocatable      :: Nii
  real(8),dimension(:),allocatable          :: Tzii,Szii,N1ii,N2ii
  complex(8),dimension(:,:,:,:),allocatable :: Gf
  integer                                   :: ilat,iorb,ispin,io,ii,id
  logical                                   :: converged,iexist,ioidum,post_processing
  integer                                   :: Iter,Nsuccess=2
  real(8),dimension(:,:),allocatable        :: params,params_prev
  real(8),dimension(:),allocatable          :: Z2data,Tzdata,Szdata
  integer,dimension(:),allocatable :: AFMv

  call init_parallel()
  !  
  !Read input:
  call parse_cmd_variable(inputFILE,"inputFILE",default="inputABHZ.conf")
  call parse_cmd_variable(idumFILE,"idumFILE",default="list_idum")
  call parse_input_variable(Nx,"Nx",inputFILE,default=10)
  call parse_input_variable(Wdis,"WDIS",inputFILE,default=0d0)
  call parse_input_variable(idum,"IDUM",inputFILE,default=1234567)
  call parse_input_variable(disorder_type,"DISORDER_TYPE",inputFILE,default=0)
  call parse_input_variable(bhz_pbc,"BHZ_PBC",inputFILE,default=.true.)
  call parse_input_variable(lambda,"LAMBDA",inputFILE,default=0.3d0)
  call parse_input_variable(Jhratio,"JhRatio",inputFILE,default=0.25d0)
  call parse_input_variable(mu,"MU",inputFILE,default=0.d0)
  call parse_input_variable(temp,"temp",inputFILE,default=0.0001d0)
  call parse_input_variable(wmix,"WMIX",inputFILE,default=0.5d0)
  call parse_input_variable(p_field,"P_FIELD",inputFILE,default=0.01d0)
  call parse_input_variable(it_error,"IT_ERROR",inputFILE,default=1d-5)
  call parse_input_variable(maxiter,"MAXITER",inputFILE,default=100)
  call parse_input_variable(Nblock,"NBLOCK",inputFILE,default=4)
  if(MPImaster)call save_input_file(inputFILE)
  !
  Mh = M
  Jh = JhRatio*Uloc
  !
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


  !Set the basis vectors square lattice
  call TB_set_ei([1d0,0d0],[0d0,1d0])
  call TB_set_bk([pi2,0d0],[0d0,pi2])

  !SOLVE THE HOMOGENOUS PROBLEM:  
  if(MPImaster)write(*,*) "Solve homogeneous model with Nk="//str(Nk)
  allocate(Hk(Nso,Nso,Nk))
  call TB_build_model(Hk,hk_model,Nso,[Nx,Ny])
  z2 = hk_to_spin_Chern(Hk,[Nx,Ny])

  allocate(AFMv(Nlat))
  do ilat=1,Nlat
     AFMv(ilat) = afm_sign(ilat)
  enddo
  !
  call read_list_idum(str(idumFILE),Nidum)

  call getcwd(here)


  !####################################
  do ii=1,Nidum
     idum  = list_idum(ii)
     dir   = "IDUM_"//str(idum)
     !
     !Create IDUM directory and ENTER it     
     if(MpiMaster)call system("mkdir -p "//str(dir))
     call chdir(str(dir))
     !
     !< Build up disorder:
     call setup_Abhz()
     !
     !>MF Solution:
     !read params
     allocate(params(Nlat,2),params_prev(Nlat,2))
     params_prev=0d0
     do ilat=1,Nlat
        params(ilat,:)= [p_field,p_field*afm_sign(ilat)]
     enddo
     !
     !> Mean-Field cycle with linear mixing
     converged=.false. ; iter=0
     do while(.not.converged.AND.iter<maxiter)
        iter=iter+1
        call solve_MF_Abhz(iter,params)           
        if(iter>1)params = wmix*params + (1d0-wmix)*params_prev
        params_prev = params
        converged = check_convergence_local(params,it_error,1,maxiter) 
     end do
     if(MPImaster)call save_array("params.restart",params)
     !
     !Get topological info:
     call push_Bloch(H,Ev)
     sp_chern(1) = single_point_spin_chern(spin=1)
     sp_chern(2) = single_point_spin_chern(spin=2)
     z2          = (sp_chern(1)-sp_chern(2))/2d0
     if(MPImaster)then
        write(*,*)iter,get_mean(Tzii),get_AFMmean(Szii),z2
        call save_array("z2.dat", z2)
        call save_array("Ebhz.dat",Ev)
        call save_array("tz.dat",Tzii)
        call save_array("sz.dat",Szii)

        call add_to(Tzdata,Tzii)
        call add_to(Szdata,Szii*Afmv)
        call add_to(Z2data,[Z2])
     endif
     !

     !< Free memory:
     call chdir(str(here))
     call free_Abhz()
     deallocate(params,params_prev)
  enddo
  !####################################


  !P = [E(Tz),E(Sz),E(Z2)]
  ![!=0,0,1]: Res=0 Topological phase (Z2=1, Sz=0, Tz>0)
  ![!=0,0,0]: Res=1 Trivial phase (Z2=0, Sz=0, Tz>0)
  ![!=0,!=0,0]: Res=4  phase 4 (Z2=0, Sz=1, Tz!=0)
  ![!=0,!=0,1]: Res=5  phase 5 (Z2=1, Sz=1, Tz!=0)
  !
  ![0,!=0,0]: Res=2 AFM phase (Z2=0, Sz=1, Tz=0)
  ![0,!=0,1]: Res=3 Top-AFM phase (Z2=1, Sz=1, Tz=0)
  ![0,0,0]: Res=6 AFM phase (Z2=0, Sz=1, Tz=0)
  ![0,0,1]: Res=7 Top-AFM phase (Z2=1, Sz=1, Tz=0)
  threshold=1d-2
  Ez2 = get_mean(Z2data)
  Etz = abs(get_mean(Tzdata))
  Esz = get_mean(Szdata)
  if(Etz>=threshold)then         !Tz > 0
     if(Esz>=threshold)then      ! Sz > 0
        if(Ez2>=threshold)then   !  Z2 = 1
           Res=5
        else
           Res=4                !  Z2 = 0
        endif
     else                       ! Sz = 0
        if(Ez2>=threshold)then   !  Z2 = 1 TI
           Res=0
        else 
           Res=1                !  Z2 = 0 BI
        endif
     endif
  else                          !Tz = 0
     if(Esz>=threshold)then      ! Sz > 0
        if(Ez2>=threshold)then   !  Z2 = 1 AFM-TI
           Res=3
        else
           Res=2                !  Z2 = 0 AFM
        endif
     else                       ! Sz = 0
        if(Ez2>=threshold)then   !  Z2 = 1 [0,0,1] 
           Res=6
        else 
           Res=7                !  Z2 = 0 [0,0,0]
        endif
     endif
  endif


  call end_parallel()


contains



  subroutine solve_MF_Abhz(iter,a)
    integer                                 :: iter
    real(8),dimension(Nlat,2),intent(inout) :: a
    complex(8),dimension(Nlso,NLso)         :: fE
    complex(8),dimension(Nlso,NLso)         :: Rho
    integer                                 :: io,jo,ilat,iorb,ispin
    !
    !
    if(.not.allocated(Nii))allocate(Nii(Nlat,Nspin,Norb))
    if(.not.allocated(Szii))allocate(Szii(Nlat))
    if(.not.allocated(Tzii))allocate(Tzii(Nlat))
    if(.not.allocated(N1ii))allocate(N1ii(Nlat))
    if(.not.allocated(N2ii))allocate(N2ii(Nlat))
    if(.not.allocated(H))allocate(H(Nlso,Nlso))
    if(.not.allocated(Ev))allocate(Ev(Nlso))
    !
    H = Hij(:,:,1) + mf_Hij_correction(a)
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
    !Actual solution:
    fE  = one*diag(fermi(Ev,1d0/temp))
    Rho = ( H.mx.fE ).mx.conjg(transpose(H))
    ! Rho = matmul(H , matmul(fE, conjg(transpose(H))) )
    !
    !Get observables:
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,iorb=1:Norb)
       io = iorb+(ispin-1)*Norb+(ilat-1)*Nspin*Norb
       Nii(ilat,ispin,iorb) = Rho(io,io)
    enddo
    do ilat=1,Nlat
       N1ii(ilat) = sum(Nii(ilat,:,1))
       N2ii(ilat) = sum(Nii(ilat,:,2))
       Tzii(ilat) = 0.5d0*(N1ii(ilat) - N2ii(ilat))
       Szii(ilat) = 0.5d0*sum(Nii(ilat,1,:)) - 0.5d0*sum(Nii(ilat,2,:)) !N_up - N_dw
    enddo
    a(:,1) = Tzii
    a(:,2) = Szii
  end subroutine solve_MF_Abhz


  function mf_Hij_correction(a) result(HijMF)
    real(8),dimension(Nlat,2)               :: a
    complex(8),dimension(Nlat,Nso,Nso)      :: Htmp
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: HijMF
    integer                                 :: ilat,io,jo,i,j
    HijMF=zero
    do ilat=1,Nlat
       Htmp(ilat,:,:) = -a(ilat,1)*(Uloc-5d0*Jh)/2d0*Gamma5 -a(ilat,2)*(Uloc+Jh)/2d0*GammaS
       do io=1,Nso
          do jo=1,Nso
             i = io + (ilat-1)*Nso
             j = jo + (ilat-1)*Nso
             HijMF(i,j) = Htmp(ilat,io,jo)
          enddo
       enddo
    enddo
  end function mf_Hij_correction





  function get_Hmf(a,Hij) result(Hmat)
    real(8),dimension(Nlat,2)                :: a
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Hij
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Hmat
    Hmat = Hij + mf_Hij_correction(a)
  end function get_Hmf



  subroutine symmetrize_params(a)
    real(8),dimension(Nlat,Nso) :: a
    integer                     :: ilat
    do ilat=2,Nlat
       a(ilat,:) = a(1,:)
    enddo
  end subroutine symmetrize_params



END subroutine mf_anderson_bhz_2d






