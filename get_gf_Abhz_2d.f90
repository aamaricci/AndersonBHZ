program analysis
  USE COMMON
  USE LCM_SQUARE
#ifdef _SCALAPACK
  USE MPI
#endif
  implicit none



  integer,parameter                   :: Norb=2,Nspin=2
  real(8)                             :: bT
  !Solution
  real(8),dimension(:),allocatable    :: w
  real(8),dimension(:,:),allocatable  :: Aw
  real(8),dimension(:),allocatable    :: Egf,Egf2,Ggf,Sgf
  integer                             :: ilat,iorb,ispin,i,id,io
  complex(8),dimension(:),allocatable :: Gi
  character(len=64)                   :: suffix,iostr

  call init_parallel()  
  !  
  !Read input:
  call parse_cmd_variable(inputFILE,"inputFILE",default="used.inputABHZ.conf")
  call parse_cmd_variable(idumFILE,"idumFILE",default="list_idum")
  call parse_input_variable(Nx,"Nx",inputFILE,default=10)
  call parse_input_variable(Wdis,"WDIS",inputFILE,default=0d0)
  call parse_input_variable(idum,"IDUM",inputFILE,default=1234567)
  call parse_input_variable(disorder_type,"DISORDER_TYPE",inputFILE,default=0)
  call parse_input_variable(bhz_pbc,"BHZ_PBC",inputFILE,default=.true.)
  call parse_input_variable(mh,"MH",inputFILE,default=3.d0)
  call parse_input_variable(lambda,"LAMBDA",inputFILE,default=0.3d0)  
  call parse_input_variable(mu,"MU",inputFILE,default=0.d0)
  call parse_input_variable(temp,"temp",inputFILE,default=0.0001d0)
  call parse_input_variable(wmix,"WMIX",inputFILE,default=0.5d0)
  call parse_input_variable(Lfreq,"Lfreq",inputFILE,default=1024)
  call parse_input_variable(wmin,"WMIN",inputFILE,default=-5d0)
  call parse_input_variable(wmax,"wmax",inputFILE,default= 5d0)
  call parse_input_variable(ieta,"iETA",inputFILE,default=4.d-2)
  call parse_input_variable(p_field,"P_FIELD",inputFILE,default=0.01d0)
  call parse_input_variable(it_error,"IT_ERROR",inputFILE,default=1d-5)
  call parse_input_variable(maxiter,"MAXITER",inputFILE,default=100)
  call parse_input_variable(with_lcm,"WITH_lcm",inputFILE,default=.false.)
  call parse_input_variable(with_mats_gf,"WITH_MATS_GF",inputFILE,default=.false.)
  call parse_input_variable(with_real_gf,"WITH_REAL_GF",inputFILE,default=.false.)
  call parse_input_variable(Nblock,"NBLOCK",inputFILE,default=4)
  !
  !Save variables into DMFT_TOOLS memory pool
  bT = 1d0/temp
  call add_ctrl_var(bT,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(mu,"xmu")
  call add_ctrl_var(wmin,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(ieta,"eps")



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




  Nidum = file_length(str(idumFILE))
  allocate(list_idum(Nidum))
  open(unit=100,file=str(idumFILE))
  do i=1,Nidum
     read(100,*)id,list_idum(i)
  enddo
  close(100)


  call getcwd(here)


  allocate(w(Lfreq))
  allocate(Gi(Lfreq))
  !
  allocate(Aw(Nidum*Nlat,Lfreq))
  allocate(Egf(Lfreq),Egf2(Lfreq),Ggf(Lfreq),Sgf(Lfreq))


  do ispin=1,Nspin
     do iorb=1,Norb
        !
        iostr = "_l"//str(iorb)//"m"//str(iorb)//"_s"//str(ispin)
        !
        if(MPImaster)call start_timer(str(iostr))
        Aw = zero
        do i=1+MPIrank,Nidum,MPIsize
           idum = list_idum(i)
           !< Build up disorder:
           call setup_Abhz(verbose=.false.)
           dir="IDUM_"//str(list_idum(i))
           suffix = str(idum)//str(iostr)
           call chdir(str(dir))
           call system("tar -xjf tar_Gloc_"//str(suffix)//"_realw.tar.bz2")
           do ilat=1,Nlat
              call gf_read("Gloc_"//str(suffix)//"_realw_indx"//str(ilat,6)//".dat",w,Gi)
              io = ilat + (i-1)*Nlat
              Aw(io,:)  = abs(-dimag(Gi))/pi
           enddo
           call system("rm -rf Gloc_"//str(suffix)//"_realw_indx*.dat")
           call chdir(str(here))
           if(MPImaster)call eta(i,Nidum)
           !< Free memory:
           call free_Abhz()
        enddo
        if(MPImaster)call stop_timer()
#ifdef _SCALAPACK
        call MPI_AllReduce(MPI_IN_PLACE,Aw,size(Aw),MPI_Double_Precision,MPI_Sum,MPI_COMM_WORLD,mpi_ierr)
#endif
        !
        Egf = sum(Aw,dim=1)/Nlat/Nidum
        Egf2= sum(Aw*AW,dim=1)/Nlat/Nidum

        Ggf = sum(log(Aw),dim=1)/Nlat/Nidum
        Ggf = exp(Ggf)
        Sgf = sqrt(Egf2 - Egf**2)
        !
        call gf_plot("EAw"//str(iostr)//"_realw.dat",w,Egf,Sgf)
        call gf_plot("GAw"//str(iostr)//"_realw.dat",w,Ggf)
        !
     enddo
  enddo



  call end_parallel()


contains


  subroutine gf_read(pname,X,Y1)
    integer                       :: i,Np,unit
    character(len=*)              :: pname
    real(8),dimension(:)          :: X
    complex(8),dimension(size(X)) :: Y1
    real(8),dimension(size(X))    :: reY,imY
    open(free_unit(unit),file=reg(pname))
    Np=size(X)
    do i=1,Np
       read(unit,*)X(i),imY(i),reY(i)
    enddo
    Y1=dcmplx(reY,imY)
    close(unit)
  end subroutine gf_read


  subroutine gf_plot(pname,X,Y,S)
    integer                             :: i,Np,unit
    character(len=*)                    :: pname
    real(8),dimension(:)                :: X
    real(8),dimension(size(X))          :: Y
    real(8),dimension(size(X)),optional :: S
    open(free_unit(unit),file=reg(pname))
    Np=size(X)
    if(present(S))then
       do i=1,Np
          write(unit,*)X(i),Y(i),S(i)
       enddo
    else
       do i=1,Np
          write(unit,*)X(i),Y(i)
       enddo
    endif
    close(unit)
  end subroutine gf_plot


end program analysis
