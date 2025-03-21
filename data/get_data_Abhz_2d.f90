program analysis
  USE COMMON
  implicit none



  integer,parameter                     :: Norb=2,Nspin=2
  real(8)                               :: bT
  !Solution
  real(8),dimension(:),allocatable      :: Edata,Tzdata,Szdata,Gapdata,PSzPdata,Z2data
  real(8),dimension(:),allocatable      :: Tzii,Szii,Z2ii
  integer                               :: io,i,id
  type(pdf_kernel)                      :: pdf_E,pdf_Tz,pdf_Sz,pdf_Gap,pdf_PSzP,pdf_Z2

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



  allocate(E(Nlso),Epsp(Nocc))
  allocate(Tzii(Nlat))
  allocate(Szii(Nlat))
  allocate(Z2ii(1))

  Nidum = file_length(str(idumFILE))
  allocate(list_idum(Nidum))
  open(unit=100,file=str(idumFILE))
  do i=1,Nidum
     read(100,*)id,list_idum(i)
  enddo
  close(100)


  call getcwd(here)
  call system("rm -fv pdf_{E,PSzP,Tz,Sz,Egap}.dat")

  do i=1,Nidum
     idum = list_idum(i)
     dir="IDUM_"//str(list_idum(i))
     print*,dir
     call chdir(str(dir))

     call read_array("E.dat",E)
     call read_array("Epszp.dat",Epsp)
     call read_array("tz.dat",Tzii)
     call read_array("sz.dat",Szii)
     call read_array("z2.dat",z2ii)
     !
     call save_array("Egap.dat",E(Nocc+1)-E(Nocc))
     call save_array("PSzPgap.dat",Epsp(Nocc/2+1)-Epsp(Nocc/2))

     ! Szii = Szii*afm_vec()
     call add_to(Edata,E)
     call add_to(PSzPdata,Epsp)
     call add_to(Tzdata,Tzii)
     call add_to(Szdata,Szii)
     call add_to(Z2data,z2ii)
     call add_to(Gapdata,[E(Nocc+1)-E(Nocc)])

     call chdir(str(here))
  enddo


  ! call save_array("Szdata.dat",Szdata)
  call save_array("Etz.dat",get_mean(abs(Tzdata)))
  call save_array("Esz.dat",get_mean(abs(Szdata)))
  call save_array("Ez2.dat",get_mean(Z2data))

  pdf_E = get_pdf_from_data(Edata,500,ieta,ini=1d-5)
  call pdf_print(pdf_E,"pdf_E.dat")

  pdf_PSzP = get_pdf_from_data(PSzPdata,500,ieta)
  call pdf_print(pdf_PSzP,"pdf_PSzP.dat")


  pdf_Tz = get_pdf_from_data(-Tzdata,500,ini=1d-5,end=1d0)
  call pdf_print(pdf_Tz,"pdf_Tz.dat")  

  pdf_Sz = get_pdf_from_data(Szdata,500,ini=-2d0,end=2d0)
  call pdf_print(pdf_Sz,"pdf_Sz.dat")  

  pdf_Gap = get_pdf_from_data(Gapdata,500,ini=1d-5)
  call pdf_print(pdf_Gap,"pdf_Egap.dat")  


  pdf_Z2 = get_pdf_from_data(Z2data,500,ini=1d-5)
  call pdf_print(pdf_Z2,"pdf_Z2.dat")  


contains


  function get_pdf_from_data(data,N,sdev,ini,end) result(pdf)
    real(8),dimension(:),intent(in) :: data
    integer,intent(in)              :: N
    real(8),optional                :: sdev,ini,end  
    type(pdf_kernel)                :: pdf
    real(8)                         :: a,b
    real(8)                         :: sigma
    !
    a = minval(data) ;a = a - 0.1d0*abs(a)
    b = maxval(data) ;b = b + 0.1d0*abs(b)
    if(present(ini))a = ini
    if(present(end))b = end
    !
    call pdf_allocate(pdf,N)
    call pdf_set_range(pdf,a,b)
    call pdf_sigma(pdf,data,sigma)
    if(present(sdev))sigma=sdev
    print*,"PDF a,b,Sigma=",a,b,sigma
    call pdf_push_sigma(pdf,sigma)
    call pdf_accumulate(pdf,data)
  end function get_pdf_from_data




end program analysis
