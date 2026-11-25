module BHZ_LOCALIZER
  USE COMMON
  USE, INTRINSIC :: IEEE_ARITHMETIC, only: ieee_value,ieee_quiet_nan
  implicit none
  private

  public :: periodic_localizer_chern
  public :: localizer_scale_bound

contains

  ! Return the maximum absolute row sum of H-MU for square CSR SPH.
  function localizer_scale_bound(spH) result(bound)
    type(sparse_zmatrix_csr),intent(in) :: spH
    real(8)                            :: bound,row_sum
    complex(8)                         :: diagonal
    integer                            :: i
    !
    if(.not.spH%status.or.spH%Nrow<1.or.spH%Nrow/=spH%Ncol) &
      error stop "localizer_scale_bound: invalid SPH"
    bound = 0d0
    do i=1,spH%Nrow
      diagonal = spH%get(i,i)
      row_sum = sum(abs(spH%row(i)%vals))-abs(diagonal)+abs(diagonal-mu)
      bound = max(bound,row_sum)
    enddo
  end function localizer_scale_bound

  ! Return the periodic localizer index of conserved spin ISPIN=1/2.
  ! SPH is the site-major CSR Hamiltonian with 2 or 4 states per site;
  ! NINTERNAL selects the layout and SCALE is the positive eta.
  ! Serial row-slice elimination uses
  ! pivoted LAPACK blocks. An unsafe slice pivot triggers a dense
  ! fallback; an unresolved final inertia returns NaN.
  function periodic_localizer_chern(spH,ninternal,ispin,scale) result(chern)
    type(sparse_zmatrix_csr),intent(in) :: spH
    integer,intent(in)                 :: ninternal,ispin
    real(8),intent(in)                 :: scale
    real(8)                            :: chern
    complex(8),allocatable             :: onsite(:,:,:),link(:,:,:)
    complex(8)                         :: hopping,dirac
    real(8)                            :: theta_x,theta_y,mass
    integer                            :: i,j,k,x,y,yp,a,b,offset,ns,sig
    integer                            :: ix,iy,orb,other
    logical                            :: ok
    !
    if(ninternal/=2.and.ninternal/=4) &
      error stop "periodic_localizer_chern: invalid NINTERNAL"
    if(ispin<1.or.ispin>ninternal/2) &
      error stop "periodic_localizer_chern: invalid ISPIN"
    if(.not.spH%status.or.spH%Nrow/=ninternal*nx*ny.or. &
      spH%Ncol/=spH%Nrow.or.nx<2.or.ny<2) &
      error stop "periodic_localizer_chern: invalid SPH or lattice"
    !
    if(.not.(scale>0d0))error stop "periodic_localizer_chern: invalid SCALE"
    call start_timer("periodic_localizer_chern")
    offset = 2*(ispin-1)
    !
    ! Site-major torus blocks, Eq. (1) of Doll, Loring and
    ! Schulz-Baldes, doi:10.1007/s11040-025-09508-0.
    ns = 4*nx
    allocate(onsite(ns,ns,ny),link(ns,ns,ny),source=(0d0,0d0))
    do y=1,ny
      yp = modulo(y,ny)+1
      theta_y = pi2*(y-1)/ny
      do x=1,nx
        theta_x = pi2*(x-1)/nx
        mass = 2d0-cos(theta_x)-cos(theta_y)+mu/scale
        dirac = cmplx(sin(theta_x),-sin(theta_y),8)
        do orb=1,2
          a = 2*(x-1)+orb
          i = ninternal*((y-1)*nx+x-1)+offset+orb
          onsite(a,a,y) = onsite(a,a,y)+mass
          onsite(a+2*nx,a+2*nx,y) = onsite(a+2*nx,a+2*nx,y)-mass
          onsite(a,a+2*nx,y) = dirac
          onsite(a+2*nx,a,y) = conjg(dirac)
          do k=1,spH%row(i)%Size
            j = spH%row(i)%cols(k)
            hopping = spH%row(i)%vals(k)/scale
            if(hopping==(0d0,0d0))cycle
            other = modulo(j-1,ninternal)-offset+1
            if(other<1.or.other>2) &
              error stop "periodic_localizer_chern: SPH must conserve Sz"
            ix = modulo((j-1)/ninternal,nx)+1
            iy = (j-1)/(ninternal*nx)+1
            b = 2*(ix-1)+other
            if(iy==y)then
              onsite(a,b,y) = onsite(a,b,y)-hopping
              onsite(a+2*nx,b+2*nx,y) = onsite(a+2*nx,b+2*nx,y)+hopping
            elseif(iy==yp)then
              if(ny==2.and.y==2)cycle
              link(a,b,y) = link(a,b,y)-hopping
              link(a+2*nx,b+2*nx,y) = link(a+2*nx,b+2*nx,y)+hopping
            elseif(iy/=modulo(y-2,ny)+1)then
              error stop "periodic_localizer_chern: invalid row-slice hopping"
            endif
          enddo
        enddo
      enddo
    enddo
    !
    call slice_inertia(onsite,link,sig,ok)
    if(.not.ok)then
      write(*,*) "Periodic localizer: unsafe slice pivot, using dense inertia"
      call dense_inertia(onsite,link,sig,ok)
      if(.not.ok)then
        chern = ieee_value(0d0,ieee_quiet_nan)
        write(*,*) "Periodic localizer: unresolved inertia at scale =",scale
        call stop_timer()
        return
      endif
    endif
    ! Match the orientation of the BHZ single-point convention.
    chern = -0.5d0*sig
    write(*,*) "Periodic localizer: scale, signature, BHZ Chern =", &
      scale,sig,chern
    call stop_timer()
  end function periodic_localizer_chern


  ! Eliminate slices 2,...,NY-1 by Hermitian congruences. Keep slice 1
  ! as a border so each Schur update couples only two remaining rows.
  subroutine slice_inertia(onsite,link,sig,ok)
    complex(8),intent(in) :: onsite(:,:,:),link(:,:,:)
    integer,intent(out)   :: sig
    logical,intent(out)   :: ok
    complex(8),allocatable :: first(:,:),current(:,:),edge(:,:)
    complex(8),allocatable :: rhs(:,:),tail(:,:)
    integer,allocatable    :: piv(:)
    integer                :: b,m,y,part,info
    !
    b = size(onsite,1)
    m = size(onsite,3)
    allocate(rhs(b,2*b),tail(2*b,2*b),piv(b))
    first = onsite(:,:,1)
    current = onsite(:,:,2)
    edge = link(:,:,1)
    sig = 0
    do y=2,m-1
      call factor_inertia(current,piv,part,ok)
      if(.not.ok)return
      sig = sig+part
      rhs(:,:b) = conjg(transpose(edge))
      rhs(:,b+1:) = link(:,:,y)
      call zhetrs('L',b,2*b,current,b,piv,rhs,b,info)
      if(info/=0)error stop "slice_inertia: ZHETRS failed"
      first = first-matmul(edge,rhs(:,:b))
      current = onsite(:,:,y+1)- &
        matmul(conjg(transpose(link(:,:,y))),rhs(:,b+1:))
      edge = -matmul(edge,rhs(:,b+1:))
    enddo
    if(m>2)edge = edge+conjg(transpose(link(:,:,m)))
    tail(:b,:b) = first
    tail(b+1:,b+1:) = current
    tail(:b,b+1:) = edge
    tail(b+1:,:b) = conjg(transpose(edge))
    deallocate(piv)
    allocate(piv(2*b))
    call factor_inertia(tail,piv,part,ok)
    sig = sig+part
  end subroutine slice_inertia


  ! Factor the complete finite matrix only when slice pivoting fails.
  subroutine dense_inertia(onsite,link,sig,ok)
    complex(8),intent(in) :: onsite(:,:,:),link(:,:,:)
    integer,intent(out)   :: sig
    complex(8),allocatable :: dense(:,:)
    integer,allocatable    :: piv(:)
    integer                :: b,m,y,yp,i,j
    logical,intent(out)    :: ok
    !
    b = size(onsite,1)
    m = size(onsite,3)
    allocate(dense(b*m,b*m),source=(0d0,0d0))
    allocate(piv(b*m))
    do y=1,m
      i = b*(y-1)
      dense(i+1:i+b,i+1:i+b) = onsite(:,:,y)
      if(m==2.and.y==2)cycle
      yp = modulo(y,m)+1
      j = b*(yp-1)
      dense(i+1:i+b,j+1:j+b) = link(:,:,y)
      dense(j+1:j+b,i+1:i+b) = conjg(transpose(link(:,:,y)))
    enddo
    call factor_inertia(dense,piv,sig,ok)
  end subroutine dense_inertia


  ! Factor one Hermitian block in place and count signs of its
  ! Bunch-Kaufman pivots. Reject numerically unsafe inverses using
  ! a reciprocal condition estimate, not a localizer gap estimate.
  subroutine factor_inertia(a,piv,sig,ok)
    complex(8),intent(inout) :: a(:,:)
    integer,intent(out)      :: piv(:),sig
    logical,intent(out)      :: ok
    complex(8),allocatable   :: work(:)
    complex(8)               :: query(1)
    real(8)                  :: anorm,rcond,aa,bb,ab,s,det
    integer                  :: n,i,info,lwork
    !
    n = size(piv)
    anorm = maxval(sum(abs(a),dim=1))
    sig = 0
    ok = .false.
    call zhetrf('L',n,a,n,piv,query,-1,info)
    if(info/=0)error stop "factor_inertia: ZHETRF workspace query failed"
    lwork = max(2*n,int(real(query(1),8)))
    allocate(work(lwork))
    call zhetrf('L',n,a,n,piv,work,lwork,info)
    if(info<0)error stop "factor_inertia: ZHETRF failed"
    if(info>0)return
    call zhecon('L',n,a,n,piv,anorm,rcond,work,info)
    if(info/=0)error stop "factor_inertia: ZHECON failed"
    if(.not.(rcond>1d-12))return
    !
    i = 1
    do while(i<=n)
      aa = real(a(i,i),8)
      if(piv(i)>0)then
        sig = sig+merge(1,-1,aa>0d0)
        i = i+1
      else
        bb = real(a(i+1,i+1),8)
        ab = abs(a(i+1,i))
        s = max(abs(aa),abs(bb),ab)
        if(s==0d0)return
        det = (aa/s)*(bb/s)-(ab/s)**2
        if(det==0d0)return
        if(det>0d0)sig = sig+merge(2,-2,aa>0d0)
        i = i+2
      endif
    enddo
    ok = .true.
  end subroutine factor_inertia

end module BHZ_LOCALIZER
