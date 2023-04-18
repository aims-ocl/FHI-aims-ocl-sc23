
 subroutine x2c_fdirac_solver( dim_mat, vr,vi, wr,wi, tr,ti, sr,si, ssr,ssi, eigenvector, eigenval)
  use constants, only : light_speed
  implicit none
  integer,intent(in) :: dim_mat
  real*8,intent(in) :: sr(dim_mat*(2*dim_mat+1)),  si(dim_mat*(2*dim_mat+1)), &
                      ssr(dim_mat*(2*dim_mat+1)), ssi(dim_mat*(2*dim_mat+1)), &
                       tr(dim_mat*(2*dim_mat+1)),  ti(dim_mat*(2*dim_mat+1)), &
                       vr(dim_mat*(2*dim_mat+1)),  vi(dim_mat*(2*dim_mat+1)), &
                       wr(dim_mat*(2*dim_mat+1)),  wi(dim_mat*(2*dim_mat+1))
  real*8,intent(out) :: eigenvector(16*dim_mat*dim_mat,2), eigenval(4*dim_mat)
  integer :: i,j,k,l,m,n
  real*8,allocatable :: s_square(:,:,:),ss_square(:,:,:),t_square(:,:,:),v_square(:,:,:),w_square(:,:,:)
  real*8,allocatable :: ham_square(:,:,:), ovlp_square(:,:,:)
  real*8,allocatable :: ham(:,:), ovlp(:,:)

  allocate( s_square(2*dim_mat,2*dim_mat,2), ss_square(2*dim_mat,2*dim_mat,2), t_square(2*dim_mat,2*dim_mat,2), &
            v_square(2*dim_mat,2*dim_mat,2), w_square(2*dim_mat,2*dim_mat,2) )

  allocate( ham_square(4*dim_mat,4*dim_mat,2), ovlp_square(4*dim_mat,4*dim_mat,2) )
  allocate( ham(2*dim_mat*(4*dim_mat+1),2), ovlp(2*dim_mat*(4*dim_mat+1),2) )

  ham_square=0.d0; ovlp_square=0.d0; ham=0.d0; ovlp=0.d0

  ! Expand the LOWER triangular matrices to square matrice:
  m=0
  do i=1, 2*dim_mat
    do j=1, i
       m=m+1  
       !     Rel part:                   Imag part:
       s_square(j,i,1)  = sr(m);    s_square(j,i,2) = -si(m)
       s_square(i,j,1)  = sr(m);    s_square(i,j,2) =  si(m)
       ss_square(j,i,1) = ssr(m);  ss_square(j,i,2) = -ssi(m)
       ss_square(i,j,1) = ssr(m);  ss_square(i,j,2) =  ssi(m)
       t_square(j,i,1)  = tr(m);    t_square(j,i,2) = -ti(m)
       t_square(i,j,1)  = tr(m);    t_square(i,j,2) =  ti(m) 
       v_square(j,i,1)  = vr(m);    v_square(j,i,2) = -vi(m) 
       v_square(i,j,1)  = vr(m);    v_square(i,j,2) =  vi(m) 
       w_square(j,i,1)  = wr(m);    w_square(j,i,2) = -wi(m) 
       w_square(i,j,1)  = wr(m);    w_square(i,j,2) =  wi(m)
    enddo
  enddo

  ! Generate 4c DKS Fock matrix and overlap matrix:
  do i=1, 2*dim_mat
     do j=1, 2*dim_mat
        ham_square (j,i,:) = v_square(j,i,:)
        ovlp_square(j,i,:) = s_square(j,i,:)
     enddo
  enddo
  do i=2*dim_mat+1, 4*dim_mat
     do j=2*dim_mat+1, 4*dim_mat
        ham_square (j,i,:) = w_square ( j-2*dim_mat, i-2*dim_mat, : ) - 2*light_speed**2* ss_square( j-2*dim_mat, i-2*dim_mat, : )
        ovlp_square(j,i,:) = ss_square( j-2*dim_mat, i-2*dim_mat, : )
     enddo
  enddo
  do i=2*dim_mat+1, 4*dim_mat
     do j=1, 2*dim_mat
        ham_square(j,i,:) = t_square( j, i-2*dim_mat, : )
     enddo
  enddo
  do i=1, 2*dim_mat
     do j=2*dim_mat+1, 4*dim_mat
        ham_square(j,i,:) = t_square( j-2*dim_mat, i, : )
     enddo
  enddo

  m = 0
  do i=1, 4*dim_mat
  do j=1, i
    m=m+1
    ham(m,:) = ham_square(j,i,:)
    ovlp(m,:)=ovlp_square(j,i,:)
  enddo
  enddo

  call atom_diagonalization(2*dim_mat,ovlp,ham,eigenvector,eigenval)

 end subroutine x2c_fdirac_solver

subroutine x2c_calc_xmat_generator(ifrel,nspn,nsym,jsyml,norb,cr,ci,xr,xi)
!> calculate X matrix in X2C method: X == BA^{-1}
  implicit none
  integer,intent(in)  :: ifrel,nspn
  integer,intent(in)  :: nsym,jsyml(nsym),norb(nsym)
  real*8, intent(in)  :: cr(*),ci(*)  ! MO coef. (both negative and positive states inputed)
  real*8, intent(out) :: xr(*),xi(*)  ! symm. block square form
  integer :: q,p,nsos,isym,nn

  q   = 1
  p   = 1
  do isym=1,nsym
    nsos= norb(isym)
    if(jsyml(isym).eq.0.or.nsos.le.0) cycle ! skip diagonalization for equivalent representations
    nn  = nsos*nsos
    q   = q+nn+nn ! for MO Coef., first nsos*nsos*2 are negative states
    call x2c_bainv2xmat_onesym(ifrel,nsos,cr(q),ci(q),xr(p),xi(p))
    q   = q+nn+nn
    p   = p+nn
  enddo
end subroutine x2c_calc_xmat_generator

subroutine x2c_bainv2xmat_onesym(ifrel,n,cr,ci,xr,xi)
!> from gotten positive state MO coefficients:  A (large component) and B (small component)
!! to get X2C X matrix : AX = B
  implicit none
  integer,intent(in)  :: n,ifrel
  real*8, intent(in)  :: cr(n,2,n),ci(n,2,n)
  real*8, intent(out) :: xr(n,n),xi(n,n)
  complex*16, parameter   :: alph=dcmplx(1.d0,0.d0),beta=dcmplx(0.d0,0.d0)
  real*8,     allocatable :: ar(:,:),br(:,:)
  complex*16, allocatable :: az(:,:),bz(:,:),xz(:,:)
  integer :: i,j  

  if( ifrel.eq.0 )then  ! X1C case
   !allocate( ar(n,n),br(n,n) )
   !do i=1,n
   !  ar(:,i) = cr(:,1,i)
   !  br(:,i) = cr(:,2,i)
   !enddo
   !call dmatinv(ar,n)
   !call dgemm('N','N',n,n,n,1.d0,br,n,ar,n, 0.d0,xr,n)
   !deallocate( ar,br )
  else
    allocate( az(n,n),bz(n,n),xz(n,n) )
    do i=1,n
      az(:,i) = dcmplx( cr(:,1,i), ci(:,1,i) )
      bz(:,i) = dcmplx( cr(:,2,i), ci(:,2,i) )
    enddo
    call zmatinv(az,n)
    call zgemm('N','N',n,n,n,alph,bz,n,az,n, beta,xz,n)
    xr = dreal( xz )
    xi = dimag( xz )
    deallocate( az,bz,xz )
  endif
end subroutine x2c_bainv2xmat_onesym


subroutine x2c_gen_sescu_nesct_fw_mat(iesc,fwcr,ifrel,nsym,jsyml,norb,clight, &
  sr,si,tr,ti,xr,xi,fwr,fwi,tnr,tni,ssr,ssi,ur,ui)
!> generate FW trans. matrix and SESC U (\tilde{S}S^{-1}) for SESC method
!! and also NESC form kinetic matrix: X^{\dag}T + TX - X^{\dag}TX
  implicit none
  integer,intent(in)  :: iesc,fwcr,ifrel
  integer,intent(in)  :: nsym,jsyml(nsym),norb(nsym)
  real*8, intent(in)  :: sr(*),si(*),ssr(*),ssi(*),tr(*),ti(*),xr(*),xi(*),clight
  real*8, intent(out) :: fwr(*),fwi(*),ur(*),ui(*),tnr(*),tni(*)
  complex*16, parameter   :: alph = dcmplx(1.d0,0.d0), beta = dcmplx(0.d0,0.d0)
  real*8,     allocatable :: ds1(:),ds2(:),dt1(:),dtx(:)
  complex*16, allocatable :: zs1(:),zs2(:),zt1(:),ztx(:),zx(:)
  integer :: i,n,k,l,q,nn
  real*8  :: c2

  c2= 0.5d0/(clight*clight)
  k = 1
  l = 1
  do i=1,nsym
    n = norb(i)
    if( n.le.0.or.jsyml(i).le.0 ) cycle
    nn= n*n

    if(ifrel.eq.0)then
    ! scalar relativistic: X1C
    else
    ! 2c
      q = k+nn-1
      allocate( zs1(nn), zs2(nn), ztx(nn), zt1(nn), zx(nn) )
      call cvtzmx(n,sr(l),si(l),zs1)
     !call cvtzmx(n,ssr(l),ssi(l),zss1) ! aims
      call cvtzmx(n,tr(l),ti(l),zt1)
      zx(:) = dcmplx( xr(k:q),xi(k:q) )
      call zgemm('N','N',n,n,n,alph,zt1,n,zx,n,beta,ztx,n)    ! TX
      call zgemm('C','N',n,n,n,alph,zx,n,ztx,n,beta,zt1,n)    !!! X^{\dag}TX for BDF STO scheme
      zs2 = zs1 + c2*zt1        !!! \tilde{S} = S + X^{\dag}TX/2c^2 , this is the BDF STO scheme
   !!!call zgemm('N','N',n,n,n,alph,zss1,n,zx,n,beta,zssx,n)  ! (S_s)X
   !!!call zgemm('C','N',n,n,n,alph,zx,n,zssx,n,beta,zt1,n)   ! X^{\dag}(S_s)X
   !!!zs2 = zs1 + zt1        ! \tilde{S} = S + X^{\dag}(Ss)TX

   !!!zt1=zt1*2*clight**2 ! for aim NAO scheme
      call x2c_form_nesct_mat_onesym(n,ztx,zt1,tnr(l),tni(l)) ! nesc form kinetic matrix
      call x2c_fwmat_gen(fwcr,n,zs1,zs2,fwr(k),fwi(k))        ! FW trans. matrix

      if( iesc.eq.2 )then       ! SESC U matrix
        call zmatinv(zs1,n)     ! S^{-1}
        call zgemm('N','N',n,n,n,alph,zs2,n,zs1,n,beta,zx,n)
        ur(k:q) = dreal( zx(:) )
        ui(k:q) = dimag( zx(:) )
      endif

      deallocate( zs1,zs2,ztx,zt1,zx )
    endif

    l = l+(nn+n)/2
    k = k+nn
  enddo
end subroutine x2c_gen_sescu_nesct_fw_mat

!subroutine x2c_form_nesct_mat_onesym(n,tx,xssx,tnescr,tnesci) !!! aims NAO
subroutine x2c_form_nesct_mat_onesym(n,tx,xtx,tnescr,tnesci) !!! BDF
  implicit none
  integer,    intent(in)  :: n
 !complex*16, intent(in)  :: tx(n,n),xssx(n,n) !!! aims
  complex*16, intent(in)  :: tx(n,n),xtx(n,n) !!! BDF
  real*8,     intent(out) :: tnescr(n*(n+1)/2),tnesci(n*(n+1)/2)
  complex*16  :: aa
  integer     :: i,j,k
  k = 1
  do i=1,n
  do j=1,i
   !aa = tx(j,i) + dconjg(tx(i,j)) - xssx(j,i) !!! aims
    aa = tx(j,i) + dconjg(tx(i,j)) - xtx(j,i)  !!! BDF
    tnescr(k) = dreal(aa)
    tnesci(k) =-dimag(aa)
    k = k+1
  enddo
  enddo
end subroutine x2c_form_nesct_mat_onesym


subroutine x2c_gen_sescu_nesct_fw_mat_atom(iesc,fwcr,ifrel,nsym,jsyml,norb,clight, &
  sr,si,tr,ti,tsr,tsi,xr,xi,fwr,fwi,tnr,tni,ssr,ssi,ur,ui)
!> generate FW trans. matrix and SESC U (\tilde{S}S^{-1}) for SESC method
!! and also NESC form kinetic matrix: X^{\dag}T + TX - X^{\dag}TX
  implicit none
  integer,intent(in)  :: iesc,fwcr,ifrel
  integer,intent(in)  :: nsym,jsyml(nsym),norb(nsym)
  real*8, intent(in)  :: sr(*),si(*),ssr(*),ssi(*),tr(*),ti(*),tsr(*),tsi(*),xr(*),xi(*),clight
  real*8, intent(out) :: fwr(*),fwi(*),ur(*),ui(*),tnr(*),tni(*)
  complex*16, parameter   :: alph = dcmplx(1.d0,0.d0), beta = dcmplx(0.d0,0.d0)
  real*8,     allocatable :: ds1(:),ds2(:),dt1(:),dtx(:)
  complex*16, allocatable :: zs1(:),zss1(:),zs2(:),zt1(:),zt2(:),ztx(:),zxt(:),zssx(:),zx(:)
  integer :: i,n,k,l,q,nn
  real*8  :: c2

  c2= 0.5d0/(clight*clight)
  k = 1
  l = 1
  do i=1,nsym
    n = norb(i)
    if( n.le.0.or.jsyml(i).le.0 ) cycle
    nn= n*n

    if(ifrel.eq.0)then
    ! X1C can be added here
    else
    ! 2c
      q = k+nn-1
      allocate( zs1(nn), zss1(nn), zs2(nn), ztx(nn), zxt(nn), zssx(nn), zt1(nn), zt2(nn), zx(nn) )
      call cvtzmx(n,sr(l),si(l),zs1)
      call cvtzmx(n,ssr(l),ssi(l),zss1)
      call cvtzmx(n,tr(l),ti(l),zt1)
      call cvtzmx(n,tsr(l),tsi(l),zt2)
      zx(:) = dcmplx( xr(k:q),xi(k:q) )
      call zgemm('N','N',n,n,n,alph,zt1,n,zx,n,beta,ztx,n)    ! TX
     !call zgemm('C','N',n,n,n,alph,zx,n,ztx,n,beta,zt1,n)    ! X^{\dag}TX
     !zs2 = zs1 + c2*zt1        ! \tilde{S} = S + X^{\dag}TX/2c^2
      call zgemm('C','N',n,n,n,alph,zx,n,zt2,n,beta,zxt,n)    ! X^{\dag}T_s
      call zgemm('N','N',n,n,n,alph,zss1,n,zx,n,beta,zssx,n)  ! (S_s)X
      call zgemm('C','N',n,n,n,alph,zx,n,zssx,n,beta,zt1,n)   ! X^{\dag}(S_s)X
      zs2 = zs1 + zt1        ! \tilde{S} = S + X^{\dag}(Ss)TX

      zt1=zt1*2*clight**2
      call x2c_form_nesct_mat_atom(n,ztx,zxt,zt1,tnr(l),tni(l)) ! nesc form kinetic matrix
      call x2c_fwmat_gen(fwcr,n,zs1,zs2,fwr(k),fwi(k))         ! FW trans. matrix

      if( iesc.eq.2 )then       ! SESC U matrix
        call zmatinv(zs1,n)     ! S^{-1}
        call zgemm('N','N',n,n,n,alph,zs2,n,zs1,n,beta,zx,n)
        ur(k:q) = dreal( zx(:) )
        ui(k:q) = dimag( zx(:) )
      endif

      deallocate( zs1,zss1,zs2,ztx,zxt,zssx,zt1,zt2,zx )
    endif

    l = l+(nn+n)/2
    k = k+nn
  enddo
end subroutine x2c_gen_sescu_nesct_fw_mat_atom

subroutine x2c_form_nesct_mat_atom(n,tx,xt,xssx,tnescr,tnesci)
  implicit none
  integer,    intent(in)  :: n
  complex*16, intent(in)  :: tx(n,n),xt(n,n),xssx(n,n)
  real*8,     intent(out) :: tnescr(n*(n+1)/2),tnesci(n*(n+1)/2)
  complex*16  :: aa
  integer     :: i,j,k
  real*8      :: tsqr(n,n),tsqi(n,n)
  k = 1
  do i=1,n
  do j=1,n
    aa = tx(j,i) + xt(j,i) - xssx(j,i)
    tsqr(j,i) = dreal(aa)
    tsqi(j,i) = dimag(aa)
    k = k+1
  enddo
  enddo

  k = 1
  do i=1,n
  do j=1,i
    tnescr(k) = tsqr(i,j)
    tnesci(k) = tsqi(i,j)
    k = k+1
  enddo
  enddo

end subroutine x2c_form_nesct_mat_atom

subroutine x2c_fwmat_gen(iop,n,si,si2,fwmr,fwmi)
!> generate FW transformation matrix for X1c method
!! iop: control paramter for FW transformation matrix
!! =0 fwmat=\tilde{s}^{-1/2}s^{1/2}
!! =1 fwmat=s^{-1/2} (s^{-1/2}\tilde{s}s^{-1/2})^{-1/2} s^{1/2}
!! <0 fwmat= I (not really do FW transformation)
  use localorb_io, only: use_unit
  implicit none
  integer,    intent(in)  :: n,iop
  complex*16, intent(in)  :: si2(n,n),si(n,n)     ! \tilde{S}, S matrix
  real*8,     intent(out) :: fwmr(n,n),fwmi(n,n)  ! FW matrix: real/imag
  complex*16, allocatable :: s(:,:),s2(:,:),t(:,:)
  complex*16, parameter   :: alph=dcmplx(1.d0,0.d0),beta=dcmplx(0.d0,0.d0)
  integer :: i

  if( iop.lt.0 )then
  ! identity
    fwmr= 0.d0
    fwmi= 0.d0
    do i=1,n
      fwmr(i,i) = 1.d0
    enddo
    return
  endif

  allocate( t(n,n),s(n,n),s2(n,n) )
  s = si
  if( iop.eq.0 )then
  ! =0 fwmat=\tilde{s}^{-1/2}s^{1/2}
    s2 = si2
    call zmatsqrtz(s, n,1)   ! s^{1/2}
    call zmatsqrtz(s2,n,0)   ! \tilde{s}^{-1/2}
    call zgemm('N','N',n,n,n,alph,s2,n,s,n,beta,t,n)
  elseif( iop.eq.1.or.iop.eq.2 )then
  ! =1 fwmat=s^{-1/2} (s^{-1/2}\tilde{s}s^{-1/2})^{-1/2} s^{1/2}
    t = si2
    call zmatsqrtz(s,n,0)  ! s^{-1/2}
    call zgemm('N','N',n,n,n,alph,s, n,t, n,beta,s2,n)
    call zgemm('N','N',n,n,n,alph,s2,n,s, n,beta,t, n)
    call zmatsqrtz(t, n,0)
    call zgemm('N','N',n,n,n,alph,s, n,t, n,beta,s2,n)
    s = si
    call zmatsqrtz(s,n,1)  ! s^{ 1/2}
    call zgemm('N','N',n,n,n,alph,s2,n,s, n,beta,t, n)
  else
    write(use_unit,'(a,i9)') " Wrong input control of FW matrix: iop = ",iop
    stop
  endif

  fwmr = dreal(t)
  fwmi = dimag(t)
  deallocate( t,s,s2 )
end subroutine x2c_fwmat_gen


subroutine x2c_w_to_xwx_nao(n,xr,xi,wqr,wqi)
  implicit none
  integer,intent(in)  :: n
  real*8, intent(in)  :: xr(n*n), xi(n*n)
  real*8, intent(inout) :: wqr(n*n), wqi(n*n)
  complex*16, parameter   :: alph=dcmplx(1.d0,0.d0),beta=dcmplx(0.d0,0.d0)
  complex*16, allocatable :: zt(:),zw(:),zx(:)
  integer :: i,j,nn

  nn = n*n

  allocate( zw(nn),zt(nn),zx(nn) )
  zw(1:nn) = dcmplx( wqr(1:nn),wqi(1:nn) )
  zx(1:nn) = dcmplx( xr(1:nn),xi(1:nn) )
  call zgemm('N','N',n,n,n,alph,zw,n,zx,n,beta,zt,n) ! WX
  call zgemm('C','N',n,n,n,alph,zx,n,zt,n,beta,zw,n) ! X^{\dag}WX
  wqr(1:nn) = dreal(zw(1:nn))
  wqi(1:nn) = dimag(zw(1:nn))
  deallocate( zw,zt,zx )
end subroutine x2c_w_to_xwx_nao

subroutine x2c_w_to_xwx(ifrel,nspn,nsym,jsyml,norb,xr,xi,wr,wi)
  implicit none
  integer,intent(in)    :: ifrel,nspn
  integer,intent(in)    :: nsym,jsyml(nsym),norb(nsym)
  real*8, intent(in)    :: xr(*),xi(*)
  real*8, intent(inout) :: wr(*),wi(*)
  complex*16, parameter   :: alph=dcmplx(1.d0,0.d0),beta=dcmplx(0.d0,0.d0)
  real*8,     allocatable :: tt(:),tw(:)
  complex*16, allocatable :: zt(:),zw(:),zx(:)
  integer :: k,i,n,nn,q,l,j

  l = 1
  do j=1,nspn
    k = 1
    do i=1,nsym
      n = norb(i)
      if( n.le.0.or.jsyml(i).le.0 ) cycle

      nn  = n*n
      if(ifrel.eq.0)then
    ! ! scalar relativistic: X1C
    !   allocate( tt(nn),tw(nn) )
    !   call cvtsmx(n,wr(l),tw )
    !   call dgemm('N','N',n,n,n,1.d0,tw,n,xr(k),n,0.d0,tt,n) ! WX
    !   call dgemm('C','N',n,n,n,1.d0,xr(k),n,tt,n,0.d0,tw,n) ! X^{\dag}WX
    !   call cvtsmx2(n,wr(l),tw)
    !   deallocate( tt,tw )
      else
      ! 2c
        q = k+nn-1
        allocate( zw(nn),zt(nn),zx(nn) )
        call cvtzmx(n,wr(l),wi(l),zw)
        zx(1:nn)= dcmplx( xr(k:q),xi(k:q) )
        call zgemm('N','N',n,n,n,alph,zw,n,zx,n,beta,zt,n) ! WX
        call zgemm('C','N',n,n,n,alph,zx,n,zt,n,beta,zw,n) ! X^{\dag}WX
        call cvtrmx(n,wr(l),wi(l),zw)
        deallocate( zw,zt,zx )
      endif
      k = k+nn
      l = l+(nn+n)/2
    enddo
  enddo
end subroutine x2c_w_to_xwx



subroutine x2c_sesc_gen_total_fockmat(ifrel,nspn,nsym,jsyml,norb,xr,xi,tr,ti,vr,vi,ur,ui,fr,fi)
  implicit none
  integer,intent(in)    :: ifrel,nspn
  integer,intent(in)    :: nsym,jsyml(nsym),norb(nsym)
  real*8, intent(in)    :: xr(*),xi(*),tr(*),ti(*),vr(*),vi(*),ur(*),ui(*)
  real*8, intent(out)   :: fr(*),fi(*)
  real*8,     allocatable :: dt(:),df(:)
  complex*16, parameter   :: alph=dcmplx(1.d0,0.d0),beta=dcmplx(0.d0,0.d0)
  complex*16, allocatable :: zt(:),zf(:),zx(:)
  integer :: j,i,l,k,q,n,nn

  l = 1
  do j=1,nspn
    k = 1
    do i=1,nsym
      n = norb(i)
      if( n.le.0.or.jsyml(i).le.0 ) cycle

      nn  = n*n
      if(ifrel.eq.0)then
      ! scalar relativistic: SESC1C
        allocate( dt(nn),df(nn) )
        call cvtsmx(n,tr(l),dt)
        call dgemm('N','N',n,n,n,1.d0,dt,n,xr(k),n,0.d0,df,n) ! TX
        call cvtsmx(n,vr(l),dt)
        df = df+dt  ! TX +V
        call dgemm('N','N',n,n,n,1.d0,ur(k),n,df,n,0.d0,dt,n) ! U(TX+V)
        call cvtsmx2(n,fr(l),dt) ! U(TX+V)+c.c. ==> triangular form
        deallocate(dt,df)
      else
      ! 2c
        q = k+nn-1
        allocate( zt(nn),zf(nn),zx(nn) )
        call cvtzmx(n,tr(l),ti(l),zt)
        zx(:) = dcmplx( xr(k:q),xi(k:q) )
        call zgemm('N','N',n,n,n,alph,zt,n,zx,n,beta,zf,n) ! TX
        call cvtzmx(n,vr(l),vi(l),zt)
        zf = zf+zt  ! TX +V
        zx(:) = dcmplx( ur(k:q),ui(k:q) )
        call zgemm('N','N',n,n,n,alph,zx,n,zf,n,beta,zt,n) ! U(TX+V)
        call cvtrmx(n,fr(l),fi(l),zt) ! U(TX+V)+c.c. ==> triangular form
        deallocate( zt,zf,zx )
      endif
      k = k+nn
      l = l+(nn+n)/2
    enddo
  enddo
end subroutine x2c_sesc_gen_total_fockmat

