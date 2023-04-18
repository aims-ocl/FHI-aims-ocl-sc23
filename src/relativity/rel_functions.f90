
! This file contains functions/subroutines that are used in fully-relativistic integrations.
! -- Rundong Zhao, Oct. 2018

 subroutine tran_mat2kspace( n, mat_real, mat_complex )
  use pbc_lists, only: n_cells_in_hamiltonian, index_hamiltonian, column_index_hamiltonian, k_phase
  use dimensions, only: n_basis
  use scalapack_wrapper, only: my_k_point
  implicit none
  integer,intent(in) :: n
  real*8,intent(in) :: mat_real(n*n)
  complex*16,intent(out) :: mat_complex(n,n)

  integer:: i_cell, i_col, i_row, idx

  mat_complex = (0.d0, 0.d0)

  do i_cell = 1, n_cells_in_hamiltonian-1
     do i_col = 1, n_basis
        if(index_hamiltonian(1,i_cell,i_col) > 0) then
        do idx = index_hamiltonian(1,i_cell,i_col),index_hamiltonian(2,i_cell,i_col)
           i_row = column_index_hamiltonian(idx)
           mat_complex(i_row,i_col) = mat_complex(i_row,i_col) + k_phase(i_cell,my_k_point) * mat_real(idx)
        end do
        end if
     end do
  end do ! i_cell

 end subroutine tran_mat2kspace


      subroutine apbtoc(n,dx,incx,dy,incy,dz,incz)
!
!     dz = dy + dx
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
      double precision dx(*),dy(*),dz(*)
      integer i,incx,incy,incz,ix,iy,iz,m,mp1,n
 
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1.and.incz.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      iz = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      if(incz.lt.0)iz = (-n+1)*incz + 1
      do 10 i = 1,n
        dz(iz) = dy(iy) + dx(ix)
        ix = ix + incx
        iy = iy + incy
        iz = iz + incz
   10 continue
      return
!
!        code for both increments equal to 1
!
!        clean-up loop
!
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dz(i) = dy(i) + dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dz(i    ) = dy(i    ) + dx(i    )
        dz(i + 1) = dy(i + 1) + dx(i + 1)
        dz(i + 2) = dy(i + 2) + dx(i + 2)
        dz(i + 3) = dy(i + 3) + dx(i + 3)
   50 continue
      return
      end



!-----------------------------------------------------------------------------       
      subroutine cvtzmx(n,ar,ai,az)
!
! Combine two real triangle matrices to a square complex matrix
!      
      implicit double precision (a-h,o-z)
      dimension ar(*),ai(*)
      complex*16 az(*)      
      k=0
      do i=1,n
      do j=1,i
      k=k+1
      l=(j-1)*n+i
      az(l)=dcmplx(ar(k), ai(k))
      l=(i-1)*n+j
      az(l)=dcmplx(ar(k),-ai(k))
      enddo
      enddo      
      return
      end
!----------------------------------------------------------------------------- 



!-----------------------------------------------------------------------------   
      subroutine cvtsmx2(n,ar,az)
!
! Transform a real square matrices to a real triangle matrix
!        
      implicit double precision (a-h,o-z)
      dimension ar(*),az(*)     
      k=0
      do i=1,n
      do j=1,i
      k=k+1
      l1=(j-1)*n+i
      l2=(i-1)*n+j
      ar(k)=(az(l1)+az(l2))*0.5d0
      enddo
      enddo      
      return
      end        
!-----------------------------------------------------------------------------



!-----------------------------------------------------------------------------       
      subroutine cvtsmx(n,ar,az)
!
! Transform a real triangle matrices to a square real matrix
!        
      implicit double precision (a-h,o-z)
      dimension ar(*),az(*)     
      k=0
      do i=1,n
      do j=1,i
      k=k+1
      l=(j-1)*n+i
      az(l)=ar(k)
      l=(i-1)*n+j
      az(l)=ar(k)
      enddo
      enddo      
      return
      end        
!----------------------------------------------------------------------------- 



!----------------------------------------------------------------------------- 
      subroutine zmatinv(a,n)
!
! Invert a complex square matrix
!      
      use localorb_io, only: use_unit
      implicit none
      integer,    intent(in)    :: n
      complex*16, intent(inout) :: a(n*n)
      integer,    allocatable   :: ipiv(:)
      complex*16, allocatable   :: wrk(:)
      integer :: info1,info2
      info1 = 0
      info2 = 0
      allocate( wrk(n),ipiv(n) )
      call zgetrf(n,n,a(1),n,ipiv,info1)
      if( info1.ne.0 )then
        write(use_unit,'(a,i10)') ' zgetrf failed in zmatinv: info=',info1
        stop
      endif
      call zgetri(n  ,a(1),n,ipiv,wrk,n,info2)
      if( info2.ne.0 )then
        write(use_unit,'(a,i10)') ' zgetri failed in zmatinv: info=',info2
        stop
      endif
      deallocate( wrk,ipiv )
      end subroutine zmatinv
!-----------------------------------------------------------------------------



!-----------------------------------------------------------------------------      
      subroutine cvtrmx(n,ar,ai,az)      
!
! Convert a square complex matrix to two real triangle matrices 
!       
      implicit double precision (a-h,o-z)
      dimension ar(*),ai(*)
      complex*16 az(*)      
      k=0
      do i=1,n
      do j=1,i
      k=k+1
      l1=(j-1)*n+i
      l2=(i-1)*n+j
      ar(k)=(dreal(az(l1))+dreal(az(l2)))*0.5d0
      ai(k)=(dimag(az(l1))-dimag(az(l2)))*0.5d0
      enddo
      enddo      
      return
      end    
!----------------------------------------------------------------------------- 




subroutine x2c_calc_xmat_get_fdirac(nsos,vr,vi,tr,ti,wr,wi,ssr,ssi,sh,th,c2,f1,s1)
!> get full Dirac fock and root square of overlap matrix

!(Rundong) The following comments are for BDF codes:
!!  f = ( V  T ),  s = ( S   0    )^{-1/2}
!!      ( T W-T)       ( 0  T/2c^2)
  implicit none
  integer,    intent(in)  :: nsos
  real*8,     intent(in)  :: vr(*),vi(*),tr(*),ti(*),wr(*),wi(*),ssr(*),ssi(*),c2
  complex*16, intent(in)  :: sh(nsos,nsos),th(nsos,nsos)  ! S^{-1/2}, T^{-1/2}
  complex*16, intent(out) :: f1(nsos,2,nsos,2),s1(nsos,2,nsos,2)
  complex*16, parameter :: z0 = dcmplx(0.d0,0.d0)
  integer :: k,i,j
  real*8 :: mc2

  mc2=c2*c2 ! 2*c^2

  k = 0
  do i=1,nsos
  do j=1,i
    k = k+1
    f1(j,1,i,1) = dcmplx( vr(k), -vi(k) )
    f1(i,1,j,1) = dconjg( f1(j,1,i,1) )
    f1(j,2,i,1) = dcmplx( tr(k), -ti(k) )
    f1(i,2,j,1) = dconjg( f1(j,2,i,1) )
    f1(i,1,j,2) = f1(i,2,j,1)
    f1(j,1,i,2) = f1(j,2,i,1)
   !f1(j,2,i,2) = dcmplx( wr(k)-tr(k), ti(k)-wi(k) )           ! BDF-STO
    f1(j,2,i,2) = dcmplx( wr(k)-tr(k), ti(k)-wi(k) )           ! aims-STO scheme
   !f1(j,2,i,2) = dcmplx( wr(k)-mc2*ssr(k), mc2*ssi(k)-wi(k) ) ! aims NAO scheme
    f1(i,2,j,2) = dconjg( f1(j,2,i,2) )
    s1(j,1,i,1) = sh(j,i)
    s1(i,1,j,1) = sh(i,j)
    s1(j,2,i,1) = z0
    s1(i,1,j,2) = z0
    s1(j,1,i,2) = z0
    s1(i,2,j,1) = z0
    s1(j,2,i,2) = th(j,i)*c2  ! STO scheme
   !s1(j,2,i,2) = th(j,i)     ! NAO scheme
    s1(i,2,j,2) = dconjg( s1(j,2,i,2) )
  enddo
  enddo
end subroutine x2c_calc_xmat_get_fdirac



subroutine x2c_calc_xmat_get_fdirac_nao(nsos,vr,vi,tr,ti,tsr,tsi,wr,wi,ssr,ssi,sh,th,c2,f1,s1)
!> get full Dirac fock and root square of overlap matrix

!(Rundong) The following comments are for BDF codes:
!!  f = ( V  T ),  s = ( S   0    )^{-1/2}
!!      ( T W-T)       ( 0  T/2c^2)
!(Rundong) For aims, we modify it to:
!!  f = ( V      T    ),  s = ( S   0 )^{-1/2}
!!      ( Ts W-2c^2*Ss)       ( 0  Ss )
  implicit none
  integer,    intent(in)  :: nsos
  real*8,     intent(in)  :: vr(*),vi(*),tr(*),ti(*),tsr(*),tsi(*),wr(*),wi(*),ssr(*),ssi(*),c2
  complex*16, intent(in)  :: sh(nsos,nsos),th(nsos,nsos)  ! S^{-1/2}, Ss^{-1/2}
  complex*16, intent(out) :: f1(nsos,2,nsos,2),s1(nsos,2,nsos,2)
  complex*16, parameter :: z0 = dcmplx(0.d0,0.d0)
  integer :: k,i,j
  real*8 :: mc2

  mc2=c2*c2 ! 2*c^2

  k = 0
  do i=1,nsos
  do j=1,i
    k = k+1
    f1(j,1,i,1) = dcmplx( vr(k), -vi(k) )
    f1(i,1,j,1) = dconjg( f1(j,1,i,1) )
    f1(j,2,i,1) = dcmplx( tsr(k), -tsi(k) )
    f1(i,2,j,1) = dconjg( f1(j,2,i,1) )
    f1(j,1,i,2) = dcmplx( tr(k), -ti(k) )
    f1(i,1,j,2) = dconjg( f1(j,1,i,2) )
   !f1(j,2,i,2) = dcmplx( wr(k)-tr(k), ti(k)-wi(k) )           ! BDF
    f1(j,2,i,2) = dcmplx( wr(k)-mc2*ssr(k), mc2*ssi(k)-wi(k) ) ! aims
    f1(i,2,j,2) = dconjg( f1(j,2,i,2) )
    s1(j,1,i,1) = sh(j,i)
    s1(i,1,j,1) = sh(i,j)
    s1(j,2,i,1) = z0
    s1(i,1,j,2) = z0
    s1(j,1,i,2) = z0
    s1(i,2,j,1) = z0
   !s1(j,2,i,2) = th(j,i)*c2  ! BDF
    s1(j,2,i,2) = th(j,i)     ! aims
    s1(i,2,j,2) = dconjg( s1(j,2,i,2) )
  enddo
  enddo
end subroutine x2c_calc_xmat_get_fdirac_nao







subroutine zmatsqrtz(a,n,icase)
!> Compute a square root of a complex Hermitian matrix
!! icase : =0, a^{-1/2}; =1 a^{1/2}; else a^{-1}
  implicit none
  integer,    intent(in)    :: n,icase
  complex*16, intent(inout) :: a(n,n)
  complex*16, parameter     :: alph = dcmplx(1.d0,0.d0), beta = dcmplx(0.d0,0.d0)
  real*8,     parameter     :: dm = 1.d-30
  real*8,     allocatable   :: e(:),rwork(:)
  complex*16, allocatable   :: c(:,:),zwork(:)
  integer :: info,i,lwork,j

  if(icase.ne.1.and.icase.ne.0)then
    call zmatinv(a,n)
    return
  endif

  lwork = 8*n
  allocate( e(n),c(n,n),zwork(lwork),rwork(n*3) )
  c = a
  info  = 0
  call zheev( 'V', 'U', n, c, n, e, zwork, lwork, rwork, info)

  if(icase.eq.1)then
    e = dsqrt(dsqrt(e))
    do i=1,n
      c(:,i)= c(:,i)*e(i)
    enddo
  elseif(icase.eq.0)then
    do i=1,n
      if( e(i).le.dm )then
        c(:,i)= beta
      else
        e(i)  = 1.d0/dsqrt(dsqrt(e(i)))
        c(:,i)= c(:,i)*e(i)
      endif
    enddo
  endif
 ! call zherk('U','N',n,n,alph,c,n,beta,a,n) ! CC^{\dag} ==> A
 ! do i=2,n
 ! do j=1,i-1
 !   a(i,j)= a(j,i)
 ! enddo
 ! enddo
  call zgemm('N','C',n,n,n,alph,c,n,c,n, beta,a,n)
  deallocate( e,c,zwork,rwork )
end subroutine zmatsqrtz


 subroutine pbc_trans_fabc_to_fflt(n,k,fabc,fflt)
  implicit none
  integer,    intent(in)  :: n,k
  complex*16, intent(in)  :: fabc(n,n,3,k)
  real*8,     intent(out) :: fflt(n*(n+n+1),k,2)
  complex*16  :: aa
  integer     :: i,l,u,v

  do i = 1,k
    call pbc_expand_mat_2_lr(n,fabc(:,:,:,i),fflt(:,i,1),fflt(:,i,2),1)
  enddo

  return

  do i = 1,k
    l = 0
    do u = 1,n
    do v = 1,u
      l = l+1
      aa = (fabc(v,u,1,i)+dconjg(fabc(u,v,1,i)))*0.5d0
      fflt(l,i,1) = dreal(aa)
      fflt(l,i,2) = dimag(aa)
    enddo
    enddo

    do u = 1,n
      do v = 1,n
        l = l+1
        aa = fabc(v,u,2,i)
        fflt(l,i,1) = dreal(aa)
        fflt(l,i,2) = dimag(aa)
      enddo

      do v = 1,u
        l = l+1
        aa = (fabc(v,u,3,i)+dconjg(fabc(u,v,3,i)))*0.5d0
        fflt(l,i,1) = dreal(aa)
        fflt(l,i,2) = dimag(aa)
      enddo
    enddo

  enddo
 end subroutine pbc_trans_fabc_to_fflt

 subroutine pbc_expand_mat_2_lr(n,mq,mlr,mli,iop)
  implicit none
  integer,    intent(in)  :: iop,n
  complex*16, intent(in)  :: mq(n*n*3)
  real*8,     intent(out) :: mlr(*),mli(*)
  real*8, allocatable :: tt(:,:,:)
  integer :: i,j,m,n2
  n2 = n+n
  allocate( tt(n2,n2,2) )
  call pbc_expand_wmat(n,mq,tt(1,1,1),tt(1,1,2) )
  if( iop.ne.0 ) call pbc_x2c_scale_kinetic(n2,tt,tt(1,1,2))  ! scale kinetic matrix
  m = 1
  do i = 1,n2
  do j = 1,i
    mlr(m) = tt(i,j,1)
    mli(m) = tt(i,j,2)
    m=m+1
  enddo
  enddo
  deallocate( tt )
 end subroutine pbc_expand_mat_2_lr

 subroutine pbc_expand_wmat(nao,mat,wr,wi)
  implicit none
  integer,intent(in) :: nao
  complex*16,intent(in) :: mat(nao,nao,3)  ! 1 A, 2 B, 3 C
  real*8,intent(out) :: wr(2*nao,2*nao),wi(2*nao,2*nao)

  integer i,j,k,m,n
  complex*16,allocatable :: bdag(:,:)

  allocate( bdag(nao,nao) )
  
  do i=1,nao
    do j=1,nao
      bdag(j,i)=dconjg(mat(i,j,2)) ! B^+
    enddo
  enddo

 ! now, piece together A, B, B^+, C for output:
  do i=1,nao
    do j=1,nao
      wr(j,i)=dreal(mat(j,i,1))
      wi(j,i)=dimag(mat(j,i,1))
    enddo
  enddo

  do i=1,nao
    m=1
    do j=nao+1,2*nao
      wr(i,j)=dble(mat(i,m,2))
      wi(i,j)=aimag(mat(i,m,2))
      m=m+1
    enddo
  enddo

  m=1
  do i=nao+1,2*nao
    do j=1,nao
      wr(i,j)=dble(bdag(m,j))
      wi(i,j)=aimag(bdag(m,j))
    enddo
    m=m+1
  enddo

  m=1
  do i=nao+1,2*nao
    n=1
    do j=nao+1,2*nao
      wr(i,j)=dble(mat(m,n,3))
      wi(i,j)=aimag(mat(m,n,3))
      n=n+1
    enddo
    m=m+1
  enddo

  deallocate( bdag )
 end subroutine pbc_expand_wmat

 subroutine pbc_x2c_scale_kinetic(n2,tr,ti)
  implicit none
  integer,intent(in) :: n2
  real*8,intent(inout) :: tr(n2,n2),ti(n2,n2)

  integer :: k,i
  real*8,allocatable :: krs(:,:),kis(:,:)

  allocate( krs(n2,n2),kis(n2,n2) )

  call conjugate_symmetrize_r(n2,tr,krs)
  call conjugate_symmetrize_i(n2,ti,kis)
  tr=krs
  ti=kis

  deallocate( krs,kis )

 end subroutine

 subroutine conjugate_symmetrize_r(n,kr,krs)
   implicit none
   integer,intent(in)::n
   real*8,intent(in) ::kr(n,n)
   real*8,intent(out)::krs(n,n)
   integer :: i,j
   do i=1,n
     do j=1,i
       krs(i,j)=0.5*(kr(i,j)+kr(j,i))
       krs(j,i)=0.5*(kr(i,j)+kr(j,i))
     enddo
   enddo
 end subroutine conjugate_symmetrize_r

 subroutine conjugate_symmetrize_i(n,ki,kis)
    implicit none
   integer,intent(in)::n
   real*8,intent(in) ::ki(n,n)
   real*8,intent(out)::kis(n,n)
   integer :: i,j
   do i=1,n
     do j=1,i
       kis(i,j)=0.5*(ki(i,j)-ki(j,i))
       kis(j,i)=0.5*(ki(j,i)-ki(i,j))
     enddo
   enddo 
 end subroutine conjugate_symmetrize_i


 subroutine reorder_ascending_eigenvec(n,m,e,cr,ci)
 ! Reorder the elements in e: on output, they are in an ascending order.
 ! So to cr and ci.
  implicit none
  integer,intent(in)   :: n,m
  real*8,intent(inout) :: e(m),cr(n,m),ci(n,m)
  integer :: i,j,k
  real*8  :: tmpe,tmpcr(n),tmpci(n)

  do i=1, m
  do j=i+1, m
     if( e(j).lt.e(i) )then
        tmpe     = e(j)
        tmpcr(:) = cr(:,j)
        tmpci(:) = ci(:,j)
        e(j)     = e(i)
        e(i)     = tmpe
        cr(:,j)  = cr(:,i)
        cr(:,i)  = tmpcr(:)
        ci(:,j)  = ci(:,i)
        ci(:,i)  = tmpci(:)
     endif
  enddo 
  enddo

 end subroutine reorder_ascending_eigenvec




 subroutine fn_dftseq_u(r,z,c2,l,alpha,u)
  use physics, only: finite_nuclear_radius
  implicit none
  integer,intent(in) :: l
  real*8,intent(in) :: r,z,alpha
  real*8,intent(in) :: c2  ! 1/c^2
  real*8,intent(out) :: u

  real*8 :: fn_r, gamma
  real*8 :: term1, term2
  ! PGI bug requires this intrinsic function be defined
  real*8 :: derf

  fn_r = finite_nuclear_radius
  if(fn_r.le.0.d0) return

  call fn_dftseq_gamma(r,z,c2,l,fn_r,gamma)

  term1 = z*derf(r/fn_r)*r

  term2 = 0.5d0*(z*r*derf(r/fn_r))**2

  u = r**gamma*(1.d0-term1+term2)

 end subroutine

 subroutine fn_dftseq_up(r,z,c2,l,alpha,up)
  use physics, only: finite_nuclear_radius
  use constants, only: pi
  implicit none
  integer,intent(in) :: l
  real*8,intent(in) :: z,c2,alpha,r
  real*8,intent(out) :: up

  real*8 :: fn_r
  real*8 :: gamma, gammap ! gamma and gamma_prime
  real*8 :: ggp, ggp1, ggp2 ! (gamma^gamma)', (gamma^(gamma+1))', (gamma^(gamma+2))'
  real*8 :: term2, term3, term4, term5
  ! PGI bug requires this intrinsic function be defined
  real*8 :: derf

  fn_r = finite_nuclear_radius
  if(fn_r.le.0.d0) return

  call fn_dftseq_gamma(r,z,c2,l,fn_r,gamma)
  call fn_dftseq_gammap(r,z,c2,l,fn_r,gammap)

  ggp = gamma**gamma * (gammap*dlog(r) + gamma/r)
  ggp1 = gamma**(gamma+1) * (gammap*dlog(r) + (gamma+1)/r)
  ggp2 = gamma**(gamma+2) * (gammap*dlog(r) + (gamma+2)/r)

  term2 = 2*z/fn_r/dsqrt(pi) * dexp(-(r/fn_r)**2) * r**(gamma+1)

  term3 = -z*derf(r/fn_r)*ggp1

  term4 = 2*z**2/fn_r/dsqrt(pi) * derf(r/fn_r) * dexp(-(r/fn_r)**2) * r**(gamma+2)

  term5 = z**2/2*derf(r/fn_r)**2 * ggp2

  up = alpha*r * (ggp - term2 - term3 + term4 + term5)

 end subroutine

 subroutine fn_dftseq_gamma(r,z,c2,l,fn_r,gamma)
  implicit none
  integer,intent(in) :: l
  real*8,intent(in) :: r,z
  real*8,intent(in) :: c2  ! 1/c^2
  real*8,intent(in) :: fn_r ! finite nucleus radius
  real*8,intent(out) :: gamma
  real*8 factor, numerator, denominator
  ! PGI bug requires this intrinsic function be defined
  real*8 :: derf

  denominator = 2*l+1

  factor = z**2 * c2 * derf(r/fn_r)**2
  numerator = l*dsqrt(l**2-factor) + (l+1)*dsqrt((l+1)**2-factor)

  gamma = numerator/denominator

 end subroutine fn_dftseq_gamma

 subroutine fn_dftseq_gammap(r,z,c2,l,fn_r,gammap)
  use constants, only: pi
  implicit none
  integer,intent(in) :: l
  real*8,intent(in) :: r,z,c2,fn_r
  real*8,intent(out) :: gammap
  real*8 constant, factor, bracket, b1,b2, suffix
  ! PGI bug requires this intrinsic function be defined
  real*8 :: derf

  constant = -2*z**2*c2/(2*l+1)/fn_r/dsqrt(pi)

  factor = z**2 * c2 * derf(r/fn_r)
  b1 = l/dsqrt(l**2-factor)
  b2 = (l+1)/dsqrt((l+1)**2-factor)
  bracket = b1 + b2

  suffix = dexp(-(r/fn_r)**2) * derf(r/fn_r)

  gammap = constant * bracket * suffix

 end subroutine fn_dftseq_gammap


 subroutine expand_diagonal(n,mat)
  use localorb_io, only: use_unit
  implicit none
  integer n,m,i,j,k
  real*8 mat(n)
  real*8,allocatable :: mat_sq(:,:)

  m= int(dsqrt(8*n+1.d0)-1.d0+0.1d0)/2
  write(use_unit,*)'n=',n,'   m=',m
  allocate(mat_sq(m,m))

  k=1
  do i=1, m
    do j=1, i
       mat_sq(j,i) = mat(k)
       mat_sq(i,j) = mat(k)
       k=k+1
    enddo
  enddo

  do i=1, m
    write(use_unit,"(20(f11.5))")mat_sq(1:m,i)
  enddo

 end subroutine
 subroutine expand_diagonal_complex(n,mat)
  use localorb_io, only: use_unit
  implicit none
  integer n,m,i,j,k
  complex*16 :: mat(n)
  complex*16,allocatable :: mat_sq(:,:)

  m= int(dsqrt(8*n+1.d0)-1.d0+0.1d0)/2
  write(use_unit,*)'n=',n,'   m=',m
  allocate(mat_sq(m,m))

  k=1
  do i=1, m
    do j=1, i
       mat_sq(j,i) = mat(k)
       mat_sq(i,j) = mat(k)
       k=k+1
    enddo
  enddo

  do i=1, m
    write(use_unit,"(10(f11.5,f11.5,2x))")mat_sq(1:m,i)
  enddo

 end subroutine


