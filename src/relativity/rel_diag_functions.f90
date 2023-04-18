
! Functions for diagonalization of the H matrix.
! Lowdin orthogonalization is involved here.
! For NAO in aims, I found that the 4C and X2C Hamiltonian matrices are not always hermitian,
! thus we should use zgeev instead of zheev in LAPACK.
! -- Rundong Zhao, March. 2019 @ Durham, NC
!
 subroutine pbc_diagmat_SX(sq,Sdim,cut,Sr,Si,Es,Xr,Xi,xdim)
  !------------------------------------------------------------------
  !> Given a symmetric matrix S, find a unitaray matrix that can 
  ! remove linear dependence, -- if there exists.
  ! 
  ! Input:  S - Symmetric Matrix
  ! Output: X
  !
  !    X=S^(-1/2)=Us^(-1/2)U^+, if no linear dependence exists.
  !    X=Us^(-1/2), with U n by m, n>=m. --- if linear dependence exists.
  !   (s means the diagonalized S matrix, which satisfies U^+SU=s)
  !!
  !! the output X satisfies: X^{\dag}SX=I
  !------------------------------------------------------------------
  use localorb_io, only: use_unit
  implicit none
  logical,intent(in) :: sq ! if sq=.true., the input S matrix is an Sdim by Sdim square matrix;
                           ! or it is a triangular matrix
  integer,intent(in) :: Sdim ! dim of input S matrix
  real*8, intent(in) :: cut  ! cutoff value for linear dependence judgement
  real*8, intent(in) :: Sr(*),Si(*) ! real and imag part of S matrix
  real*8, intent(out) :: Es(Sdim),Xr(Sdim,Sdim),Xi(Sdim,Sdim)  ! Eigen value of S matrix; X matrix used for FC=SCE later on
  integer,intent(out) :: xdim ! output dimemsion of X matrix, xdim<=Sdim

  integer :: i,j,k,lwork,info
  real*8 :: dia
  real*8, allocatable :: rwork(:)
  complex*16,allocatable :: zwork(:),U(:,:),tmp(:,:),etmp(:,:),xtmp(:,:)
  complex*16,parameter :: alph = dcmplx(1.d0,0.d0), beta = dcmplx(0.d0,0.d0)
 
  info=0
  lwork=8*Sdim
  allocate(U(Sdim,Sdim),zwork(lwork),rwork(3*Sdim))

  if(sq)then
    k=0
    do i=1,Sdim
    do j=1,Sdim
      k=k+1
      U(j,i)=dcmplx( Sr(k),Si(k) )
    enddo
    enddo
  else
    k=0
    do i=1,Sdim
    do j=1,i
      k=k+1
      U(j,i)=dcmplx( Sr(k),-Si(k) )
      U(i,j)=dconjg( U(j,i) )
    enddo
    enddo
  endif

  call zheev('V','U',Sdim,U,Sdim,Es,zwork,lwork,rwork,info)
  deallocate(zwork,rwork)
  if( info.ne.0 )then
    write(use_unit,'(a,i5)') 'Error in pbc_diagmat_SX: zheev failed, info=',info
    stop
  endif

 ! check linear dependent
  if(cut.gt.0.d0) then
    xdim =0
    do i=1,Sdim
      if(Es(i).ge.cut) then
        xdim=xdim+1
      elseif(Es(i).lt.0.d0)then
        write(use_unit,"(a,e20.12)")'Warning, Eigen value of S matrix < 0 : ',Es(i)
       !stop
      endif
    enddo
  else
    ! if cut.lt.0.d0, keep "xdim" orbitals
    xdim=Sdim
  endif 

  ! Get X matrix:
  if(xdim.eq.Sdim) then
   ! when no linear dependence exists:
    allocate(tmp(Sdim,Sdim),xtmp(Sdim,Sdim))
    call zcopy(sdim*sdim,U,1,tmp,1)
    do i=1,Sdim
      dia=1.d0/dsqrt(dsqrt(Es(i)))
      call zdscal(sdim,dia,tmp(1,i),1)
    enddo
    call zgemm('N','C',Sdim,Sdim,Sdim,alph,tmp,Sdim,tmp,Sdim,beta,xtmp,Sdim)
    do i = 1,sdim
    do j = 1,sdim
      Xr(j,i) = dreal(xtmp(j,i))
      Xi(j,i) = dimag(xtmp(j,i))
    enddo
    enddo
    deallocate(tmp,xtmp)
  else
   ! for linear dependence case:
    allocate(etmp(xdim,xdim),xtmp(Sdim,xdim))
    call zcopy(xdim*sdim,U(1,sdim-xdim+1),1,xtmp,1)
    j=0
    do i=Sdim-xdim+1,Sdim
      j=j+1
      dia=1.d0/dsqrt(dsqrt(Es(i)))
      call zdscal(sdim,dia,xtmp(1,j),1)
    enddo
    do i=1,xdim
      do j=1,Sdim
        Xr(j,i)=dreal(xtmp(j,i))
        Xi(j,i)=dimag(xtmp(j,i))
      enddo
    enddo
    do i=xdim+1,Sdim
      Xr(:,i)=0.d0
      Xi(:,i)=0.d0
    enddo
    deallocate(xtmp)
  endif
  deallocate(U)
 end subroutine pbc_diagmat_SX


 subroutine pbc_diagmat_F(sq,n,m,fr,fi,xr,xi,cr,ci,e)
 !----------------------------------------------------------------------
 !> solve FC=SCE:
 !  FC=SCE ==> X^{\dag}FX C'= C' E
 !  C'=X^{-1}C ==> C=XC'
 ! Non-rel case now.
 !
 ! Input:  F - lower triangular, n by n
 !         X - n by m, m<=n. Here the m+1~n column are 0.
 ! Output: C - Eigen vector, n by m, m<=n. Here, set the m+1~n column to 0.
 !         E - Eigen value of F, the valid dimension is m, but here we fix the
 !               m+1~n orbital with a large-enough value, making sure it will
 !               not be occupied in the following code.
 !----------------------------------------------------------------------
  use localorb_io, only: use_unit
  implicit none
  logical,intent(in)  :: sq ! if sq=.true., the input F matrix is an n by n square NON-hermitian matrix;
                            ! or it is a triangular hermitian matrix
  integer,intent(in)  :: n,m
  real*8, intent(in)  :: fr(*),fi(*),xr(n,*),xi(n,*)
  real*8, intent(out) :: cr(n,n),ci(n,n),e(n)
  complex*16, parameter   :: alph = dcmplx(1.d0,0.d0), beta = dcmplx(0.d0,0.d0)
  complex*16, allocatable :: cz(:,:),xz(:,:),tz(:,:),zvl(:,:),zvr(:,:),ee(:),zwork(:)
  real*8,     allocatable :: rwork(:)
  integer :: lwork,i,j,k,info

  lwork = 8*n
  info  = 0
  allocate( cz(n,n),xz(n,m),tz(n,n),zwork(lwork),rwork(3*n) )
  allocate( zvl(m,m),zvr(m,m),ee(m) )

  if(sq)then
    k=0
    do i=1,n
    do j=1,n
      k=k+1
      cz(j,i) = dcmplx( fr(k),fi(k) )
    enddo
    enddo
  else
    k=0
    do i=1,n
    do j=1,i
      k=k+1
      cz(j,i) = dcmplx( fr(k),-fi(k) )
      cz(i,j) = dconjg( cz(j,i) )
    enddo
    enddo
  endif

  do i=1,m
  do j=1,n
    xz(j,i) = dcmplx( xr(j,i),xi(j,i) ) ! xz, n by m
  enddo
  enddo

 ! F'= X^{\dag}FX :
  call zgemm('C','N',m,n,n,alph,xz,n,cz,n,beta,tz,m) ! tz, m by n now
  call zgemm('N','N',m,m,n,alph,tz,m,xz,n,beta,cz,m) ! cz, m by m now

 ! F'C'=C'E :
  if(sq)then
     call zgeev( 'N', 'V', m, cz, m, ee, zvl,m, zvr,m, zwork, lwork, rwork, info)
  else
     call zheev( 'V', 'U', m, cz, m, e, zwork, lwork, rwork, info)
  endif
  if( info.ne.0 )then
    write(use_unit,'(a,i5)')' diagonalization failed in subroutine pbc_diagmat_F, info=',info
    stop
  endif
 
 ! C=XC' :
  if(sq)then
     call zgemm('N','N',n,m,m,alph,xz,n,zvr,m,beta,tz,n) ! tz, n by m now
  else
     call zgemm('N','N',n,m,m,alph,xz,n,cz,m,beta,tz,n) ! tz, n by m now
  endif
  do i=1,m
  do j=1,n
    cr(j,i) = dreal( tz(j,i) )
    ci(j,i) = dimag( tz(j,i) )
  enddo
  enddo
  do i=m+1,n
    cr(:,i)=0.d0
    ci(:,i)=0.d0
  enddo

  if(sq)then
    e(1:m) = dreal(ee(1:m))
    if(sum(dimag(ee)).gt.1.d-4)&
    write(use_unit,"('Warning: eigenvalues from X2C are complex! Sum of the imaginary part:',f10.5)")sum(dimag(ee))
    call reorder_ascending_eigenvec(n,m,e,cr,ci)
  endif

 ! Set the m+1~n orbital:
  do i=m+1,n
    e(i)=1.d30
  enddo

  deallocate(zwork,rwork,cz,xz,tz,zvl,zvr)
 end subroutine pbc_diagmat_F


 subroutine diagonalization_zgeev(n,s,v,eigenvector,eigenval)
  use localorb_io, only: use_unit
  implicit none
  integer,intent(in) :: n
  real*8,intent(in) :: s(n,n,2), v(n*n,2)
  real*8,intent(out) :: eigenvector(n*n,2), eigenval(n)

  integer :: xdim
  real*8  :: eps
  real*8,allocatable :: smihalf(:,:),e(:)
  integer :: k,l

  eps = 1.d-12
  allocate( smihalf(n*n,2), e(n) )

  call pbc_diagmat_SX(.true., n,eps, s,s(1,1,2), e, smihalf,smihalf(1,2), xdim )
  write(use_unit,"(' Xdim after Lowdin orthorgonalization:',i5)")xdim
  ! solve FC = SCE
  call pbc_diagmat_F(.true., n,xdim, v,v(1,2), smihalf,smihalf(1,2), eigenvector(:,1),eigenvector(:,2),eigenval)

 end subroutine diagonalization_zgeev


