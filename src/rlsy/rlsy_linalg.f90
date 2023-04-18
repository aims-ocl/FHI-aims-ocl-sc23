module rlsy_linalg
!! Basic geometric linear algebra functions neatly wrapped.
!! Not meant to be fast or scalable, only convenient.
use rlsy_constants, only: r8, rl_exitcode_blaslapack, rl_exitcode_param, rl_status
use rlsy_helpers, only: tochar,rl_chop
use rlsy_mpi_helper, only: rl_stop_gracefully

implicit none
private

public :: rl_general_real_eigenvalues_eigenvectors
public :: rl_invert_real_matrix
public :: rl_triplegemm
public :: rl_real_singular_value_decomposition
public :: rl_real_pseudoinverse
public :: rl_symmetric_eigensystem_3x3matrix

contains

!> get left and right eigenvectors + eigenvalues of a general matrix
subroutine rl_general_real_eigenvalues_eigenvectors(A,eigenvalues,vec_left,vec_right)
    !> matrix to check
    real(r8), dimension(:,:), intent(in) :: A
    !> eigenvalues
    complex(r8), dimension(:), intent(out) :: eigenvalues
    !> left and right eigenvectors. left are (:,i), right are also (:,i)
    real(r8), dimension(:,:), intent(out) :: vec_left,vec_right

    integer :: n

    !solve: block
        real(r8), dimension(:,:), allocatable :: rA
        real(r8), dimension(:), allocatable :: rwork,wr,wi
        real(r8), dimension(1) :: drwork
        integer :: lwork

    ! Set some stuff
    !init: block
        ! Set size of problem
        if ( size(A,1) .ne. size(A,2) ) then
            call rl_stop_gracefully(['Can only get eigenvectors of square matrices'],&
            rl_exitcode_param)
        else
            n=size(A,1)
        endif
        vec_left=0.0_r8
        vec_right=0.0_r8
        eigenvalues=0.0_r8
    !end block init

    ! First do the actual solution
    !solve: block

        ! size of problem and space for stuff
        ! solve original eigenproblem, temporary space
        allocate(rA(n,n))
        allocate(wr(n))
        allocate(wi(n))
        rA=0.0_r8
        wr=0.0_r8
        wi=0.0_r8
        lwork=-1
        call dgeev('V','V', n, rA, n, wr, wi, vec_left, n, vec_right, n, drwork, lwork, rl_status)
        if ( rl_status .ne. 0 ) then
            call rl_stop_gracefully(['dgeev exit code '//tochar(rl_status)],&
            rl_exitcode_blaslapack)
        endif
        lwork=int(drwork(1))
        allocate(rwork(lwork))
        rwork=0.0_r8
        rA=A
        call dgeev('V','V', n, rA, n, wr, wi, vec_left, n, vec_right, n, rwork, lwork, rl_status)
        if ( rl_status .ne. 0 ) then
            call rl_stop_gracefully(['dgeev exit code '//tochar(rl_status)],&
            rl_exitcode_blaslapack)
        endif
        eigenvalues=cmplx(wr,wi,r8)
        deallocate(wr,wi,rA,rwork)
    !end block solve
end subroutine

!> Matrix inversion.
subroutine rl_invert_real_matrix(A,iA)
    !> matrix to invert
    real(r8), dimension(:,:), intent(in) :: A
    !> inverse of matrix
    real(r8), dimension(:,:), intent(out) :: iA

    real(r8), dimension(:), allocatable :: work,ipiv
    integer :: n

    n=size(A,1)
    ! Some workspace
    allocate(work(n))
    allocate(ipiv(n))
    iA=A
    work=0.0_r8
    ipiv=0
    call dgetrf(n,n,iA,n,ipiv,rl_status)
    if ( rl_status .ne. 0 ) then
       call rl_stop_gracefully(['dgetrf exit code '//tochar(rl_status)],&
       rl_exitcode_blaslapack)
    endif
    call dgetri(n,iA,n,ipiv,work,n,rl_status)
    if ( rl_status .ne. 0 ) then
       call rl_stop_gracefully(['dgetri exit code '//tochar(rl_status)],&
       rl_exitcode_blaslapack)
    endif
    deallocate(work)
    deallocate(ipiv)
end subroutine

!> calculates D = A*B*C in one go
subroutine rl_triplegemm(A,B,C,D,transa,transb,transc)
    !> matrices
    real(r8), dimension(:,:), intent(in) :: A,B,C
    !> output
    real(r8), dimension(:,:), intent(out) :: D
    !> transpose stuff?
    character(len=1), intent(in), optional :: transa,transb,transc

    real(r8), dimension(:,:), allocatable :: BC
    character(len=1) :: ta,tb,tc

    integer :: lda,ldb,ldc,ldd,ldbc
    integer :: ai,aj,bi,bj,ci,cj,di,dj

    ! Sort out all the ways things can be transposed
    if ( present(transa) ) then
        ta=transa
    else
        ta='N'
    endif
    if ( present(transb) ) then
        tb=transb
    else
        tb='N'
    endif
    if ( present(transc) ) then
        tc=transc
    else
        tc='N'
    endif

    ! Then all the funny dimensions, depending on how stuff got transposed
    if ( ta .eq. 'T' ) then
        ai=size(A,2)
        aj=size(A,1)
    else
        ai=size(A,1)
        aj=size(A,2)
    endif
    if ( tb .eq. 'T' ) then
        bi=size(B,2)
        bj=size(B,1)
    else
        bi=size(B,1)
        bj=size(B,2)
    endif
    if ( tc .eq. 'T' ) then
        ci=size(C,2)
        cj=size(C,1)
    else
        ci=size(C,1)
        cj=size(C,2)
    endif
    di=size(D,1)
    dj=size(D,2)

    ! Now I know enough to make temporary space
    allocate(BC(bi,cj))
    BC=0.0_r8
    lda=size(A,1)
    ldb=size(B,1)
    ldc=size(C,1)
    ldbc=size(BC,1)
    ldd=size(D,1)

    ! Begin operation actual operation
    call dgemm(tb,tc,bi,cj,bj,1.0_r8,B,ldb,C,ldc,0.0_r8,BC,ldbc)
    call dgemm(ta,'N',di,dj,aj,1.0_r8,A,lda,BC,ldbc,0.0_r8,D,ldd)
    deallocate(BC)
end subroutine

!> calculate the pseudoinverse via SVD
subroutine rl_real_pseudoinverse(A,B,tolerance)
    !> matrix to pseudoinvert, n x m
    real(r8), dimension(:,:), intent(in) :: A
    !> pseudoinverse, m x n
    real(r8), dimension(:,:), intent(out) :: B
    !> optional tolerance
    real(r8), intent(in), optional :: tolerance

    real(r8), dimension(:,:), allocatable :: U,V,dm
    real(r8), dimension(:), allocatable :: S
    real(r8) :: tol
    integer :: i,n,m

    ! No idea how I chose the default here.
    if ( present(tolerance) ) then
        tol=tolerance
    else
        tol=1E-12_r8
    endif
    n=size(A,1)
    m=size(A,2)
    allocate(dm(m,n))
    dm=0.0_r8
    call rl_real_singular_value_decomposition(A,S,U,V)
    ! Not sure about tolerance here. Think it's fine.
    do i=1,size(S)
        if ( abs(S(i)) .gt. tol ) then
            S(i)=1.0_r8/S(i)
        else
            S(i)=0.0_r8
        endif
    enddo
    do i=1,size(S)
        dm(i,:)=S(i)*U(:,i)
    enddo
    call dgemm('T','N',m,n,m,1.0_r8,V,m,DM,m,0.0_r8,B,m)
    B=rl_chop( B, 1E-13_r8 )
    ! Check that it actually worked? Nah. Cleanup instead.
    deallocate(dm)
    deallocate(S)
    deallocate(U)
    deallocate(V)
end subroutine

!> wrapper around gesvd. Will add things to return the singular vectors and that stuff when I need it.
subroutine rl_real_singular_value_decomposition(A,S,U,V)
    !> matrix A
    real(r8), dimension(:,:), intent(in) :: A
    !> singular values
    real(r8), dimension(:), allocatable, intent(out) :: S
    !> left singular vectors
    real(r8), dimension(:,:), allocatable, intent(out), optional :: U
    !> right singular vectors
    real(r8), dimension(:,:), allocatable, intent(out), optional :: V

    real(r8), dimension(:,:), allocatable :: Sigma
    real(r8), dimension(:), allocatable :: wrk
    real(r8), dimension(1,1) :: dV,dU
    real(r8), dimension(1) :: dwrk
    integer :: nc,nr,ns,lwork

    nc=size(A,1)
    nr=size(A,2)
    ns=min(nc,nr)
    allocate(Sigma(nc,nr))
    allocate(S(ns))
    Sigma=A
    S=0.0_r8
    lwork=-1

    if ( present(U) .and. present(V) ) then
        ! Return left and right vectors
        allocate(U(nc,nc))
        allocate(V(nr,nr))
        U=0.0_r8
        V=0.0_r8
        call dgesvd('A','A',nc,nr,Sigma,nc,S,U,nc,V,nr,dwrk,lwork,rl_status)
        if ( rl_status .ne. 0 ) call rl_stop_gracefully(['dgesvd exit code '//tochar(rl_status)],rl_exitcode_blaslapack)
        lwork=int(anint(dwrk(1)))
        allocate(wrk(lwork))
        call dgesvd('A','A',nc,nr,Sigma,nc,S,U,nc,V,nr,wrk,lwork,rl_status)
        if ( rl_status .ne. 0 ) call rl_stop_gracefully(['dgesvd exit code '//tochar(rl_status)],rl_exitcode_blaslapack)
    else
        ! No vectors returned
        call dgesvd('N','N',nc,nr,Sigma,nc,S,dU,nc,dV,1,dwrk,lwork,rl_status)
        if ( rl_status .ne. 0 ) call rl_stop_gracefully(['dgesvd exit code '//tochar(rl_status)],rl_exitcode_blaslapack)
        lwork=int(anint(dwrk(1)))
        allocate(wrk(lwork))
        call dgesvd('N','N',nc,nr,Sigma,nc,S,dU,nc,dV,1,wrk,lwork,rl_status)
        if ( rl_status .ne. 0 ) call rl_stop_gracefully(['dgesvd exit code '//tochar(rl_status)],rl_exitcode_blaslapack)
    endif
    deallocate(wrk)
    deallocate(Sigma)
end subroutine

!> Solve 3x3 symmetric eigenproblem. Quite fast.
subroutine rl_symmetric_eigensystem_3x3matrix(matrix,eigenvalues,eigenvectors)
    !> matrix
    real(r8), dimension(3,3), intent(in) :: matrix
    !> eigenvalues
    real(r8), dimension(3), intent(out) :: eigenvalues
    !> eigenvectors
    real(r8), dimension(3,3), intent(out) :: eigenvectors

    real(r8), dimension(3) :: e

    !tridiag: block
        real(r8), dimension(3) :: u,vp
        real(r8) :: omega,f,k,h,g
        integer :: i,j
    !ql: block
        !real(r8) :: g, r, p, f, b, s, c, t
        real(r8) :: r, p, b, s, c, t
        !integer :: l, m, i, j, k
        integer :: l, m, ik

    !tridiag: block
        ! initialize q to the identitity matrix
        eigenvectors=reshape([1.0_r8,0.0_r8,0.0_r8,0.0_r8,1.0_r8,0.0_r8,0.0_r8,0.0_r8,1.0_r8],[3,3])
        ! bring first row and column to the desired form
        h = matrix(1,2)**2 + matrix(1,3)**2
        if ( matrix(1,2) .gt. 0.0_r8 ) then
            g = -sqrt(h)
        else
            g = sqrt(h)
        endif
        e(1)  = g
        f     = g * matrix(1,2)
        u(2)  = matrix(1,2) - g
        u(3)  = matrix(1,3)
        omega = h - f
        if ( omega .gt. 0.0_r8 ) then
            omega = 1.0_r8/omega
            k     = 0.0d0
            do i=2,3
                f = matrix(2,i)*u(2) + matrix(i,3)*u(3)
                vp(i) = omega * f
                k    = k + u(i) * f
            enddo
            k=0.5_r8*k*omega**2
            do i=2,3
                vp(i)=vp(i)-k*u(i)
            enddo
            eigenvalues(1) = matrix(1,1)
            eigenvalues(2) = matrix(2,2) - 2.0_r8 * vp(2) * u(2)
            eigenvalues(3) = matrix(3,3) - 2.0_r8 * vp(3) * u(3)
            ! store inverse householder transformation in q
            do j = 2, 3
                f=omega*u(j)
                do i = 2, 3
                    eigenvectors(i,j) = eigenvectors(i,j) - f * u(i)
                enddo
            enddo
            ! calculated updated a(2, 3) and store it in e(2)
            e(2) = matrix(2,3) - vp(2)*u(3) - u(2)*vp(3)
        else
            do i=1,3
                eigenvalues(i) = matrix(i,i)
                e(2) = matrix(2,3)
            enddo
        endif
    !end block tridiag

    ! Calculate eigensystem of the remaining real symmetric tridiagonal
    ! matrix with the QL method
    !ql: block

        ! Loop over all off-diagonal elements
        outloop: do l=1,2
        inloop: do i=1,50
            ! Check for convergence and exit iteration loop if off-diagonal
            ! element E(L) is zero
            do m=l,2
                g = abs(eigenvalues(m)) + abs(eigenvalues(m+1))
                if ( abs(e(m)) + g .eq. g ) exit
            enddo
            if ( m .eq. l ) cycle outloop
            ! Calculate G = D(M) - K
            g = (eigenvalues(l+1) - eigenvalues(l)) / (2.0_r8 * e(l))
            r = sqrt(1.0_r8 + g**2)
            if ( g .ge. 0.0_r8 ) then
                g = eigenvalues(m) - eigenvalues(l) + e(l)/(g + r)
            else
                g = eigenvalues(m) - eigenvalues(l) + e(l)/(g - r)
            end if
            S = 1.0_r8
            C = 1.0_r8
            P = 0.0_r8
            do j=m-1,l,-1
                f = s*e(j)
                b = c*e(j)
                if ( abs(f) .gt. abs(g) ) then
                    c      = g / f
                    r      = sqrt(1.0_r8 + c**2)
                    e(j+1) = f * r
                    s      = 1.0_r8/r
                    c      = c * s
                else
                    s      = f / g
                    r      = sqrt(1.0_r8 + s**2)
                    e(j+1) = g * r
                    c      = 1.0_r8 / r
                    s      = s * c
                end if
                g      = eigenvalues(j+1) - p
                r      = (eigenvalues(j) - g) * s + 2.0_r8 * c * b
                p      = s * r
                eigenvalues(j+1) = g + p
                g      = c * r - b
                ! Form eigenvectors
                do ik=1,3
                    t        = eigenvectors(ik,j+1)
                    eigenvectors(ik,j+1) = s*eigenvectors(ik,j)+c*t
                    eigenvectors(ik,j)   = c*eigenvectors(ik,j)-s*t
                enddo
            enddo
            eigenvalues(l) = eigenvalues(l) - p
            e(l) = g
            e(m) = 0.0_r8
        enddo inloop
        enddo outloop
    !end block ql
end subroutine

end module
