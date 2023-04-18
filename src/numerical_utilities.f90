!****h* FHI-aims/numerical_utilities
!  NAME
!    numerical_utilities
!  SYNOPSIS

module numerical_utilities

  !  PURPOSE
  !    Mainly provide some convenience wrappers around lapack routines
  !    which cope with all the allocation stuff.
  !  USES

  implicit none

  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  character*200, private              :: info_str

contains

  !----------------------------------------------------------------------------
  !****s* numerical_utilities/diagonalize_rmatrix
  !  NAME
  !    diagonalize_rmatrix
  !  SYNOPSIS

  subroutine diagonalize_rmatrix(n, M, D, need_vectors)

    !  PURPOSE
    !
    !    Diagonalize M.
    !
    !    If (need_vectors), A being the input M, then
    !
    !        A M = M D
    !
    !    where D is the diagonal matrix with D(i) as i-th entry.
    !
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    integer, intent(IN)                 :: n
    real*8, intent(INOUT)               :: M(n, n)
    real*8, intent(OUT)                 :: D(n)
    logical, intent(IN)                 :: need_vectors

    !  INPUTS
    !    o n -- dimension
    !    o M -- matrix to diagonalize
    !    o need_vectors -- flag if eigenvectors are needed
    !  OUTPUTS
    !    o M -- eigenvectors in case of (need_vectors), garbage otherwise
    !    o D -- eigenvalues
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character(1)                        :: JOBZ
    character(1), parameter             :: UPLO = 'U'
    integer                             :: LDA
    integer                             :: LWORK
    real*8,    allocatable              :: WORK(:)
    integer                             :: INFO
    real*8                              :: WORK_tmp(1)
    external                               DSYEV
    character(*), parameter             :: func = 'diagonalize_rmatrix'

    LDA = n
    if (need_vectors) then
       JOBZ = 'V'
    else
       JOBZ = 'N'
    end if
    
    LWORK = -1
    call DSYEV(JOBZ, UPLO, n, M, LDA, D, WORK_tmp, LWORK, INFO)
    LWORK = nint(WORK_tmp(1))
    allocate(WORK(LWORK))

    ! DSYEV(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)
    call DSYEV(JOBZ, UPLO, n, M, LDA, D, WORK, LWORK, INFO)

    if (INFO /= 0) then
       write(info_str, "('Diagonalization failed; info =',I6)") info
       call aims_stop(info_str, func)
    end if

    deallocate(WORK)

  end subroutine diagonalize_rmatrix
  !******
  !----------------------------------------------------------------------------
  !****s* numerical_utilities/solve_LEQ
  !  NAME
  !    solve_LEQ
  !  SYNOPSIS

  subroutine solve_LEQ(caller, N, A, b, x)

    !  PURPOSE
    !    Solve A x == b, keeping A intact.
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    character(*), intent(IN)            :: caller
    integer, intent(IN)                 :: N
    real*8, intent(IN)                  :: A(N,N), b(N)
    real*8, intent(OUT)                 :: x(N)

    !  INPUTS
    !    o N -- matrix dimension
    !    o A -- N x M input matrix
    !    o b -- right hand side
    !  OUTPUTS
    !    o x -- solution
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer                             :: NRHS, LDA, LDx, INFO
    integer, allocatable                :: IPIV(:)
    real*8, allocatable                 :: A_tmp(:,:)
    external                               DGESV
    character(*), parameter             :: func = 'solve_LEQ'

    NRHS = 1
    LDA = N
    LDx = N
    allocate(IPIV(N), A_tmp(N,N))

    x = b
    A_tmp = A
    call DGESV(N, NRHS, A_tmp, LDA, IPIV, x, LDx, INFO)
    if (INFO /= 0) then
       write(info_str, "('For ',A,': DGESV info =',I5)") trim(caller), INFO
       call aims_stop(info_str)
    end if

    deallocate(IPIV, A_tmp)

  end subroutine solve_LEQ
  !******
  !----------------------------------------------------------------------------
  !****s* numerical_utilities/solve_LSQ
  !  NAME
  !    solve_LSQ
  !  SYNOPSIS

  subroutine solve_LSQ(caller, M, N, A, b, x, rank, sing_thres)

    !  PURPOSE
    !    Solves A.x = b in a least squares sense (argmin_x |Ax-b|^2).
    !      => x(i): how often do I need A(:,i) in b(:)?
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS
    character(*), intent(IN)            :: caller
    integer, intent(IN)                 :: M, N
    real*8, intent(IN)                  :: A(M,N)
    real*8, intent(IN)                  :: b(M)
    real*8, intent(OUT)                 :: x(N)
    integer, intent(OUT), optional      :: rank
    real*8, intent(IN), optional        :: sing_thres

    !  INPUTS
    !    o caller -- calling procedure (for error messages)
    !    o M, N -- matrix dimensions
    !    o A -- M x N matrix
    !    o b -- right hand side
    !    o sing_thres (optional)
    !        -- Assume singularity for singular values below sing_thres.
    !           Default: 1d-10.
    !  OUTPUTS
    !    o x -- least squares solution: argmin_x |Ax-b|^2
    !    o rank (optional) -- number of singular values larger than sing_thres
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer                             :: NRHS, LDA, LDB
    real*8, allocatable                 :: A_tmp(:,:), B_tmp(:)
    integer                             :: INFO
    integer, allocatable                :: JPVT(:)
    integer                             :: myRANK
    real*8                              :: RCOND
    real*8, allocatable                 :: WORK(:)
    real*8                              :: WORK_tmp(1)
    integer                             :: LWORK
    character(*), parameter             :: func = 'solve_LSQ'
    external                               DGELSY

    allocate(A_tmp(M, N), B_tmp(max(M, N)))
    A_tmp = A
    B_tmp = 0.d0
    B_tmp(1:M) = b
    NRHS = 1
    LDA = M
    LDB = max(M, N)
    allocate(JPVT(N))
    JPVT = 0
    RCOND = 1d-10; if (present(sing_thres)) RCOND = sing_thres

    ! --- work space query

    LWORK = -1
    call DGELSY(M, N, NRHS, A_tmp, LDA, B_tmp, LDB, JPVT, RCOND, myRANK, &
    &           WORK_tmp, LWORK, INFO)
    if (INFO /= 0) then
       write(info_str, "('For ',A,': DGELSY workspace info =',I5)") &
       & trim(caller), INFO
       call aims_stop(info_str)
    end if
    LWORK = nint(WORK_tmp(1))
    allocate(WORK(LWORK))


    ! --- call

    call DGELSY(M, N, NRHS, A_tmp, LDA, B_tmp, LDB, JPVT, RCOND, myRANK, &
    &           WORK, LWORK, INFO)
    if (INFO /= 0) then
       write(info_str, "('For ',A,': DGELSY info =',I5)") trim(caller), INFO
       call aims_stop(info_str)
    end if

    ! --- test

    x = B_tmp(1:N)
    deallocate(A_tmp, B_tmp, JPVT, WORK)
    if (present(rank)) rank = myRANK

  end subroutine solve_LSQ
  !******
  !----------------------------------------------------------------------------
  !****s* numerical_utilities/pseudo_inverse
  !  NAME
  !    pseudo_inverse
  !  SYNOPSIS

  subroutine pseudo_inverse(caller, M, N, A, Ainv, eps, rank)
    !  PURPOSE
    !    Calculate the pseudoinverse A^+ of A with the properties
    !      A A^+ A == A, A^+ A A^+ == A, (A A^+)^* == A A^+, (A^+ A)^* == A A^+
    !    as described in
    !    http://en.wikipedia.org/wiki/Moore-Penrose_pseudoinverse
    !    The Pseudoinverse has the property that
    !      x := A^+ b
    !    gives the 'least squares' solution of the problem
    !      A b == x
    !    in the sense that |Ax - b|_2 and |x|_2 are both minimal.
    !  USES
    use localorb_io, only: use_unit
    use mpi_tasks, only: myid, aims_stop
    implicit none
    character(*), intent(IN)            :: caller
    real*8, intent(IN)                  :: A(M,N)
    real*8, intent(OUT)                 :: Ainv(N,M)
    real*8, intent(IN)                  :: eps
    integer, intent(OUT), optional      :: rank

    !  INPUTS
    !    o caller -- calling procedure (for error messages)
    !    o M, N -- matrix dimensions
    !    o A -- M x N matrix
    !    o eps -- threshold for small eigenvalues
    !  OUTPUTS
    !    o Ainv -- right hand side
    !    o rank (optional) -- number of nonzero singular values
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer                             :: M, N, MNmax, MNmin, LDA, LDU, LDVT
    real*8, allocatable                 :: A_tmp(:,:)
    real*8, allocatable                 :: U(:,:), VT(:,:), S(:)
    integer                             :: INFO
    real*8, allocatable                 :: WORK(:)
    real*8                              :: WORK_tmp(1)
    integer                             :: LWORK
    character(*), parameter             :: JOBU = 'A' ! All columns of U
    character(*), parameter             :: JOBVT = 'A' ! All rows of VT
    character(*), parameter             :: func = 'pseudo_inverse'
    character*150                       :: info_str
    integer                             :: i, j
    external                               dgesvd, dgemm

    ! --- prepare

    MNmax = max(M, N)
    MNmin = min(M, N)
    allocate(A_tmp(M, N), U(M, M), VT(N, N), S(MNmin))
    A_tmp = A
    LDA = M
    LDU = M
    LDVT = N

    ! --- work space query

    LWORK = -1
    call dgesvd(JOBU, JOBVT, M, N, A_tmp, LDA, S, U, LDU, VT, LDVT, &
    &           WORK_tmp, LWORK, INFO)
    if (INFO /= 0) then
       write(info_str, "('For ',A,': DGESVD workspace info =',I5)") &
       & trim(caller), INFO
       call aims_stop(info_str)
    end if
    LWORK = nint(WORK_tmp(1))
    allocate(WORK(LWORK))

    ! --- call

    call dgesvd(JOBU, JOBVT, M, N, A_tmp, LDA, S, U, LDU, VT, LDVT, &
    &           WORK, LWORK, INFO)
    if (INFO /= 0) then
       write(info_str, "('For ',A,': DGESVD workspace info =',I5)") &
       & trim(caller), INFO
       call aims_stop(info_str)
    end if

    ! --- invert

    if (present(rank)) rank = 0
    do i = 1, MNmin
       if (S(i) > eps .and. present(rank)) rank = rank + 1
       ! S(i) = S(i) / (S(i)**2 + eps**2)
       if (S(i) > eps) then
          S(i) = 1.d0 / S(i)
       else
          S(i) = 0.d0
       end if
    end do

    ! --- calculate Ainv

    ! A_tmp = ( V.S^+ )^T = S^+ . VT
    do j = 1, N
       do i = 1, MNmin
          A_tmp(i, j) = S(i) * VT(i, j)
       end do
       do i = MNmin+1, M
          A_tmp(i, j) = 0.d0
       end do
    end do

    ! Ainv = A_tmp^T . U^T
    ! Ainv = matmul(transpose(A_tmp), transpose(U))
    call dgemm('T', 'T', N, M, M, 1.d0, A_tmp, M, U, M, 0.d0, Ainv, N)

    ! --- test

    A_tmp = matmul(A, matmul(Ainv, A))
    if (maxval(A_tmp - A) > 2*MNmax*eps) then
       write(info_str, &
       &     "('*** ',A,', for ',A,': error of',ES10.2,'; eps:',ES10.2)") &
       & trim(func), trim(caller), maxval(A_tmp - A), eps
       write(use_unit,"(2X,'Task',I8,': ',A)") myid, trim(info_str)
       ! call aims_stop(info_str)
    end if

  end subroutine pseudo_inverse
  !******
  !----------------------------------------------------------------------------
  !****s* numerical_utilities/orthonormalize_by_QR
  !  NAME
  !    orthonormalize_by_QR
  !  SYNOPSIS

  subroutine orthonormalize_by_QR(caller, M, N, A)

    !  PURPOSE
    !    Orthonormalize columns of A by QR (without column pivoting).
    !    BEWARE:  This routine is not thouroughly tested!
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    character(*), intent(IN)            :: caller
    integer, intent(IN)                 :: M, N
    real*8, intent(INOUT)               :: A(M,N)

    !  INPUTS
    !    o N, M -- matrix dimensions, M >= N
    !    o A -- N x M input matrix
    !  OUTPUTS
    !    o A -- matrix with orthonormal columns
    !           A_out(:,i) is a linear combination of input A_in(:,1:i).
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8, allocatable                 :: TAU(:), WORK(:)
    integer                             :: LDA, INFO, LWORK, K
    real*8                              :: WORK_tmp(1)
    external                               DGEQR2, DORGQR
    character(*), parameter             :: func = 'orthonormalize_by_QR'

    ! --- preparation

    LDA = M
    K = N ! == min(N, M)
    allocate(TAU(K))

    ! --- decomposition

    allocate(WORK(N))
    call DGEQR2(M, N, A, LDA, TAU, WORK, INFO)
    deallocate(WORK)
    if (INFO /= 0) then
       write(info_str, "('For ',A,': dgeqr2 info =',I5)") trim(caller), INFO
       call aims_stop(info_str)
    end if

    ! --- Q construction

    ! workspace query
    LWORK = -1
    call DORGQR(M, N, K, A, LDA, TAU, WORK_tmp, LWORK, INFO)
    if (INFO /= 0) then
       write(info_str, "('For ',A,': dorgqr workspace info =',I5)") &
       & trim(caller), INFO
       call aims_stop(info_str)
    end if
    LWORK = nint(WORK_tmp(1))
    ! actual call
    allocate(WORK(LWORK))
    call DORGQR(M, N, K, A, LDA, TAU, WORK, LWORK, INFO)
    if (INFO /= 0) then
       write(info_str, "('For ',A,': dorgqr info =',I5)") trim(caller), INFO
       call aims_stop(info_str)
    end if
    deallocate(WORK)

    ! --- tidy up

    deallocate(TAU)

  end subroutine orthonormalize_by_QR
  !******
  !----------------------------------------------------------------------------
  !****s* numerical_utilities/twonorm
  !  NAME
  !    twonorm
  !  SYNOPSIS
  real*8 function twonorm(v)
    !  PURPOSE
    !    Calculate the two-norm sqrt(sum(v**2)) of vector v.
    !  USES
    implicit none
    !  ARGUMENTS
    real*8, intent(IN)                  :: v(:)
    !  INPUTS
    !    o v -- vector
    !  OUTPUTS
    !    o twonorm -- ||v||_2
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE
    twonorm = sqrt(dot_product(v, v))
  end function twonorm
  !******
end module numerical_utilities
!******
