!****h* FHI-aims/linalg
!  NAME
!    linalg
!  SYNOPSIS

module linalg

  !  PURPOSE
  !    Provide generic methods for typical linear algebra tasks.    
  !
  implicit none
  !  AUTHOR
  !    Christoph Schober
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications (2008), submitted.
  !  HISTORY
  !    First version (2015) 
  !  SOURCE

  interface simple_matmul
      module procedure simple_matmul_real_lapack
      module procedure simple_matmul_real_scalapack
      module procedure simple_matmul_cmplx_lapack
      module procedure simple_matmul_cmplx_scalapack
  end interface

private

public :: simple_matmul
public :: lowdin_lapack
public :: diagonalize_rmatrix_elpa
public :: lowdin_elpa
public :: get_rinverse_lapack
public :: get_rinverse_elpa


contains

  !****s* linalg/lowdin_lapack
  !  NAME
  !    lowdin_lapack
  !  SYNOPSIS
  subroutine lowdin_lapack(M, X)
    !  PURPOSE
    !  Do a symmetric Lowdin-Orthogonalisation of matrix M
    !     U**T*(s**-0.5)*U = S**-0.5 = X
    !     1) s, U = diag(S)
    !     2) s**-0.5 = diag_mat(s**-0.5)
    !     3) us = dgemm(U, s**-0.5)
    !     4) X = S**-0.5 = dgemm(us, U)

    !  USES
    use numerical_utilities, only: diagonalize_rmatrix
    implicit none
    !  ARGUMENTS
    real*8, dimension(:,:), intent(IN) :: M 
    real*8, dimension(:,:), intent(INOUT) :: X
    !  INPUTS
    !    o M -- matrix
    !  OUTPUTS
    !    o X -- Lowdin transformation matrix 
    !  AUTHOR
    !    Christoph Schober 
    !  HISTORY
    !    initial version, 2015

    !  SOURCE 
    
    integer :: n_size, i_size
    real*8, dimension(:), allocatable :: s
    real*8, dimension(:,:), allocatable :: U, s_diag

    
    real*8, dimension(:,:), allocatable :: tmp
    
    n_size = size(M, 1)

    allocate(s(n_size))
    allocate(U(n_size, n_size))

    ! need to have a copy of M since diagonalize_rmatrix overwrites
    ! M with the eigenvector
    U = M
    call diagonalize_rmatrix(n_size, U, s, .true.)
    
    ! get s**-0.5-matrix
    s = s**(-0.5)     

    allocate(s_diag(n_size, n_size))
    s_diag = 0.d0
    do i_size = 1, n_size
        s_diag(i_size,i_size) = s(i_size)
    end do

    ! get X_mat = U*s**-0.5*U**T
    allocate(tmp(n_size, n_size))
    call dgemm('N','N', n_size,n_size,n_size, 1.d0, &
           U,ubound(U,1), &
           s_diag,ubound(s_diag,1), &
           0.d0, tmp,ubound(tmp,1))
    !dgemm function for U*s_diag*U**T
    call dgemm('N','T', n_size,n_size,n_size, 1.d0, &
           tmp,ubound(tmp,1), &
           U,ubound(U,1), &
           0.d0, X, ubound(X,1))

end subroutine lowdin_lapack

  !****s* linalg/diagonalize_rmatrix_elpa
  !  NAME
  !    diagonalize_rmatrix_elpa
  !  SYNOPSIS
  subroutine diagonalize_rmatrix_elpa(n_size, M, D)
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
    !  The diagonalization is done using the ELPA solve_evp_real routine.
    !  USES
    use elpa2_2013, only: solve_evp_real_2stage_2013
    use scalapack_wrapper, only: mxld, mxcol, nb, mpi_comm_cols, mpi_comm_rows, set_full_matrix_real
    use mpi_tasks, only: check_allocation, mpi_comm_global

    implicit none

    !  ARGUMENTS
    integer, intent(IN)                 :: n_size
    real*8, intent(INOUT)               :: M(mxld, mxcol)
    real*8, intent(OUT)                 :: D(n_size)

    !  INPUTS
    !    o n_size -- dimension of (undistributed) matrix (e.g. ovlp = n_basis)
    !    o M -- matrix to diagonalize
    !  OUTPUTS
    !    o M -- eigenvectors
    !    o D -- eigenvalues
    !  AUTHOR
    !    Christoph Schober
    !  HISTORY
    !    First version, 2015 
    !  SOURCE
    integer, dimension(:), allocatable :: iwork

    real*8, dimension(:,:), allocatable :: tmp, scalapack_matrix, matrix_eigenvec

    integer :: info, i_row, i_col
    real*8 :: ev_sqrt

    allocate(tmp(mxld, mxcol),stat=info)
    call check_allocation(info, 'tmp')

    tmp(:,:) = -M(:,:)
    call set_full_matrix_real(tmp) 
    call solve_evp_real_2stage_2013(n_size, n_size, tmp, mxld, D, M, mxld, &
         nb, mpi_comm_rows, mpi_comm_cols, mpi_comm_global)

    D(1:n_size) = -D(1:n_size)

end subroutine diagonalize_rmatrix_elpa

  !****s* linalg/lowdin_elpa
  !  NAME
  !    lowdin_elpa
  !  SYNOPSIS
  subroutine lowdin_elpa(M, X, n_size)
    !  PURPOSE
    !  Do a symmetric Lowdin-Orthogonalisation of scalapack matrix M 
    !     U**T*(s**-0.5)*U = S**-0.5 = X
    !     1) s, U = diag(S)
    !     2) s**-0.5 = diag_mat(s**-0.5)
    !     3) us = dgemm(U, s**-0.5)
    !     4) X = S**-0.5 = dgemm(us, U)
    !  
    !  USES
    use scalapack_wrapper, only: mxld, mxcol, nb, mpi_comm_cols, mpi_comm_rows, set_full_matrix_real, sc_desc, setup_scalapack_rmatrix
    use mpi_tasks, only: check_allocation
    implicit none
    !  ARGUMENTS
    real*8, dimension(:,:), intent(IN) :: M 
    real*8, dimension(:,:), intent(INOUT) :: X
    integer, intent(IN) :: n_size 
    !  INPUTS
    !    o M -- matrix
    !    o n_size -- dimension of (undistributed) matrix (e.g. ovlp = n_basis)
    !  OUTPUTS
    !    o X -- Lowdin transformation matrix 
    !  AUTHOR
    !    Christoph Schober 
    !  HISTORY
    !    initial version, 2015

    !  SOURCE 
    integer :: i_size, n
    real*8, dimension(:), allocatable :: s
    real*8, dimension(:,:), allocatable :: U, s_diag, s_sca
    
    real*8, dimension(:,:), allocatable :: tmp

    allocate(U(mxld, mxcol))
    allocate(s(n_size))

    ! need to have a copy of M since diagonalize_rmatrix overwrites
    ! M with the eigenvector
    U = M

    call diagonalize_rmatrix_elpa(n_size, U, s)
    
    ! get s**-0.5-matrix
    ! missing: scaLAPACKify with more intelligent loop to directly set s**-0.5 in distributed storage
    allocate(s_diag(n_size, n_size))
    
    s = s**(-0.5)
    
    s_diag = 0.d0
    do i_size = 1, n_size
        s_diag(i_size,i_size) = s(i_size)
    end do
        
    deallocate(s)

    allocate(s_sca(mxld, mxcol))
    s_sca = 0.d0
    call setup_scalapack_rmatrix(s_diag, s_sca)
    deallocate(s_diag)

    ! get X_mat = U*s**-0.5*U**T

    allocate(tmp(mxld, mxcol))

    call pdgemm('N', 'N', n_size, n_size, n_size, 1.d0, U, 1, 1, sc_desc,&
            & s_sca, 1, 1, sc_desc, 0.d0, tmp, 1, 1, sc_desc)
    call pdgemm('N', 'T', n_size, n_size, n_size, 1.d0, tmp, 1, 1, sc_desc,&
            & U, 1, 1, sc_desc, 0.d0, X, 1, 1, sc_desc)

end subroutine lowdin_elpa



  !****s* linalg/get_inverse_real
  !  NAME
  !    get_inverse_real
  !  SYNOPSIS
  subroutine get_rinverse_lapack(M, n)
    !  PURPOSE
    !    Calculate the inverse of a real matrix, using the lapack
    !    functions DGETRF and DGETRI
    !    
    !    1.) LU factorization A = P * L * U
    !    2.) Calculate inverse with above factorization
    !           inv(A)*L = inv(U) for inv(A).
    !  USES
    use localorb_io, only: use_unit
    implicit none
    !  ARGUMENTS
    real*8, dimension(:,:), intent(INOUT) :: M 
    integer, intent(IN) :: n
    !  INPUTS
    !    o M -- matrix
    !    o n_size -- dimension of n x n matrix 
    !  OUTPUTS
    !    o M -- inverse of matrix M 
    !  AUTHOR
    !    Christoph Schober 
    !  HISTORY
    !    initial version, 2015

    !  SOURCE 
    integer,    dimension(n, n) :: ipiv_ev
    integer                     :: info
    real*8,     dimension(n) :: work

    ipiv_ev = 0 
    work = 0.d0

    call DGETRF(n, n, M, n, ipiv_ev, info)
    if(info /= 0) then
        write(use_unit,*) 'ERROR: Matrix is singular!'
        stop
    end if

    call DGETRI(n, M, n, ipiv_ev, work, n, info)
    if(info /= 0) then
        write(use_unit,*) 'ERROR: Matrix is singular!'
        stop
    end if


end subroutine get_rinverse_lapack

  !****s* linalg/get_rinverse_elpa
  !  NAME
  !    get_rinverse_elpa
  !  SYNOPSIS
  subroutine get_rinverse_elpa(M, n_size)
    !  PURPOSE
    !    Calculate the inverse of a real (scaLAPACK distributed) matrix,
    !    using the Cholesky decomposition.
    ! 
    !  USES
    use scalapack_wrapper, only: mxld, mxcol, nb, mpi_comm_cols, mpi_comm_rows, set_full_matrix_real, sc_desc
    use mpi_tasks, only: check_allocation
    use elpa1_2013, only: cholesky_real_2013, invert_trm_real_2013

    implicit none
    !  ARGUMENTS

    real*8, dimension(:,:), intent(INOUT) :: M 
    integer, intent(IN) :: n_size
    !  INPUTS
    !    o M -- matrix
    !    o n_size -- dimension of n x n matrix
    !                (un-distributed, e.g. for ovlp = n_basis)
    !  OUTPUTS
    !    o M -- inverse of matrix M 
    !  AUTHOR
    !    Christoph Schober 
    !  HISTORY
    !    initial version, 2015
    !
    !  SOURCE 
    real*8,     dimension(mxld, mxcol) :: sca_work

    call cholesky_real_2013(n_size, M, mxld, nb, mpi_comm_rows, mpi_comm_cols)
    call invert_trm_real_2013(n_size, M, mxld, nb, mpi_comm_rows, mpi_comm_cols)

    sca_work = M
    call pdgemm('N','T', n_size, n_size, n_size,1.d0,M,1,1,sc_desc, &
        M,1,1,sc_desc,0.d0,sca_work,1,1,sc_desc)
    M = sca_work

end subroutine get_rinverse_elpa

  !****s* linalg/simple_matmul_real_scalapack
  !  NAME
  !    simple_matmul_real_scalapack
  !  SYNOPSIS
  subroutine simple_matmul_real_scalapack(TRANSA, TRANSB, A, B, C, M, N, K)
    !  PURPOSE
    !    Do a matrix matrix multiplication using the scaLAPACK pdgemm
    !    for the simplest possible case of two 2D matrices.
    !    
    !    C = A*B
    !    
    !    The dimensions of the arrays M, N, K need to be specified.
    !
    !  USES
    use scalapack_wrapper, only: sc_desc
    implicit none
    !  ARGUMENTS
    character*1, intent(IN) :: TRANSA, TRANSB 
    real*8, dimension(:,:), intent(IN) :: A, B 
    real*8, dimension(:,:), intent(INOUT) :: C 
    !  INPUTS
    !    o TRANSA -- form of matrix ['N': M=M, 'T': M=transpose(M)]
    !    o TRANSB -- form of matrix ['N': M=M, 'T': M=transpose(M)]
    !       (see dgemm or pdgemm reference for more information)
    !    o A -- matrix A
    !    o B -- matrix B 
    !    o M -- num of rows (A + C)
    !    o N -- num of cols (B + C)
    !    o K -- num of cols (A) or number of rows (B) 
    !  OUTPUTS
    !    o C -- matrix product of A*B 
    !  AUTHOR
    !    Christoph Schober 
    !  HISTORY
    !    initial version, 2015

    !  SOURCE 
    real*8 :: ALPHA, BETA
    integer :: M, N, K, LDA, LDB, LDC

    ALPHA = 1.d0
    BETA = 0.d0
   
    call pdgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, 1, 1, sc_desc, &
                B, 1, 1, sc_desc, BETA, C, 1, 1, sc_desc)

    end subroutine simple_matmul_real_scalapack

  !****s* linalg/simple_matmul_real_lapack
  !  NAME
  !    simple_matmul_real
  !  SYNOPSIS
  subroutine simple_matmul_real_lapack(TRANSA, TRANSB, A, B, C)
    !  PURPOSE
    !    Do a matrix matrix multiplication using the LAPACK dgemm
    !    for the simplest possible case of two 2D matrices.
    !    
    !    C = A*B
    !    
    !    All dimensions and bounds are determined automatically from the
    !    input arrays.
    !
    !  USES
    implicit none
    !  ARGUMENTS
    character*1, intent(IN) :: TRANSA, TRANSB 
    real*8, dimension(:,:), intent(IN) :: A, B 
    real*8, dimension(:,:), intent(INOUT) :: C 
    !  INPUTS
    !    o TRANSA -- form of matrix ['N': M=M, 'T': M=transpose(M)]
    !    o TRANSB -- form of matrix ['N': M=M, 'T': M=transpose(M)]
    !       (see dgemm or pdgemm reference for more information)
    !    o A -- matrix A
    !    o B -- matrix B 
    !  OUTPUTS
    !    o C -- matrix product of A*B 
    !  AUTHOR
    !    Christoph Schober 
    !  HISTORY
    !    initial version, 2015

    !  SOURCE 
    real*8 :: ALPHA, BETA
    integer :: M, N, K, LDA, LDB, LDC

    ALPHA = 1.d0
    BETA = 0.d0
   
    ! Automatic determination of dimensions
    ! 
    ! M = num of rows (A + C)
    ! N = num of cols (B + C)
    ! K = num of cols (A), number of rows (B) 

    M = size(A, 1)
    N = size(B, 2)
    K = size(A, 2)

    ! If input array is transposed ('T'), adjust dimensions accordingly
    if (TRANSA.eq.'T') then
        M = size(A, 2)
        K = size(A, 1)
    elseif (TRANSB.eq.'T') then
        N = size(B, 1)
        K = size(B, 2)
    endif

    LDA = ubound(A,1)
    LDB = ubound(B,1)
    LDC = ubound(C,1)

    call dgemm(TRANSA, TRANSB, M, N, K, ALPHA, &
               A, LDA, B, LDB, BETA, C, LDC)

    end subroutine simple_matmul_real_lapack

  !****s* linalg/simple_matmul_cmplx_scalapack
  !  NAME
  !    simple_matmul_cmplx_scalapack
  !  SYNOPSIS
  subroutine simple_matmul_cmplx_scalapack(TRANSA, TRANSB, A, B, C, M, N, K)
    !  PURPOSE
    !    Do a matrix matrix multiplication using the scaLAPACK pdgemm
    !    for the simplest possible case of two 2D matrices.
    !    
    !    C = A*B
    !    
    !    The dimensions of the arrays M, N, K need to be specified.
    !
    !  USES
    use scalapack_wrapper, only: sc_desc
    implicit none
    !  ARGUMENTS
    character*1, intent(IN) :: TRANSA, TRANSB 
    complex*16, dimension(:,:), intent(IN) :: A, B 
    complex*16, dimension(:,:), intent(INOUT) :: C 
    !  INPUTS
    !    o TRANSA -- form of matrix ['N': M=M, 'T': M=transpose(M)]
    !    o TRANSB -- form of matrix ['N': M=M, 'T': M=transpose(M)]
    !       (see dgemm or pdgemm reference for more information)
    !    o A -- matrix A
    !    o B -- matrix B 
    !    o M -- num of rows (A + C)
    !    o N -- num of cols (B + C)
    !    o K -- num of cols (A) or number of rows (B) 
    !  OUTPUTS
    !    o C -- matrix product of A*B 
    !  AUTHOR
    !    Christoph Schober 
    !  HISTORY
    !    initial version, 2015

    !  SOURCE 
    complex*16 :: ALPHA, BETA
    integer :: M, N, K, LDA, LDB, LDC

    ALPHA = (1.d0,0.d0)
    BETA = (0.d0,0.d0)
   
    call pzgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, 1, 1, sc_desc, &
                B, 1, 1, sc_desc, BETA, C, 1, 1, sc_desc)

    end subroutine simple_matmul_cmplx_scalapack

  !****s* linalg/simple_matmul_cmplx_lapack
  !  NAME
  !    simple_matmul_cmplx
  !  SYNOPSIS
  subroutine simple_matmul_cmplx_lapack(TRANSA, TRANSB, A, B, C)
    !  PURPOSE
    !    Do a matrix matrix multiplication using the LAPACK dgemm
    !    for the simplest possible case of two 2D matrices.
    !    
    !    C = A*B
    !    
    !    All dimensions and bounds are determined automatically from the
    !    input arrays.
    !
    !  USES
    implicit none
    !  ARGUMENTS
    character*1, intent(IN) :: TRANSA, TRANSB 
    complex*16, dimension(:,:), intent(IN) :: A, B 
    complex*16, dimension(:,:), intent(INOUT) :: C 
    !  INPUTS
    !    o TRANSA -- form of matrix ['N': M=M, 'T': M=transpose(M)]
    !    o TRANSB -- form of matrix ['N': M=M, 'T': M=transpose(M)]
    !       (see dgemm or pdgemm reference for more information)
    !    o A -- matrix A
    !    o B -- matrix B 
    !  OUTPUTS
    !    o C -- matrix product of A*B 
    !  AUTHOR
    !    Christoph Schober 
    !  HISTORY
    !    initial version, 2015

    !  SOURCE 
    complex*16 :: ALPHA, BETA
    integer :: M, N, K, LDA, LDB, LDC

    ALPHA = (1.d0,0.d0)
    BETA = (0.d0,0.d0)
   
    ! Automatic determination of dimensions
    ! 
    ! M = num of rows (A + C)
    ! N = num of cols (B + C)
    ! K = num of cols (A), number of rows (B) 

    M = size(A, 1)
    N = size(B, 2)
    K = size(A, 2)

    ! If input array is transposed ('T'), adjust dimensions accordingly
    if (TRANSA.eq.'T'.or.TRANSA.eq.'C') then
        M = size(A, 2)
        K = size(A, 1)
    elseif (TRANSB.eq.'T'.or.TRANSB.eq.'C') then
        N = size(B, 1)
        K = size(B, 2)
    endif

    LDA = ubound(A,1)
    LDB = ubound(B,1)
    LDC = ubound(C,1)

    call zgemm(TRANSA, TRANSB, M, N, K, ALPHA, &
               A, LDA, B, LDB, BETA, C, LDC)

    end subroutine simple_matmul_cmplx_lapack

end module linalg
!******
